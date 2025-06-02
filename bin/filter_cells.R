#!/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)
library(DropletUtils)
library(ggplot2)

source('/usr/local/bin/malat1_function.R')  # Path in the container

mito_filter <- function(sobj, mt_threshold, outdir) {
    # Calculate the percentage of mitochondrial genes
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

    # Open a PNG device to save the plot
    png(filename = paste0(outdir, "/mito_filter.png"), width = 800, height = 600)

    # Create the Violin Plot for visualization of percent.mt and the threshold
    p = VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
    p[[3]] = p[[3]] + 
        geom_hline(yintercept = mt_threshold, color = "red")

    # Show the plot and close the device to write the file
    print(p)
    dev.off()

    # Create a df which cells pass the threshold
    df = data.frame(
        cell_id = colnames(sobj),
        mito_percent = sobj[["percent.mt"]],
        mito_filter = ifelse(sobj[["percent.mt"]] < mt_threshold, 'PASS', 'FAIL')
    )
    colnames(df) = c('cell_id', 'mito_percent', 'mito_filter')
    row.names(df) <- NULL

    return(df)
}

malat1_filter <- function(sobj, outdir) {
    # Create a density plot of MALAT1 expression
    gene = "MALAT1"
    gene_expression = FetchData(sobj, vars = gene)[[gene]]
    malat1_threshold = define_malat1_threshold(counts = gene_expression)  # From the sourced file

    png(filename = paste0(outdir, '/malat1_filter.png'), width = 800, height = 600)

    hist(gene_expression, breaks = 100, main = gene, xlab = "Ln(counts_per_10k + 1)")
    abline(v = malat1_threshold, col = "red", lwd = 2)

    dev.off()

    # Creata a df which cells pass the threshold
    df = data.frame(
        cell_id = colnames(sobj),
        norm_malat1_expression = gene_expression,
        malat1_filter = ifelse(gene_expression > malat1_threshold, 'PASS', 'FAIL')
    )

    return(df)
}

save_data <- function(filtered_mtx_outdir, name, obj) {
    # Save the file as ...
    saveRDS(obj, paste0(filtered_mtx_outdir, "/", name, ".rds"))  # RDS
    SaveH5Seurat(obj, paste0(filtered_mtx_outdir, "/", name, ".h5seurat"), overwrite = TRUE)  # h5Seurat (needed as intermediary for h5ad)
    Convert(paste0(filtered_mtx_outdir, "/", name, ".h5seurat"), paste0(filtered_mtx_outdir, "/", name, ".h5ad"), overwrite = TRUE)  # h5ad
}

save_df <- function(df, outpath) {
    write.table(
        df, 
        file = outpath, 
        sep = '\t', 
        quote = FALSE, 
        row.names = FALSE
    )
}

filtered_cells_10X_dir = 'outs/filtered_feature_bc_matrix'
cli_args = commandArgs(trailingOnly = TRUE)
mt_threshold = as.numeric(cli_args[1])  # Percentage; mitochondrial content must be below this threshold
filtered_mtx_outdir = cli_args[2]
filter_summary_outpath = cli_args[3]
sample_id = cli_args[4]
n_cells_outpath = cli_args[5]
outdir = getwd()  # Output to the NF work dir

filtered_cells.data = Read10X(data.dir = filtered_cells_10X_dir, gene.column = 2)
filtered_cells = CreateSeuratObject(counts = filtered_cells.data, project = 'filtered_cells')

# Also read the gene_id column from the features.tsv file and add it to the Seurat object
filtered_cells[["RNA"]]@meta.features$gene_id = read.table(paste0(filtered_cells_10X_dir, '/features.tsv.gz'), header = FALSE, sep = '\t')[, 1]

mito_filtered = mito_filter(filtered_cells, mt_threshold, outdir)
head(mito_filtered)

# Normalize: scale to 10,000 reads per cell, then take log (base e) with pseudocounts
filtered_cells_norm = NormalizeData(filtered_cells, normalization.method = "LogNormalize", scale.factor = 1e4)

malat1_filtered = malat1_filter(filtered_cells_norm, outdir)
head(malat1_filtered)

my_filter_results = merge(mito_filtered, malat1_filtered, by = 'cell_id')
my_filter_results['is_high_quality'] = ifelse( (my_filter_results$mito_filter == 'PASS') & (my_filter_results$malat1_filter == 'PASS') , TRUE, FALSE)
head(my_filter_results)

# Assert that the number of cells in the filter results matches the number of cells in the Seurat object
if (nrow(my_filter_results) != ncol(filtered_cells_norm)) {
    stop("The number of cells in the filtering summary dataframe does not match the input Seurat object.")
}

save_df(my_filter_results, paste0(outdir, '/', filter_summary_outpath))

my_filtered_cells = subset(filtered_cells, cells = my_filter_results$cell_id[my_filter_results$is_high_quality])
my_filtered_cells
tail(my_filtered_cells[["RNA"]][[]])

# Number of droplets in the whole pool
num_lines = length(readLines(gzfile('outs/raw_feature_bc_matrix/barcodes.tsv.gz')))

save_df(
    df = data.frame(
        sample = sample_id, 
        n_droplets = num_lines, 
        n_cells = ncol(filtered_cells), 
        n_cells_HQ = ncol(my_filtered_cells)
    ), 
    outpath = paste0(outdir, '/', n_cells_outpath)
)

my_filtered_cells_norm = subset(filtered_cells_norm, cells = my_filter_results$cell_id[my_filter_results$is_high_quality])

# Overwrite the files in the folder that was copied to the NF work dir so that
# subsequent processes use the filtered countmtx but original bam files.
write10xCounts(
    path = filtered_cells_10X_dir, 
    x = GetAssayData(my_filtered_cells, slot = "counts"), 
    barcodes = colnames(my_filtered_cells), 
    gene.id = my_filtered_cells[["RNA"]]@meta.features$gene_id, 
    gene.symbol = rownames(my_filtered_cells), 
    overwrite = TRUE, 
    type = 'sparse', 
    version = '3'
)

save_data(filtered_mtx_outdir, "filtered_cells", my_filtered_cells)
save_data(filtered_mtx_outdir, "filtered_cells.norm", my_filtered_cells_norm)
