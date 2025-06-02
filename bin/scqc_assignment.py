#!/usr/bin/env python

import pandas as pd
import anndata as ad
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('--sample', type=str, required=True)
argparser.add_argument('--genotype_tools', type=str, required=True)
argparser.add_argument('--transcription_tools', type=str, required=True)
argparser.add_argument('--combined_results', type=str, required=True)
argparser.add_argument('--cell_meta', type=str, required=True)
argparser.add_argument('--n_singlets', type=str, required=True)
argparser.add_argument('--mtx_conversions', type=str, required=True)
argparser.add_argument('--outdir', type=str, required=True)

args = argparser.parse_args()

sample = args.sample
genotype_tools = args.genotype_tools.split(',')
transcription_tools = args.transcription_tools.split(',')
combined_results = args.combined_results
cell_meta = args.cell_meta
n_singlets = args.n_singlets
mtx_conversions = args.mtx_conversions
outdir = args.outdir

# To make a confident decision about a cell, we require a strict majority (>50%) of the tools to agree on the decision.
N_genotype_majority = len(genotype_tools) // 2 + 1
N_transcription_majority = len(transcription_tools) // 2 + 1

genotype_ind_assignment_colnames = [f'{genotype_tool}_Individual_Assignment' for genotype_tool in genotype_tools]
genotype_droplet_type_colnames = [f'{genotype_tool}_DropletType' for genotype_tool in genotype_tools]
transcription_droplet_type_colnames = [f'{transcription_tool}_DropletType' for transcription_tool in transcription_tools]

combined_results = pd.read_csv(combined_results, sep='\t', index_col='Barcode')

# The consensus droplettype is 'singlet' if the majority of the tools agree on 'singlet'; otherwise 'doublet'
combined_results['genetic_consensus_DropletType'] = combined_results[genotype_droplet_type_colnames].apply(lambda x: 'singlet' if (x == 'singlet').sum() >= N_genotype_majority else 'doublet', axis=1)
combined_results['transcriptional_consensus_DropletType'] = combined_results[transcription_droplet_type_colnames].apply(lambda x: 'singlet' if (x == 'singlet').sum() >= N_transcription_majority else 'doublet', axis=1)

# The droplet type is 'singlet' if both the genetic and transcriptional consensus are 'singlet'; otherwise 'doublet'
combined_results['consensus_DropletType'] = combined_results[['genetic_consensus_DropletType', 'transcriptional_consensus_DropletType']].apply(lambda x: 'singlet' if x[0] == 'singlet' and x[1] == 'singlet' else 'doublet', axis=1)

# The donor assignment is the donor ID if the majority of the genetic tools agree on the same donor (that is not 'doublet'); otherwise it is unassigned
combined_results['consensus_Individual_Assignment'] = combined_results[genotype_ind_assignment_colnames].apply(lambda x: x.mode()[0] if x[x != 'doublet'].value_counts().max() >= N_genotype_majority else 'unassigned', axis=1)

# Write the results to a file
combined_results.to_csv(f'{outdir}/combined_results_w_scqc_assignments.tsv', sep='\t')

# The true singlets are the cells that have 'singlet' as their consensus droplet type and are assigned to a donor, not 'unassigned'
n_total_cells = combined_results.shape[0]
true_HQ_singlets = combined_results[['consensus_DropletType', 'consensus_Individual_Assignment']].query('consensus_DropletType == "singlet" and consensus_Individual_Assignment != "unassigned"')
n_true_HQ_singlets = true_HQ_singlets.shape[0]
print(f'Number of successfully assigned singlets is {n_true_HQ_singlets} out of {n_total_cells} cells.')

# Write the true singlets barcodes to a file
true_HQ_singlets.index.to_series().to_csv(f'{outdir}/singlets_barcodes.tsv', index=False, header=False)

# Add the individual labels to the mtx files
for suffix in ['', '.norm']:
    adata = ad.read_h5ad(f'{mtx_conversions}/filtered_cells{suffix}.h5ad')
    adata = adata[true_HQ_singlets.index]
    adata.obs['individual'] = true_HQ_singlets['consensus_Individual_Assignment']
    adata.write_h5ad(f'{outdir}/labelled_singlets{suffix}.h5ad', compression='gzip')

# How many singlets are shared with the cell metadata and also share the same donor label?
cell_meta = pd.read_csv(cell_meta, sep='\t', index_col=0)
true_HQ_singlets = true_HQ_singlets.merge(cell_meta, how='inner', left_index=True, right_index=True, validate='1:1')
n_shared_singlets = true_HQ_singlets.shape[0]
n_shared_singlets_same_donor = (true_HQ_singlets['consensus_Individual_Assignment'] == true_HQ_singlets['individual']).sum()
print(f'Number of singlets shared with the cell metadata is {n_shared_singlets} and {n_shared_singlets_same_donor} of them also share the same donor label.')

# And save the number of singlets to a file
pd.DataFrame({
    'sample': [sample], 
    'n_meta_cells': [cell_meta.shape[0]],
    'n_singlets': [n_true_HQ_singlets], 
    'n_shared': [n_shared_singlets], 
    'n_shared_donor': [n_shared_singlets_same_donor]
}).to_csv(n_singlets, sep='\t', index=False)
