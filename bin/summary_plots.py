#!/usr/bin/env python

import sys
import pandas as pd
import matplotlib.pyplot as plt

n_cells = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
n_singlets = pd.read_csv(sys.argv[2], sep='\t', index_col=0)

df = n_cells.join(n_singlets, how='inner', validate='1:1').sort_index(ascending=False)

# Summary plot
df1 = df[['n_singlets', 'n_cells_HQ', 'n_cells', 'n_droplets']]
df1[['n_droplets', 'n_cells', 'n_cells_HQ', 'n_singlets']].to_csv('filtering_summary.tsv', sep='\t')
fig, ax = plt.subplots(figsize=(8, max(3, df1.shape[0])))
df1[['n_singlets', 'n_cells_HQ', 'n_cells']].plot(kind='barh', stacked=False, ax=ax)
ax.set_xlabel('# cells')
ax.set_ylabel(None)
ax.set_yticklabels([f'{item.get_text()}\n{df1.loc[item.get_text(), "n_droplets"]:,} droplets' for item in ax.get_yticklabels()])
ax.set_title('Filtering Summary')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[::-1], labels=['Non-Empty', 'HQ Cells', 'HQ Singlets'], loc='best')
plt.tight_layout()
plt.savefig('filtering_summary.png')
plt.close()

# Comparison to metadata plot
df2 = df[['n_shared_donor', 'n_shared', 'n_singlets', 'n_meta_cells']]
df2[['n_meta_cells', 'n_singlets', 'n_shared', 'n_shared_donor']].to_csv('comparison_to_metadata.tsv', sep='\t')
fig, ax = plt.subplots(figsize=(8, max(3, df2.shape[0])))
df2.plot(kind='barh', stacked=False, ax=ax)
ax.set_xlabel('# cells')
ax.set_ylabel(None)
ax.set_title('Comparison to Metadata')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[::-1], labels=['Metadata', 'My HQ Singlets', 'Shared', 'Shared + same donor'], loc='best')
plt.tight_layout()
plt.savefig('comparison_to_metadata.png')
plt.close()
