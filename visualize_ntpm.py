import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the trace data
print("Loading trace data...")
trace_data = pd.read_csv('Trace/Data tracing.csv')

# Filter genes with liver nTPM values and sort by value
liver_genes = trace_data[trace_data['liver_nTPM_value'].notna()].copy()
liver_genes = liver_genes.sort_values('liver_nTPM_value', ascending=False)

# Check liver protein_2 group (no nTPM values expected)
liver_protein_2_count = len(trace_data[trace_data['Classification'] == 'liver protein_2'])

print(f"Total genes with liver nTPM values: {len(liver_genes)}")
print(f"Liver protein_2 genes (cluster only, no nTPM): {liver_protein_2_count}")
print(f"\nTop 10 genes with highest liver nTPM:")
print(liver_genes[['Gene', 'liver_nTPM_value', 'Classification']].head(10))

# 1. Create horizontal bar plot for top genes
fig, axes = plt.subplots(2, 1, figsize=(14, 16))

# Top 50 genes
top_n = 50
top_genes = liver_genes.head(top_n)

# Color by classification - improved contrast
color_map = {
    'liver protein (all 3)': '#FF6B6B',  # Bright red for all 3 - HIGHEST PRIORITY
    'liver protein (nTPM + cluster)': '#4ECDC4',  # Teal for nTPM + cluster
    'liver protein (nTPM + enrichment)': '#FFE66D',  # Bright yellow for nTPM + enrichment
    'liver protein (cluster + enrichment)': '#A8E6CF',  # Mint green for cluster + enrichment
    'liver protein_1': '#FF8C42',  # Bright orange for nTPM only
    'liver protein_2': '#95E1D3',  # Light teal for cluster only
    'liver protein_3': '#C7CEEA',  # Light periwinkle for enrichment only
    'non-liver protein': '#E8E8E8'  # Light gray for non-liver
}

colors = []
for classification in top_genes['Classification']:
    colors.append(color_map.get(classification, '#CCCCCC'))

# Plot top 50
ax1 = axes[0]
bars = ax1.barh(range(len(top_genes)), top_genes['liver_nTPM_value'], color=colors, edgecolor='black', linewidth=0.5)
ax1.set_yticks(range(len(top_genes)))
ax1.set_yticklabels(top_genes['Gene'], fontsize=9)
ax1.set_xlabel('Liver nTPM Value', fontsize=12, fontweight='bold')
ax1.set_title(f'Top {top_n} Genes by Liver nTPM Value', fontsize=14, fontweight='bold', pad=15)
ax1.invert_yaxis()
ax1.grid(axis='x', alpha=0.3, linestyle='--')

# Add value labels on bars
for i, (idx, row) in enumerate(top_genes.iterrows()):
    value = row['liver_nTPM_value']
    ax1.text(value + max(top_genes['liver_nTPM_value']) * 0.01, i, f'{value:.1f}',
             va='center', fontsize=8, fontweight='bold')

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#FF6B6B', edgecolor='black', label='All 3 criteria', linewidth=1.5),
    Patch(facecolor='#4ECDC4', edgecolor='black', label='nTPM + cluster', linewidth=1.5),
    Patch(facecolor='#FFE66D', edgecolor='black', label='nTPM + enrichment', linewidth=1.5),
    Patch(facecolor='#A8E6CF', edgecolor='black', label='cluster + enrichment', linewidth=1.5),
    Patch(facecolor='#FF8C42', edgecolor='black', label='nTPM only', linewidth=1.5),
    Patch(facecolor='#95E1D3', edgecolor='black', label='cluster only', linewidth=1.5),
    Patch(facecolor='#C7CEEA', edgecolor='black', label='enrichment only', linewidth=1.5)
]
ax1.legend(handles=legend_elements, loc='lower right', fontsize=9, ncol=2)

# 2. Distribution plot (log scale)
ax2 = axes[1]
liver_genes_sorted = liver_genes.sort_values('liver_nTPM_value', ascending=True)
x_positions = range(len(liver_genes_sorted))
values = liver_genes_sorted['liver_nTPM_value'].values

# Color by classification
scatter_colors = []
for classification in liver_genes_sorted['Classification']:
    scatter_colors.append(color_map.get(classification, '#CCCCCC'))

ax2.scatter(x_positions, values, c=scatter_colors, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
ax2.set_xlabel('Gene Index (sorted by nTPM value)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Liver nTPM Value (log scale)', fontsize=12, fontweight='bold')
ax2.set_yscale('log')
ax2.set_title('Distribution of Liver nTPM Values Across All Genes', fontsize=14, fontweight='bold', pad=15)
ax2.grid(True, alpha=0.3, linestyle='--')

# Create legend with updated classifications
legend_elements_scatter = [
    Patch(facecolor='#FF6B6B', edgecolor='black', label='All 3 criteria', alpha=0.7, linewidth=1.5),
    Patch(facecolor='#4ECDC4', edgecolor='black', label='nTPM + cluster', alpha=0.7, linewidth=1.5),
    Patch(facecolor='#FFE66D', edgecolor='black', label='nTPM + enrichment', alpha=0.7, linewidth=1.5),
    Patch(facecolor='#FF8C42', edgecolor='black', label='nTPM only', alpha=0.7, linewidth=1.5)
]
ax2.legend(handles=legend_elements_scatter, loc='upper left', fontsize=10)

# Add statistics text
# Count each classification
all_3 = len(liver_genes[liver_genes["Classification"] == "liver protein (all 3)"])
ntpm_cluster = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM + cluster)"])
ntpm_enrichment = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM + enrichment)"])
ntpm_only = len(liver_genes[liver_genes["Classification"] == "liver protein_1"])

# Count those without nTPM values
liver_protein_2_count = len(trace_data[trace_data['Classification'] == 'liver protein_2'])
liver_protein_3_count = len(trace_data[trace_data['Classification'] == 'liver protein_3'])
cluster_enrichment = len(trace_data[trace_data['Classification'] == 'liver protein (cluster + enrichment)'])

stats_text = f'Genes with nTPM values: {len(liver_genes)}\n'
stats_text += f'  All 3 criteria: {all_3}\n'
stats_text += f'  nTPM + cluster: {ntpm_cluster}\n'
stats_text += f'  nTPM + enrichment: {ntpm_enrichment}\n'
stats_text += f'  nTPM only: {ntpm_only}\n'
stats_text += f'\nGenes without nTPM:\n'
stats_text += f'  cluster only: {liver_protein_2_count}\n'
stats_text += f'  enrichment only: {liver_protein_3_count}\n'
stats_text += f'  cluster + enrichment: {cluster_enrichment}\n'
stats_text += f'\nStats:\n'
stats_text += f'  Mean: {liver_genes["liver_nTPM_value"].mean():.2f}\n'
stats_text += f'  Median: {liver_genes["liver_nTPM_value"].median():.2f}\n'
stats_text += f'  Max: {liver_genes["liver_nTPM_value"].max():.2f}\n'
stats_text += f'  Min: {liver_genes["liver_nTPM_value"].min():.2f}'
ax2.text(0.98, 0.05, stats_text, transform=ax2.transAxes,
         fontsize=8, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('output/liver_ntpm_distribution.png', dpi=300, bbox_inches='tight')
print(f"\nDistribution plot saved to: output/liver_ntpm_distribution.png")
plt.close()

# 3. Create a detailed table with all genes sorted by nTPM value
print("\nCreating detailed CSV file sorted by nTPM value...")
liver_genes_detailed = liver_genes[['Gene', 'liver_nTPM_value', 'Classification',
                                     'RNA tissue specific nTPM', 'Tissue expression cluster']].copy()
liver_genes_detailed = liver_genes_detailed.sort_values('liver_nTPM_value', ascending=False)
liver_genes_detailed.to_csv('Trace/liver_genes_sorted_by_ntpm.csv', index=False)
print(f"Detailed sorted gene list saved to: Trace/liver_genes_sorted_by_ntpm.csv")

print("\nAll visualizations completed successfully!")
print(f"\nGenerated files:")
print(f"1. output/liver_ntpm_distribution.png - Top 50 bar chart and distribution")
print(f"2. Trace/liver_genes_sorted_by_ntpm.csv - Complete sorted gene list")
