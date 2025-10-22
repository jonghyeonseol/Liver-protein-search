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

# Color by classification - improved contrast with new 4-category system
color_map = {
    # High nTPM (>=100) categories
    'liver protein (nTPM≥100 + cluster + enrichment)': '#FF0000',  # Bright red - highest priority
    'liver protein (nTPM≥100 + cluster)': '#FF6B6B',  # Light red
    'liver protein (nTPM≥100 + enrichment)': '#FF8C42',  # Orange-red
    'liver protein (nTPM≥100 only)': '#FFB6B6',  # Very light red

    # Low nTPM (<100) categories
    'liver protein (nTPM<100 + cluster + enrichment)': '#4169E1',  # Royal blue
    'liver protein (nTPM<100 + cluster)': '#87CEEB',  # Sky blue
    'liver protein (nTPM<100 + enrichment)': '#ADD8E6',  # Light blue
    'liver protein (nTPM<100 only)': '#E0F0FF',  # Very light blue

    # Both nTPM levels
    'liver protein (both nTPM levels)': '#9370DB',  # Medium purple
    'liver protein (both nTPM + cluster)': '#8B008B',  # Dark magenta
    'liver protein (both nTPM + enrichment)': '#BA55D3',  # Medium orchid
    'liver protein (all 4)': '#4B0082',  # Indigo - ultimate priority

    # No nTPM categories
    'liver protein (cluster + enrichment)': '#A8E6CF',  # Mint green
    'liver protein (cluster only)': '#95E1D3',  # Light teal
    'liver protein (enrichment only)': '#C7CEEA',  # Light periwinkle
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

# Add legend - show only categories that appear in top genes
from matplotlib.patches import Patch
unique_classifications = top_genes['Classification'].unique()
legend_elements = []

# High nTPM categories
if 'liver protein (nTPM≥100 + cluster + enrichment)' in unique_classifications:
    legend_elements.append(Patch(facecolor='#FF0000', edgecolor='black', label='nTPM≥100 + C + E', linewidth=1.5))
if 'liver protein (nTPM≥100 + cluster)' in unique_classifications:
    legend_elements.append(Patch(facecolor='#FF6B6B', edgecolor='black', label='nTPM≥100 + C', linewidth=1.5))
if 'liver protein (nTPM≥100 + enrichment)' in unique_classifications:
    legend_elements.append(Patch(facecolor='#FF8C42', edgecolor='black', label='nTPM≥100 + E', linewidth=1.5))
if 'liver protein (nTPM≥100 only)' in unique_classifications:
    legend_elements.append(Patch(facecolor='#FFB6B6', edgecolor='black', label='nTPM≥100 only', linewidth=1.5))

# Low nTPM categories
if 'liver protein (nTPM<100 + cluster + enrichment)' in unique_classifications:
    legend_elements.append(Patch(facecolor='#4169E1', edgecolor='black', label='nTPM<100 + C + E', linewidth=1.5))
if 'liver protein (nTPM<100 + cluster)' in unique_classifications:
    legend_elements.append(Patch(facecolor='#87CEEB', edgecolor='black', label='nTPM<100 + C', linewidth=1.5))
if 'liver protein (nTPM<100 + enrichment)' in unique_classifications:
    legend_elements.append(Patch(facecolor='#ADD8E6', edgecolor='black', label='nTPM<100 + E', linewidth=1.5))
if 'liver protein (nTPM<100 only)' in unique_classifications:
    legend_elements.append(Patch(facecolor='#E0F0FF', edgecolor='black', label='nTPM<100 only', linewidth=1.5))

ax1.legend(handles=legend_elements, loc='lower right', fontsize=8, ncol=2)

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
unique_scatter_classifications = liver_genes_sorted['Classification'].unique()
legend_elements_scatter = []

# High nTPM categories
if 'liver protein (nTPM≥100 + cluster + enrichment)' in unique_scatter_classifications:
    legend_elements_scatter.append(Patch(facecolor='#FF0000', edgecolor='black', label='nTPM≥100 + C + E', alpha=0.7, linewidth=1.5))
if 'liver protein (nTPM≥100 + cluster)' in unique_scatter_classifications:
    legend_elements_scatter.append(Patch(facecolor='#FF6B6B', edgecolor='black', label='nTPM≥100 + C', alpha=0.7, linewidth=1.5))
if 'liver protein (nTPM≥100 + enrichment)' in unique_scatter_classifications:
    legend_elements_scatter.append(Patch(facecolor='#FF8C42', edgecolor='black', label='nTPM≥100 + E', alpha=0.7, linewidth=1.5))
if 'liver protein (nTPM≥100 only)' in unique_scatter_classifications:
    legend_elements_scatter.append(Patch(facecolor='#FFB6B6', edgecolor='black', label='nTPM≥100 only', alpha=0.7, linewidth=1.5))

# Low nTPM categories
if 'liver protein (nTPM<100 + cluster + enrichment)' in unique_scatter_classifications:
    legend_elements_scatter.append(Patch(facecolor='#4169E1', edgecolor='black', label='nTPM<100 + C + E', alpha=0.7, linewidth=1.5))
if 'liver protein (nTPM<100 + cluster)' in unique_scatter_classifications:
    legend_elements_scatter.append(Patch(facecolor='#87CEEB', edgecolor='black', label='nTPM<100 + C', alpha=0.7, linewidth=1.5))
if 'liver protein (nTPM<100 + enrichment)' in unique_scatter_classifications:
    legend_elements_scatter.append(Patch(facecolor='#ADD8E6', edgecolor='black', label='nTPM<100 + E', alpha=0.7, linewidth=1.5))
if 'liver protein (nTPM<100 only)' in unique_scatter_classifications:
    legend_elements_scatter.append(Patch(facecolor='#E0F0FF', edgecolor='black', label='nTPM<100 only', alpha=0.7, linewidth=1.5))

ax2.legend(handles=legend_elements_scatter, loc='upper left', fontsize=9, ncol=2)

# Add statistics text
# Count each classification with nTPM values (high and low)
ntpm_high_ce = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM≥100 + cluster + enrichment)"])
ntpm_high_c = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM≥100 + cluster)"])
ntpm_high_e = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM≥100 + enrichment)"])
ntpm_high_only = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM≥100 only)"])

ntpm_low_ce = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM<100 + cluster + enrichment)"])
ntpm_low_c = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM<100 + cluster)"])
ntpm_low_e = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM<100 + enrichment)"])
ntpm_low_only = len(liver_genes[liver_genes["Classification"] == "liver protein (nTPM<100 only)"])

# Count those without nTPM values
cluster_only = len(trace_data[trace_data['Classification'] == 'liver protein (cluster only)'])
enrichment_only = len(trace_data[trace_data['Classification'] == 'liver protein (enrichment only)'])
cluster_enrichment = len(trace_data[trace_data['Classification'] == 'liver protein (cluster + enrichment)'])

stats_text = f'nTPM values: {len(liver_genes)}\n'
stats_text += f'nTPM≥100:\n'
stats_text += f'  +C+E: {ntpm_high_ce}\n'
stats_text += f'  +C: {ntpm_high_c}\n'
stats_text += f'  +E: {ntpm_high_e}\n'
stats_text += f'  only: {ntpm_high_only}\n'
stats_text += f'nTPM<100:\n'
stats_text += f'  +C+E: {ntpm_low_ce}\n'
stats_text += f'  +C: {ntpm_low_c}\n'
stats_text += f'  +E: {ntpm_low_e}\n'
stats_text += f'  only: {ntpm_low_only}\n'
stats_text += f'\nNo nTPM:\n'
stats_text += f'  C only: {cluster_only}\n'
stats_text += f'  E only: {enrichment_only}\n'
stats_text += f'  C+E: {cluster_enrichment}\n'
stats_text += f'\nStats:\n'
stats_text += f'Mean: {liver_genes["liver_nTPM_value"].mean():.2f}\n'
stats_text += f'Median: {liver_genes["liver_nTPM_value"].median():.2f}\n'
stats_text += f'Max: {liver_genes["liver_nTPM_value"].max():.2f}\n'
stats_text += f'Min: {liver_genes["liver_nTPM_value"].min():.2f}'
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
