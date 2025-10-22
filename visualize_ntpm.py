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

# Color by classification
colors = []
for classification in top_genes['Classification']:
    if classification == 'liver protein':
        colors.append('#4CAF50')  # Green for both
    elif classification == 'liver protein_1':
        colors.append('#FF9999')  # Pink for nTPM only
    elif classification == 'liver protein_2':
        colors.append('#99CCFF')  # Blue for cluster only
    else:
        colors.append('#CCCCCC')  # Gray for non-liver

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
    Patch(facecolor='#4CAF50', edgecolor='black', label='Liver protein (both)'),
    Patch(facecolor='#FF9999', edgecolor='black', label='Liver protein_1 (nTPM only)'),
    Patch(facecolor='#99CCFF', edgecolor='black', label='Liver protein_2 (cluster only)')
]
ax1.legend(handles=legend_elements, loc='lower right', fontsize=10)

# 2. Distribution plot (log scale)
ax2 = axes[1]
liver_genes_sorted = liver_genes.sort_values('liver_nTPM_value', ascending=True)
x_positions = range(len(liver_genes_sorted))
values = liver_genes_sorted['liver_nTPM_value'].values

# Color by classification
scatter_colors = []
for classification in liver_genes_sorted['Classification']:
    if classification == 'liver protein':
        scatter_colors.append('#4CAF50')
    elif classification == 'liver protein_1':
        scatter_colors.append('#FF9999')
    elif classification == 'liver protein_2':
        scatter_colors.append('#99CCFF')
    else:
        scatter_colors.append('#CCCCCC')

ax2.scatter(x_positions, values, c=scatter_colors, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
ax2.set_xlabel('Gene Index (sorted by nTPM value)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Liver nTPM Value (log scale)', fontsize=12, fontweight='bold')
ax2.set_yscale('log')
ax2.set_title('Distribution of Liver nTPM Values Across All Genes', fontsize=14, fontweight='bold', pad=15)
ax2.grid(True, alpha=0.3, linestyle='--')
ax2.legend(handles=legend_elements, loc='upper left', fontsize=10)

# Add statistics text
stats_text = f'Genes with nTPM values: {len(liver_genes)}\n'
stats_text += f'liver protein (both): {len(liver_genes[liver_genes["Classification"] == "liver protein"])}\n'
stats_text += f'liver protein_1 (nTPM only): {len(liver_genes[liver_genes["Classification"] == "liver protein_1"])}\n'
stats_text += f'liver protein_2 (cluster only): {liver_protein_2_count} (no nTPM)\n'
stats_text += f'\nMean: {liver_genes["liver_nTPM_value"].mean():.2f}\n'
stats_text += f'Median: {liver_genes["liver_nTPM_value"].median():.2f}\n'
stats_text += f'Max: {liver_genes["liver_nTPM_value"].max():.2f}\n'
stats_text += f'Min: {liver_genes["liver_nTPM_value"].min():.2f}'
ax2.text(0.98, 0.05, stats_text, transform=ax2.transAxes,
         fontsize=9, verticalalignment='bottom', horizontalalignment='right',
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
