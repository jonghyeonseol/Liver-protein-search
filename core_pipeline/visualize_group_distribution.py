import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Set style for better aesthetics
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# Read the trace data
print("Loading trace data...")
trace_data = pd.read_csv('Trace/Data tracing.csv')

# Filter genes with liver nTPM values
liver_genes = trace_data[trace_data['liver_nTPM_value'].notna()].copy()

print(f"Total genes with liver nTPM values: {len(liver_genes)}")

# Group by classification
groups = liver_genes.groupby('Classification')['liver_nTPM_value'].apply(list).to_dict()

# Color mapping for new 4-category system
color_map = {
    # High nTPM (>=100) categories
    'liver protein (nTPM≥100 + cluster + enrichment)': '#FF0000',
    'liver protein (nTPM≥100 + cluster)': '#FF6B6B',
    'liver protein (nTPM≥100 + enrichment)': '#FF8C42',
    'liver protein (nTPM≥100 only)': '#FFB6B6',

    # Low nTPM (<100) categories
    'liver protein (nTPM<100 + cluster + enrichment)': '#4169E1',
    'liver protein (nTPM<100 + cluster)': '#87CEEB',
    'liver protein (nTPM<100 + enrichment)': '#ADD8E6',
    'liver protein (nTPM<100 only)': '#E0F0FF',

    # Both nTPM levels
    'liver protein (both nTPM levels)': '#9370DB',
    'liver protein (both nTPM + cluster)': '#8B008B',
    'liver protein (both nTPM + enrichment)': '#BA55D3',
    'liver protein (all 4)': '#4B0082',

    # No nTPM categories
    'liver protein (cluster + enrichment)': '#A8E6CF',
    'liver protein (cluster only)': '#95E1D3',
    'liver protein (enrichment only)': '#C7CEEA'
}

# Order of classifications for display (only those with nTPM values)
classification_order = [
    'liver protein (nTPM≥100 + cluster + enrichment)',
    'liver protein (nTPM≥100 + cluster)',
    'liver protein (nTPM≥100 + enrichment)',
    'liver protein (nTPM≥100 only)',
    'liver protein (nTPM<100 + cluster + enrichment)',
    'liver protein (nTPM<100 + cluster)',
    'liver protein (nTPM<100 + enrichment)',
    'liver protein (nTPM<100 only)'
]

# Filter groups that have nTPM values
groups_with_values = {k: v for k, v in groups.items() if k in classification_order}

# Create figure with multiple subplots - improved layout
fig = plt.figure(figsize=(22, 16), facecolor='white')
gs = fig.add_gridspec(3, 2, hspace=0.40, wspace=0.35, top=0.94, bottom=0.06, left=0.08, right=0.96)

# Add overall title
fig.suptitle('Liver Protein nTPM Distribution Analysis by Classification Group',
             fontsize=20, fontweight='bold', y=0.98)

# 1. Box Plot (Top Left)
ax1 = fig.add_subplot(gs[0, 0])
data_for_box = [groups_with_values[cat] for cat in classification_order if cat in groups_with_values]
# Shorten labels for better display
labels_for_box = []
for cat in classification_order:
    if cat in groups_with_values:
        label = cat.replace('liver protein ', '').replace('(', '').replace(')', '')
        label = label.replace('nTPM≥100', '≥100').replace('nTPM<100', '<100')
        label = label.replace('cluster', 'C').replace('enrichment', 'E')
        label = label.replace(' + ', '+')
        labels_for_box.append(label)

bp = ax1.boxplot(data_for_box, labels=labels_for_box, patch_artist=True,
                 showmeans=True, meanline=True,
                 boxprops=dict(linewidth=2),
                 medianprops=dict(linewidth=2.5, color='red'),
                 meanprops=dict(linewidth=2.5, color='blue', linestyle='--'),
                 whiskerprops=dict(linewidth=1.5),
                 capprops=dict(linewidth=1.5),
                 flierprops=dict(marker='o', markersize=4, alpha=0.5))

# Color the boxes
for patch, cat in zip(bp['boxes'], [c for c in classification_order if c in groups_with_values]):
    patch.set_facecolor(color_map[cat])
    patch.set_alpha(0.7)

ax1.set_yscale('log')
ax1.set_ylabel('Liver nTPM Value (log scale)', fontsize=13, fontweight='bold')
ax1.set_title('Box Plot: nTPM Distribution by Group', fontsize=15, fontweight='bold', pad=20)
ax1.grid(axis='y', alpha=0.25, linestyle='--', linewidth=1)
ax1.set_facecolor('#F8F9FA')
plt.setp(ax1.xaxis.get_majorticklabels(), rotation=20, ha='right', fontsize=11)

# Add legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='red', linewidth=2.5, label='Median'),
    Line2D([0], [0], color='blue', linewidth=2.5, linestyle='--', label='Mean')
]
ax1.legend(handles=legend_elements, loc='upper right', fontsize=9)

# 2. Violin Plot (Top Right)
ax2 = fig.add_subplot(gs[0, 1])

# Prepare data for violin plot
violin_data = []
violin_labels = []
violin_colors = []

for cat in classification_order:
    if cat in groups_with_values:
        violin_data.append(groups_with_values[cat])
        label = cat.replace('liver protein ', '').replace('(', '').replace(')', '')
        label = label.replace('nTPM≥100', '≥100').replace('nTPM<100', '<100')
        label = label.replace('cluster', 'C').replace('enrichment', 'E')
        label = label.replace(' + ', '+')
        violin_labels.append(label)
        violin_colors.append(color_map[cat])

parts = ax2.violinplot(violin_data, positions=range(len(violin_data)),
                       showmeans=True, showmedians=True, widths=0.7)

# Color the violins
for pc, color in zip(parts['bodies'], violin_colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)
    pc.set_edgecolor('black')
    pc.set_linewidth(1.5)

# Style the components
parts['cmeans'].set_color('blue')
parts['cmeans'].set_linewidth(2)
parts['cmedians'].set_color('red')
parts['cmedians'].set_linewidth(2)

ax2.set_xticks(range(len(violin_labels)))
ax2.set_xticklabels(violin_labels, rotation=20, ha='right', fontsize=11)
ax2.set_yscale('log')
ax2.set_ylabel('Liver nTPM Value (log scale)', fontsize=13, fontweight='bold')
ax2.set_title('Violin Plot: Distribution Shape & Density', fontsize=15, fontweight='bold', pad=20)
ax2.grid(axis='y', alpha=0.25, linestyle='--', linewidth=1)
ax2.set_facecolor('#F8F9FA')

# 3. Strip Plot with Box (Middle Left)
ax3 = fig.add_subplot(gs[1, 0])

positions = []
all_values = []
all_colors = []

for i, cat in enumerate(classification_order):
    if cat in groups_with_values:
        values = groups_with_values[cat]
        positions.extend([i] * len(values))
        all_values.extend(values)
        all_colors.extend([color_map[cat]] * len(values))

# Add jitter to x positions
np.random.seed(42)
jittered_positions = np.array(positions) + np.random.normal(0, 0.1, len(positions))

ax3.scatter(jittered_positions, all_values, c=all_colors, alpha=0.5, s=30, edgecolors='black', linewidth=0.5)

# Overlay box plot
bp3 = ax3.boxplot(data_for_box, positions=range(len(data_for_box)),
                  widths=0.3, patch_artist=False,
                  boxprops=dict(linewidth=2, color='black'),
                  medianprops=dict(linewidth=2.5, color='red'),
                  whiskerprops=dict(linewidth=1.5, color='black'),
                  capprops=dict(linewidth=1.5, color='black'),
                  showfliers=False)

ax3.set_xticks(range(len(labels_for_box)))
ax3.set_xticklabels(labels_for_box, rotation=20, ha='right', fontsize=11)
ax3.set_yscale('log')
ax3.set_ylabel('Liver nTPM Value (log scale)', fontsize=13, fontweight='bold')
ax3.set_title('Strip Plot: Individual Gene Values', fontsize=15, fontweight='bold', pad=20)
ax3.grid(axis='y', alpha=0.25, linestyle='--', linewidth=1)
ax3.set_facecolor('#F8F9FA')

# 4. Histogram/Density Plot (Middle Right)
ax4 = fig.add_subplot(gs[1, 1])

for cat in classification_order:
    if cat in groups_with_values:
        values = np.array(groups_with_values[cat])
        log_values = np.log10(values)
        label = cat.replace('liver protein ', '').replace('(', '').replace(')', '')
        label = label.replace('nTPM≥100', '≥100').replace('nTPM<100', '<100')
        label = label.replace('cluster', 'C').replace('enrichment', 'E')
        label = label.replace(' + ', '+')
        ax4.hist(log_values, bins=20, alpha=0.5, label=label,
                color=color_map[cat], edgecolor='black', linewidth=1)

ax4.set_xlabel('Log10(Liver nTPM Value)', fontsize=13, fontweight='bold')
ax4.set_ylabel('Frequency', fontsize=13, fontweight='bold')
ax4.set_title('Histogram: nTPM Distribution Overlap', fontsize=15, fontweight='bold', pad=20)
ax4.legend(loc='upper right', fontsize=10.5, framealpha=0.95, edgecolor='black')
ax4.grid(axis='y', alpha=0.25, linestyle='--', linewidth=1)
ax4.set_facecolor('#F8F9FA')

# 5. Cumulative Distribution (Bottom Left)
ax5 = fig.add_subplot(gs[2, 0])

for cat in classification_order:
    if cat in groups_with_values:
        values = sorted(groups_with_values[cat])
        cumulative = np.arange(1, len(values) + 1) / len(values) * 100
        label = cat.replace('liver protein ', '').replace('(', '').replace(')', '')
        label = label.replace('nTPM≥100', '≥100').replace('nTPM<100', '<100')
        label = label.replace('cluster', 'C').replace('enrichment', 'E')
        label = label.replace(' + ', '+')
        ax5.plot(values, cumulative, linewidth=2.5, label=label,
                color=color_map[cat], marker='o', markersize=3, markevery=max(1, len(values)//20))

ax5.set_xscale('log')
ax5.set_xlabel('Liver nTPM Value (log scale)', fontsize=13, fontweight='bold')
ax5.set_ylabel('Cumulative Percentage (%)', fontsize=13, fontweight='bold')
ax5.set_title('Cumulative Distribution (CDF)', fontsize=15, fontweight='bold', pad=20)
ax5.legend(loc='lower right', fontsize=10.5, framealpha=0.95, edgecolor='black')
ax5.grid(True, alpha=0.25, linestyle='--', linewidth=1)
ax5.set_facecolor('#F8F9FA')
ax5.set_ylim(0, 105)

# 6. Summary Statistics Table (Bottom Right)
ax6 = fig.add_subplot(gs[2, 1])
ax6.axis('off')

# Calculate statistics
stats_data = []
for cat in classification_order:
    if cat in groups_with_values:
        values = np.array(groups_with_values[cat])
        label = cat.replace('liver protein ', '').replace('(', '').replace(')', '')
        label = label.replace('nTPM≥100', '≥100').replace('nTPM<100', '<100')
        label = label.replace('cluster', 'C').replace('enrichment', 'E')
        label = label.replace(' + ', '+')
        stats_data.append([
            label,
            len(values),
            f"{np.mean(values):.2f}",
            f"{np.median(values):.2f}",
            f"{np.std(values):.2f}",
            f"{np.min(values):.2f}",
            f"{np.max(values):.2f}",
            f"{np.percentile(values, 25):.2f}",
            f"{np.percentile(values, 75):.2f}"
        ])

# Create table
column_labels = ['Group', 'N', 'Mean', 'Median', 'Std Dev', 'Min', 'Max', 'Q1', 'Q3']
table = ax6.table(cellText=stats_data, colLabels=column_labels,
                  cellLoc='center', loc='center',
                  colWidths=[0.19, 0.09, 0.11, 0.11, 0.11, 0.10, 0.10, 0.09, 0.09])

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2.2)

# Style header
for i in range(len(column_labels)):
    cell = table[(0, i)]
    cell.set_facecolor('#4A90E2')
    cell.set_text_props(weight='bold', color='white', fontsize=11)

# Color rows by group
for i, cat in enumerate([c for c in classification_order if c in groups_with_values]):
    for j in range(len(column_labels)):
        cell = table[(i+1, j)]
        cell.set_facecolor(color_map[cat])
        cell.set_alpha(0.3)
        cell.set_text_props(fontsize=10)

ax6.set_title('Summary Statistics by Group', fontsize=15, fontweight='bold', pad=25)

# Save the comprehensive figure
plt.savefig('output/liver_ntpm_group_distribution.png', dpi=300, bbox_inches='tight')
print(f"\nComprehensive group distribution saved to: output/liver_ntpm_group_distribution.png")

# Print summary
print("\n=== Summary Statistics ===")
for cat in classification_order:
    if cat in groups_with_values:
        values = np.array(groups_with_values[cat])
        print(f"\n{cat}:")
        print(f"  N: {len(values)}")
        print(f"  Mean: {np.mean(values):.2f}")
        print(f"  Median: {np.median(values):.2f}")
        print(f"  Std Dev: {np.std(values):.2f}")
        print(f"  Range: {np.min(values):.2f} - {np.max(values):.2f}")
        print(f"  IQR: {np.percentile(values, 25):.2f} - {np.percentile(values, 75):.2f}")

plt.close()
print("\nVisualization complete!")
