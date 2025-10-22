import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
import re

# Read the input files
print("Loading data files...")
gene_data = pd.read_csv('input/log2_center_normalization_t_test_Human.csv')
protein_atlas = pd.read_csv('input/proteinatlas.tsv', sep='\t')

# Extract Gene column from gene_data (it's in the "Gene" column)
print(f"Total genes in input file: {len(gene_data)}")
print(f"Total genes in protein atlas: {len(protein_atlas)}")

# Merge the two dataframes on Gene column
merged_data = gene_data.merge(
    protein_atlas[['Gene', 'RNA tissue specific nTPM', 'Tissue expression cluster']],
    on='Gene',
    how='left'
)

print(f"Total merged records: {len(merged_data)}")

# Function to extract liver nTPM value
def extract_liver_ntpm(text):
    if pd.isna(text) or text == '':
        return None
    # Search for pattern "liver: number" or "liver:number"
    pattern = r'liver:\s*(\d+\.?\d*)'
    match = re.search(pattern, str(text).lower())
    if match:
        return float(match.group(1))
    return None

# Classification based on the criteria
# 1. liver protein_1: "liver" in "RNA tissue specific nTPM" column
# 2. liver protein_2: "liver" in "Tissue expression cluster" column
# 3. liver protein: "liver" in both columns

# Extract liver nTPM values
merged_data['liver_nTPM_value'] = merged_data['RNA tissue specific nTPM'].apply(extract_liver_ntpm)

merged_data['has_liver_nTPM'] = merged_data['RNA tissue specific nTPM'].fillna('').str.lower().str.contains('liver')
merged_data['has_liver_cluster'] = merged_data['Tissue expression cluster'].fillna('').str.lower().str.contains('liver')

# Create classification column
def classify_liver_protein(row):
    if row['has_liver_nTPM'] and row['has_liver_cluster']:
        return 'liver protein'
    elif row['has_liver_nTPM']:
        return 'liver protein_1'
    elif row['has_liver_cluster']:
        return 'liver protein_2'
    else:
        return 'non-liver protein'

merged_data['Classification'] = merged_data.apply(classify_liver_protein, axis=1)

# Count each category
classification_counts = merged_data['Classification'].value_counts()
print("\nClassification Results:")
print(classification_counts)

# Print statistics about extracted liver nTPM values
liver_ntpm_stats = merged_data[merged_data['liver_nTPM_value'].notna()]['liver_nTPM_value']
if len(liver_ntpm_stats) > 0:
    print(f"\nLiver nTPM Value Statistics:")
    print(f"Total genes with liver nTPM values: {len(liver_ntpm_stats)}")
    print(f"Mean: {liver_ntpm_stats.mean():.2f}")
    print(f"Median: {liver_ntpm_stats.median():.2f}")
    print(f"Min: {liver_ntpm_stats.min():.2f}")
    print(f"Max: {liver_ntpm_stats.max():.2f}")

# Get gene lists for Venn diagram
all_genes = set(merged_data['Gene'].tolist())
liver_protein_1_genes = set(merged_data[merged_data['has_liver_nTPM']]['Gene'].tolist())
liver_protein_2_genes = set(merged_data[merged_data['has_liver_cluster']]['Gene'].tolist())
non_liver_genes = all_genes - liver_protein_1_genes - liver_protein_2_genes
liver_protein_both = liver_protein_1_genes & liver_protein_2_genes

total_genes = len(all_genes)
liver_candidates = len(liver_protein_1_genes | liver_protein_2_genes)

print(f"\nVenn Diagram Statistics:")
print(f"Total genes: {total_genes}")
print(f"Liver protein candidates (union): {liver_candidates}")
print(f"Non-liver proteins: {len(non_liver_genes)}")
print(f"Genes with liver in RNA tissue specific nTPM (liver protein_1 total): {len(liver_protein_1_genes)}")
print(f"Genes with liver in Tissue expression cluster (liver protein_2 total): {len(liver_protein_2_genes)}")
print(f"Genes with liver in both (liver protein): {len(liver_protein_both)}")
print(f"Genes only in RNA tissue specific nTPM: {len(liver_protein_1_genes - liver_protein_2_genes)}")
print(f"Genes only in Tissue expression cluster: {len(liver_protein_2_genes - liver_protein_1_genes)}")

# Create improved Venn diagram
fig = plt.figure(figsize=(14, 12))
ax = fig.add_subplot(111)

venn = venn3(
    [liver_protein_1_genes, liver_protein_2_genes, non_liver_genes],
    set_labels=('', '', ''),  # We'll add custom labels
    set_colors=('#FFB6C1', '#87CEEB', '#E8E8E8'),
    alpha=0.65,
    ax=ax
)

# Add title with statistics
title_text = f'Liver Protein Classification'
subtitle_text = f'Total: {total_genes} proteins | Liver candidates: {liver_candidates} ({liver_candidates/total_genes*100:.1f}%)'
plt.title(title_text, fontsize=18, fontweight='bold', pad=25, y=0.98)
plt.text(0.5, 0.93, subtitle_text, ha='center', fontsize=13, transform=fig.transFigure)

# Customize the subset labels with better positioning
if venn.get_label_by_id('100'):
    venn.get_label_by_id('100').set_text(f"{len(liver_protein_1_genes - liver_protein_2_genes)}")
    venn.get_label_by_id('100').set_fontsize(14)
    venn.get_label_by_id('100').set_fontweight('bold')
    venn.get_label_by_id('100').set_color('#8B0000')

if venn.get_label_by_id('010'):
    venn.get_label_by_id('010').set_text(f"{len(liver_protein_2_genes - liver_protein_1_genes)}")
    venn.get_label_by_id('010').set_fontsize(14)
    venn.get_label_by_id('010').set_fontweight('bold')
    venn.get_label_by_id('010').set_color('#00008B')

if venn.get_label_by_id('110'):
    venn.get_label_by_id('110').set_text(f"{len(liver_protein_both)}")
    venn.get_label_by_id('110').set_fontsize(14)
    venn.get_label_by_id('110').set_fontweight('bold')
    venn.get_label_by_id('110').set_color('#006400')

if venn.get_label_by_id('001'):
    venn.get_label_by_id('001').set_text(f"{len(non_liver_genes)}")
    venn.get_label_by_id('001').set_fontsize(14)
    venn.get_label_by_id('001').set_fontweight('bold')
    venn.get_label_by_id('001').set_color('#404040')

# Add custom set labels with better positioning
# RNA tissue specific nTPM label (left circle)
ax.text(-0.55, 0.35, 'RNA tissue\nspecific nTPM', ha='center', va='center',
        fontsize=12, fontweight='bold', color='#8B0000',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='#FFB6C1', linewidth=2))

# Tissue expression cluster label (right circle)
ax.text(0.55, 0.35, 'Tissue expression\ncluster', ha='center', va='center',
        fontsize=12, fontweight='bold', color='#00008B',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='#87CEEB', linewidth=2))

# Non-liver protein label (bottom circle)
ax.text(0, -0.75, 'Non-liver protein', ha='center', va='center',
        fontsize=12, fontweight='bold', color='#404040',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='#E8E8E8', linewidth=2))

# Add legend with detailed information
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#FFB6C1', edgecolor='#8B0000', linewidth=2,
          label=f'RNA tissue specific nTPM only: {len(liver_protein_1_genes - liver_protein_2_genes)}', alpha=0.65),
    Patch(facecolor='#87CEEB', edgecolor='#00008B', linewidth=2,
          label=f'Tissue expression cluster only: {len(liver_protein_2_genes - liver_protein_1_genes)}', alpha=0.65),
    Patch(facecolor='#98FB98', edgecolor='#006400', linewidth=2,
          label=f'Both conditions met: {len(liver_protein_both)}', alpha=0.65),
    Patch(facecolor='#E8E8E8', edgecolor='#404040', linewidth=2,
          label=f'Non-liver protein: {len(non_liver_genes)}', alpha=0.65)
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=11, frameon=True,
          fancybox=True, shadow=True, bbox_to_anchor=(0.02, 0.98))

# Save the Venn diagram
output_dir = 'output'
os.makedirs(output_dir, exist_ok=True)
plt.savefig(f'{output_dir}/liver_protein_venn_diagram.png', dpi=300, bbox_inches='tight')
print(f"\nVenn diagram saved to: {output_dir}/liver_protein_venn_diagram.png")

# Create Trace folder and save tracing data
trace_dir = 'Trace'
os.makedirs(trace_dir, exist_ok=True)

# Prepare tracing data with all relevant columns
trace_data = merged_data[[
    'Gene',
    'RNA tissue specific nTPM',
    'liver_nTPM_value',
    'Tissue expression cluster',
    'has_liver_nTPM',
    'has_liver_cluster',
    'Classification'
]].copy()

# Save tracing data
trace_file = f'{trace_dir}/Data tracing.csv'
trace_data.to_csv(trace_file, index=False)
print(f"Tracing data saved to: {trace_file}")

# Save summary statistics to Trace folder
summary_file = f'{trace_dir}/classification_summary.csv'
classification_counts.to_csv(summary_file, header=['Count'])
print(f"Summary statistics saved to: {summary_file}")

# Save classified gene lists with liver nTPM values to Trace folder
for classification in ['liver protein', 'liver protein_1', 'liver protein_2']:
    subset = merged_data[merged_data['Classification'] == classification]
    if len(subset) > 0:
        genes_file = f'{trace_dir}/{classification.replace(" ", "_")}_genes.csv'
        subset[['Gene', 'liver_nTPM_value']].to_csv(genes_file, index=False)
        print(f"{classification} genes saved to: {genes_file}")

print("\nClassification complete!")
