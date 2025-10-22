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
    protein_atlas[['Gene', 'RNA tissue specific nTPM', 'Tissue expression cluster', 'RNA tissue cell type enrichment']],
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
# 3. liver protein_3: "liver" in "RNA tissue cell type enrichment" column
# 4. liver protein: "liver" in all three columns (or combinations)

# Extract liver nTPM values
merged_data['liver_nTPM_value'] = merged_data['RNA tissue specific nTPM'].apply(extract_liver_ntpm)

merged_data['has_liver_nTPM'] = merged_data['RNA tissue specific nTPM'].fillna('').str.lower().str.contains('liver')
merged_data['has_liver_cluster'] = merged_data['Tissue expression cluster'].fillna('').str.lower().str.contains('liver')
merged_data['has_liver_enrichment'] = merged_data['RNA tissue cell type enrichment'].fillna('').str.lower().str.contains('liver')

# Create classification column with more detailed logic
def classify_liver_protein(row):
    has_ntpm = row['has_liver_nTPM']
    has_cluster = row['has_liver_cluster']
    has_enrichment = row['has_liver_enrichment']

    # Count how many conditions are met
    count = sum([has_ntpm, has_cluster, has_enrichment])

    if count == 0:
        return 'non-liver protein'
    elif count == 3:
        return 'liver protein (all 3)'
    elif count == 2:
        if has_ntpm and has_cluster:
            return 'liver protein (nTPM + cluster)'
        elif has_ntpm and has_enrichment:
            return 'liver protein (nTPM + enrichment)'
        elif has_cluster and has_enrichment:
            return 'liver protein (cluster + enrichment)'
    else:  # count == 1
        if has_ntpm:
            return 'liver protein_1'
        elif has_cluster:
            return 'liver protein_2'
        elif has_enrichment:
            return 'liver protein_3'

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
liver_protein_3_genes = set(merged_data[merged_data['has_liver_enrichment']]['Gene'].tolist())

total_genes = len(all_genes)
liver_candidates_set = liver_protein_1_genes | liver_protein_2_genes | liver_protein_3_genes
liver_candidates = len(liver_candidates_set)
non_liver_genes = all_genes - liver_candidates_set

print(f"\nVenn Diagram Statistics:")
print(f"Total genes: {total_genes}")
print(f"Liver protein candidates (union): {liver_candidates}")
print(f"Non-liver proteins: {len(non_liver_genes)}")
print(f"Genes with liver in RNA tissue specific nTPM: {len(liver_protein_1_genes)}")
print(f"Genes with liver in Tissue expression cluster: {len(liver_protein_2_genes)}")
print(f"Genes with liver in RNA tissue cell type enrichment: {len(liver_protein_3_genes)}")
print(f"Genes only in nTPM: {len(liver_protein_1_genes - liver_protein_2_genes - liver_protein_3_genes)}")
print(f"Genes only in cluster: {len(liver_protein_2_genes - liver_protein_1_genes - liver_protein_3_genes)}")
print(f"Genes only in enrichment: {len(liver_protein_3_genes - liver_protein_1_genes - liver_protein_2_genes)}")
print(f"Genes in all 3: {len(liver_protein_1_genes & liver_protein_2_genes & liver_protein_3_genes)}")

# Create improved 3-way Venn diagram
fig = plt.figure(figsize=(14, 12))
ax = fig.add_subplot(111)

venn = venn3(
    [liver_protein_1_genes, liver_protein_2_genes, liver_protein_3_genes],
    set_labels=('', '', ''),  # We'll add custom labels
    set_colors=('#FFB6C1', '#87CEEB', '#98FB98'),
    alpha=0.65,
    ax=ax
)

# Add title with statistics
title_text = f'Liver Protein Classification (3 Criteria)'
subtitle_text = f'Total: {total_genes} proteins | Liver candidates: {liver_candidates} ({liver_candidates/total_genes*100:.1f}%) | Non-liver: {len(non_liver_genes)}'
plt.title(title_text, fontsize=18, fontweight='bold', pad=25, y=0.98)
plt.text(0.5, 0.93, subtitle_text, ha='center', fontsize=13, transform=fig.transFigure)

# Customize all subset labels
only_ntpm = len(liver_protein_1_genes - liver_protein_2_genes - liver_protein_3_genes)
only_cluster = len(liver_protein_2_genes - liver_protein_1_genes - liver_protein_3_genes)
only_enrichment = len(liver_protein_3_genes - liver_protein_1_genes - liver_protein_2_genes)
ntpm_cluster = len((liver_protein_1_genes & liver_protein_2_genes) - liver_protein_3_genes)
ntpm_enrichment = len((liver_protein_1_genes & liver_protein_3_genes) - liver_protein_2_genes)
cluster_enrichment = len((liver_protein_2_genes & liver_protein_3_genes) - liver_protein_1_genes)
all_three = len(liver_protein_1_genes & liver_protein_2_genes & liver_protein_3_genes)

if venn.get_label_by_id('100'):
    venn.get_label_by_id('100').set_text(f"{only_ntpm}")
    venn.get_label_by_id('100').set_fontsize(13)
    venn.get_label_by_id('100').set_fontweight('bold')

if venn.get_label_by_id('010'):
    venn.get_label_by_id('010').set_text(f"{only_cluster}")
    venn.get_label_by_id('010').set_fontsize(13)
    venn.get_label_by_id('010').set_fontweight('bold')

if venn.get_label_by_id('001'):
    venn.get_label_by_id('001').set_text(f"{only_enrichment}")
    venn.get_label_by_id('001').set_fontsize(13)
    venn.get_label_by_id('001').set_fontweight('bold')

if venn.get_label_by_id('110'):
    venn.get_label_by_id('110').set_text(f"{ntpm_cluster}")
    venn.get_label_by_id('110').set_fontsize(13)
    venn.get_label_by_id('110').set_fontweight('bold')

if venn.get_label_by_id('101'):
    venn.get_label_by_id('101').set_text(f"{ntpm_enrichment}")
    venn.get_label_by_id('101').set_fontsize(13)
    venn.get_label_by_id('101').set_fontweight('bold')

if venn.get_label_by_id('011'):
    venn.get_label_by_id('011').set_text(f"{cluster_enrichment}")
    venn.get_label_by_id('011').set_fontsize(13)
    venn.get_label_by_id('011').set_fontweight('bold')

if venn.get_label_by_id('111'):
    venn.get_label_by_id('111').set_text(f"{all_three}")
    venn.get_label_by_id('111').set_fontsize(14)
    venn.get_label_by_id('111').set_fontweight('bold')
    venn.get_label_by_id('111').set_color('#006400')

# Add custom set labels with better positioning
ax.text(-0.65, 0.4, 'RNA tissue\nspecific nTPM', ha='center', va='center',
        fontsize=11, fontweight='bold', color='#8B0000',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#FFB6C1', linewidth=2))

ax.text(0.65, 0.4, 'Tissue\nexpression cluster', ha='center', va='center',
        fontsize=11, fontweight='bold', color='#00008B',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#87CEEB', linewidth=2))

ax.text(0, -0.72, 'RNA tissue cell\ntype enrichment', ha='center', va='center',
        fontsize=11, fontweight='bold', color='#006400',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#98FB98', linewidth=2))

# Add legend with detailed information
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#FFB6C1', edgecolor='#8B0000', linewidth=2,
          label=f'Only nTPM: {only_ntpm}', alpha=0.65),
    Patch(facecolor='#87CEEB', edgecolor='#00008B', linewidth=2,
          label=f'Only cluster: {only_cluster}', alpha=0.65),
    Patch(facecolor='#98FB98', edgecolor='#006400', linewidth=2,
          label=f'Only enrichment: {only_enrichment}', alpha=0.65),
    Patch(facecolor='#DDA0DD', edgecolor='#4B0082', linewidth=2,
          label=f'nTPM + cluster: {ntpm_cluster}', alpha=0.65),
    Patch(facecolor='#F0E68C', edgecolor='#DAA520', linewidth=2,
          label=f'nTPM + enrichment: {ntpm_enrichment}', alpha=0.65),
    Patch(facecolor='#87CEFA', edgecolor='#4682B4', linewidth=2,
          label=f'cluster + enrichment: {cluster_enrichment}', alpha=0.65),
    Patch(facecolor='#FFD700', edgecolor='#FF8C00', linewidth=2,
          label=f'All 3 criteria: {all_three}', alpha=0.8)
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=10, frameon=True,
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
    'RNA tissue cell type enrichment',
    'has_liver_nTPM',
    'has_liver_cluster',
    'has_liver_enrichment',
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
# Get all unique classifications
all_classifications = merged_data['Classification'].unique()
liver_classifications = [c for c in all_classifications if c.startswith('liver protein')]

for classification in liver_classifications:
    subset = merged_data[merged_data['Classification'] == classification]
    if len(subset) > 0:
        genes_file = f'{trace_dir}/{classification.replace(" ", "_").replace("(", "").replace(")", "").replace("+", "and")}_genes.csv'
        subset[['Gene', 'liver_nTPM_value']].to_csv(genes_file, index=False)
        print(f"{classification} genes ({len(subset)}): saved to {genes_file}")

print("\nClassification complete!")
