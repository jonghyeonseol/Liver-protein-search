import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
import re

# Read the input files
print("Loading data files...")

# Find CSV file in input folder
import glob
csv_files = glob.glob('input/*.csv')
if not csv_files:
    raise FileNotFoundError("No CSV file found in input folder")
if len(csv_files) > 1:
    print(f"Warning: Multiple CSV files found. Using: {csv_files[0]}")
    print(f"Other files: {csv_files[1:]}")

gene_data_file = csv_files[0]
print(f"Loading gene data from: {gene_data_file}")
gene_data = pd.read_csv(gene_data_file)

# Load protein atlas TSV file
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
# 1. liver protein_1: "liver" in "RNA tissue specific nTPM" column OR nTPM value >= 100
# 2. liver protein_2: "liver" in "Tissue expression cluster" column
# 3. liver protein_3: "liver" in "RNA tissue cell type enrichment" column
# 4. liver protein: "liver" in all three columns (or combinations)

# Extract liver nTPM values
merged_data['liver_nTPM_value'] = merged_data['RNA tissue specific nTPM'].apply(extract_liver_ntpm)

# Split nTPM into two categories: high (>=100) and low (<100)
merged_data['has_liver_in_ntpm'] = merged_data['RNA tissue specific nTPM'].fillna('').str.lower().str.contains('liver')
merged_data['has_liver_nTPM_high'] = (
    merged_data['has_liver_in_ntpm'] &
    (merged_data['liver_nTPM_value'] >= 100)
)
merged_data['has_liver_nTPM_low'] = (
    merged_data['has_liver_in_ntpm'] &
    (merged_data['liver_nTPM_value'] < 100)
)
merged_data['has_liver_cluster'] = merged_data['Tissue expression cluster'].fillna('').str.lower().str.contains('liver')
merged_data['has_liver_enrichment'] = merged_data['RNA tissue cell type enrichment'].fillna('').str.lower().str.contains('liver')

# Create classification column with more detailed logic
def classify_liver_protein(row):
    has_ntpm_high = row['has_liver_nTPM_high']
    has_ntpm_low = row['has_liver_nTPM_low']
    has_cluster = row['has_liver_cluster']
    has_enrichment = row['has_liver_enrichment']

    # Count how many conditions are met (treating high and low nTPM as separate)
    conditions = [has_ntpm_high, has_ntpm_low, has_cluster, has_enrichment]
    count = sum(conditions)

    if count == 0:
        return 'non-liver protein'
    elif count == 4:
        return 'liver protein (all 4)'
    elif count == 3:
        if not has_ntpm_high:
            return 'liver protein (nTPM<100 + cluster + enrichment)'
        elif not has_ntpm_low:
            return 'liver protein (nTPM≥100 + cluster + enrichment)'
        elif not has_cluster:
            return 'liver protein (both nTPM + enrichment)'
        elif not has_enrichment:
            return 'liver protein (both nTPM + cluster)'
    elif count == 2:
        if has_ntpm_high and has_ntpm_low:
            return 'liver protein (both nTPM levels)'
        elif has_ntpm_high and has_cluster:
            return 'liver protein (nTPM≥100 + cluster)'
        elif has_ntpm_high and has_enrichment:
            return 'liver protein (nTPM≥100 + enrichment)'
        elif has_ntpm_low and has_cluster:
            return 'liver protein (nTPM<100 + cluster)'
        elif has_ntpm_low and has_enrichment:
            return 'liver protein (nTPM<100 + enrichment)'
        elif has_cluster and has_enrichment:
            return 'liver protein (cluster + enrichment)'
    else:  # count == 1
        if has_ntpm_high:
            return 'liver protein (nTPM≥100 only)'
        elif has_ntpm_low:
            return 'liver protein (nTPM<100 only)'
        elif has_cluster:
            return 'liver protein (cluster only)'
        elif has_enrichment:
            return 'liver protein (enrichment only)'

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

# Get gene lists for Venn diagram (4 categories)
all_genes = set(merged_data['Gene'].tolist())
liver_ntpm_high_genes = set(merged_data[merged_data['has_liver_nTPM_high']]['Gene'].tolist())
liver_ntpm_low_genes = set(merged_data[merged_data['has_liver_nTPM_low']]['Gene'].tolist())
liver_cluster_genes = set(merged_data[merged_data['has_liver_cluster']]['Gene'].tolist())
liver_enrichment_genes = set(merged_data[merged_data['has_liver_enrichment']]['Gene'].tolist())

total_genes = len(all_genes)
liver_candidates_set = liver_ntpm_high_genes | liver_ntpm_low_genes | liver_cluster_genes | liver_enrichment_genes
liver_candidates = len(liver_candidates_set)
non_liver_genes = all_genes - liver_candidates_set

print(f"\nVenn Diagram Statistics (4 categories):")
print(f"Total genes: {total_genes}")
print(f"Liver protein candidates (union): {liver_candidates}")
print(f"Non-liver proteins: {len(non_liver_genes)}")
print(f"Genes with liver nTPM ≥100: {len(liver_ntpm_high_genes)}")
print(f"Genes with liver nTPM <100: {len(liver_ntpm_low_genes)}")
print(f"Genes with liver in Tissue expression cluster: {len(liver_cluster_genes)}")
print(f"Genes with liver in RNA tissue cell type enrichment: {len(liver_enrichment_genes)}")
print(f"Genes in all 4 categories: {len(liver_ntpm_high_genes & liver_ntpm_low_genes & liver_cluster_genes & liver_enrichment_genes)}")

# Create 4-way Venn diagram using two 3-way diagrams
fig = plt.figure(figsize=(20, 10))

# Add overall title with statistics
title_text = f'Liver Protein Classification (4 Criteria: nTPM Split by Threshold 100)'
subtitle_text = f'Total: {total_genes} proteins | Liver candidates: {liver_candidates} ({liver_candidates/total_genes*100:.1f}%) | Non-liver: {len(non_liver_genes)}'
fig.suptitle(title_text, fontsize=18, fontweight='bold', y=0.98)
plt.text(0.5, 0.93, subtitle_text, ha='center', fontsize=13, transform=fig.transFigure)

# Left subplot: nTPM >= 100 with cluster and enrichment
ax1 = fig.add_subplot(121)
venn1 = venn3(
    [liver_ntpm_high_genes, liver_cluster_genes, liver_enrichment_genes],
    set_labels=('', '', ''),
    set_colors=('#FF6B6B', '#87CEEB', '#98FB98'),
    alpha=0.65,
    ax=ax1
)
ax1.set_title('nTPM ≥ 100 (liver mentioned)', fontsize=14, fontweight='bold', pad=15)

# Customize subset labels for left diagram (nTPM >= 100)
only_ntpm_high = len(liver_ntpm_high_genes - liver_cluster_genes - liver_enrichment_genes)
only_cluster_left = len(liver_cluster_genes - liver_ntpm_high_genes - liver_enrichment_genes)
only_enrichment_left = len(liver_enrichment_genes - liver_ntpm_high_genes - liver_cluster_genes)
ntpm_high_cluster = len((liver_ntpm_high_genes & liver_cluster_genes) - liver_enrichment_genes)
ntpm_high_enrichment = len((liver_ntpm_high_genes & liver_enrichment_genes) - liver_cluster_genes)
cluster_enrichment_left = len((liver_cluster_genes & liver_enrichment_genes) - liver_ntpm_high_genes)
all_three_high = len(liver_ntpm_high_genes & liver_cluster_genes & liver_enrichment_genes)

if venn1.get_label_by_id('100'):
    venn1.get_label_by_id('100').set_text(f"{only_ntpm_high}")
    venn1.get_label_by_id('100').set_fontsize(12)
    venn1.get_label_by_id('100').set_fontweight('bold')
if venn1.get_label_by_id('010'):
    venn1.get_label_by_id('010').set_text(f"{only_cluster_left}")
    venn1.get_label_by_id('010').set_fontsize(12)
    venn1.get_label_by_id('010').set_fontweight('bold')
if venn1.get_label_by_id('001'):
    venn1.get_label_by_id('001').set_text(f"{only_enrichment_left}")
    venn1.get_label_by_id('001').set_fontsize(12)
    venn1.get_label_by_id('001').set_fontweight('bold')
if venn1.get_label_by_id('110'):
    venn1.get_label_by_id('110').set_text(f"{ntpm_high_cluster}")
    venn1.get_label_by_id('110').set_fontsize(12)
    venn1.get_label_by_id('110').set_fontweight('bold')
if venn1.get_label_by_id('101'):
    venn1.get_label_by_id('101').set_text(f"{ntpm_high_enrichment}")
    venn1.get_label_by_id('101').set_fontsize(12)
    venn1.get_label_by_id('101').set_fontweight('bold')
if venn1.get_label_by_id('011'):
    venn1.get_label_by_id('011').set_text(f"{cluster_enrichment_left}")
    venn1.get_label_by_id('011').set_fontsize(12)
    venn1.get_label_by_id('011').set_fontweight('bold')
if venn1.get_label_by_id('111'):
    venn1.get_label_by_id('111').set_text(f"{all_three_high}")
    venn1.get_label_by_id('111').set_fontsize(13)
    venn1.get_label_by_id('111').set_fontweight('bold')
    venn1.get_label_by_id('111').set_color('#8B0000')

# Add custom set labels for left diagram
ax1.text(-0.65, 0.4, 'nTPM ≥ 100\n(liver)', ha='center', va='center',
        fontsize=10, fontweight='bold', color='#8B0000',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#FF6B6B', linewidth=2))
ax1.text(0.65, 0.4, 'Cluster', ha='center', va='center',
        fontsize=10, fontweight='bold', color='#00008B',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#87CEEB', linewidth=2))
ax1.text(0, -0.72, 'Enrichment', ha='center', va='center',
        fontsize=10, fontweight='bold', color='#006400',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#98FB98', linewidth=2))

# Right subplot: nTPM < 100 with cluster and enrichment
ax2 = fig.add_subplot(122)
venn2 = venn3(
    [liver_ntpm_low_genes, liver_cluster_genes, liver_enrichment_genes],
    set_labels=('', '', ''),
    set_colors=('#FFB6C1', '#87CEEB', '#98FB98'),
    alpha=0.65,
    ax=ax2
)
ax2.set_title('nTPM < 100 (liver mentioned)', fontsize=14, fontweight='bold', pad=15)

# Customize subset labels for right diagram (nTPM < 100)
only_ntpm_low = len(liver_ntpm_low_genes - liver_cluster_genes - liver_enrichment_genes)
only_cluster_right = len(liver_cluster_genes - liver_ntpm_low_genes - liver_enrichment_genes)
only_enrichment_right = len(liver_enrichment_genes - liver_ntpm_low_genes - liver_cluster_genes)
ntpm_low_cluster = len((liver_ntpm_low_genes & liver_cluster_genes) - liver_enrichment_genes)
ntpm_low_enrichment = len((liver_ntpm_low_genes & liver_enrichment_genes) - liver_cluster_genes)
cluster_enrichment_right = len((liver_cluster_genes & liver_enrichment_genes) - liver_ntpm_low_genes)
all_three_low = len(liver_ntpm_low_genes & liver_cluster_genes & liver_enrichment_genes)

if venn2.get_label_by_id('100'):
    venn2.get_label_by_id('100').set_text(f"{only_ntpm_low}")
    venn2.get_label_by_id('100').set_fontsize(12)
    venn2.get_label_by_id('100').set_fontweight('bold')
if venn2.get_label_by_id('010'):
    venn2.get_label_by_id('010').set_text(f"{only_cluster_right}")
    venn2.get_label_by_id('010').set_fontsize(12)
    venn2.get_label_by_id('010').set_fontweight('bold')
if venn2.get_label_by_id('001'):
    venn2.get_label_by_id('001').set_text(f"{only_enrichment_right}")
    venn2.get_label_by_id('001').set_fontsize(12)
    venn2.get_label_by_id('001').set_fontweight('bold')
if venn2.get_label_by_id('110'):
    venn2.get_label_by_id('110').set_text(f"{ntpm_low_cluster}")
    venn2.get_label_by_id('110').set_fontsize(12)
    venn2.get_label_by_id('110').set_fontweight('bold')
if venn2.get_label_by_id('101'):
    venn2.get_label_by_id('101').set_text(f"{ntpm_low_enrichment}")
    venn2.get_label_by_id('101').set_fontsize(12)
    venn2.get_label_by_id('101').set_fontweight('bold')
if venn2.get_label_by_id('011'):
    venn2.get_label_by_id('011').set_text(f"{cluster_enrichment_right}")
    venn2.get_label_by_id('011').set_fontsize(12)
    venn2.get_label_by_id('011').set_fontweight('bold')
if venn2.get_label_by_id('111'):
    venn2.get_label_by_id('111').set_text(f"{all_three_low}")
    venn2.get_label_by_id('111').set_fontsize(13)
    venn2.get_label_by_id('111').set_fontweight('bold')
    venn2.get_label_by_id('111').set_color('#8B0000')

# Add custom set labels for right diagram
ax2.text(-0.65, 0.4, 'nTPM < 100\n(liver)', ha='center', va='center',
        fontsize=10, fontweight='bold', color='#8B0000',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#FFB6C1', linewidth=2))
ax2.text(0.65, 0.4, 'Cluster', ha='center', va='center',
        fontsize=10, fontweight='bold', color='#00008B',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#87CEEB', linewidth=2))
ax2.text(0, -0.72, 'Enrichment', ha='center', va='center',
        fontsize=10, fontweight='bold', color='#006400',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#98FB98', linewidth=2))

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
    'has_liver_nTPM_high',
    'has_liver_nTPM_low',
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
