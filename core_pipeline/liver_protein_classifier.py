# Standard library imports
import re
import sys

# Third-party imports
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Local imports
from config import (
    INPUT_DIR, OUTPUT_DIR, TRACE_DIR,
    PROTEIN_ATLAS_FILE, CSV_PATTERN,
    NTPM_HIGH_THRESHOLD, DPI, FIGURE_SIZE_VENN,
    TRACE_DATA_FILE, CLASSIFICATION_SUMMARY_FILE,
    VENN_DIAGRAM_OUTPUT
)

# =============================================================================
# INPUT VALIDATION FUNCTIONS
# =============================================================================

def validate_input_files():
    """
    Validate that all required input files exist and are readable.

    Returns:
        bool: True if all files are valid, False otherwise
    """
    errors = []

    # Check if input directory exists
    if not INPUT_DIR.exists():
        errors.append(f"Input directory not found: {INPUT_DIR}")
        return False, errors

    # Check if protein atlas file exists
    if not PROTEIN_ATLAS_FILE.exists():
        errors.append(f"Protein atlas file not found: {PROTEIN_ATLAS_FILE}")

    # Check if at least one CSV file exists
    csv_files = list(INPUT_DIR.glob(CSV_PATTERN))
    if not csv_files:
        errors.append(f"No CSV files found in {INPUT_DIR} matching pattern: {CSV_PATTERN}")

    # Check if output and trace directories exist, create if not
    for directory in [OUTPUT_DIR, TRACE_DIR]:
        if not directory.exists():
            try:
                directory.mkdir(parents=True, exist_ok=True)
                print(f"Created directory: {directory}")
            except Exception as e:
                errors.append(f"Failed to create directory {directory}: {str(e)}")

    if errors:
        return False, errors
    return True, []


def validate_dataframe(df, name, required_columns):
    """
    Validate that a dataframe contains required columns and is not empty.

    Args:
        df: pandas DataFrame to validate
        name: Name of the dataframe (for error messages)
        required_columns: List of column names that must be present

    Returns:
        tuple: (is_valid, errors, warnings)
            is_valid (bool): True if no critical errors
            errors (list): List of critical error messages
            warnings (list): List of warning messages
    """
    errors = []
    warnings = []

    # Check if dataframe is empty
    if df.empty:
        errors.append(f"{name} is empty")
        return False, errors, warnings

    # Check for required columns
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        errors.append(f"{name} missing required columns: {missing_columns}")
        errors.append(f"Available columns: {list(df.columns)}")

    # Check for duplicate gene names (if Gene column exists) - WARNING, not error
    if 'Gene' in df.columns:
        duplicates = df[df['Gene'].duplicated(keep=False)]['Gene'].dropna().unique()
        if len(duplicates) > 0:
            warnings.append(f"{name} contains {len(duplicates)} duplicate gene names (first 5): {list(duplicates[:5])}")

    if errors:
        return False, errors, warnings
    return True, [], warnings

# =============================================================================
# MAIN EXECUTION
# =============================================================================

print("="*80)
print("LIVER PROTEIN CLASSIFICATION PIPELINE")
print("="*80)

# Step 1: Validate input files
print("\n[Step 1/5] Validating input files...")
valid, errors = validate_input_files()
if not valid:
    print("\n❌ Input validation failed:")
    for error in errors:
        print(f"  - {error}")
    sys.exit(1)
print("✓ All input files validated successfully")

# Step 2: Load data files
print("\n[Step 2/5] Loading data files...")

# Find CSV file in input folder
csv_files = list(INPUT_DIR.glob(CSV_PATTERN))
if len(csv_files) > 1:
    print(f"⚠ Warning: Multiple CSV files found. Using: {csv_files[0]}")
    print(f"  Other files: {csv_files[1:]}")

gene_data_file = csv_files[0]
print(f"  Loading gene data from: {gene_data_file}")
try:
    gene_data = pd.read_csv(gene_data_file)
except Exception as e:
    print(f"❌ Failed to read {gene_data_file}: {str(e)}")
    sys.exit(1)

# Load protein atlas TSV file
print(f"  Loading protein atlas from: {PROTEIN_ATLAS_FILE}")
try:
    protein_atlas = pd.read_csv(PROTEIN_ATLAS_FILE, sep='\t')
except Exception as e:
    print(f"❌ Failed to read {PROTEIN_ATLAS_FILE}: {str(e)}")
    sys.exit(1)

print("✓ Data files loaded successfully")

# Step 3: Validate data structure
print("\n[Step 3/5] Validating data structure...")

# Validate gene_data
valid, errors, warnings = validate_dataframe(gene_data, "Gene data", ['Gene'])
if not valid:
    print("\n❌ Gene data validation failed:")
    for error in errors:
        print(f"  - {error}")
    sys.exit(1)

# Print warnings if any
for warning in warnings:
    print(f"  ⚠ Warning: {warning}")

# Validate protein_atlas
required_protein_atlas_columns = [
    'Gene',
    'RNA tissue specific nTPM',
    'Tissue expression cluster',
    'RNA tissue cell type enrichment'
]
valid, errors, warnings_pa = validate_dataframe(protein_atlas, "Protein atlas", required_protein_atlas_columns)
if not valid:
    print("\n❌ Protein atlas validation failed:")
    for error in errors:
        print(f"  - {error}")
    sys.exit(1)

# Print warnings if any
for warning in warnings_pa:
    print(f"  ⚠ Warning: {warning}")

print(f"  Gene data: {len(gene_data)} genes")
print(f"  Protein atlas: {len(protein_atlas)} genes")
print("✓ Data structure validated successfully")

# Step 4: Merge datasets
print("\n[Step 4/5] Merging datasets...")

merged_data = gene_data.merge(
    protein_atlas[['Gene', 'RNA tissue specific nTPM', 'Tissue expression cluster', 'RNA tissue cell type enrichment']],
    on='Gene',
    how='left'
)

print(f"  Merged {len(merged_data)} records")
if len(merged_data) < len(gene_data):
    print(f"  ⚠ Warning: {len(gene_data) - len(merged_data)} genes lost during merge")

print("✓ Datasets merged successfully")

# =============================================================================
# CLASSIFICATION FUNCTIONS
# =============================================================================

def extract_liver_ntpm(text):
    """
    Extract numerical liver nTPM value from text.

    Args:
        text: String containing liver nTPM data

    Returns:
        float: Liver nTPM value if found, None otherwise
    """
    if pd.isna(text) or text == '':
        return None
    # Search for pattern "liver: number" or "liver:number"
    pattern = r'liver:\s*(\d+\.?\d*)'
    match = re.search(pattern, str(text).lower())
    if match:
        return float(match.group(1))
    return None

# Step 5: Classify liver proteins
print("\n[Step 5/5] Classifying liver proteins...")

# Extract liver nTPM values
print("  Extracting liver nTPM values...")
merged_data['liver_nTPM_value'] = merged_data['RNA tissue specific nTPM'].apply(extract_liver_ntpm)

# Split nTPM into two categories: high (>=100) and low (<100)
merged_data['has_liver_in_ntpm'] = merged_data['RNA tissue specific nTPM'].fillna('').str.lower().str.contains('liver')
merged_data['has_liver_nTPM_high'] = (
    merged_data['has_liver_in_ntpm'] &
    (merged_data['liver_nTPM_value'] >= NTPM_HIGH_THRESHOLD)
)
merged_data['has_liver_nTPM_low'] = (
    merged_data['has_liver_in_ntpm'] &
    (merged_data['liver_nTPM_value'] < NTPM_HIGH_THRESHOLD)
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

print("  Applying classification rules...")
merged_data['Classification'] = merged_data.apply(classify_liver_protein, axis=1)

# Count each category
classification_counts = merged_data['Classification'].value_counts()
print("\n" + "="*80)
print("CLASSIFICATION RESULTS")
print("="*80)
print(classification_counts)
print("="*80)

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
fig = plt.figure(figsize=FIGURE_SIZE_VENN)

# Add overall title with statistics (improved layout)
title_text = f'Liver Protein Classification'
subtitle_text = f'nTPM Threshold: {NTPM_HIGH_THRESHOLD} | Total: {total_genes:,} genes | Liver: {liver_candidates:,} ({liver_candidates/total_genes*100:.1f}%) | Non-liver: {len(non_liver_genes):,} ({len(non_liver_genes)/total_genes*100:.1f}%)'
fig.suptitle(title_text, fontsize=22, fontweight='bold', y=0.97)
plt.text(0.5, 0.93, subtitle_text, ha='center', fontsize=12, transform=fig.transFigure,
         bbox=dict(boxstyle='round,pad=0.5', facecolor='#F0F0F0', edgecolor='#CCCCCC', linewidth=1.5, alpha=0.8))

# Left subplot: nTPM >= threshold with cluster and enrichment
ax1 = fig.add_subplot(121)
venn1 = venn3(
    [liver_ntpm_high_genes, liver_cluster_genes, liver_enrichment_genes],
    set_labels=('', '', ''),
    set_colors=('#FF6B6B', '#87CEEB', '#98FB98'),
    alpha=0.7,
    ax=ax1
)
ax1.set_title(f'High Expression (nTPM ≥ {NTPM_HIGH_THRESHOLD})', fontsize=15, fontweight='bold', pad=20,
              bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFE6E6', edgecolor='#FF6B6B', linewidth=2))

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
    venn1.get_label_by_id('100').set_fontsize(14)
    venn1.get_label_by_id('100').set_fontweight('bold')
if venn1.get_label_by_id('010'):
    venn1.get_label_by_id('010').set_text(f"{only_cluster_left}")
    venn1.get_label_by_id('010').set_fontsize(14)
    venn1.get_label_by_id('010').set_fontweight('bold')
if venn1.get_label_by_id('001'):
    venn1.get_label_by_id('001').set_text(f"{only_enrichment_left}")
    venn1.get_label_by_id('001').set_fontsize(14)
    venn1.get_label_by_id('001').set_fontweight('bold')
if venn1.get_label_by_id('110'):
    venn1.get_label_by_id('110').set_text(f"{ntpm_high_cluster}")
    venn1.get_label_by_id('110').set_fontsize(14)
    venn1.get_label_by_id('110').set_fontweight('bold')
if venn1.get_label_by_id('101'):
    venn1.get_label_by_id('101').set_text(f"{ntpm_high_enrichment}")
    venn1.get_label_by_id('101').set_fontsize(14)
    venn1.get_label_by_id('101').set_fontweight('bold')
if venn1.get_label_by_id('011'):
    venn1.get_label_by_id('011').set_text(f"{cluster_enrichment_left}")
    venn1.get_label_by_id('011').set_fontsize(14)
    venn1.get_label_by_id('011').set_fontweight('bold')
if venn1.get_label_by_id('111'):
    venn1.get_label_by_id('111').set_text(f"{all_three_high}")
    venn1.get_label_by_id('111').set_fontsize(16)
    venn1.get_label_by_id('111').set_fontweight('bold')
    venn1.get_label_by_id('111').set_color('#FFFFFF')
    venn1.get_label_by_id('111').set_bbox(dict(boxstyle='round,pad=0.5', facecolor='#8B0000', edgecolor='white', linewidth=2))

# Add custom set labels for left diagram (improved positioning and styling)
ax1.text(-0.7, 0.45, f'nTPM ≥ {NTPM_HIGH_THRESHOLD}', ha='center', va='center',
        fontsize=11, fontweight='bold', color='white',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FF6B6B', edgecolor='#8B0000', linewidth=2.5, alpha=0.9))
ax1.text(0.7, 0.45, 'Tissue\nCluster', ha='center', va='center',
        fontsize=11, fontweight='bold', color='white',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#5DADE2', edgecolor='#1F618D', linewidth=2.5, alpha=0.9))
ax1.text(0, -0.75, 'Cell Type\nEnrichment', ha='center', va='center',
        fontsize=11, fontweight='bold', color='white',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#58D68D', edgecolor='#1E8449', linewidth=2.5, alpha=0.9))

# Right subplot: nTPM < threshold with cluster and enrichment
ax2 = fig.add_subplot(122)
venn2 = venn3(
    [liver_ntpm_low_genes, liver_cluster_genes, liver_enrichment_genes],
    set_labels=('', '', ''),
    set_colors=('#FFB6C1', '#87CEEB', '#98FB98'),
    alpha=0.7,
    ax=ax2
)
ax2.set_title(f'Low Expression (nTPM < {NTPM_HIGH_THRESHOLD})', fontsize=15, fontweight='bold', pad=20,
              bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFE6F0', edgecolor='#FFB6C1', linewidth=2))

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
    venn2.get_label_by_id('100').set_fontsize(14)
    venn2.get_label_by_id('100').set_fontweight('bold')
if venn2.get_label_by_id('010'):
    venn2.get_label_by_id('010').set_text(f"{only_cluster_right}")
    venn2.get_label_by_id('010').set_fontsize(14)
    venn2.get_label_by_id('010').set_fontweight('bold')
if venn2.get_label_by_id('001'):
    venn2.get_label_by_id('001').set_text(f"{only_enrichment_right}")
    venn2.get_label_by_id('001').set_fontsize(14)
    venn2.get_label_by_id('001').set_fontweight('bold')
if venn2.get_label_by_id('110'):
    venn2.get_label_by_id('110').set_text(f"{ntpm_low_cluster}")
    venn2.get_label_by_id('110').set_fontsize(14)
    venn2.get_label_by_id('110').set_fontweight('bold')
if venn2.get_label_by_id('101'):
    venn2.get_label_by_id('101').set_text(f"{ntpm_low_enrichment}")
    venn2.get_label_by_id('101').set_fontsize(14)
    venn2.get_label_by_id('101').set_fontweight('bold')
if venn2.get_label_by_id('011'):
    venn2.get_label_by_id('011').set_text(f"{cluster_enrichment_right}")
    venn2.get_label_by_id('011').set_fontsize(14)
    venn2.get_label_by_id('011').set_fontweight('bold')
if venn2.get_label_by_id('111'):
    venn2.get_label_by_id('111').set_text(f"{all_three_low}")
    venn2.get_label_by_id('111').set_fontsize(16)
    venn2.get_label_by_id('111').set_fontweight('bold')
    venn2.get_label_by_id('111').set_color('#FFFFFF')
    venn2.get_label_by_id('111').set_bbox(dict(boxstyle='round,pad=0.5', facecolor='#8B0000', edgecolor='white', linewidth=2))

# Add custom set labels for right diagram (improved positioning and styling)
ax2.text(-0.7, 0.45, f'nTPM < {NTPM_HIGH_THRESHOLD}', ha='center', va='center',
        fontsize=11, fontweight='bold', color='white',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFB6C1', edgecolor='#C70039', linewidth=2.5, alpha=0.9))
ax2.text(0.7, 0.45, 'Tissue\nCluster', ha='center', va='center',
        fontsize=11, fontweight='bold', color='white',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#5DADE2', edgecolor='#1F618D', linewidth=2.5, alpha=0.9))
ax2.text(0, -0.75, 'Cell Type\nEnrichment', ha='center', va='center',
        fontsize=11, fontweight='bold', color='white',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#58D68D', edgecolor='#1E8449', linewidth=2.5, alpha=0.9))

# Save the Venn diagram
OUTPUT_DIR.mkdir(exist_ok=True)
plt.savefig(VENN_DIAGRAM_OUTPUT, dpi=DPI, bbox_inches='tight')
print(f"\nVenn diagram saved to: {VENN_DIAGRAM_OUTPUT}")
plt.close()  # Free memory

# Create Trace folder and save tracing data
TRACE_DIR.mkdir(exist_ok=True)

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
trace_data.to_csv(TRACE_DATA_FILE, index=False)
print(f"Tracing data saved to: {TRACE_DATA_FILE}")

# Save summary statistics to Trace folder
classification_counts.to_csv(CLASSIFICATION_SUMMARY_FILE, header=['Count'])
print(f"Summary statistics saved to: {CLASSIFICATION_SUMMARY_FILE}")

# Save classified gene lists with liver nTPM values to Trace folder
# Get all unique classifications
all_classifications = merged_data['Classification'].unique()
liver_classifications = [c for c in all_classifications if c.startswith('liver protein')]

for classification in liver_classifications:
    subset = merged_data[merged_data['Classification'] == classification]
    if len(subset) > 0:
        genes_file = TRACE_DIR / f'{classification.replace(" ", "_").replace("(", "").replace(")", "").replace("+", "and")}_genes.csv'
        subset[['Gene', 'liver_nTPM_value']].to_csv(genes_file, index=False)
        print(f"{classification} genes ({len(subset)}): saved to {genes_file}")

print("\n" + "="*80)
print("✓ CLASSIFICATION PIPELINE COMPLETED SUCCESSFULLY")
print("="*80)
