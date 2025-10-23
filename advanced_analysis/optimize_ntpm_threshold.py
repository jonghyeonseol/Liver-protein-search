"""
nTPM Threshold Optimization Analysis
Compare classification results across different nTPM threshold values.
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent / 'core_pipeline'))
from config import (
    INPUT_DIR, OUTPUT_DIR, TRACE_DIR, PROTEIN_ATLAS_FILE,
    CSV_PATTERN, KNOWN_LIVER_MARKERS, DPI
)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Threshold values to test
THRESHOLDS_TO_TEST = [10, 20, 50, 100, 150, 200, 300]

# Output file
THRESHOLD_ANALYSIS_OUTPUT = OUTPUT_DIR / 'ntpm_threshold_optimization.png'
THRESHOLD_RESULTS_CSV = TRACE_DIR / 'threshold_optimization_results.csv'

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def extract_liver_ntpm(text):
    """Extract numerical liver nTPM value from text."""
    if pd.isna(text) or text == '':
        return None
    pattern = r'liver:\s*(\d+\.?\d*)'
    match = re.search(pattern, str(text).lower())
    if match:
        return float(match.group(1))
    return None

def classify_with_threshold(merged_data, threshold):
    """
    Classify liver proteins using a specific nTPM threshold.

    Args:
        merged_data: DataFrame with merged gene and protein atlas data
        threshold: nTPM threshold value to use

    Returns:
        DataFrame with classification results
    """
    df = merged_data.copy()

    # Apply threshold
    df['has_liver_nTPM_high'] = (
        df['has_liver_in_ntpm'] &
        (df['liver_nTPM_value'] >= threshold)
    )
    df['has_liver_nTPM_low'] = (
        df['has_liver_in_ntpm'] &
        (df['liver_nTPM_value'] < threshold)
    )

    # Classification function
    def classify_liver_protein(row):
        has_ntpm_high = row['has_liver_nTPM_high']
        has_ntpm_low = row['has_liver_nTPM_low']
        has_cluster = row['has_liver_cluster']
        has_enrichment = row['has_liver_enrichment']

        conditions = [has_ntpm_high, has_ntpm_low, has_cluster, has_enrichment]
        count = sum(conditions)

        if count == 0:
            return 'non-liver protein'
        elif count >= 3:
            if has_ntpm_high:
                return f'liver protein (nTPM≥{threshold} + cluster + enrichment)'
            else:
                return f'liver protein (nTPM<{threshold} + cluster + enrichment)'
        elif count == 2:
            if has_ntpm_high and has_cluster:
                return f'liver protein (nTPM≥{threshold} + cluster)'
            elif has_ntpm_high and has_enrichment:
                return f'liver protein (nTPM≥{threshold} + enrichment)'
            elif has_ntpm_low and has_cluster:
                return f'liver protein (nTPM<{threshold} + cluster)'
            elif has_ntpm_low and has_enrichment:
                return f'liver protein (nTPM<{threshold} + enrichment)'
            elif has_cluster and has_enrichment:
                return 'liver protein (cluster + enrichment)'
        else:  # count == 1
            if has_ntpm_high:
                return f'liver protein (nTPM≥{threshold} only)'
            elif has_ntpm_low:
                return f'liver protein (nTPM<{threshold} only)'
            elif has_cluster:
                return 'liver protein (cluster only)'
            elif has_enrichment:
                return 'liver protein (enrichment only)'

        return 'non-liver protein'

    df['Classification'] = df.apply(classify_liver_protein, axis=1)
    return df

def calculate_metrics(df, threshold):
    """
    Calculate performance metrics for a specific threshold.

    Returns:
        dict with various metrics
    """
    metrics = {
        'threshold': threshold,
    }

    # Count genes in different categories
    all_genes = len(df)
    liver_genes = len(df[df['Classification'] != 'non-liver protein'])

    # High confidence categories (all 3 criteria met)
    high_conf_high = len(df[df['Classification'].str.contains(f'nTPM≥{threshold}.*cluster.*enrichment', na=False)])
    high_conf_low = len(df[df['Classification'].str.contains(f'nTPM<{threshold}.*cluster.*enrichment', na=False)])

    # Medium confidence (2 criteria)
    med_conf_high = len(df[df['Classification'].str.contains(f'nTPM≥{threshold}', na=False)]) - high_conf_high
    med_conf_low = len(df[df['Classification'].str.contains(f'nTPM<{threshold}', na=False)]) - high_conf_low

    # Low confidence (1 criterion only)
    low_conf = len(df[df['Classification'].str.contains('only', na=False)])

    # nTPM distribution
    ntpm_high_count = len(df[(df['liver_nTPM_value'] >= threshold) & (df['liver_nTPM_value'].notna())])
    ntpm_low_count = len(df[(df['liver_nTPM_value'] < threshold) & (df['liver_nTPM_value'].notna())])

    # Known liver markers coverage
    known_markers_found = len(df[(df['Gene'].isin(KNOWN_LIVER_MARKERS)) &
                                  (df['Classification'] != 'non-liver protein')])
    known_markers_high_conf = len(df[(df['Gene'].isin(KNOWN_LIVER_MARKERS)) &
                                      (df['Classification'].str.contains('cluster.*enrichment', na=False))])

    metrics.update({
        'total_genes': all_genes,
        'liver_genes': liver_genes,
        'liver_percentage': liver_genes / all_genes * 100,
        'high_conf_ntpm_high': high_conf_high,
        'high_conf_ntpm_low': high_conf_low,
        'high_conf_total': high_conf_high + high_conf_low,
        'med_conf_ntpm_high': med_conf_high,
        'med_conf_ntpm_low': med_conf_low,
        'med_conf_total': med_conf_high + med_conf_low,
        'low_conf_total': low_conf,
        'ntpm_high_count': ntpm_high_count,
        'ntpm_low_count': ntpm_low_count,
        'ntpm_balance_ratio': ntpm_high_count / ntpm_low_count if ntpm_low_count > 0 else float('inf'),
        'known_markers_total': len(KNOWN_LIVER_MARKERS),
        'known_markers_found': known_markers_found,
        'known_markers_coverage': known_markers_found / len(KNOWN_LIVER_MARKERS) * 100,
        'known_markers_high_conf': known_markers_high_conf,
        'known_markers_high_conf_rate': known_markers_high_conf / len(KNOWN_LIVER_MARKERS) * 100,
    })

    return metrics

# =============================================================================
# MAIN EXECUTION
# =============================================================================

print("="*80)
print("nTPM THRESHOLD OPTIMIZATION ANALYSIS")
print("="*80)

# Load data
print("\n[1/4] Loading data...")
csv_files = list(INPUT_DIR.glob(CSV_PATTERN))
if not csv_files:
    print(f"❌ No CSV files found in {INPUT_DIR}")
    sys.exit(1)

gene_data = pd.read_csv(csv_files[0])
protein_atlas = pd.read_csv(PROTEIN_ATLAS_FILE, sep='\t')

print(f"  Gene data: {len(gene_data)} genes")
print(f"  Protein atlas: {len(protein_atlas)} genes")

# Merge datasets
print("\n[2/4] Merging datasets...")
merged_data = gene_data.merge(
    protein_atlas[['Gene', 'RNA tissue specific nTPM', 'Tissue expression cluster', 'RNA tissue cell type enrichment']],
    on='Gene',
    how='left'
)

# Extract liver nTPM values
merged_data['liver_nTPM_value'] = merged_data['RNA tissue specific nTPM'].apply(extract_liver_ntpm)
merged_data['has_liver_in_ntpm'] = merged_data['RNA tissue specific nTPM'].fillna('').str.lower().str.contains('liver')
merged_data['has_liver_cluster'] = merged_data['Tissue expression cluster'].fillna('').str.lower().str.contains('liver')
merged_data['has_liver_enrichment'] = merged_data['RNA tissue cell type enrichment'].fillna('').str.lower().str.contains('liver')

print(f"  Merged: {len(merged_data)} records")
print(f"  Genes with liver nTPM values: {merged_data['liver_nTPM_value'].notna().sum()}")

# Test different thresholds
print(f"\n[3/4] Testing {len(THRESHOLDS_TO_TEST)} different thresholds...")
results = []

for threshold in THRESHOLDS_TO_TEST:
    print(f"\n  Testing threshold: {threshold}")

    # Classify with this threshold
    classified_data = classify_with_threshold(merged_data, threshold)

    # Calculate metrics
    metrics = calculate_metrics(classified_data, threshold)
    results.append(metrics)

    print(f"    Liver genes: {metrics['liver_genes']} ({metrics['liver_percentage']:.1f}%)")
    print(f"    High confidence: {metrics['high_conf_total']} (≥{threshold}: {metrics['high_conf_ntpm_high']}, <{threshold}: {metrics['high_conf_ntpm_low']})")
    print(f"    Known markers coverage: {metrics['known_markers_found']}/{metrics['known_markers_total']} ({metrics['known_markers_coverage']:.1f}%)")
    print(f"    Known markers in high conf: {metrics['known_markers_high_conf']} ({metrics['known_markers_high_conf_rate']:.1f}%)")

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Save results
results_df.to_csv(THRESHOLD_RESULTS_CSV, index=False)
print(f"\n✓ Results saved to: {THRESHOLD_RESULTS_CSV}")

# Create visualization
print("\n[4/4] Creating visualization...")

fig = plt.figure(figsize=(20, 12))
fig.suptitle('nTPM Threshold Optimization Analysis', fontsize=20, fontweight='bold', y=0.98)

# Plot 1: Total liver genes identified
ax1 = plt.subplot(2, 3, 1)
ax1.plot(results_df['threshold'], results_df['liver_genes'], 'o-', linewidth=2, markersize=8, color='#2E86AB')
ax1.set_xlabel('nTPM Threshold', fontsize=12, fontweight='bold')
ax1.set_ylabel('Number of Liver Genes', fontsize=12, fontweight='bold')
ax1.set_title('Total Liver Genes Identified', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
for i, row in results_df.iterrows():
    ax1.text(row['threshold'], row['liver_genes'], f"{row['liver_genes']}",
             ha='center', va='bottom', fontsize=9)

# Plot 2: High confidence genes
ax2 = plt.subplot(2, 3, 2)
ax2.plot(results_df['threshold'], results_df['high_conf_total'], 'o-', linewidth=2, markersize=8, color='#A23B72', label='Total')
ax2.plot(results_df['threshold'], results_df['high_conf_ntpm_high'], 's--', linewidth=2, markersize=6, color='#F18F01', label='nTPM ≥ threshold')
ax2.plot(results_df['threshold'], results_df['high_conf_ntpm_low'], '^--', linewidth=2, markersize=6, color='#C73E1D', label='nTPM < threshold')
ax2.set_xlabel('nTPM Threshold', fontsize=12, fontweight='bold')
ax2.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
ax2.set_title('High Confidence Genes (3 Criteria)', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

# Plot 3: Known liver markers coverage
ax3 = plt.subplot(2, 3, 3)
ax3.plot(results_df['threshold'], results_df['known_markers_coverage'], 'o-', linewidth=2, markersize=8, color='#16DB93', label='Any category')
ax3.plot(results_df['threshold'], results_df['known_markers_high_conf_rate'], 's-', linewidth=2, markersize=8, color='#048A81', label='High confidence')
ax3.set_xlabel('nTPM Threshold', fontsize=12, fontweight='bold')
ax3.set_ylabel('Coverage (%)', fontsize=12, fontweight='bold')
ax3.set_title(f'Known Liver Markers Coverage (n={len(KNOWN_LIVER_MARKERS)})', fontsize=14, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.axhline(y=100, color='gray', linestyle=':', alpha=0.5)

# Plot 4: nTPM distribution balance
ax4 = plt.subplot(2, 3, 4)
width = 0.35
x = np.arange(len(results_df))
ax4.bar(x - width/2, results_df['ntpm_high_count'], width, label='nTPM ≥ threshold', color='#FF6B6B', alpha=0.8)
ax4.bar(x + width/2, results_df['ntpm_low_count'], width, label='nTPM < threshold', color='#4169E1', alpha=0.8)
ax4.set_xlabel('nTPM Threshold', fontsize=12, fontweight='bold')
ax4.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
ax4.set_title('nTPM Value Distribution by Threshold', fontsize=14, fontweight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels(results_df['threshold'])
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3, axis='y')

# Plot 5: Confidence distribution
ax5 = plt.subplot(2, 3, 5)
width = 0.25
x = np.arange(len(results_df))
ax5.bar(x - width, results_df['high_conf_total'], width, label='High (3 criteria)', color='#27AE60', alpha=0.8)
ax5.bar(x, results_df['med_conf_total'], width, label='Medium (2 criteria)', color='#F39C12', alpha=0.8)
ax5.bar(x + width, results_df['low_conf_total'], width, label='Low (1 criterion)', color='#E74C3C', alpha=0.8)
ax5.set_xlabel('nTPM Threshold', fontsize=12, fontweight='bold')
ax5.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
ax5.set_title('Confidence Distribution', fontsize=14, fontweight='bold')
ax5.set_xticks(x)
ax5.set_xticklabels(results_df['threshold'])
ax5.legend(fontsize=10)
ax5.grid(True, alpha=0.3, axis='y')

# Plot 6: Summary metrics table
ax6 = plt.subplot(2, 3, 6)
ax6.axis('off')

# Create summary table
summary_data = []
for _, row in results_df.iterrows():
    summary_data.append([
        f"{row['threshold']}",
        f"{row['liver_genes']}",
        f"{row['high_conf_total']}",
        f"{row['known_markers_found']}/{row['known_markers_total']}",
        f"{row['known_markers_high_conf']}",
    ])

table = ax6.table(cellText=summary_data,
                  colLabels=['Threshold', 'Liver\nGenes', 'High\nConf', 'Known\nMarkers', 'Markers\nin High'],
                  cellLoc='center',
                  loc='center',
                  bbox=[0, 0, 1, 1])

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)

# Color header
for i in range(5):
    table[(0, i)].set_facecolor('#34495E')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Color rows alternately
for i in range(1, len(summary_data) + 1):
    for j in range(5):
        if i % 2 == 0:
            table[(i, j)].set_facecolor('#ECF0F1')

ax6.set_title('Summary Table', fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig(THRESHOLD_ANALYSIS_OUTPUT, dpi=DPI, bbox_inches='tight')
print(f"✓ Visualization saved to: {THRESHOLD_ANALYSIS_OUTPUT}")

# Print recommendation
print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)

# Find best threshold based on known markers coverage
best_by_markers = results_df.loc[results_df['known_markers_coverage'].idxmax()]
best_by_high_conf = results_df.loc[results_df['known_markers_high_conf_rate'].idxmax()]

print(f"\nBest threshold by known markers coverage: {best_by_markers['threshold']}")
print(f"  - Coverage: {best_by_markers['known_markers_coverage']:.1f}%")
print(f"  - High confidence markers: {best_by_markers['known_markers_high_conf']}")
print(f"  - Total liver genes: {best_by_markers['liver_genes']}")

print(f"\nBest threshold by known markers in high confidence: {best_by_high_conf['threshold']}")
print(f"  - High confidence markers: {best_by_high_conf['known_markers_high_conf']} ({best_by_high_conf['known_markers_high_conf_rate']:.1f}%)")
print(f"  - Total coverage: {best_by_high_conf['known_markers_coverage']:.1f}%")
print(f"  - Total liver genes: {best_by_high_conf['liver_genes']}")

print("\n" + "="*80)
print("✓ THRESHOLD OPTIMIZATION ANALYSIS COMPLETED")
print("="*80)
