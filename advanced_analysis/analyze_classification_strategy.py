import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# Read the trace data
print("Loading trace data...")
trace_data = pd.read_csv('Trace/Data tracing.csv')

# Create combined nTPM column for compatibility
trace_data['has_liver_in_ntpm'] = (trace_data['has_liver_nTPM_high'] == True) | (trace_data['has_liver_nTPM_low'] == True)

print("\n" + "="*80)
print("LIVER PROTEIN CLASSIFICATION STRATEGY ANALYSIS")
print("="*80)

# ============================================================================
# STRATEGY 1: Current Multi-Criteria Approach
# ============================================================================
print("\n### STRATEGY 1: Multi-Criteria Overlap ###")
print("-" * 80)

current_stats = trace_data['Classification'].value_counts()
print("\nCurrent Classification:")
print(current_stats)

liver_candidates = len(trace_data[trace_data['Classification'].str.startswith('liver protein')])
print(f"\nTotal Liver Candidates: {liver_candidates}")
print(f"Non-liver: {len(trace_data[trace_data['Classification'] == 'non-liver protein'])}")

# ============================================================================
# STRATEGY 2: nTPM Threshold-Based Approach
# ============================================================================
print("\n\n### STRATEGY 2: nTPM Threshold-Based Classification ###")
print("-" * 80)

liver_genes = trace_data[trace_data['liver_nTPM_value'].notna()].copy()
print(f"\nTotal genes with nTPM values: {len(liver_genes)}")

# Analyze nTPM value distribution
ntpm_values = liver_genes['liver_nTPM_value'].values
percentiles = [25, 50, 75, 90, 95, 99]
print("\nnTPM Value Percentiles:")
for p in percentiles:
    print(f"  {p}th percentile: {np.percentile(ntpm_values, p):.2f}")

# Test different thresholds
thresholds = [10, 50, 100, 200, 500, 1000]
print("\nGenes above different nTPM thresholds:")
for threshold in thresholds:
    count = len(liver_genes[liver_genes['liver_nTPM_value'] >= threshold])
    percentage = (count / len(liver_genes)) * 100
    print(f"  nTPM >= {threshold}: {count} genes ({percentage:.1f}%)")

# ============================================================================
# STRATEGY 3: Weighted Scoring System
# ============================================================================
print("\n\n### STRATEGY 3: Weighted Scoring System ###")
print("-" * 80)

def calculate_confidence_score(row):
    """
    Calculate confidence score based on:
    - nTPM value (0-40 points)
    - Number of criteria met (0-30 points)
    - Presence in cluster (0-15 points)
    - Presence in enrichment (0-15 points)
    """
    score = 0

    # nTPM-based scoring (0-40 points)
    if pd.notna(row['liver_nTPM_value']):
        ntpm = row['liver_nTPM_value']
        if ntpm >= 1000:
            score += 40
        elif ntpm >= 500:
            score += 35
        elif ntpm >= 200:
            score += 30
        elif ntpm >= 100:
            score += 25
        elif ntpm >= 50:
            score += 20
        elif ntpm >= 10:
            score += 15
        else:
            score += 10

    # Criteria count (0-30 points)
    # Check if gene has liver in nTPM (either high or low threshold)
    has_any_ntpm = row.get('has_liver_nTPM_high', False) or row.get('has_liver_nTPM_low', False)
    criteria_count = sum([has_any_ntpm, row['has_liver_cluster'], row['has_liver_enrichment']])
    if criteria_count == 3:
        score += 30
    elif criteria_count == 2:
        score += 20
    elif criteria_count == 1:
        score += 10

    # Tissue expression cluster bonus (0-15 points)
    if row['has_liver_cluster']:
        score += 15

    # Cell type enrichment bonus (0-15 points)
    if row['has_liver_enrichment']:
        score += 15

    return score

trace_data['confidence_score'] = trace_data.apply(calculate_confidence_score, axis=1)

# Classify based on confidence scores
def classify_by_confidence(score):
    if score >= 80:
        return 'High Confidence'
    elif score >= 60:
        return 'Medium-High Confidence'
    elif score >= 40:
        return 'Medium Confidence'
    elif score >= 20:
        return 'Low-Medium Confidence'
    else:
        return 'Low Confidence / Non-liver'

trace_data['confidence_class'] = trace_data['confidence_score'].apply(classify_by_confidence)

print("\nConfidence-based Classification:")
print(trace_data['confidence_class'].value_counts().sort_index(ascending=False))

print("\nTop 10 genes by confidence score:")
top_confidence = trace_data.nlargest(10, 'confidence_score')[['Gene', 'liver_nTPM_value', 'Classification', 'confidence_score']]
print(top_confidence.to_string(index=False))

# ============================================================================
# STRATEGY 4: Restrictive High-Stringency Approach
# ============================================================================
print("\n\n### STRATEGY 4: High-Stringency Filtering ###")
print("-" * 80)

# Define high-stringency criteria
high_stringency = trace_data[
    (trace_data['has_liver_in_ntpm'] == True) &
    (trace_data['has_liver_cluster'] == True) &
    (trace_data['liver_nTPM_value'] >= 100)
]

print(f"\nHigh-Stringency Criteria:")
print(f"  - Must have liver in nTPM")
print(f"  - Must have liver in tissue cluster")
print(f"  - nTPM value >= 100")
print(f"\nResult: {len(high_stringency)} genes")

very_high_stringency = trace_data[
    (trace_data['has_liver_in_ntpm'] == True) &
    (trace_data['has_liver_cluster'] == True) &
    (trace_data['has_liver_enrichment'] == True) &
    (trace_data['liver_nTPM_value'] >= 200)
]

print(f"\nVery High-Stringency Criteria:")
print(f"  - All 3 criteria must be met")
print(f"  - nTPM value >= 200")
print(f"\nResult: {len(very_high_stringency)} genes")

# ============================================================================
# STRATEGY 5: Literature/Known Liver Marker Validation
# ============================================================================
print("\n\n### STRATEGY 5: Known Liver Marker Validation ###")
print("-" * 80)

# Known liver-specific markers
known_liver_markers = [
    'ALB', 'HP', 'APOA1', 'APOA2', 'APOC3', 'SERPINA1', 'FGA', 'FGB', 'FGG',
    'CYP3A4', 'CYP2E1', 'FABP1', 'ALDOB', 'ORM1', 'ORM2', 'GC', 'TF',
    'APOH', 'APOE', 'CRP', 'SAA1', 'SAA2', 'RBP4', 'AMBP', 'SELENOP'
]

print(f"\nValidation against {len(known_liver_markers)} known liver markers:")

found_markers = trace_data[trace_data['Gene'].isin(known_liver_markers)]
print(f"  Found in dataset: {len(found_markers)}/{len(known_liver_markers)}")

classified_as_liver = found_markers[found_markers['Classification'].str.startswith('liver protein')]
print(f"  Classified as liver protein: {len(classified_as_liver)}/{len(found_markers)}")

if len(found_markers) > 0:
    print(f"  Capture rate: {(len(classified_as_liver)/len(found_markers))*100:.1f}%")

    print("\nClassification of known markers:")
    marker_class = found_markers.groupby('Classification').size()
    print(marker_class)

    # Show markers NOT classified as liver
    not_classified = found_markers[~found_markers['Classification'].str.startswith('liver protein')]
    if len(not_classified) > 0:
        print(f"\nWarning: {len(not_classified)} known markers NOT classified as liver:")
        print(not_classified[['Gene', 'Classification', 'liver_nTPM_value']].to_string(index=False))

# ============================================================================
# COMPARISON VISUALIZATION
# ============================================================================
print("\n\n### Creating Comparison Visualization ###")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# 1. Current Classification Distribution
ax1 = axes[0, 0]
current_liver = trace_data[trace_data['Classification'].str.startswith('liver protein')]['Classification'].value_counts()
colors1 = ['#FF6B6B', '#4ECDC4', '#FFE66D', '#A8E6CF', '#FF8C42', '#95E1D3', '#C7CEEA']
current_liver.plot(kind='barh', ax=ax1, color=colors1[:len(current_liver)])
ax1.set_xlabel('Number of Genes', fontsize=11, fontweight='bold')
ax1.set_title('Strategy 1: Multi-Criteria Classification', fontsize=12, fontweight='bold')
ax1.grid(axis='x', alpha=0.3)

# 2. Confidence Score Distribution
ax2 = axes[0, 1]
confidence_dist = trace_data['confidence_class'].value_counts().sort_index(ascending=False)
colors2 = ['#D32F2F', '#FF6B6B', '#FFA726', '#FFE66D', '#E8E8E8']
confidence_dist.plot(kind='barh', ax=ax2, color=colors2)
ax2.set_xlabel('Number of Genes', fontsize=11, fontweight='bold')
ax2.set_title('Strategy 3: Confidence Score Classification', fontsize=12, fontweight='bold')
ax2.grid(axis='x', alpha=0.3)

# 3. nTPM Threshold Comparison
ax3 = axes[1, 0]
threshold_data = []
threshold_labels = []
for threshold in [0, 10, 50, 100, 200, 500, 1000]:
    count = len(liver_genes[liver_genes['liver_nTPM_value'] >= threshold])
    threshold_data.append(count)
    threshold_labels.append(f'≥{threshold}')

bars = ax3.barh(threshold_labels, threshold_data, color='#4ECDC4', edgecolor='black', linewidth=1.5)
ax3.set_xlabel('Number of Genes', fontsize=11, fontweight='bold')
ax3.set_title('Strategy 2: nTPM Threshold Impact', fontsize=12, fontweight='bold')
ax3.grid(axis='x', alpha=0.3)

# Add value labels
for i, (bar, val) in enumerate(zip(bars, threshold_data)):
    ax3.text(val + max(threshold_data)*0.01, i, f'{val}', va='center', fontweight='bold')

# 4. High-Stringency Tiers
ax4 = axes[1, 1]
stringency_data = {
    'All liver\ncandidates\n(any criterion)': liver_candidates,
    'nTPM only\n(any value)': len(trace_data[trace_data['has_liver_in_ntpm'] == True]),
    'nTPM + cluster\n(any value)': len(trace_data[(trace_data['has_liver_in_ntpm'] == True) &
                                                    (trace_data['has_liver_cluster'] == True)]),
    'nTPM + cluster\n(≥100)': len(high_stringency),
    'All 3 criteria\n(≥200)': len(very_high_stringency)
}

y_pos = np.arange(len(stringency_data))
colors4 = ['#E8E8E8', '#C7CEEA', '#95E1D3', '#4ECDC4', '#FF6B6B']
bars4 = ax4.barh(y_pos, list(stringency_data.values()), color=colors4, edgecolor='black', linewidth=1.5)
ax4.set_yticks(y_pos)
ax4.set_yticklabels(list(stringency_data.keys()), fontsize=9)
ax4.set_xlabel('Number of Genes', fontsize=11, fontweight='bold')
ax4.set_title('Strategy 4: Progressive Stringency Filtering', fontsize=12, fontweight='bold')
ax4.grid(axis='x', alpha=0.3)

# Add value labels
for i, (bar, val) in enumerate(zip(bars4, stringency_data.values())):
    ax4.text(val + max(stringency_data.values())*0.01, i, f'{val}', va='center', fontweight='bold')

plt.tight_layout()
plt.savefig('output/classification_strategy_comparison.png', dpi=300, bbox_inches='tight')
print("Visualization saved to: output/classification_strategy_comparison.png")
plt.close()  # Free memory

# ============================================================================
# RECOMMENDATIONS
# ============================================================================
print("\n\n" + "="*80)
print("RECOMMENDATIONS FOR HIGH-CONFIDENCE LIVER PROTEIN CLASSIFICATION")
print("="*80)

print("\n### TIER-BASED RECOMMENDATION ###\n")

# Tier 1: Highest confidence
tier1 = trace_data[
    (trace_data['has_liver_in_ntpm'] == True) &
    (trace_data['has_liver_cluster'] == True) &
    (trace_data['liver_nTPM_value'] >= 200)
]
print(f"TIER 1 (Highest Confidence): {len(tier1)} genes")
print(f"  Criteria: nTPM present + Cluster present + nTPM ≥ 200")
print(f"  Recommended for: Primary liver-specific protein candidates")

# Tier 2: High confidence
tier2 = trace_data[
    (trace_data['has_liver_in_ntpm'] == True) &
    (trace_data['has_liver_cluster'] == True) &
    (trace_data['liver_nTPM_value'] >= 50) &
    (~trace_data.index.isin(tier1.index))
]
print(f"\nTIER 2 (High Confidence): {len(tier2)} genes")
print(f"  Criteria: nTPM present + Cluster present + nTPM ≥ 50")
print(f"  Recommended for: Secondary liver-specific candidates")

# Tier 3: Medium confidence
tier3 = trace_data[
    (trace_data['has_liver_in_ntpm'] == True) &
    (trace_data['has_liver_cluster'] == True) &
    (~trace_data.index.isin(tier1.index)) &
    (~trace_data.index.isin(tier2.index))
]
print(f"\nTIER 3 (Medium Confidence): {len(tier3)} genes")
print(f"  Criteria: nTPM present + Cluster present (any nTPM value)")
print(f"  Recommended for: Potential liver-associated proteins")

# Tier 4: Exploratory
tier4 = trace_data[
    (trace_data['Classification'] == 'liver protein (all 3)') &
    (~trace_data.index.isin(tier1.index)) &
    (~trace_data.index.isin(tier2.index)) &
    (~trace_data.index.isin(tier3.index))
]
print(f"\nTIER 4 (Exploratory): {len(tier4)} genes")
print(f"  Criteria: All 3 criteria met but low nTPM OR enrichment-based")
print(f"  Recommended for: Exploratory validation studies")

# Save tier assignments
trace_data['recommended_tier'] = 'Not recommended'
trace_data.loc[tier1.index, 'recommended_tier'] = 'Tier 1 (Highest)'
trace_data.loc[tier2.index, 'recommended_tier'] = 'Tier 2 (High)'
trace_data.loc[tier3.index, 'recommended_tier'] = 'Tier 3 (Medium)'
trace_data.loc[tier4.index, 'recommended_tier'] = 'Tier 4 (Exploratory)'

# Save results
output_file = 'Trace/classification_strategy_analysis.csv'
trace_data.to_csv(output_file, index=False)
print(f"\n\nDetailed analysis saved to: {output_file}")

# Save tier-specific gene lists
for tier_name, tier_df in [('tier1', tier1), ('tier2', tier2), ('tier3', tier3), ('tier4', tier4)]:
    if len(tier_df) > 0:
        tier_file = f'Trace/{tier_name}_high_confidence_genes.csv'
        tier_df[['Gene', 'liver_nTPM_value', 'Classification', 'confidence_score']].to_csv(tier_file, index=False)
        print(f"  {tier_name} genes saved to: {tier_file}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
