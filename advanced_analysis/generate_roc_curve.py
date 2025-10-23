import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.metrics import confusion_matrix, classification_report

# Read the analysis data
print("Loading classification data...")
trace_data = pd.read_csv('Trace/classification_strategy_analysis.csv')

# Use has_liver_in_ntpm column (already exists in file)
# Fallback: create combined nTPM column if not exists
if 'has_liver_in_ntpm' in trace_data.columns:
    trace_data['has_liver_nTPM'] = trace_data['has_liver_in_ntpm']
else:
    trace_data['has_liver_nTPM'] = (trace_data.get('has_liver_nTPM_high', pd.Series([False]*len(trace_data))) == True) | (trace_data.get('has_liver_nTPM_low', pd.Series([False]*len(trace_data))) == True)

# Define known liver markers as ground truth
known_liver_markers = [
    'ALB', 'HP', 'APOA1', 'APOA2', 'APOC3', 'SERPINA1', 'FGA', 'FGB', 'FGG',
    'CYP3A4', 'CYP2E1', 'FABP1', 'ALDOB', 'ORM1', 'ORM2', 'GC', 'TF',
    'APOH', 'APOE', 'CRP', 'SAA1', 'SAA2', 'RBP4', 'AMBP', 'SELENOP',
    'HRG', 'F2', 'FTL', 'AHSG', 'C3', 'HPX', 'VTN', 'ITIH4', 'ITIH3',
    'FETUB', 'KNG1', 'SERPINC1', 'SERPIND1', 'F9', 'F10', 'F11',
    'CYP1A2', 'CYP2A6', 'CYP2C8', 'CYP2C9', 'CYP2D6',
    'UGT1A1', 'UGT2B7', 'GSTA1', 'GSTA2', 'G6PC', 'PCK1'
]

# Create ground truth labels
trace_data['is_known_liver'] = trace_data['Gene'].isin(known_liver_markers).astype(int)

known_count = trace_data['is_known_liver'].sum()
print(f"\nKnown liver markers in dataset: {known_count}")
print(f"Total genes: {len(trace_data)}")

# Create binary predictions based on different criteria
# We'll test multiple classification strategies

# Strategy 1: Confidence score threshold
trace_data['pred_confidence'] = (trace_data['confidence_score'] >= 60).astype(int)

# Strategy 2: nTPM thresholds (different values)
trace_data['pred_ntpm_0'] = (trace_data['liver_nTPM_value'] >= 0).astype(int)
trace_data['pred_ntpm_0'] = trace_data['pred_ntpm_0'].fillna(0)

trace_data['pred_ntpm_10'] = (trace_data['liver_nTPM_value'] >= 10).astype(int)
trace_data['pred_ntpm_10'] = trace_data['pred_ntpm_10'].fillna(0)

# Strategy 3: Multi-criteria (2 or more)
trace_data['criteria_count'] = (
    trace_data['has_liver_nTPM'].fillna(False).astype(int) +
    trace_data['has_liver_cluster'].fillna(False).astype(int) +
    trace_data['has_liver_enrichment'].fillna(False).astype(int)
)
trace_data['pred_multicriteria'] = (trace_data['criteria_count'] >= 2).astype(int)

# Strategy 4: Tier 1 classification
trace_data['pred_tier1'] = (trace_data['recommended_tier'] == 'Tier 1 (Highest)').astype(int)

# Strategy 5: nTPM + cluster
trace_data['pred_ntpm_cluster'] = (
    (trace_data['has_liver_nTPM'] == True) &
    (trace_data['has_liver_cluster'] == True)
).astype(int)

# Create figure with multiple subplots
fig = plt.figure(figsize=(20, 13))

# ============================================================================
# 1. ROC Curves Comparison (Top Left)
# ============================================================================
ax1 = plt.subplot(2, 3, 1)

strategies = {
    'Confidence Score ≥60': ('confidence_score', trace_data['confidence_score'].fillna(0)),
    'nTPM Threshold ≥0 (any nTPM)': ('pred_ntpm_0', trace_data['liver_nTPM_value'].fillna(-1).apply(lambda x: max(0, x))),
    'nTPM Threshold ≥10': ('pred_ntpm_10', trace_data['liver_nTPM_value'].fillna(0)),
    'Multi-Criteria (≥2)': ('criteria_count', trace_data['criteria_count']),
    'nTPM+Cluster (any nTPM)': ('pred_ntpm_cluster', trace_data['pred_ntpm_cluster']),
    'nTPM+Cluster (≥200)': ('pred_tier1', trace_data['pred_tier1'])
}

colors = ['#FF6B6B', '#4ECDC4', '#A569BD', '#FFE66D', '#95E1D3', '#FF8C42']
y_true = trace_data['is_known_liver']

roc_results = {}

for (label, (col, scores)), color in zip(strategies.items(), colors):
    fpr, tpr, thresholds = roc_curve(y_true, scores)
    roc_auc = auc(fpr, tpr)
    roc_results[label] = {'fpr': fpr, 'tpr': tpr, 'auc': roc_auc}

    ax1.plot(fpr, tpr, color=color, linewidth=2.5,
             label=f'{label} (AUC = {roc_auc:.3f})')

ax1.plot([0, 1], [0, 1], 'k--', linewidth=2, label='Random Classifier')
ax1.set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
ax1.set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
ax1.set_title('ROC Curves: Strategy Comparison', fontsize=14, fontweight='bold', pad=15)
ax1.legend(loc='lower right', fontsize=9.5, framealpha=0.95)
ax1.grid(True, alpha=0.3)

# ============================================================================
# 2. Precision-Recall Curves (Top Middle)
# ============================================================================
ax2 = plt.subplot(2, 3, 2)

for (label, (col, scores)), color in zip(strategies.items(), colors):
    precision, recall, thresholds = precision_recall_curve(y_true, scores)
    avg_precision = average_precision_score(y_true, scores)

    ax2.plot(recall, precision, color=color, linewidth=2.5,
             label=f'{label} (AP = {avg_precision:.3f})')

ax2.set_xlabel('Recall', fontsize=12, fontweight='bold')
ax2.set_ylabel('Precision', fontsize=12, fontweight='bold')
ax2.set_title('Precision-Recall Curves', fontsize=14, fontweight='bold', pad=15)
ax2.legend(loc='upper right', fontsize=8.5, framealpha=0.95)
ax2.grid(True, alpha=0.3)

# ============================================================================
# 3. Confusion Matrices (Top Right - nTPM Threshold ≥10)
# ============================================================================
ax3 = plt.subplot(2, 3, 3)

y_pred_ntpm = trace_data['pred_ntpm_10']
cm = confusion_matrix(y_true, y_pred_ntpm)

im = ax3.imshow(cm, cmap='Blues', alpha=0.8)
ax3.set_xticks([0, 1])
ax3.set_yticks([0, 1])
ax3.set_xticklabels(['Not Liver', 'Liver'], fontsize=11)
ax3.set_yticklabels(['Not Liver', 'Liver'], fontsize=11)
ax3.set_xlabel('Predicted', fontsize=12, fontweight='bold')
ax3.set_ylabel('True (Known Markers)', fontsize=12, fontweight='bold')
ax3.set_title('Confusion Matrix: nTPM Threshold ≥10\n(Best Overall Performance)', fontsize=13, fontweight='bold', pad=15)

# Add text annotations
for i in range(2):
    for j in range(2):
        text = ax3.text(j, i, cm[i, j], ha="center", va="center",
                       fontsize=20, fontweight='bold',
                       color="white" if cm[i, j] > cm.max()/2 else "black")

# Add colorbar
cbar = plt.colorbar(im, ax=ax3)
cbar.set_label('Count', fontsize=10)

# ============================================================================
# 4. Performance Metrics Comparison (Bottom Left)
# ============================================================================
ax4 = plt.subplot(2, 3, 4)

metrics_data = []
for label, (col, scores) in strategies.items():
    if label == 'nTPM Threshold ≥0 (any nTPM)':
        y_pred = trace_data['pred_ntpm_0']
    elif label == 'nTPM Threshold ≥10':
        y_pred = trace_data['pred_ntpm_10']
    elif label == 'Multi-Criteria (≥2)':
        y_pred = trace_data['pred_multicriteria']
    elif label == 'nTPM+Cluster (≥200)':
        y_pred = trace_data['pred_tier1']
    elif label == 'nTPM+Cluster (any nTPM)':
        y_pred = trace_data['pred_ntpm_cluster']
    else:
        y_pred = trace_data['pred_confidence']

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1 = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0

    metrics_data.append({
        'Strategy': label.replace(' ', '\n'),
        'Sensitivity': sensitivity,
        'Specificity': specificity,
        'Precision': precision,
        'F1-Score': f1
    })

metrics_df = pd.DataFrame(metrics_data)

x = np.arange(len(metrics_df))
width = 0.2

bars1 = ax4.bar(x - 1.5*width, metrics_df['Sensitivity'], width, label='Sensitivity', color='#FF6B6B', alpha=0.8)
bars2 = ax4.bar(x - 0.5*width, metrics_df['Specificity'], width, label='Specificity', color='#4ECDC4', alpha=0.8)
bars3 = ax4.bar(x + 0.5*width, metrics_df['Precision'], width, label='Precision', color='#FFE66D', alpha=0.8)
bars4 = ax4.bar(x + 1.5*width, metrics_df['F1-Score'], width, label='F1-Score', color='#95E1D3', alpha=0.8)

ax4.set_ylabel('Score', fontsize=12, fontweight='bold')
ax4.set_title('Performance Metrics Comparison', fontsize=14, fontweight='bold', pad=15)
ax4.set_xticks(x)
ax4.set_xticklabels(metrics_df['Strategy'], fontsize=7.5, rotation=0, ha='center')
ax4.legend(loc='upper left', fontsize=9.5, framealpha=0.95)
ax4.set_ylim(0, 1.1)
ax4.grid(axis='y', alpha=0.3)

# ============================================================================
# 5. Score Distribution (Bottom Middle)
# ============================================================================
ax5 = plt.subplot(2, 3, 5)

known_liver = trace_data[trace_data['is_known_liver'] == 1]['confidence_score']
non_liver = trace_data[trace_data['is_known_liver'] == 0]['confidence_score']

ax5.hist(non_liver, bins=30, alpha=0.6, label='Non-liver proteins', color='#E8E8E8', edgecolor='black')
ax5.hist(known_liver, bins=15, alpha=0.8, label='Known liver markers', color='#FF6B6B', edgecolor='black')

ax5.axvline(60, color='red', linestyle='--', linewidth=2.5, label='Threshold = 60')
ax5.set_xlabel('Confidence Score', fontsize=12, fontweight='bold')
ax5.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax5.set_title('Confidence Score Distribution', fontsize=14, fontweight='bold', pad=15)
ax5.legend(loc='upper right', fontsize=9.5, framealpha=0.95)
ax5.grid(axis='y', alpha=0.3)

# ============================================================================
# 6. Performance Summary Table (Bottom Right)
# ============================================================================
ax6 = plt.subplot(2, 3, 6)
ax6.axis('off')

# Create summary table
summary_data = []
for label, (col, scores) in strategies.items():
    roc_auc = roc_results[label]['auc']
    avg_precision = average_precision_score(y_true, scores)

    if label == 'nTPM Threshold ≥0 (any nTPM)':
        y_pred = trace_data['pred_ntpm_0']
    elif label == 'nTPM Threshold ≥10':
        y_pred = trace_data['pred_ntpm_10']
    elif label == 'Multi-Criteria (≥2)':
        y_pred = trace_data['pred_multicriteria']
    elif label == 'nTPM+Cluster (≥200)':
        y_pred = trace_data['pred_tier1']
    elif label == 'nTPM+Cluster (any nTPM)':
        y_pred = trace_data['pred_ntpm_cluster']
    else:
        y_pred = trace_data['pred_confidence']

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    summary_data.append([
        label,
        f"{roc_auc:.3f}",
        f"{avg_precision:.3f}",
        f"{tp}/{tp+fn}",
        f"{tn}/{tn+fp}"
    ])

column_labels = ['Strategy', 'AUC', 'AP', 'TP/P', 'TN/N']
table = ax6.table(cellText=summary_data, colLabels=column_labels,
                  cellLoc='center', loc='center',
                  colWidths=[0.38, 0.14, 0.14, 0.16, 0.18])

table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 2.8)

# Style header
for i in range(len(column_labels)):
    cell = table[(0, i)]
    cell.set_facecolor('#4A90E2')
    cell.set_text_props(weight='bold', color='white', fontsize=10)

# Color rows
row_colors = colors
for i in range(len(summary_data)):
    for j in range(len(column_labels)):
        cell = table[(i+1, j)]
        cell.set_facecolor(row_colors[i])
        cell.set_alpha(0.3)
        cell.set_text_props(fontsize=9)

ax6.set_title('Performance Summary', fontsize=14, fontweight='bold', pad=25)

plt.tight_layout()
plt.savefig('output/roc_auc_analysis.png', dpi=300, bbox_inches='tight')
print("\nROC/AUC analysis saved to: output/roc_auc_analysis.png")

# ============================================================================
# Print detailed metrics
# ============================================================================
print("\n" + "="*80)
print("DETAILED PERFORMANCE METRICS")
print("="*80)

for i, (label, (col, scores)) in enumerate(strategies.items()):
    print(f"\n### {label} ###")

    if label == 'nTPM Threshold ≥0 (any nTPM)':
        y_pred = trace_data['pred_ntpm_0']
    elif label == 'nTPM Threshold ≥10':
        y_pred = trace_data['pred_ntpm_10']
    elif label == 'Multi-Criteria (≥2)':
        y_pred = trace_data['pred_multicriteria']
    elif label == 'nTPM+Cluster (≥200)':
        y_pred = trace_data['pred_tier1']
    elif label == 'nTPM+Cluster (any nTPM)':
        y_pred = trace_data['pred_ntpm_cluster']
    else:
        y_pred = trace_data['pred_confidence']

    print(classification_report(y_true, y_pred,
                                target_names=['Non-liver', 'Known Liver Marker'],
                                digits=3))

    print(f"AUC-ROC: {roc_results[label]['auc']:.3f}")
    print(f"Average Precision: {average_precision_score(y_true, scores):.3f}")

print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)

# Calculate AP scores for all strategies
ap_scores = {}
for label, (col, scores) in strategies.items():
    ap_scores[label] = average_precision_score(y_true, scores)

# Find best strategies
best_auc = max(roc_results.items(), key=lambda x: x[1]['auc'])
best_ap = max(ap_scores.items(), key=lambda x: x[1])

print(f"\nBest by AUC-ROC: {best_auc[0]}")
print(f"  AUC-ROC: {best_auc[1]['auc']:.3f}")

print(f"\nBest by AP (Average Precision): {best_ap[0]}")
print(f"  AP: {best_ap[1]:.3f}")

# Combined recommendation
print("\n" + "-"*80)
print("STRATEGY RECOMMENDATIONS BY USE CASE:")
print("-"*80)
print("\nFor liver protein classification, we recommend:")
print("  1. nTPM+Cluster (≥200) - For highest confidence/specificity")
print("  2. nTPM Threshold ≥10 - For best overall performance (AUC & AP)")
print("  3. Confidence Score ≥60 - For balanced sensitivity/specificity")
print("  4. Multi-Criteria (≥2) - For broadest coverage")

print("\nKey Insight:")
if best_auc[0] == best_ap[0]:
    print(f"  • {best_auc[0]} excels in both AUC-ROC and AP metrics")
else:
    print(f"  • AUC-ROC favors: {best_auc[0]} (better overall discrimination)")
    print(f"  • AP favors: {best_ap[0]} (better for imbalanced data)")

plt.close()
print("\nAnalysis complete!")
