# Liver Protein Classification Strategy Report

## Executive Summary

This report analyzes **5 different strategies** for classifying liver-specific proteins with varying levels of stringency and reliability. Based on comprehensive analysis including validation against 24 known liver markers (100% capture rate), we provide **tier-based recommendations** for selecting high-confidence liver protein candidates.

---

## Strategy Comparison

### Strategy 1: Multi-Criteria Overlap (Current Approach)
**Total Liver Candidates: 1,593 genes**

**Approach**: Classify based on presence of "liver" in three independent criteria
- RNA tissue specific nTPM
- Tissue expression cluster
- RNA tissue cell type enrichment

**Results**:
- liver protein (all 3): 391 genes
- liver protein (nTPM + cluster): 184 genes
- liver protein (nTPM + enrichment): 71 genes
- liver protein (cluster + enrichment): 57 genes
- liver protein_1 (nTPM only): 61 genes
- liver protein_2 (cluster only): 157 genes
- liver protein_3 (enrichment only): 672 genes

**Pros**:
- Comprehensive coverage
- Multiple independent lines of evidence
- Captures diverse liver-associated proteins

**Cons**:
- Includes many low-confidence candidates
- No quantitative threshold
- High false positive risk for single-criterion matches

---

### Strategy 2: nTPM Threshold-Based Classification
**Focus**: Expression level as primary filter

**nTPM Distribution Analysis** (707 genes with nTPM values):
- 25th percentile: 88.30
- Median: 207.40
- 75th percentile: 571.65
- 90th percentile: 1,807.10
- 95th percentile: 5,086.98

**Threshold Impact**:
| Threshold | Genes | Percentage |
|-----------|-------|------------|
| ≥10       | 701   | 99.2%      |
| ≥50       | 614   | 86.8%      |
| ≥100      | 509   | 72.0%      |
| ≥200      | 356   | 50.4%      |
| ≥500      | 199   | 28.1%      |
| ≥1000     | 114   | 16.1%      |

**Pros**:
- Quantitative, objective cutoff
- Directly measures expression level
- Easy to adjust stringency

**Cons**:
- Ignores qualitative evidence (cluster, enrichment)
- May miss important low-abundance liver proteins
- Threshold selection is somewhat arbitrary

---

### Strategy 3: Weighted Confidence Scoring System
**Approach**: Composite score (0-100 points) based on:
- nTPM value (0-40 points)
- Number of criteria met (0-30 points)
- Tissue cluster presence (0-15 points)
- Cell type enrichment presence (0-15 points)

**Results**:
- High Confidence (≥80): 358 genes
- Medium-High Confidence (60-79): 203 genes
- Medium Confidence (40-59): 151 genes
- Low-Medium Confidence (20-39): 881 genes
- Low Confidence (<20): 4,930 genes

**Top Scoring Genes** (score = 100):
- HSD17B6, CYB5A, ADH1B, ADH1C, ALDH1A1, PAH, SOD1, CP, F2, C1R

**Pros**:
- Integrates multiple factors
- Provides graded confidence levels
- Flexible weighting system

**Cons**:
- Somewhat complex
- Weighting choices affect results
- Less intuitive than simple thresholds

---

### Strategy 4: High-Stringency Filtering
**Approach**: Progressive filtering with increasingly strict criteria

**Stringency Levels**:

1. **Basic**: nTPM + Cluster (any value) → **575 genes**
2. **High**: nTPM + Cluster (≥100) → **436 genes**
3. **Very High**: All 3 criteria + nTPM ≥200 → **239 genes**

**Pros**:
- Clear, interpretable criteria
- High specificity
- Low false positive rate

**Cons**:
- May miss some true liver proteins
- Rigid cutoffs
- Lower sensitivity

---

### Strategy 5: Known Liver Marker Validation
**Validation Set**: 25 well-established liver markers

**Results**:
- Found in dataset: 24/25 (96%)
- Correctly classified: 24/24 (100% capture rate)
- Classification breakdown:
  - liver protein (all 3): 21 markers
  - liver protein (nTPM + cluster): 3 markers

**Known Markers Validated**:
ALB, HP, APOA1, APOA2, APOC3, SERPINA1, FGA, FGB, FGG, CYP3A4, CYP2E1, FABP1, ALDOB, ORM1, ORM2, GC, TF, APOH, APOE, CRP, SAA2, RBP4, AMBP, SELENOP

**Key Finding**:
✅ **All known liver markers are captured**, with most (87.5%) meeting the strictest criteria (all 3)

**Conclusion**:
The classification system successfully identifies established liver proteins, validating the approach.

---

## Recommended Tier-Based Classification System

Based on comprehensive analysis, we recommend a **4-tier system** combining multiple criteria:

### 🔴 TIER 1: Highest Confidence (319 genes)
**Criteria**:
- ✅ nTPM present
- ✅ Tissue cluster present
- ✅ nTPM ≥ 200

**Recommended Use**: Primary liver-specific protein candidates for functional studies

**Representative Genes**: ALB (205,952), HP (54,575), ORM1 (41,892), APOA2 (34,742), SERPINA1 (33,156)

**File**: `Trace/tier1_high_confidence_genes.csv`

---

### 🟠 TIER 2: High Confidence (193 genes)
**Criteria**:
- ✅ nTPM present
- ✅ Tissue cluster present
- ✅ nTPM ≥ 50

**Recommended Use**: Secondary liver-specific candidates, validation studies

**File**: `Trace/tier2_high_confidence_genes.csv`

---

### 🟡 TIER 3: Medium Confidence (63 genes)
**Criteria**:
- ✅ nTPM present
- ✅ Tissue cluster present
- ⚠️ Any nTPM value

**Recommended Use**: Potential liver-associated proteins requiring additional validation

**File**: `Trace/tier3_high_confidence_genes.csv`

---

### ⚪ TIER 4: Exploratory
**Criteria**:
- Single criterion OR
- Low nTPM values

**Recommended Use**: Exploratory research only, high validation burden

---

## Key Insights from Analysis

### 1. Criteria Importance Ranking
Based on known marker validation and nTPM distribution:

1. **RNA tissue specific nTPM** (Primary) - Direct expression evidence
2. **Tissue expression cluster** (Secondary) - Functional context
3. **RNA tissue cell type enrichment** (Supporting) - Cell-type specificity

### 2. Optimal nTPM Threshold
- **Median of liver proteins**: 207.40
- **Recommended threshold**: ≥200 for high confidence
- **Alternative threshold**: ≥50 for broader coverage

### 3. Multi-Criteria Advantage
Genes meeting multiple criteria show:
- Higher average nTPM values
- Better overlap with known markers
- Greater consistency across data types

### 4. False Positive Risk
- Single criterion: HIGH risk (especially enrichment-only)
- Two criteria: MODERATE risk
- Three criteria + high nTPM: LOW risk

---

## Practical Recommendations

### For Biomarker Discovery
✅ **Use TIER 1** (319 genes)
- Highest specificity
- Best supported by evidence
- Suitable for clinical applications

### For Functional Studies
✅ **Use TIER 1 + TIER 2** (512 genes)
- Good balance of coverage and confidence
- Includes well-characterized and novel candidates

### For Comprehensive Analysis
✅ **Use All Multi-Criteria Groups** (646 genes)
- Genes meeting ≥2 criteria
- Reasonable confidence level
- Broader biological coverage

### For Exploratory Research
✅ **Consider Single-Criterion Groups**
- Requires extensive validation
- May include tissue-resident immune cells
- Useful for hypothesis generation

---

## Statistical Summary

| Category | Total Genes | % of Dataset |
|----------|-------------|--------------|
| Total analyzed | 6,523 | 100% |
| Any liver criterion | 1,593 | 24.4% |
| ≥2 criteria met | 646 | 9.9% |
| All 3 criteria | 391 | 6.0% |
| **TIER 1 (Recommended)** | **319** | **4.9%** |
| Known markers captured | 24/24 | 100% |

---

## Conclusion

The **TIER-based classification system** provides a robust framework for selecting liver protein candidates with appropriate confidence levels for different research applications.

**Key Recommendations**:

1. ✅ **Prioritize TIER 1** for high-confidence liver-specific proteins
2. ✅ **Combine nTPM threshold (≥200) with multi-criteria approach** for optimal specificity
3. ✅ **Validate single-criterion candidates** before downstream applications
4. ✅ **Use confidence scores** for ranking within tiers

The 100% capture rate of known liver markers validates the classification approach and supports its use for identifying novel liver-specific proteins.

---

## Files Generated

- `output/classification_strategy_comparison.png` - Visual comparison of all strategies
- `output/liver_ntpm_group_distribution.png` - Group-wise nTPM distribution analysis
- `Trace/classification_strategy_analysis.csv` - Complete analysis with confidence scores
- `Trace/tier1_high_confidence_genes.csv` - Top priority candidates (319 genes)
- `Trace/tier2_high_confidence_genes.csv` - Secondary candidates (193 genes)
- `Trace/tier3_high_confidence_genes.csv` - Medium confidence (63 genes)

---

**Analysis Date**: 2025-10-22
**Total Genes Analyzed**: 6,523
**Recommended High-Confidence Candidates**: 319 (TIER 1)
