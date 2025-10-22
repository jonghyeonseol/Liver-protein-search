# Liver Protein Classification Tool

This tool classifies proteins based on liver-specific expression data from the Human Protein Atlas, with advanced nTPM threshold analysis.

## Features

- **4-Criteria Classification**: Classifies proteins based on four liver-related criteria
  - RNA tissue specific nTPM ≥100 (liver mentioned)
  - RNA tissue specific nTPM <100 (liver mentioned)
  - Tissue expression cluster
  - RNA tissue cell type enrichment

- **Automated Pipeline**: One-command execution of entire analysis workflow
- **Comprehensive Visualization**: Dual Venn diagrams, distribution plots, and group statistics
- **Data Tracing**: Complete tracking of all classification decisions
- **Advanced Analysis**: ROC curves and classification strategy evaluation

## Requirements

```bash
pip install pandas matplotlib matplotlib-venn seaborn numpy
```

## Quick Start

### Option 1: Run Full Pipeline (Recommended)

```bash
# Full analysis with advanced metrics
python3 run_pipeline.py

# Quick mode (classification + basic visualization only)
python3 run_pipeline.py --quick

# Visualization only (requires existing data)
python3 run_pipeline.py --viz-only
```

### Option 2: Run Individual Steps

```bash
# Step 1: Classification
python3 liver_protein_classifier.py

# Step 2: All visualizations
python3 run_visualizations.py

# Or run visualizations individually
python3 visualize_ntpm.py
python3 visualize_group_distribution.py
```

## Input Files

Place the following files in the `input/` folder:

1. **Any CSV file** containing gene expression data with a "Gene" column
   - Example: `log2_center_normalization_t_test_Human.csv`
   - If multiple CSV files exist, the first one will be used

2. **proteinatlas.tsv** - Human Protein Atlas data
   - Must contain columns: Gene, RNA tissue specific nTPM, Tissue expression cluster, RNA tissue cell type enrichment

## Pipeline Structure

```
┌─────────────────────────────────────────────────────────────┐
│                   LIVER PROTEIN ANALYSIS PIPELINE            │
└─────────────────────────────────────────────────────────────┘

┌──────────────────────┐
│  1. CLASSIFICATION   │  (liver_protein_classifier.py)
│  ───────────────────│
│  • Load input data   │
│  • Extract nTPM      │
│  • Split by nTPM≥100 │
│  • Generate 4-way    │
│    Venn diagram      │
└──────┬───────────────┘
       │
       ▼
┌──────────────────────┐
│  2. VISUALIZATION    │  (run_visualizations.py)
│  ───────────────────│
│  A. nTPM Analysis    │  (visualize_ntpm.py)
│     • Top 50 genes   │
│     • Distribution   │
│                      │
│  B. Group Analysis   │  (visualize_group_distribution.py)
│     • Box plots      │
│     • Violin plots   │
│     • Strip plots    │
│     • Histograms     │
│     • CDF curves     │
│     • Statistics     │
└──────┬───────────────┘
       │
       ▼
┌──────────────────────┐
│  3. ADVANCED         │  (Optional)
│  ───────────────────│
│  • ROC Curves        │  (generate_roc_curve.py)
│  • Strategy Analysis │  (analyze_classification_strategy.py)
└──────────────────────┘
```

## Output Files

### output/ folder (Visualizations)
- `liver_protein_venn_diagram.png` - Dual 3-way Venn diagrams (nTPM≥100 vs <100)
- `liver_ntpm_distribution.png` - nTPM distribution and top 50 genes
- `liver_ntpm_group_distribution.png` - Comprehensive 6-subplot group analysis
- `liver_protein_roc_curve.png` - ROC curve analysis (if advanced analysis enabled)

### Trace/ folder (Data)
- `Data tracing.csv` - Complete dataset with all classifications
- `classification_summary.csv` - Summary statistics
- `liver_genes_sorted_by_ntpm.csv` - All genes sorted by nTPM value

**Classification-specific gene lists:**
- `liver_protein_nTPM≥100_and_cluster_and_enrichment_genes.csv` - 309 genes (highest confidence)
- `liver_protein_nTPM≥100_and_cluster_genes.csv` - 127 genes
- `liver_protein_nTPM≥100_and_enrichment_genes.csv` - 43 genes
- `liver_protein_nTPM≥100_only_genes.csv` - 30 genes
- `liver_protein_nTPM<100_and_cluster_and_enrichment_genes.csv` - 82 genes
- `liver_protein_nTPM<100_and_cluster_genes.csv` - 57 genes
- `liver_protein_nTPM<100_and_enrichment_genes.csv` - 28 genes
- `liver_protein_nTPM<100_only_genes.csv` - 31 genes
- `liver_protein_cluster_only_genes.csv` - 157 genes
- `liver_protein_enrichment_only_genes.csv` - 672 genes
- `liver_protein_cluster_and_enrichment_genes.csv` - 57 genes

## Classification Logic

### 4-Criteria System:

The classification splits liver nTPM data by a threshold of 100:

1. **High nTPM Categories (nTPM ≥100)**:
   - `nTPM≥100 + cluster + enrichment` - 309 genes (HIGHEST CONFIDENCE)
   - `nTPM≥100 + cluster` - 127 genes
   - `nTPM≥100 + enrichment` - 43 genes
   - `nTPM≥100 only` - 30 genes

2. **Low nTPM Categories (nTPM <100)**:
   - `nTPM<100 + cluster + enrichment` - 82 genes
   - `nTPM<100 + cluster` - 57 genes
   - `nTPM<100 + enrichment` - 28 genes
   - `nTPM<100 only` - 31 genes

3. **No nTPM Categories**:
   - `cluster + enrichment` - 57 genes
   - `cluster only` - 157 genes
   - `enrichment only` - 672 genes

4. **Non-liver protein**: Does not meet any criteria

### Statistics (Example Dataset):

- **Total proteins**: 6,517
- **Liver candidates**: 1,592 (24.4%)
- **Non-liver proteins**: 4,925 (75.6%)
- **With nTPM ≥100**: 509 genes
- **With nTPM <100**: 198 genes

## Key Results

### Top Liver Proteins (by nTPM):
1. ALB (205,952.1 nTPM) - nTPM≥100 + cluster + enrichment
2. HP (54,575.6 nTPM) - nTPM≥100 + cluster + enrichment
3. ORM1 (41,892.9 nTPM) - nTPM≥100 + cluster + enrichment
4. SAA2 (35,358.8 nTPM) - nTPM≥100 + cluster
5. APOA2 (34,742.6 nTPM) - nTPM≥100 + cluster + enrichment

### Recommended Priority:

For most reliable liver protein candidates, use genes in:
1. `liver_protein_nTPM≥100_and_cluster_and_enrichment_genes.csv` (309 genes) - Highest confidence
2. `liver_protein_nTPM≥100_and_cluster_genes.csv` (127 genes) - High confidence
3. Other high nTPM combinations as needed

## File Structure

```
Liver-protein-search/
├── input/
│   ├── *.csv (any gene expression data)
│   └── proteinatlas.tsv
│
├── output/
│   ├── liver_protein_venn_diagram.png
│   ├── liver_ntpm_distribution.png
│   ├── liver_ntpm_group_distribution.png
│   └── liver_protein_roc_curve.png (optional)
│
├── Trace/
│   ├── Data tracing.csv
│   ├── classification_summary.csv
│   ├── liver_genes_sorted_by_ntpm.csv
│   └── [11 classification-specific CSV files]
│
├── Core Pipeline:
│   ├── liver_protein_classifier.py      # Main classification
│   ├── visualize_ntpm.py                 # nTPM distribution plots
│   ├── visualize_group_distribution.py   # Group analysis plots
│   ├── run_visualizations.py             # Visualization runner
│   └── run_pipeline.py                   # Main pipeline runner
│
├── Advanced Analysis:
│   ├── generate_roc_curve.py             # ROC curve analysis
│   └── analyze_classification_strategy.py # Strategy evaluation
│
├── Documentation:
│   ├── README.md
│   └── CLASSIFICATION_STRATEGY_REPORT.md
│
└── Configuration:
    └── .vscode/
```

## Pipeline Commands

### Full Pipeline
```bash
# Complete analysis with all steps
python3 run_pipeline.py
```

### Quick Mode
```bash
# Classification + basic visualization only
python3 run_pipeline.py --quick
```

### Visualization Only
```bash
# Re-run visualizations without re-classifying
python3 run_pipeline.py --viz-only
```

### Individual Components
```bash
# Run classification only
python3 liver_protein_classifier.py

# Run all visualizations
python3 run_visualizations.py

# Run specific visualization
python3 visualize_ntpm.py
python3 visualize_group_distribution.py

# Run advanced analysis
python3 generate_roc_curve.py
python3 analyze_classification_strategy.py
```

## Visualization Outputs

### 1. Venn Diagram (liver_protein_venn_diagram.png)
- Left: nTPM ≥100 overlap with cluster and enrichment
- Right: nTPM <100 overlap with cluster and enrichment
- Color-coded by category

### 2. nTPM Distribution (liver_ntpm_distribution.png)
- Top 50 genes horizontal bar chart
- Distribution scatter plot (log scale)
- Statistics summary box

### 3. Group Distribution (liver_ntpm_group_distribution.png)
Six comprehensive subplots:
- Box plot with outliers
- Violin plot showing density
- Strip plot with individual points
- Histogram with overlap
- Cumulative distribution function (CDF)
- Summary statistics table

## Notes

- The tool automatically extracts liver nTPM values from the "RNA tissue specific nTPM" column
- **New**: nTPM values are split at threshold 100 for refined classification
- Genes are searched for the keyword "liver" in all criteria columns
- All CSV outputs include the extracted nTPM values where available
- Dual Venn diagrams provide clearer visualization of nTPM threshold effects

## Troubleshooting

**No CSV file found**: Ensure your CSV file is in the `input/` folder

**Multiple CSV files**: The tool will use the first file alphabetically and show a warning

**Missing columns**: Verify that proteinatlas.tsv contains all required columns

**Visualization fails**: Ensure classification has been run first (creates Trace/Data tracing.csv)

**Pipeline errors**: Check that all required Python packages are installed

## Performance

- Classification: ~5-10 seconds for 6,500 genes
- Visualization: ~15-20 seconds for all plots
- Full pipeline: ~30-40 seconds (quick mode) / ~60-90 seconds (full mode)

## License

This tool is provided as-is for research purposes.
