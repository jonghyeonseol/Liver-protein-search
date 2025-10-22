# Liver Protein Classification Tool

This tool classifies proteins based on liver-specific expression data from the Human Protein Atlas.

## Features

- **3-Criteria Classification**: Classifies proteins based on three liver-related criteria
  - RNA tissue specific nTPM
  - Tissue expression cluster
  - RNA tissue cell type enrichment

- **Automatic File Detection**: Automatically detects CSV files in the input folder
- **Comprehensive Visualization**: Generates Venn diagrams and distribution plots
- **Data Tracing**: Complete tracking of all classification decisions

## Requirements

```bash
pip install pandas matplotlib matplotlib-venn
```

## Input Files

Place the following files in the `input/` folder:

1. **Any CSV file** containing gene expression data with a "Gene" column
   - Example: `log2_center_normalization_t_test_Human.csv`
   - If multiple CSV files exist, the first one will be used

2. **proteinatlas.tsv** - Human Protein Atlas data
   - Must contain columns: Gene, RNA tissue specific nTPM, Tissue expression cluster, RNA tissue cell type enrichment

## Usage

### 1. Run Classification

```bash
python3 liver_protein_classifier.py
```

This will:
- Automatically detect CSV file in input folder
- Classify proteins based on 3 criteria
- Generate Venn diagram and save all results

### 2. Generate Visualizations

```bash
python3 visualize_ntpm.py
```

This will:
- Create distribution plots
- Generate top 50 genes bar chart
- Save sorted gene lists

## Output Files

### output/ folder (Visualizations)
- `liver_protein_venn_diagram.png` - 3-way Venn diagram showing overlap between criteria
- `liver_ntpm_distribution.png` - nTPM value distribution and top 50 genes

### Trace/ folder (Data)
- `Data tracing.csv` - Complete dataset with all classifications
- `classification_summary.csv` - Summary statistics
- `liver_protein_all_3_genes.csv` - 391 genes meeting all 3 criteria (highest confidence)
- `liver_protein_nTPM_and_cluster_genes.csv` - 184 genes
- `liver_protein_nTPM_and_enrichment_genes.csv` - 71 genes
- `liver_protein_cluster_and_enrichment_genes.csv` - 57 genes
- `liver_protein_1_genes.csv` - 61 genes (nTPM only)
- `liver_protein_2_genes.csv` - 157 genes (cluster only)
- `liver_protein_3_genes.csv` - 672 genes (enrichment only)
- `liver_genes_sorted_by_ntpm.csv` - All genes sorted by nTPM value

## Classification Logic

### 7 Categories:

1. **liver protein (all 3)**: Meets all 3 criteria - HIGHEST CONFIDENCE
2. **liver protein (nTPM + cluster)**: nTPM and tissue cluster
3. **liver protein (nTPM + enrichment)**: nTPM and cell type enrichment
4. **liver protein (cluster + enrichment)**: Tissue cluster and cell type enrichment
5. **liver protein_1**: RNA tissue specific nTPM only
6. **liver protein_2**: Tissue expression cluster only
7. **liver protein_3**: RNA tissue cell type enrichment only
8. **non-liver protein**: Does not meet any criteria

### Statistics (Example Dataset):

- **Total proteins**: 6,517
- **Liver candidates**: 1,592 (24.4%)
- **Non-liver proteins**: 4,925 (75.6%)

## Key Results

### Top Liver Proteins (by nTPM):
1. ALB (205,952.1 nTPM) - All 3 criteria
2. HP (54,575.6 nTPM) - All 3 criteria
3. ORM1 (41,892.9 nTPM) - All 3 criteria
4. SAA2 (35,358.8 nTPM) - nTPM + cluster
5. APOA2 (34,742.6 nTPM) - All 3 criteria

### Recommended Priority:

For most reliable liver protein candidates, use genes in:
1. `liver_protein_all_3_genes.csv` (391 genes) - Highest confidence
2. `liver_protein_nTPM_and_cluster_genes.csv` (184 genes) - High confidence
3. Other combinations as needed

## File Structure

```
Liver-protein-search/
├── input/
│   ├── *.csv (any gene expression data)
│   └── proteinatlas.tsv
├── output/
│   ├── liver_protein_venn_diagram.png
│   └── liver_ntpm_distribution.png
├── Trace/
│   ├── Data tracing.csv
│   ├── classification_summary.csv
│   └── [7 classification CSV files]
├── liver_protein_classifier.py
├── visualize_ntpm.py
└── README.md
```

## Notes

- The tool automatically extracts liver nTPM values from the "RNA tissue specific nTPM" column
- Genes are searched for the keyword "liver" in all three criteria columns
- All CSV outputs include the extracted nTPM values where available
- Venn diagram shows all 7 possible combinations of the 3 criteria

## Troubleshooting

**No CSV file found**: Ensure your CSV file is in the `input/` folder

**Multiple CSV files**: The tool will use the first file alphabetically and show a warning

**Missing columns**: Verify that proteinatlas.tsv contains all required columns

## License

This tool is provided as-is for research purposes.
