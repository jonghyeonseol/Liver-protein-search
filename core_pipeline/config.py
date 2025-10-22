"""
Configuration module for Liver Protein Analysis.
Centralizes all file paths, constants, and configuration values.
"""
from pathlib import Path

# =============================================================================
# DIRECTORY PATHS
# =============================================================================

# Base directories
PROJECT_ROOT = Path(__file__).parent.parent
INPUT_DIR = PROJECT_ROOT / 'input'
OUTPUT_DIR = PROJECT_ROOT / 'output'
TRACE_DIR = PROJECT_ROOT / 'Trace'
LOG_DIR = PROJECT_ROOT / 'logs'

# Ensure directories exist
for directory in [OUTPUT_DIR, TRACE_DIR, LOG_DIR]:
    directory.mkdir(exist_ok=True)

# =============================================================================
# INPUT FILES
# =============================================================================

PROTEIN_ATLAS_FILE = INPUT_DIR / 'proteinatlas.tsv'
CSV_PATTERN = '*.csv'

# =============================================================================
# OUTPUT FILES
# =============================================================================

# Visualization outputs
VENN_DIAGRAM_OUTPUT = OUTPUT_DIR / 'liver_protein_venn_diagram.png'
NTPM_DISTRIBUTION_OUTPUT = OUTPUT_DIR / 'liver_ntpm_distribution.png'
GROUP_DISTRIBUTION_OUTPUT = OUTPUT_DIR / 'liver_ntpm_group_distribution.png'
ROC_CURVE_OUTPUT = OUTPUT_DIR / 'roc_auc_analysis.png'
STRATEGY_COMPARISON_OUTPUT = OUTPUT_DIR / 'classification_strategy_comparison.png'

# Trace/data outputs
TRACE_DATA_FILE = TRACE_DIR / 'Data tracing.csv'
CLASSIFICATION_SUMMARY_FILE = TRACE_DIR / 'classification_summary.csv'
SORTED_GENES_FILE = TRACE_DIR / 'liver_genes_sorted_by_ntpm.csv'
STRATEGY_ANALYSIS_FILE = TRACE_DIR / 'classification_strategy_analysis.csv'

# =============================================================================
# CLASSIFICATION THRESHOLDS
# =============================================================================

# nTPM thresholds
NTPM_HIGH_THRESHOLD = 100
NTPM_VERY_HIGH_THRESHOLD = 200
NTPM_LOW_THRESHOLD = 0

# Confidence scoring thresholds
CONFIDENCE_HIGH = 80
CONFIDENCE_MEDIUM_HIGH = 60
CONFIDENCE_MEDIUM = 40
CONFIDENCE_LOW_MEDIUM = 20

# =============================================================================
# VISUALIZATION SETTINGS
# =============================================================================

# General settings
DPI = 300
TOP_N_GENES = 50

# Figure sizes (width, height)
FIGURE_SIZE_LARGE = (20, 13)
FIGURE_SIZE_MEDIUM = (14, 16)
FIGURE_SIZE_SMALL = (10, 8)
FIGURE_SIZE_VENN = (20, 10)

# =============================================================================
# KNOWN LIVER MARKERS (for validation)
# =============================================================================

KNOWN_LIVER_MARKERS = [
    'ALB', 'HP', 'APOA1', 'APOA2', 'APOC3', 'SERPINA1',
    'FGA', 'FGB', 'FGG', 'CYP3A4', 'CYP2E1', 'FABP1',
    'ALDOB', 'ORM1', 'ORM2', 'GC', 'TF', 'APOH', 'APOE',
    'CRP', 'SAA1', 'SAA2', 'RBP4', 'AMBP', 'SELENOP',
    'HRG', 'F2', 'FTL', 'AHSG', 'C3', 'HPX', 'VTN',
    'ITIH4', 'ITIH3', 'FETUB', 'KNG1', 'SERPINC1',
    'SERPIND1', 'F9', 'F10', 'F11', 'CYP1A2', 'CYP2A6',
    'CYP2C8', 'CYP2C9', 'CYP2D6', 'UGT1A1', 'UGT2B7',
    'GSTA1', 'GSTA2', 'G6PC', 'PCK1'
]

# =============================================================================
# REQUIRED COLUMNS
# =============================================================================

# Columns required in protein atlas data
PROTEIN_ATLAS_REQUIRED_COLUMNS = [
    'Gene',
    'RNA tissue specific nTPM',
    'Tissue expression cluster',
    'RNA tissue cell type enrichment'
]

# Columns required in gene data
GENE_DATA_REQUIRED_COLUMNS = ['Gene']

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_output_path(filename: str) -> Path:
    """Get full path for an output file."""
    return OUTPUT_DIR / filename


def get_trace_path(filename: str) -> Path:
    """Get full path for a trace file."""
    return TRACE_DIR / filename


def get_log_path(filename: str) -> Path:
    """Get full path for a log file."""
    return LOG_DIR / filename
