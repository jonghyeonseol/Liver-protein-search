"""
Visualization configuration for liver protein analysis.
Centralized color mappings and display settings for all visualization modules.
"""

# =============================================================================
# CLASSIFICATION COLOR MAPPING
# =============================================================================

CLASSIFICATION_COLOR_MAP = {
    # High nTPM (>=100) categories - Red spectrum
    'liver protein (nTPM≥100 + cluster + enrichment)': '#FF0000',  # Bright red
    'liver protein (nTPM≥100 + cluster)': '#FF6B6B',  # Light red
    'liver protein (nTPM≥100 + enrichment)': '#FF8C42',  # Orange-red
    'liver protein (nTPM≥100 only)': '#FFB6B6',  # Very light red

    # Low nTPM (<100) categories - Blue spectrum
    'liver protein (nTPM<100 + cluster + enrichment)': '#4169E1',  # Royal blue
    'liver protein (nTPM<100 + cluster)': '#87CEEB',  # Sky blue
    'liver protein (nTPM<100 + enrichment)': '#ADD8E6',  # Light blue
    'liver protein (nTPM<100 only)': '#E0F0FF',  # Very light blue

    # Both nTPM levels - Purple spectrum
    'liver protein (both nTPM levels)': '#9370DB',  # Medium purple
    'liver protein (both nTPM + cluster)': '#8B008B',  # Dark magenta
    'liver protein (both nTPM + enrichment)': '#BA55D3',  # Medium orchid
    'liver protein (all 4)': '#4B0082',  # Indigo

    # No nTPM categories - Green spectrum
    'liver protein (cluster + enrichment)': '#A8E6CF',  # Mint green
    'liver protein (cluster only)': '#95E1D3',  # Light teal
    'liver protein (enrichment only)': '#C7CEEA',  # Light periwinkle

    # Non-liver
    'non-liver protein': '#E8E8E8'  # Light gray
}

# =============================================================================
# CLASSIFICATION DISPLAY ORDER
# =============================================================================

# Order for displaying classifications in plots and legends
CLASSIFICATION_ORDER = [
    # High nTPM first
    'liver protein (nTPM≥100 + cluster + enrichment)',
    'liver protein (nTPM≥100 + cluster)',
    'liver protein (nTPM≥100 + enrichment)',
    'liver protein (nTPM≥100 only)',
    # Low nTPM
    'liver protein (nTPM<100 + cluster + enrichment)',
    'liver protein (nTPM<100 + cluster)',
    'liver protein (nTPM<100 + enrichment)',
    'liver protein (nTPM<100 only)',
    # No nTPM
    'liver protein (cluster + enrichment)',
    'liver protein (cluster only)',
    'liver protein (enrichment only)'
]

# =============================================================================
# LABEL FORMATTING
# =============================================================================

def format_classification_label(classification: str, short: bool = True) -> str:
    """
    Format classification labels for display in plots.

    Args:
        classification: Full classification name
        short: If True, use abbreviated format (default: True)

    Returns:
        Formatted label string
    """
    if not short:
        return classification

    # Create abbreviated labels
    label = classification.replace('liver protein ', '').replace('(', '').replace(')', '')
    label = label.replace('nTPM≥100', '≥100').replace('nTPM<100', '<100')
    label = label.replace('cluster', 'C').replace('enrichment', 'E')
    label = label.replace(' + ', '+').replace(' and ', '&')

    return label


# =============================================================================
# VENN DIAGRAM COLORS
# =============================================================================

VENN_HIGH_NTPM_COLOR = '#FF6B6B'  # Light red for nTPM >= 100
VENN_LOW_NTPM_COLOR = '#FFB6C1'   # Pink for nTPM < 100
VENN_CLUSTER_COLOR = '#87CEEB'     # Sky blue for cluster
VENN_ENRICHMENT_COLOR = '#98FB98'  # Pale green for enrichment

# =============================================================================
# PLOT AESTHETICS
# =============================================================================

# Default alpha transparency
ALPHA_DEFAULT = 0.65
ALPHA_HIGHLIGHT = 0.8

# Line widths
LINE_WIDTH_THIN = 0.5
LINE_WIDTH_MEDIUM = 1.5
LINE_WIDTH_THICK = 2.0
LINE_WIDTH_VERY_THICK = 2.5

# Marker sizes
MARKER_SIZE_SMALL = 3
MARKER_SIZE_MEDIUM = 30
MARKER_SIZE_LARGE = 50

# Font sizes
FONT_SIZE_SMALL = 8
FONT_SIZE_MEDIUM = 10
FONT_SIZE_LARGE = 12
FONT_SIZE_XLARGE = 14
FONT_SIZE_TITLE = 18

# =============================================================================
# LEGEND CONFIGURATION
# =============================================================================

LEGEND_LOCATION = 'best'
LEGEND_FRAMEALPHA = 0.95
LEGEND_EDGECOLOR = 'black'
LEGEND_FANCYBOX = True
LEGEND_SHADOW = True
