#!/usr/bin/env python3
"""
Liver Protein Analysis Pipeline
Main pipeline runner that executes all analysis steps in sequence.

Pipeline Steps:
1. Classification: Classify liver proteins based on nTPM, cluster, and enrichment data
2. Visualization: Generate comprehensive visualization plots
3. (Optional) Advanced Analysis: ROC curves and classification strategy analysis

Usage:
    python run_pipeline.py              # Run full pipeline
    python run_pipeline.py --quick      # Run only classification and basic visualization
    python run_pipeline.py --viz-only   # Run only visualization (requires existing data)
"""

import sys
import subprocess
import time
import argparse
from pathlib import Path

# ANSI color codes for terminal output
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_header(message):
    """Print a formatted header message"""
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*80}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{message:^80}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*80}{Colors.ENDC}\n")

def print_step(step_num, total_steps, message):
    """Print a formatted step message"""
    print(f"\n{Colors.OKCYAN}{Colors.BOLD}[Step {step_num}/{total_steps}]{Colors.ENDC} {message}")

def print_success(message):
    """Print a success message"""
    print(f"{Colors.OKGREEN}✓ {message}{Colors.ENDC}")

def print_error(message):
    """Print an error message"""
    print(f"{Colors.FAIL}✗ {message}{Colors.ENDC}")

def print_warning(message):
    """Print a warning message"""
    print(f"{Colors.WARNING}⚠ {message}{Colors.ENDC}")

def print_info(message):
    """Print an info message"""
    print(f"{Colors.OKBLUE}ℹ {message}{Colors.ENDC}")

def run_script(script_name, description):
    """
    Run a Python script and return success status

    Args:
        script_name (str): Name of the Python script to run
        description (str): Description of what the script does

    Returns:
        bool: True if script ran successfully, False otherwise
    """
    print(f"\n{Colors.OKBLUE}Running: {script_name}{Colors.ENDC}")
    print(f"Description: {description}")

    start_time = time.time()

    try:
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=True,
            text=True,
            check=True
        )

        elapsed_time = time.time() - start_time

        # Print script output
        if result.stdout:
            print(result.stdout)

        print_success(f"Completed in {elapsed_time:.2f} seconds")
        return True

    except subprocess.CalledProcessError as e:
        elapsed_time = time.time() - start_time
        print_error(f"Failed after {elapsed_time:.2f} seconds")

        if e.stdout:
            print(f"\n{Colors.WARNING}STDOUT:{Colors.ENDC}")
            print(e.stdout)

        if e.stderr:
            print(f"\n{Colors.FAIL}STDERR:{Colors.ENDC}")
            print(e.stderr)

        return False

    except Exception as e:
        elapsed_time = time.time() - start_time
        print_error(f"Unexpected error after {elapsed_time:.2f} seconds: {str(e)}")
        return False

def check_input_files():
    """Check if required input files exist"""
    print_header("Checking Input Files")

    required_files = [
        'input/proteinatlas.tsv',
        'input/*.csv',  # At least one CSV file
    ]

    print(f"{Colors.BOLD}Required input files:{Colors.ENDC}")

    # Check TSV file
    tsv_path = Path('input/proteinatlas.tsv')
    if tsv_path.exists():
        print_success(f"Found: {tsv_path}")
    else:
        print_error(f"Missing: {tsv_path}")
        return False

    # Check CSV files
    csv_files = list(Path('input').glob('*.csv'))
    if csv_files:
        print_success(f"Found CSV file(s): {', '.join([f.name for f in csv_files])}")
    else:
        print_error("Missing: No CSV files found in input/ directory")
        return False

    return True

def check_trace_data():
    """Check if trace data exists (for visualization-only mode)"""
    trace_file = Path('Trace/Data tracing.csv')
    return trace_file.exists()

def run_full_pipeline(include_advanced=True):
    """
    Run the full analysis pipeline

    Args:
        include_advanced (bool): Whether to include advanced analysis steps
    """
    pipeline_steps = [
        {
            'name': 'Classification',
            'script': 'core_pipeline/liver_protein_classifier.py',
            'description': 'Classify liver proteins and generate Venn diagram',
            'required': True
        },
        {
            'name': 'Basic Visualization',
            'script': 'core_pipeline/run_visualizations.py',
            'description': 'Generate nTPM distribution and group analysis plots',
            'required': True
        },
    ]

    if include_advanced:
        pipeline_steps.extend([
            {
                'name': 'ROC Curve Analysis',
                'script': 'advanced_analysis/generate_roc_curve.py',
                'description': 'Generate ROC curves for classification evaluation',
                'required': False
            },
            {
                'name': 'Classification Strategy Analysis',
                'script': 'advanced_analysis/analyze_classification_strategy.py',
                'description': 'Analyze different classification strategies',
                'required': False
            },
        ])

    total_steps = len(pipeline_steps)
    successful_steps = 0
    failed_steps = []
    skipped_steps = []

    print_header("Running Analysis Pipeline")

    for idx, step in enumerate(pipeline_steps, 1):
        print_step(idx, total_steps, f"{step['name']}")

        # Check if script exists
        if not Path(step['script']).exists():
            print_warning(f"Script not found: {step['script']}")
            skipped_steps.append(step['name'])
            continue

        if run_script(step['script'], step['description']):
            successful_steps += 1
        else:
            failed_steps.append(step['name'])
            if step['required']:
                print_error(f"Required step failed: {step['name']}")
                print_error("Aborting pipeline")
                return False

    # Print summary
    print_header("Pipeline Summary")

    print(f"{Colors.BOLD}Results:{Colors.ENDC}")
    print(f"  Total steps: {total_steps}")
    print_success(f"Successful: {successful_steps}")

    if skipped_steps:
        print_warning(f"Skipped: {len(skipped_steps)}")
        print(f"\n{Colors.WARNING}Skipped steps:{Colors.ENDC}")
        for step in skipped_steps:
            print(f"  - {step}")

    if failed_steps:
        print_error(f"Failed: {len(failed_steps)}")
        print(f"\n{Colors.FAIL}Failed steps:{Colors.ENDC}")
        for step in failed_steps:
            print(f"  - {step}")
        return False

    print_success("\nAll pipeline steps completed successfully!")

    # List generated outputs
    print(f"\n{Colors.BOLD}Generated outputs:{Colors.ENDC}")

    output_categories = {
        'Venn Diagrams': ['output/liver_protein_venn_diagram.png'],
        'Distribution Plots': [
            'output/liver_ntpm_distribution.png',
            'output/liver_ntpm_group_distribution.png'
        ],
        'Data Files': [
            'Trace/Data tracing.csv',
            'Trace/classification_summary.csv',
            'Trace/liver_genes_sorted_by_ntpm.csv'
        ],
        'ROC Curves': ['output/liver_protein_roc_curve.png'] if include_advanced else [],
        'Analysis Reports': ['CLASSIFICATION_STRATEGY_REPORT.md'] if include_advanced else []
    }

    for category, files in output_categories.items():
        if not files:
            continue
        print(f"\n{Colors.BOLD}{category}:{Colors.ENDC}")
        for file_path in files:
            if Path(file_path).exists():
                print_success(f"  {file_path}")
            else:
                print_warning(f"  {file_path} (not found)")

    return True

def run_visualization_only():
    """Run only visualization scripts (requires existing trace data)"""
    print_header("Running Visualization Only")

    if not check_trace_data():
        print_error("Trace data not found: Trace/Data tracing.csv")
        print_info("Please run the full pipeline or classification step first")
        return False

    print_success("Found trace data")

    return run_script('core_pipeline/run_visualizations.py', 'Generate all visualization plots')

def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(
        description='Liver Protein Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pipeline.py              # Run full pipeline with advanced analysis
  python run_pipeline.py --quick      # Run only classification and basic visualization
  python run_pipeline.py --viz-only   # Run only visualization (requires existing data)
        """
    )

    parser.add_argument(
        '--quick',
        action='store_true',
        help='Run quick mode (skip advanced analysis)'
    )

    parser.add_argument(
        '--viz-only',
        action='store_true',
        help='Run only visualization (requires existing trace data)'
    )

    args = parser.parse_args()

    print_header("Liver Protein Analysis Pipeline")
    print(f"{Colors.BOLD}Mode:{Colors.ENDC} ", end='')

    if args.viz_only:
        print("Visualization Only")
        success = run_visualization_only()
    else:
        if args.quick:
            print("Quick (Classification + Basic Visualization)")
        else:
            print("Full (Classification + Visualization + Advanced Analysis)")

        # Check input files
        if not check_input_files():
            print_error("\nInput file check failed. Please ensure required files exist in input/ directory:")
            print("  - input/proteinatlas.tsv")
            print("  - input/*.csv (gene data file)")
            sys.exit(1)

        print_success("\nAll input files found!")

        # Run pipeline
        success = run_full_pipeline(include_advanced=not args.quick)

    # Exit with appropriate status code
    if success:
        print_header("Pipeline Completed Successfully!")
        sys.exit(0)
    else:
        print_header("Pipeline Failed")
        sys.exit(1)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print_error("\n\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print_error(f"\n\nUnexpected error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
