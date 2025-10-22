#!/usr/bin/env python3
"""
Visualization Runner Module
This module runs all visualization scripts in sequence to generate comprehensive analysis plots.
"""

# Standard library imports
import subprocess
import sys
import time
from pathlib import Path

# Local imports
from console_utils import (
    Colors, print_header, print_step, print_success,
    print_error, print_warning
)

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

def check_prerequisites():
    """Check if required files and directories exist"""
    print_header("Checking Prerequisites")

    required_files = [
        'Trace/Data tracing.csv',
    ]

    required_scripts = [
        'core_pipeline/visualize_ntpm.py',
        'core_pipeline/visualize_group_distribution.py',
    ]

    all_ok = True

    # Check required data files
    print(f"{Colors.BOLD}Required data files:{Colors.ENDC}")
    for file_path in required_files:
        if Path(file_path).exists():
            print_success(f"Found: {file_path}")
        else:
            print_error(f"Missing: {file_path}")
            all_ok = False

    # Check required scripts
    print(f"\n{Colors.BOLD}Required scripts:{Colors.ENDC}")
    for script in required_scripts:
        if Path(script).exists():
            print_success(f"Found: {script}")
        else:
            print_error(f"Missing: {script}")
            all_ok = False

    return all_ok

def main():
    """Main execution function"""
    print_header("Liver Protein Visualization Pipeline")

    print(f"{Colors.BOLD}This pipeline will generate:{Colors.ENDC}")
    print("  1. liver_ntpm_distribution.png - Top genes and distribution analysis")
    print("  2. liver_ntpm_group_distribution.png - Comprehensive group statistics")
    print("  3. liver_genes_sorted_by_ntpm.csv - Sorted gene list")

    # Check prerequisites
    if not check_prerequisites():
        print_error("\nPrerequisite check failed. Please ensure all required files exist.")
        print_warning("Hint: Run 'liver_protein_classifier.py' first to generate required data.")
        sys.exit(1)

    print_success("\nAll prerequisites satisfied!")

    # Define visualization scripts to run
    visualization_steps = [
        {
            'script': 'core_pipeline/visualize_ntpm.py',
            'description': 'Generate nTPM distribution plots (top genes and scatter plot)',
        },
        {
            'script': 'core_pipeline/visualize_group_distribution.py',
            'description': 'Generate comprehensive group distribution analysis (6 subplots)',
        },
    ]

    total_steps = len(visualization_steps)
    successful_steps = 0
    failed_steps = []

    # Run each visualization script
    print_header("Running Visualization Scripts")

    for idx, step in enumerate(visualization_steps, 1):
        print_step(idx, total_steps, step['description'])

        if run_script(step['script'], step['description']):
            successful_steps += 1
        else:
            failed_steps.append(step['script'])

    # Print summary
    print_header("Visualization Pipeline Summary")

    print(f"{Colors.BOLD}Results:{Colors.ENDC}")
    print(f"  Total steps: {total_steps}")
    print_success(f"Successful: {successful_steps}")

    if failed_steps:
        print_error(f"Failed: {len(failed_steps)}")
        print(f"\n{Colors.FAIL}Failed scripts:{Colors.ENDC}")
        for script in failed_steps:
            print(f"  - {script}")
        sys.exit(1)
    else:
        print_success("\nAll visualizations completed successfully!")

        print(f"\n{Colors.BOLD}Generated outputs:{Colors.ENDC}")
        output_files = [
            'output/liver_ntpm_distribution.png',
            'output/liver_ntpm_group_distribution.png',
            'Trace/liver_genes_sorted_by_ntpm.csv',
        ]

        for output_file in output_files:
            if Path(output_file).exists():
                print_success(f"{output_file}")
            else:
                print_warning(f"{output_file} (not found)")

        sys.exit(0)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print_error("\n\nVisualization pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print_error(f"\n\nUnexpected error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
