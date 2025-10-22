"""
Console output utilities for formatted terminal messages.
Provides consistent formatting for pipeline outputs.
"""


class Colors:
    """ANSI color codes for terminal output."""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def print_header(message: str, width: int = 80) -> None:
    """
    Print a formatted header message.

    Args:
        message: Header text to display
        width: Total width of the header (default: 80)
    """
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*width}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{message:^{width}}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*width}{Colors.ENDC}\n")


def print_step(step_num: int, total_steps: int, message: str) -> None:
    """
    Print a formatted step message.

    Args:
        step_num: Current step number
        total_steps: Total number of steps
        message: Step description
    """
    print(f"\n{Colors.OKCYAN}{Colors.BOLD}[Step {step_num}/{total_steps}]{Colors.ENDC} {message}")


def print_success(message: str) -> None:
    """
    Print a success message.

    Args:
        message: Success message to display
    """
    print(f"{Colors.OKGREEN}✓ {message}{Colors.ENDC}")


def print_error(message: str) -> None:
    """
    Print an error message.

    Args:
        message: Error message to display
    """
    print(f"{Colors.FAIL}✗ {message}{Colors.ENDC}")


def print_warning(message: str) -> None:
    """
    Print a warning message.

    Args:
        message: Warning message to display
    """
    print(f"{Colors.WARNING}⚠ {message}{Colors.ENDC}")


def print_info(message: str) -> None:
    """
    Print an info message.

    Args:
        message: Info message to display
    """
    print(f"{Colors.OKBLUE}ℹ {message}{Colors.ENDC}")


def print_section(title: str) -> None:
    """
    Print a section title.

    Args:
        title: Section title to display
    """
    print(f"\n{Colors.BOLD}{title}{Colors.ENDC}")
    print(f"{Colors.BOLD}{'-' * len(title)}{Colors.ENDC}")
