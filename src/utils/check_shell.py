import subprocess
import sys
import tempfile
from utils.misc import sprint, read, write, call_function


def can_run_command(command, q_allow_stderr=False):
    """
    Return true if the command can be run without error.
    """
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]

    if len(stdout) > 0 and (q_allow_stderr or len(stderr) == 0):
        print(" - ok")
        return True
    else:
        print(" - failed")
        return False


def can_run_specific(line, line_formatted):
    """
    Return true if the command specified in [line] can be run
    """
    print(f"Testing {line_formatted}")
    if can_run_command(line):
        return True
    else:
        print(f"ERROR: Cannot run {line_formatted}.\n    Please check {line_formatted} is installed "
              f"and that the executables are in the system path\n")
    return False


def can_run_man(package, package_formatted):
    """
    Check if the man command can be run for the given package.
    """
    return can_run_specific(f"man {package}", package_formatted)


def can_run_minus_h(package, package_formatted):
    """
    Check if the -h command can be run for the given package.
    """
    return can_run_specific(f"{package} -h", package_formatted)


def can_run_awk():
    """
    Check if awk can be run.
    """
    path_tf1 = tempfile.mktemp()
    path_tf2 = tempfile.mktemp()
    write("a", path_tf1)
    call_function("awk '{print 4 + 5}'" + " " + path_tf1 + " > " + path_tf2)
    res = read(path_tf2)
    print("Testing awk")
    if res == "9\n":
        print(f" - ok\n")
        return True
    else:
        print(f" - failed\nERROR: Cannot run awk.\n    Please check awk is installed "
              f"and that the executables are in the system path\n")
        return False


def check_shell():
    """ Run quick checks to ensure the relevant packages are installed, and that they do
        the things we need them to do (applies in particular to bedtools and to R).
        - awk           - exonerate   - bedtools
        - echo          - sed
        - cut           - grep
        - sort          -tac/cat
    """
    sprint("Checking installed programs...")
    checks = []
    for program in ["sed", "echo", "cut", "sort", "grep", "uniq", "tac", "cat", "mktemp"]:
        checks.append(can_run_man(program, program))
    for program in ["exonerate", "bedtools"]:
        checks.append(can_run_minus_h(program, program))
    checks += [can_run_awk()]
    if not all(checks):
        print("\n_some programs required to run omgene are not installed or are not callable.\n"
              "Please ensure all of the above programs are installed and in the system path.")
        sys.exit()
