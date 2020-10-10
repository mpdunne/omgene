import subprocess


##########################################
# Check programs
##########################################

def can_run_command(command, q_allow_stderr=False):
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]
    return check(len(stdout) > 0 and (q_allow_stderr or len(stderr) == 0))


def can_run_specific(line, line_formatted):
    return True if can_run_command(line) else output_bool(False, "ERROR: Cannot run " + line_formatted,
                                                          "    Please check " + line_formatted + " is installed and that the executables are in the system path\n")


def can_run_minus_h(package, package_formatted):
    return can_run_specific(package + " -h", package_formatted)


def can_run_awk():
    return returns_expected("awk '{print 4 + 5}'", "awk '{print 4 + 5}'", "a", "9\n")


def can_run_man(package, package_formatted):
    return can_run_specific("man " + package, package_formatted)


def check(bool_val, true_msg=" - ok", false_msg=" - failed", true_extra="", false_extra=""):
    return output_bool(True, true_msg, true_extra) if bool_val else output_bool(False, false_msg, false_extra)


def returns_expected(message, command, contents_in, expected_out):
    # sys.stdout.write("Test can run \""+message+"\"\n")
    path_tf1 = tempfile.mktemp()
    path_tf2 = tempfile.mktemp()
    write(contents_in, path_tf1)
    call_function(command + " " + path_tf1 + " > " + path_tf2)
    res = read(path_tf2)
    return res == expected_out


def output_bool(bool_val, msg, msg_extra):
    if msg_extra:
        print(msg_extra)
    return bool_val


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
        print(
            "\n_some programs required to run omgene are not installed or are not callable.\n_please ensure all of the above programs are installed and in the system path.")
        sys.exit()
