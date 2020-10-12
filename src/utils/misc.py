import csv
import os
import copy
import itertools
import re
import subprocess
import sys


def sprint(string):
    """
    A modifiable version of print.

    :param string: The string to print.
    """
    print(string)


def delete_if_present(path):
    """
    Delete a file or folder which may or may not exist.

    :param path: The path to delete.
    """
    try:
        os.remove(path)
    except OSError:
        pass


def read(path_file, tag="r"):
    """
    Shorthand for reading in a file.
    """
    with open(path_file, tag) as f:
        return f.read()


def write(msg, path_dest):
    """
    Shorthand for writing a file.
    """
    with open(path_dest, "w") as f:
        f.write(str(msg))


def flatten_lol(lol):
    """
    Flatten a list of lists
    """
    return [y for x in lol for y in x]


def de_dup(seq):
    """
    Get the unique items in a sequence, but preserve the order.
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def wrap_string(string, n):
    """
    Wrap a string with a maximum line length.

    :param string: (str) The string to wrap.
    :param n: (int) The line length.
    :return: (List of str) A list of chunks of the string.
    """
    return [string[start:start + n] for start in range(0, len(string), n)]


def call_function(str_function):
    """
    Call a function in the shell.
    """
    subprocess.call([str_function], shell=True)


def grab_lines(command):
    """
    Run a command in the shell and grab the output as separate lines.

    :param command: (str) The command
    :return: A list of lines as spat out by the shell.
    """
    result = subprocess.run([command], stdout=subprocess.PIPE, shell=True)
    return result.stdout.decode('ascii').split('\n')


def make_if_absent(path_dir):
    """
    Makes a directory if it doesn't exist already.
    """
    if not os.path.exists(path_dir):
        os.makedirs(path_dir)
    return path_dir


def concat_files(files, dest):
    """
    Concatenate a series of files and write the results to the destination.
    """
    with open(dest, "w") as o:
        for file in files:
            with(open(file, "r")) as b:
                o.write(b.read())


def idict(keys, defaultvalue):
    """
    A dict initialised with a default value for each key.
    """
    return {k: copy.deepcopy(defaultvalue) for k in keys}


def read_csv(path_file, ignore_blank=True, ignore_hash=False):
    """
    Read in a CSV file.

    """
    with open(path_file, 'r') as f:
        data = [*csv.reader(f, delimiter='\t')]
        if ignore_blank:
            data = [line for line in data if ''.join(line).strip()]
        if ignore_hash:
            data = [line for line in data if line[0].startswith('#')]
        return data


def write_csv(rows, path_file):
    """
    Write to a CSV file.

    """
    with open(path_file, "w") as f:
        writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE, quotechar='')
        writer.writerows(rows)


def pasync(pool, function, args):
    """
    Run asynchronously.

    """
    pool.apply_async(function, args=args)


def blank_file(path_file):
    """
    Touch a file.
    """
    open(path_file, 'a').close()
    return path_file


def interpolate(l1, l2=[]):
    """
    Fill in the values between the min and max of an integer list, or a pair of integer lists.
    """
    return range(min(l1 + l2), max(l1 + l2))


def check_file_exists(path_file):
    """
    Check whether a file exists. If it does not, exit.
    """
    if os.path.exists(path_file):
        return path_file
    else:
        sys.exit("File does not exist: " + path_file)


def del_indices(l, indices):
    """
    Delete the provided indices from a list.
    """
    to_keep = [i for i in range(len(l)) if i not in indices]
    return [l[i] for i in to_keep]


def anyperm(listin, checker):
    """
    Check whether any permutation of the inputted list is contained in the set of options.

    :param listin: (List) The list to check
    :param checker: (list of lists) A list of possible acceptable permutations.
    :return: (bool) A boolean saying whether or not some permutation of the input list is represented.
    """
    return any([(list(i) in checker) for i in itertools.permutations(listin)])


def sed_file(file_in, file_out, sub_in, sub_out):
    """
    Read in a file, make a substitution, and write out to a different file
    """
    write(re.sub(sub_in, sub_out, read(file_in)), file_out)