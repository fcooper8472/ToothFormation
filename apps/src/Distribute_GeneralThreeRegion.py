import itertools
import os
import re
import subprocess
import time

import multiprocessing as mp
import numpy as np

try:
    import svg_to_webm
except ImportError as e:
    svg_to_webm = None
    quit("Exception: " + str(e))

try:
    import dominate
    from dominate.tags import *
except ImportError as e:
    dominate = None
    quit("Exception: " + str(e))

# Globally accessible directory paths, names, and variables
chaste_build_dir = os.environ.get('CHASTE_BUILD_DIR')
chaste_test_dir = os.environ.get('CHASTE_TEST_OUTPUT')

executable = os.path.join(chaste_build_dir, 'projects/ToothFormation/apps', 'Exe_GeneralThreeRegion')
path_to_output = os.path.join(chaste_test_dir, 'tooth_formation', 'Exe_GeneralThreeRegion')
path_to_sims = os.path.join(path_to_output, 'sim')
path_to_movies = os.path.join(path_to_output, 'movies')

if not(os.path.isfile(executable)):
    quit('Py: Could not find executable: ' + executable)

# List of command line arguments for the executable, and corresponding list of parameter names
command_line_args = [' --ID ', ' --CRL ', ' --CSC ', ' --TRL ', ' --TSC ', ' --AD ', ' --DI ', ' --RM ', ' --TS ']
params_list = ['simulation_id', 'cor_rest_length', 'cor_spring_const', 'tra_rest_length', 'tra_spring_const',
               'adhesion_modifier', 'interaction_dist', 'remesh_frequency', 'num_time_steps']

# Time string when script is run, used for creating a unique archive name
today = time.strftime('%Y-%m-%dT%H%M')

# Param ranges (in lists, for itertools product)
crl = [0.25]
csc = [1e6]
trl = [0.01]
tsc = np.linspace(1e6, 2e6, num=3)
ad = np.linspace(1.4, 1.6, num=3)
di = [0.02]
rf = [2500]
ts = [100]

# An enumerated iterable containing every combination of the parameter ranges defined above
combined_iterable = enumerate(itertools.product(crl, csc, trl, tsc, ad, di, rf, ts))


def main():
    # run_simulations()
    # make_movies_parallel()
    # combine_output()
    generate_html()
    # compress_output()


# Create a list of commands and pass them to separate processes
def run_simulations():

    print("Py: Starting simulations with " + str(mp.cpu_count()) + " processes")

    # Make a list of calls to a Chaste executable
    command_list = []

    if not os.path.exists(path_to_output):
        os.makedirs(path_to_output)

    params_file = open(path_to_output + '/params_file.csv', 'w')
    params_file.write(','.join(params_list) + '\n')

    base_command = 'nice -n 19 ' + executable

    for idx, param_set in combined_iterable:

        params_file.write(str(idx) + ',' + ",".join(map(str, param_set)) + '\n')

        command = base_command
        command += ' --ID ' + str(idx)

        for arg in range(len(param_set)):
            command += command_line_args[arg+1] + str(param_set[arg])

        command_list.append(command)

    params_file.close()

    # Generate a pool of workers
    pool = mp.Pool(processes=mp.cpu_count())

    # Pass the list of bash commands to the pool and wait at most one day
    pool.map_async(execute_command, command_list).get(86400)


# Helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)


# Make a webm from each svg output
def make_movies_parallel():

    print("Py: Combining chaste output to movies with " + str(mp.cpu_count()) + " processes")

    # Validate output directories
    if not (os.path.isdir(path_to_output)):
        quit("Py: Could not find output directory: " + path_to_output)
    if not (os.path.isdir(path_to_sims)):
        quit("Py: Could not find simulation output directory: " + path_to_sims)

    # Create a set of directories containing any svg results files
    svg_pattern = re.compile('results_\d+\.svg')
    data_dirs = set()
    for root, dirs, files in os.walk(path_to_sims):
        for f in files:
            if svg_pattern.match(f) or f == 'svg_arch.tar.gz':
                data_dirs.add(root)

    for d in data_dirs:
        print (d)

    command_list = []

    idx_pattern = re.compile('sim/(\d+)/results_from_time')
    for data_dir in data_dirs:
        index_match = idx_pattern.search(data_dir)
        if not index_match:
            quit('Py: Could not determine simulation index from svg directory string: ' + data_dir)
        idx = int(index_match.group(1))

        command_list.append((data_dir, str(idx).zfill(2) + '.webm', 16.0/9, 15.0, False))

    # Generate a pool of workers
    pool = mp.Pool(processes=mp.cpu_count())

    # Wait at most one day
    pool.map_async(wrap_webm_command, command_list).get(86400)


# Helper function to wrap the movie making command so that it only takes one variable
def wrap_webm_command(args):
    try:
        svg_to_webm.svg_to_webm(*args)
    except Exception as svg_to_webm_exception:
        print("Exception: " + str(svg_to_webm_exception))


# Gather the output from all simulations and put it in the same file
def combine_output():

    print("Py: Combining output")

    if not (os.path.isdir(path_to_output)):
        quit('Py: Could not find output directory: ' + path_to_output)

    combined_results = []
    results_header = []

    for idx, param_set in combined_iterable:
        results_file = os.path.join(path_to_output, 'sim', str(idx), 'results.csv')

        if os.path.isfile(results_file):

            with open(results_file, 'r') as local_results:
                # If the results_header is empty, append the first line of the local results file
                if not results_header:
                    results_header.append(local_results.readline().strip('\n'))
                else:
                    local_results.readline()

                # Store the results (second line of the results file) in the combined_results list
                combined_results.append(local_results.readline().strip('\n'))
        else:  # results file does not exist
            combined_results.append(str(idx) + ',' + 'simulation_incomplete')

    with open(os.path.join(path_to_output, 'combined_results.csv'), 'w') as combined_results_file:
        combined_results_file.write('\n'.join(results_header + combined_results))


def generate_html():
    # Find all webm files

    # Validate output directories
    if not (os.path.isdir(path_to_output)):
        quit("Py: Could not find output directory: " + path_to_output)

    if not os.path.isdir(path_to_movies):
        os.mkdir(path_to_movies)

    # Create a set of directories containing any svg results files
    svg_pattern = re.compile('\d+\.webm')
    webm_files = []
    for root, dirs, files in os.walk(path_to_output):
        for f in files:
            if svg_pattern.match(f):
                webm_files.append(os.path.join(root, f))

    webm_files = sorted(webm_files)

    html_doc = dominate.document(title="Immersed Boundary Simulations")

    table_of_webms = table()

    table_header = thead()
    table_header += td('')
    table_header += td('phead1')
    table_header += td('phead2')
    table_header += td('rhead2')

    table_of_webms += table_header

    table_body = tbody()

    webm_pattern = re.compile('((\d+)\.webm)')
    for webm_file in webm_files:
        webm_match = webm_pattern.search(webm_file)
        if not webm_match:
            quit('Py: Could not find webm file: ' + webm_file)

        webm_name = webm_match.group(1)
        webm_idx = int(webm_match.group(2))

        row_of_table = tr()
        row_of_table += td(a(webm_name, href=webm_file))
        row_of_table += td('param_1')
        row_of_table += td('param_2')
        row_of_table += td('result_1')

        table_body += row_of_table

    table_of_webms += table_body

    html_doc += table_of_webms

    print(html_doc)
    with open(os.path.join(path_to_output, 'index.html'), 'w') as html_index:
        html_index.write(html_doc.render())


# Compress output and suffix with date run
def compress_output():

    print("Py: Compressed output to " + path_to_output + '_' + today + '.tar.gz')

    # Check that output directory exists
    if not (os.path.isdir(path_to_output)):
        raise Exception('Py: Could not find output directory: ' + path_to_output)

    # Change cwd to one above output path
    os.chdir(os.path.join(path_to_output, '..'))

    simulation_name = os.path.basename(os.path.normpath(path_to_output))

    # Compress entire output folder and append with the the simulations were started
    os.system('tar -zcf ' + simulation_name + '_' + today + '.tar.gz ' + simulation_name)


if __name__ == "__main__":
    main()
