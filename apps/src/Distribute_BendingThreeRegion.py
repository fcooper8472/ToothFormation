import itertools
import os
import re
import subprocess
import time

import multiprocessing as mp
import numpy as np

try:
    import png_to_mp4
except ImportError as e:
    png_to_mp4 = None
    quit("Exception: " + str(e))

try:
    import svg_to_png
except ImportError as e:
    svg_to_png = None
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

executable = os.path.join(chaste_build_dir, 'projects/ToothFormation/apps', 'Exe_BendingThreeRegion')
path_to_output = os.path.join(chaste_test_dir, 'tooth_formation', 'Exe_BendingThreeRegion')
path_to_sims = os.path.join(path_to_output, 'sim')
path_to_movies = os.path.join(path_to_output, 'movies')

if not(os.path.isfile(executable)):
    quit('Py: Could not find executable: ' + executable)

# List of command line arguments for the executable, and corresponding list of parameter names
command_line_args = [' --ID ', ' --CRL ', ' --CSC ', ' --TRL ', ' --TSC ', ' --KFS ', ' --ALM ', ' --DI ', ' --SM ',
                     ' --NS ', ' --RM ', ' --TS ', ' --AL ']
params_list = ['simulation_id', 'cor_rest_length', 'cor_spring_const', 'tra_rest_length', 'tra_spring_const',
               'kinematic_strength', 'apical_lam_mult', 'interaction_dist', 'stiffness_mult', 'normal_std',
               'remesh_frequency', 'num_time_steps', 'apical_lamina']

# Time string when script is run, used for creating a unique archive name
today = time.strftime('%Y-%m-%dT%H%M')

# Param ranges (in lists, for itertools product)
crl = [0.25]
csc = np.linspace(1.6 * 1e8, 2.0 * 1e8, num=1)
trl = [1.0]
tsc = np.linspace(0.5 * 1e7, 1.0 * 1e7, num=1)
kfs = np.linspace(1.0 * 1e4, 7.0 * 1e4, num=1)
alm = np.linspace(0.2, 0.2, num=1)
di = [0.02]
sm = np.linspace(0.4, 0.2, num=2)
ns = [0.00]
rf = [50]
ts = [1000]
al = [True]

# An enumerated iterable containing every combination of the parameter ranges defined above
combined_iterable = list(itertools.product(crl, csc, trl, tsc, kfs, alm, di, sm, ns, rf, ts, al))


def main():
    clean_output_dir()
    run_simulations()
    combine_output()
    make_movies_parallel()
    generate_html()
    compress_output()


# Delete output directory
def clean_output_dir():
    if os.path.exists(path_to_output):
        subprocess.call(['rm', '-rf', path_to_output])


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

    for idx, param_set in enumerate(combined_iterable):

        params_file.write(str(idx) + ',' + ",".join(map(str, param_set)) + '\n')

        this_command = base_command
        this_command += ' --ID ' + str(idx)

        for arg in range(len(param_set)):
            this_command += command_line_args[arg+1] + str(param_set[arg])

        command_list.append(this_command)

        print(this_command)

    params_file.close()

    # Generate a pool of workers
    pool = mp.Pool(processes=mp.cpu_count())

    # Pass the list of bash commands to the pool and wait at most one day
    pool.map_async(execute_command, command_list).get(86400)


# Helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)


# Make an mp4 from each svg output
def make_movies_parallel():

    print("Py: Combining chaste output to movies with " + str(mp.cpu_count()) + " processes")

    # Validate output directories
    if not (os.path.isdir(path_to_output)):
        quit("Py: Could not find output directory: " + path_to_output)
    if not (os.path.isdir(path_to_sims)):
        quit("Py: Could not find simulation output directory: " + path_to_sims)

    # Convert all svg to png using cairosvg
    svg_to_png.svg_to_png(path_to_output)

    # Create a set of directories containing any png results files
    png_pattern = re.compile('results_\d+\.png')
    data_dirs = set()
    for root, dirs, files in os.walk(path_to_sims):
        for f in files:
            if png_pattern.match(f) or f == 'svg_arch.tar.gz':
                data_dirs.add(root)

    command_list = []

    idx_pattern = re.compile('sim/(\d+)/results_from_time')
    for data_dir in data_dirs:
        index_match = idx_pattern.search(data_dir)
        if not index_match:
            quit('Py: Could not determine simulation index from directory string: ' + data_dir)
        idx = int(index_match.group(1))

        command_list.append((data_dir, str(idx).zfill(2) + '.mp4', 10.0, False))

    # Generate a pool of workers
    pool = mp.Pool(processes=mp.cpu_count())

    # Wait at most one day
    pool.map_async(wrap_mp4_command, command_list).get(86400)


# Helper function to wrap the movie making command so that it only takes one variable
def wrap_mp4_command(args):
    try:
        png_to_mp4.png_to_mp4(*args)
    except Exception as png_to_mp4_exception:
        print("Exception: " + str(png_to_mp4_exception))


# Gather the output from all simulations and put it in the same file
def combine_output():

    print("Py: Combining output")

    if not (os.path.isdir(path_to_output)):
        quit('Py: Could not find output directory: ' + path_to_output)

    combined_results = []
    results_header = []

    for idx, param_set in enumerate(combined_iterable):
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

    if os.path.getsize(os.path.join(path_to_output, 'combined_results.csv')) < 4:
        quit("Py: Combined results not generated as expected")


def generate_html():

    print("Py: Generating html")

    # Validate output directories
    if not (os.path.isdir(path_to_output)):
        quit("Py: Could not find output directory: " + path_to_output)

    if not os.path.isdir(path_to_movies):
        os.mkdir(path_to_movies)

    copy_assets()

    # Copy all mp4 files into the movies directory
    mp4_numbers = []
    mp4_pattern = re.compile('(\d+)\.mp4')
    for root, dirs, files in os.walk(path_to_sims):
        for f in files:
            mp4_match = mp4_pattern.match(f)
            if mp4_match:
                subprocess.call(['cp', os.path.join(root, f), os.path.join(path_to_movies, f), '-u'])
                mp4_numbers.append(int(mp4_match.group(1)))

    generate_movies_html(os.listdir(path_to_movies))

    # Read params and results into lists of lists
    with open(os.path.join(path_to_output, 'params_file.csv'), 'r') as params_file:
        combined_parameters = params_file.readlines()

    with open(os.path.join(path_to_output, 'combined_results.csv'), 'r') as combined_results_file:
        combined_results = combined_results_file.readlines()

    # Concatenate and format the headers
    header_row = combined_parameters[0].strip().split(',') + combined_results[0].strip().split(',')[1:]
    header_row = [header_row[x].replace('_', ' ').title() for x in range(len(header_row))]

    params_and_results = []
    for row_num in range(1, len(combined_parameters)):
        params_and_results.append(combined_parameters[row_num].strip().split(',') +
                                  combined_results[row_num].strip().split(',')[1:])

    # Determine which data columns are unchanging
    unchanging_cols = []
    for col_idx in range(len(params_and_results[0])):
        column = [row[col_idx] for row in params_and_results]
        first_value = column[0]
        if all(x == first_value for x in column):
            unchanging_cols.append(col_idx)

    # Create html document
    html_doc = dominate.document(title="Immersed Boundary Simulations")

    # Add html head containing css and js
    with html_doc.head:
        link(rel='stylesheet', href='css/theme.green.min.css')
        link(rel='stylesheet', href='css/default.css')
        script(src='js/jquery.min.js')
        script(src='js/jquery.tablesorter.min.js')
        script(src='js/sort_table.js')

    with html_doc:
        with div(_class='main').add(table(id='mp4_table', _class='tablesorter')):
            with thead():
                td('mp4 Name')
                for idx, cell in enumerate(header_row):
                    if idx in unchanging_cols:
                        td(cell, _class='cell_header cell_unchanging')
                    else:
                        td(cell, _class='cell_header')
            with tbody():
                for row in params_and_results:
                    idx = int(row[0])
                    with tr():
                        if idx in mp4_numbers:
                            mp4 = str(idx).zfill(2)
                            td(a(mp4 + '.mp4', href=os.path.join('html', mp4 + '.html'), target='_blank'))
                        else:
                            td('no mp4')
                        for idx, cell in enumerate(row):
                            if idx in unchanging_cols:
                                td(cell, _class='cell_right cell_unchanging')
                            elif idx == 0:
                                td(cell, _class='cell_centre')  # id column centred
                            else:
                                td(cell, _class='cell_right')   # other columns right-aligned

    with open(os.path.join(path_to_output, 'index.html'), 'w') as html_index:
        html_index.write(html_doc.render())


def copy_assets():

    scripts_dir = os.environ.get('SCRIPTS_DIR')
    if not scripts_dir:
        print("Py: Warning: Could not find scripts directory: css and js elements may not work")
    else:
        # create directories to copy in to
        if not os.path.isdir(os.path.join(path_to_output, 'css')):
            os.mkdir(os.path.join(path_to_output, 'css'))
        if not os.path.isdir(os.path.join(path_to_output, 'js')):
            os.mkdir(os.path.join(path_to_output, 'js'))

        # default css
        default_css_src = os.path.join(scripts_dir, 'css', 'default.css')
        default_css_dst = os.path.join(path_to_output, 'css', 'default.css')
        subprocess.call(['cp', default_css_src, default_css_dst, '-u'])

        # tablesorter theme css
        green_css_src = os.path.join(scripts_dir, 'css', 'theme.green.min.css')
        green_css_dst = os.path.join(path_to_output, 'css', 'theme.green.min.css')
        subprocess.call(['cp', green_css_src, green_css_dst, '-u'])

        # jQuery js
        jquery_js_src = os.path.join(scripts_dir, 'js', 'jquery.min.js')
        jquery_js_dst = os.path.join(path_to_output, 'js', 'jquery.min.js')
        subprocess.call(['cp', jquery_js_src, jquery_js_dst, '-u'])

        # tablesorter js
        tablesorter_js_src = os.path.join(scripts_dir, 'js', 'jquery.tablesorter.min.js')
        tablesorter_js_dst = os.path.join(path_to_output, 'js', 'jquery.tablesorter.min.js')
        subprocess.call(['cp', tablesorter_js_src, tablesorter_js_dst, '-u'])

        # sort_table js
        sort_table_js_src = os.path.join(scripts_dir, 'js', 'sort_table.js')
        sort_table_js_dst = os.path.join(path_to_output, 'js', 'sort_table.js')
        subprocess.call(['cp', sort_table_js_src, sort_table_js_dst, '-u'])


def generate_movies_html(list_of_movie_names):

    path_to_html = os.path.join(path_to_output, 'html')
    if not os.path.isdir(path_to_html):
        os.mkdir(path_to_html)

    for movie_name in list_of_movie_names:

        html_doc = dominate.document(title=movie_name)

        # Add html head containing css
        with html_doc.head:
            link(rel='stylesheet', href='../css/default.css')

        with html_doc:
            with div(_class='video'):
                video(controls='', src='../movies/'+movie_name, type='video/mp4', _class='full')

        with open(os.path.join(path_to_html, movie_name.replace('.mp4', '.html')), 'w') as mp4_html:
            mp4_html.write(html_doc.render())


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
