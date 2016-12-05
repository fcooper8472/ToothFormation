import itertools
import multiprocessing
import os
import subprocess
import time

import numpy as np
import matplotlib.pyplot as plt

try:
    import immersed_boundary
except ImportError:
    print("Py: immersed_boundary not imported - is chaste project in the python path?")

# Globally accessible directory paths, names, and variables
chaste_build_dir = os.environ.get('CHASTE_BUILD_DIR')
executable = os.path.join(chaste_build_dir, 'projects/ToothFormation/apps', 'Exe_GeneralThreeRegion')

if not(os.path.isfile(executable)):
    raise Exception('Py: Could not find executable: ' + executable)

chaste_test_dir = os.environ.get('CHASTE_TEST_OUTPUT')
path_to_output = os.path.join(chaste_test_dir, 'tooth_formation', 'Exe_GeneralThreeRegion')

command_line_args = [' --ID ', ' --CRL ', ' --CSC ', ' --TRL ', ' --TSC ', ' --AD ', ' --DI ', ' --RM ', ' --TS ']
params_list = ['simulation_id', 'cor_rest_length', 'cor_spring_const', 'tra_rest_length', 'tra_spring_const',
               'adhesion_modifier', 'interaction_dist', 'remesh_frequency', 'num_time_steps']

today = time.strftime('%Y-%m-%dT%H%M')

# Param ranges (in lists, for itertools product
crl = [0.25]
csc = [1e6]
trl = [0.01]
tsc = np.linspace(1e6, 2e6, num=3)
ad = np.linspace(1.0, 2.0, num=3)
di = np.linspace(0.01, 0.03, num=3)
rf = [500]
ts = [50000]

combined_iterable = enumerate(itertools.product(crl, csc, trl, tsc, ad, di, rf, ts))


def main():
    run_simulations()
    make_movies_parallel()
    combine_output()
    # plot_results()
    compress_output()


# Create a list of commands and pass them to separate processes
def run_simulations():

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

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Py: Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool

    # Wait at most one day
    pool.map_async(execute_command, command_list).get(86400)


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)


# Make an mp4 movie from each pvd file
def make_movies_parallel():

    print("Py: Combining chaste output to movies")

    if not (os.path.isdir(path_to_output)):
        raise Exception('Py: Could not find output directory: ' + path_to_output)

    path_to_movies = os.path.join(path_to_output, 'movies')
    if not (os.path.isdir(path_to_movies)):
        os.makedirs(path_to_movies)

    command_list = []

    for idx, param_set in combined_iterable:

        # Generate the simulation directory
        sim_dir = os.path.join(path_to_output, 'sim', str(idx))
        results_file = os.path.join(sim_dir, 'results.csv')

        # We identify whether there is valid data; if not, the paraview script will not exit gracefully
        valid_data_available = False

        possible_data_directories = []
        for directory in os.listdir(sim_dir):
            if directory.startswith('results_from_time'):
                possible_data_directories.append(os.path.join(sim_dir, directory))

        if len(possible_data_directories) > 0:
            # The last directory alphabetically will be the relevant one for visualisation
            data_directory = sorted(possible_data_directories)[-1]

            # Get the path to the pvd file
            pvd_file = os.path.join(data_directory, 'results.pvd')

            if os.path.isfile(results_file) and os.path.isfile(pvd_file) and os.path.getsize(pvd_file) > 1024:
                valid_data_available = True

        if valid_data_available:
            command_list.append((sim_dir, path_to_movies, str(idx), 'Points', 9))
        else:
            print("Py: No valid data for sim with index " + str(idx))

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Py: Creating movies from simulation output with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Wait at most one day
    pool.map_async(wrap_movie_command, command_list).get(86400)


# Helper function to wrap the movie making command so that it only takes one variable
def wrap_movie_command(args):
    immersed_boundary.pvd_to_mp4(*args)


# Gather the output from all simulations and put it in the same file
def combine_output():

    print("Py: Combining output")

    if not (os.path.isdir(path_to_output)):
        raise Exception('Py: Could not find output directory: ' + path_to_output)

    combined_results = open(os.path.join(path_to_output, 'combined_results.csv'), 'w')
    added_header = False

    for idx, param_set in combined_iterable:
        file_name = os.path.join(path_to_output, 'sim', str(idx), 'results.csv')

        if os.path.isfile(file_name):
            local_results = open(file_name, 'r')

            # Add header to combined results if not yet done - else skip the header
            if not added_header:
                combined_results.write(local_results.readline())
                added_header = True
            else:
                _ = local_results.readline()

            # Write the results to the combined results file
            combined_results.write(local_results.readline())

            local_results.close()
            combined_results.write('\n')

        else:
            combined_results.write(str(idx) + ',' + 'simulation_incomplete' + '\n')

    combined_results.close()


def plot_results():

    print("Py: Plotting results")
    #
    # bg_gray = '0.75'
    # bg_line_width = 0.5
    #
    # # Set LaTeX font rendering
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif', serif='computer modern roman')
    # plt.rc('figure', figsize=[10,6.18]) # inches
    # plt.rc('axes', linewidth=bg_line_width, edgecolor=bg_gray, axisbelow=True)
    # plt.rc('xtick', labelsize=10)
    # plt.rc('ytick', labelsize=10)
    # plt.rc('xtick.major', size=0, pad=4)
    # plt.rc('ytick.major', size=0, pad=4)
    #
    #
    # my_data = np.genfromtxt(path_to_output + exec_name + '/combined_results.dat', delimiter=',', skip_header=1)
    #
    # col_to_plot = 4 # the column to plot from the data file
    #
    # x_vals = my_data[:,1] # local const
    # y_vals = my_data[:,col_to_plot] # summary statistic
    #
    # x_min, x_max = min(x_vals), max(x_vals)
    # y_min, y_max = min(y_vals), max(y_vals)
    #
    # x_range, y_range = x_max - x_min, y_max - y_min
    # margin = 0.05;
    #
    # x_lims = [x_min - margin * x_range, x_max + margin * x_range]
    # y_lims = [y_min - margin * y_range, y_max + margin * y_range]
    #
    # num_sims_per_plot = num_local_consts * num_kicks_per_sim
    #
    # for plot in range(num_global_consts):
    #
    #     plt.clf()
    #
    #     plt.xlabel(r'Local spring constant multiple', fontsize=12, labelpad=20)
    #     plt.ylabel(r'Right asymmetry', fontsize=12, labelpad=20)
    #
    #     title = 'Global cell-cell spring constant scaled by ' + str(0.5 + 0.25 * plot)
    #
    #     plt.title(title, fontsize=14, y = 1.05)
    #
    #     highlight_x_vals = my_data[plot * num_sims_per_plot : (plot+1) * num_sims_per_plot,1]
    #     highlight_y_vals = my_data[plot * num_sims_per_plot : (plot+1) * num_sims_per_plot,col_to_plot]
    #
    #     plt.scatter(x_vals, y_vals, color = '0.5', marker = 'o', facecolors='none', edgecolors='0.5')
    #     plt.scatter(highlight_x_vals, highlight_y_vals, color = '#F39200')
    #
    #     plt.xlim(x_lims)
    #     plt.ylim(y_lims)
    #
    #
    #     ax = plt.gca()
    #     ax.grid(b=True, which='major', color=bg_gray, linestyle='dotted', dash_capstyle='round')
    #
    #     plt.savefig(path_to_output + exec_name + '/Fig_' + exec_name + str(plot) + '_pdf.pdf', bbox_inches='tight', pad_inches=0.4)
    #     plt.savefig(path_to_output + exec_name + '/Fig_' + exec_name + str(plot) + '_eps.eps', bbox_inches='tight', pad_inches=0.5)


# Compress output and suffix with date run
def compress_output():

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