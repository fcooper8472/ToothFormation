# Example CMake build configuration file for a Chaste user project.
# See also https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/CmakeBuildGuide for full documentation.

# You can set which Chaste components (or other projects) your project depends on by editing the
# find_package() call for Chaste.  E.g.
#   find_package(Chaste COMPONENTS cell_based)
# for a project just using the cell_based component (and its dependencies), or
#   find_package(Chaste COMPONENTS heart lung)
# for a project combining heart & lung simulations.
#
# Note that the order in which components are specified does not matter.

# Find eigen3
set(EIGEN3_ROOT "$ENV{EIGEN3_ROOT}")
if(EXISTS ${EIGEN3_ROOT}/signature_of_eigen3_matrix_library)
    include_directories(SYSTEM "${EIGEN3_ROOT}")
else()
    message(FATAL_ERROR "Set the environment variable EIGEN3_ROOT to the Eigen3 directory containing the file 'signature_of_eigen3_matrix_library'.")
endif()

# Find Spectra (https://github.com/yixuan/spectra)
set(SPECTRA_ROOT "$ENV{SPECTRA_ROOT}/include")
if(EXISTS ${SPECTRA_ROOT}/GenEigsComplexShiftSolver.h)
    include_directories(SYSTEM "${SPECTRA_ROOT}")
else()
    message(FATAL_ERROR "Set the environment variable SPECTRA_ROOT to the Spectra directory containing the file 'include/GenEigsComplexShiftSolver.h'.")
endif()

# Here we just depend on core components (nothing application-specific).
find_package(Chaste COMPONENTS cell_based)

# Alternatively, to specify a Chaste installation directory use a line like that below.
# This is needed if your project is not contained in the projects folder within a Chaste source tree.
#find_package(Chaste COMPONENTS heart crypt PATHS /path/to/chaste-install NO_DEFAULT_PATH)

# Change the project name in the line below to match the folder this file is in,
# i.e. the name of your project.
chaste_do_project(ToothFormation)
