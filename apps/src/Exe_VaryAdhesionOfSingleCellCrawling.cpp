/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include <cxxtest/TestSuite.h>

// Includes from trunk
#include "CellId.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ExecutableSupport.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"

// Includes from Immersed Boundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryCellCellInteractionForce.hpp"

// Program option includes for handling command line arguments
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(std::string idString, double corRestLength, double corSpringConst, double traRestLength,
                           double traSpringConst, double rhsAdhesionMod, double interactionDist,
                           unsigned numTimeSteps);
void OutputOnCompletion(std::string idString);

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StandardStartup(&argc, &argv);

    // Define command line options
    boost::program_options::options_description general_options("This is a Chaste Immersed Boundary executable.\n");
    general_options.add_options()
                    ("help", "produce help message")
                    ("ID", boost::program_options::value<std::string>(),"ID string for the simulation")
                    ("CRL", boost::program_options::value<double>()->default_value(0.0),"Cortical rest length")
                    ("CSC", boost::program_options::value<double>()->default_value(0.0),"Cortical spring constant")
                    ("TRL", boost::program_options::value<double>()->default_value(0.0),"Transmembrane rest length")
                    ("TSC", boost::program_options::value<double>()->default_value(0.0),"Transmembrane spring constant")
                    ("AD", boost::program_options::value<double>()->default_value(0.0),"RHS adhesion modifier")
                    ("DI", boost::program_options::value<double>()->default_value(0.0),"Interaction distance for cell-cell forces")
                    ("TS", boost::program_options::value<unsigned>()->default_value(1000),"Number of time steps");

    // Define parse command line into variables_map
    boost::program_options::variables_map variables_map;
    boost::program_options::store(parse_command_line(argc, argv, general_options), variables_map);

    // Print help message if wanted
    if (variables_map.count("help"))
    {
        std::cout << setprecision(3) << general_options << "\n";
        std::cout << general_options << "\n";
        return 1;
    }

    // Get ID and name from command line
    std::string id_string = variables_map["ID"].as<std::string>();
    double cor_rest_length = variables_map["CRL"].as<double>();
    double cor_spring_const = variables_map["CSC"].as<double>();
    double tra_rest_length = variables_map["TRL"].as<double>();
    double tra_spring_const = variables_map["TSC"].as<double>();
    double rhs_adhesion_mod = variables_map["AD"].as<double>();
    double interaction_dist = variables_map["DI"].as<double>();
    unsigned num_time_steps = variables_map["TS"].as<unsigned>();

    SetupSingletons();
    SetupAndRunSimulation(id_string, cor_rest_length, cor_spring_const, tra_rest_length, tra_spring_const,
                          rhs_adhesion_mod, interaction_dist, num_time_steps);
    DestroySingletons();
    OutputOnCompletion(id_string);
}

void SetupSingletons()
{
    // Set up what the test suite would do
    SimulationTime::Instance()->SetStartTime(0.0);

    // Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
    RandomNumberGenerator::Instance()->Reseed(0);
    CellPropertyRegistry::Instance()->Clear();
    CellId::ResetMaxCellId();
}

void DestroySingletons()
{
    // This is from the tearDown method of the test suite
    SimulationTime::Destroy();
    RandomNumberGenerator::Destroy();
    CellPropertyRegistry::Instance()->Clear();
}

void OutputOnCompletion(std::string idString)
{
    // Compose the message
    std::stringstream message;
    message << "Completed simulation with ID string " << idString << std::endl;

    // Send it to the console
    std::cout << message.str() << std::flush;
}
void SetupAndRunSimulation(std::string idString, double corRestLength, double corSpringConst, double traRestLength,
                           double traSpringConst, double rhsAdhesionMod, double interactionDist,
                           unsigned numTimeSteps)
{
    /*
     * 1: Num cells
     * 2: Num nodes per cell
     * 3: Superellipse exponent
     * 4: Superellipse aspect ratio
     * 5: Random y-variation
     * 6: Include membrane
     */
    ImmersedBoundaryPalisadeMeshGenerator gen(7, 128, 0.1, 2.5, 0.0, true);
    ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

    p_mesh->SetNumGridPtsXAndY(256);

    std::vector<CellPtr> cells;
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
    CellsGenerator<UniformCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

    ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
    cell_population.SetIfPopulationHasActiveSources(false);
    cell_population.SetInteractionDistance(interactionDist);

    OffLatticeSimulation<2> simulator(cell_population);

    // Add main immersed boundary simulation modifier
    MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
    simulator.AddSimulationModifier(p_main_modifier);

    // Add force law
    MAKE_PTR(ImmersedBoundaryLinearMembraneForce<2>, p_boundary_force);
    p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
    p_boundary_force->SetElementSpringConst(corSpringConst);
    p_boundary_force->SetElementRestLength(corRestLength);

    // Create and set an output directory that is different for each simulation
    std::stringstream output_directory;
    output_directory << "numerics_paper/Exe_VaryAdhesionOfSingleCellCrawling/sim/" << idString;
    simulator.SetOutputDirectory(output_directory.str());

    // Set simulation properties
    double dt = 0.075;
    simulator.SetDt(dt);
    simulator.SetSamplingTimestepMultiple(100);
    simulator.SetEndTime(100.0 * dt);
    simulator.Solve();

    // Now we have relaxed the mesh, we reset the start time to zero.  This overwrites the output, so we end up with
    // less clutter
    SimulationTime::Instance()->Destroy();
    SimulationTime::Instance()->SetStartTime(0.0);

    // Add a cell-cell interaction force with the same intrinsic strength as the membrane force
    MAKE_PTR(ImmersedBoundaryCellCellInteractionForce<2>, p_cell_cell_force);
    p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
    p_cell_cell_force->SetSpringConstant(traSpringConst);
    p_cell_cell_force->SetRestLength(traRestLength);
    p_cell_cell_force->UseMorsePotential();

    // Get the centroid of the three relevant cells before anything happens
    c_vector<double, 2> prev_centroid_start = p_mesh->GetCentroidOfElement(2);
    c_vector<double, 2> this_centroid_start = p_mesh->GetCentroidOfElement(3);
    c_vector<double, 2> next_centroid_start = p_mesh->GetCentroidOfElement(4);

    // Get average height of basement lamina
    ChasteCuboid<2> lamina_bounding_box = p_mesh->CalculateBoundingBoxOfElement(0);
    double lamina_height = 0.5 * (lamina_bounding_box.rGetLowerCorner()[1] + lamina_bounding_box.rGetUpperCorner()[1]);

    // Kick the second cell in from the left and set its E-cad level
    unsigned e_cad_location = p_cell_cell_force->rGetProteinNodeAttributeLocations()[0];
    unsigned p_cad_location = p_cell_cell_force->rGetProteinNodeAttributeLocations()[1];

    for (unsigned node_idx = 0; node_idx < p_mesh->GetElement(3)->GetNumNodes(); node_idx++)
    {
        double new_height = lamina_height + 1.05 * (p_mesh->GetElement(3)->GetNode(node_idx)->rGetLocation()[1] - lamina_height);
        p_mesh->GetElement(3)->GetNode(node_idx)->rGetModifiableLocation()[1] = new_height;

        p_mesh->GetElement(3)->GetNode(node_idx)->rGetNodeAttributes()[e_cad_location] = 1.0;
        p_mesh->GetElement(3)->GetNode(node_idx)->rGetNodeAttributes()[p_cad_location] = rhsAdhesionMod;
    }

    // In the top apical domain of the cell directly to the right, add in p_cad
    double cell_four_height = p_mesh->CalculateBoundingBoxOfElement(4).GetWidth(1);
    double cell_four_y_cent = p_mesh->GetCentroidOfElement(4)[1];
    for (unsigned node_idx = 0; node_idx < p_mesh->GetElement(4)->GetNumNodes(); node_idx++)
    {
        if (p_mesh->GetElement(4)->GetNode(node_idx)->rGetLocation()[1] - cell_four_y_cent > 0.45 * cell_four_height)
        {
            p_mesh->GetElement(4)->GetNode(node_idx)->rGetNodeAttributes()[p_cad_location] = rhsAdhesionMod;
        }
    }

    simulator.SetSamplingTimestepMultiple(100);
    simulator.SetEndTime(numTimeSteps * dt);
    simulator.Solve();

    OutputFileHandler results_handler(output_directory.str(), false);
    out_stream results_file = results_handler.OpenOutputFile("results.csv");

    // Get the centroid of the three relevant cells at end of simulation
    c_vector<double, 2> prev_centroid_end = p_mesh->GetCentroidOfElement(2);
    c_vector<double, 2> this_centroid_end = p_mesh->GetCentroidOfElement(3);
    c_vector<double, 2> next_centroid_end = p_mesh->GetCentroidOfElement(4);

    c_vector<double, 2> axis = unit_vector<double>(2,1);
    double prev_skew = p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(2, axis);
    double this_skew = p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(3, axis);
    double next_skew = p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(4, axis);

    // Output summary statistics to results file
    (*results_file) << "id" << ","
                    << "cortical_rest_length" << ","
                    << "cortical_spring_const" << ","
                    << "rhs_adhesion_mod" << ","
                    << "delta_prev_cent" << ","
                    << "delta_this_cent" << ","
                    << "delta_next_cent" << ","
                    << "delta_prev_skew" << ","
                    << "delta_this_skew" << ","
                    << "delta_next_skew" << std::endl;

    (*results_file) << idString << ","
                    << boost::lexical_cast<std::string>(corRestLength) << ","
                    << boost::lexical_cast<std::string>(corSpringConst) << ","
                    << boost::lexical_cast<std::string>(rhsAdhesionMod) << ","
                    << boost::lexical_cast<std::string>(prev_centroid_end[1] - prev_centroid_start[1]) << ","
                    << boost::lexical_cast<std::string>(this_centroid_end[1] - this_centroid_start[1]) << ","
                    << boost::lexical_cast<std::string>(next_centroid_end[1] - next_centroid_start[1]) << ","
                    << boost::lexical_cast<std::string>(prev_skew) << ","
                    << boost::lexical_cast<std::string>(this_skew) << ","
                    << boost::lexical_cast<std::string>(next_skew);

    // Tidy up
    results_file->close();
}
