/*

Copyright (c) 2005-2017, University of Oxford.
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
#include <thread>

#include "ApicalAndBasalTaggingModifier.hpp"
#include "CellId.hpp"
#include "CellRegionWriter.hpp"
#include "CellsGenerator.hpp"
#include "ChasteMakeUnique.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ContactRegionTaggingModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ExecutableSupport.hpp"
#include "FluidSource.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryMorseInteractionForce.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "ThreeRegionShearForce.hpp"
#include "ThreeRegionSvgWriter.hpp"
#include "TransitCellProliferativeType.hpp"
#include "VarAdhesionMorseMembraneForce.hpp"

#include <boost/make_shared.hpp>
#include "ForwardEulerNumericalMethod.hpp"

// Program option includes for handling command line arguments
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <mesh/src/immersed_boundary/ImmersedBoundaryPalisadeMeshGenerator.hpp>

// Make the variables map global for simplicity
namespace bpo = boost::program_options;
bpo::variables_map gVariablesMap;

/*
 * Prototype functions
 */
void SetupSingletons(unsigned seed);
void DestroySingletons();
void SetupAndRunSimulation();

int main(int argc, char* argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

    // Define command line options
    bpo::options_description general_options("This is a Chaste Immersed Boundary executable.\n");
    general_options.add_options()
        ("help", "produce help message")
        ("ID", bpo::value<std::string>(),"ID string for the simulation")
        ("CRL", bpo::value<double>()->default_value(0.0),"Cortical rest length")
        ("CSC", bpo::value<double>()->default_value(0.0),"Cortical spring constant")
        ("TRL", bpo::value<double>()->default_value(0.0),"Transmembrane rest length")
        ("TSC", bpo::value<double>()->default_value(0.0),"Transmembrane spring constant")
        ("KFS", bpo::value<double>()->default_value(0.0),"Kinematic Feedback strength")
        ("ALM", bpo::value<double>()->default_value(0.0),"Apical lamina stiffness multiplier")
        ("DI", bpo::value<double>()->default_value(0.0),"Interaction distance for cell-cell forces")
        ("SM", bpo::value<double>()->default_value(0.0),"Stiffness multiplier for membrane forces")
        ("NS", bpo::value<double>()->default_value(0.0),"Standard deviation for normal noise")
        ("RM", bpo::value<unsigned>()->default_value(1u),"ReMesh frequency")
        ("TS", bpo::value<unsigned>()->default_value(1000u),"Number of time steps")
        ("AL", bpo::value<bool>()->default_value(false),"Whether to include apical lamina");

    // Define parse command line into variables_map
    bpo::store(parse_command_line(argc, argv, general_options), gVariablesMap);

    // Print help message if wanted
    if (gVariablesMap.count("help"))
    {
        std::cout << setprecision(3) << general_options << "\n";
        std::cout << general_options << "\n";
        return 1;
    }

    // Get ID string from command line
    std::string id_string = gVariablesMap["ID"].as<std::string>();

    std::cout << "Starting simulation with ID string " << id_string << std::endl;
    SetupSingletons(0u);
    SetupAndRunSimulation();
    DestroySingletons();
    std::cout << "Completed simulation with ID string " << id_string << std::endl;
}

void SetupSingletons(unsigned seed)
{
    // Set up what the test suite would do
    SimulationTime::Instance()->SetStartTime(0.0);

    // Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
    RandomNumberGenerator::Instance()->Reseed(seed);
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

void SetupAndRunSimulation()
{
    // Get all variables from the global variables map
    std::string id_string = gVariablesMap["ID"].as<std::string>();
    double cor_rest_length = gVariablesMap["CRL"].as<double>();
    double cor_spring_const = gVariablesMap["CSC"].as<double>();
    double tra_rest_length = gVariablesMap["TRL"].as<double>();
    double tra_spring_const = gVariablesMap["TSC"].as<double>();
    double kin_feedback_str = gVariablesMap["KFS"].as<double>();
    double apical_lam_mult = gVariablesMap["ALM"].as<double>();
    double interaction_dist = gVariablesMap["DI"].as<double>();
    double stiffness_mult = gVariablesMap["SM"].as<double>();
    double normal_std = gVariablesMap["NS"].as<double>();
    unsigned remesh_freq = gVariablesMap["RM"].as<unsigned>();
    unsigned num_time_steps = gVariablesMap["TS"].as<unsigned>();
    bool apical_lamina = gVariablesMap["AL"].as<bool>();

    /*
     * 1: Num cells
     * 2: Num nodes per cell
     * 3: Superellipse exponent
     * 4: Superellipse aspect ratio
     * 5: Random y-variation
     * 6: Include basal lamina
     * 7: Include apical lamina
     * 8: Use a leaky lamina
     * 9: Num fluid mesh points: overrides nodes per cell
     */
    ImmersedBoundaryPalisadeMeshGenerator gen(15u, 128u, 0.05, 2.0, 0.0, true, apical_lamina, false, 192u);
    ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

    std::cout << p_mesh->GetSpacingRatio() << std::endl;

    std::vector<CellPtr> cells;
    CellsGenerator<NoCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

    ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
    cell_population.SetIfPopulationHasActiveSources(false);
    cell_population.SetInteractionDistance(interaction_dist);

    cell_population.SetReMeshFrequency(remesh_freq);
    cell_population.SetOutputNodeRegionToVtk(true);

    OffLatticeSimulation<2> simulator(cell_population);
    simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2, 2>>());
    simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

    // Add main immersed boundary simulation modifier
    auto p_main_modifier = boost::make_shared<ImmersedBoundarySimulationModifier<2>>();
    simulator.AddSimulationModifier(p_main_modifier);

    auto p_svg_writer = boost::make_shared<ThreeRegionSvgWriter<2>>();
    simulator.AddSimulationModifier(p_svg_writer);

    if (apical_lamina)
    {
        simulator.AddSimulationModifier(boost::make_shared<ApicalAndBasalTaggingModifier<2>>());
    }
    else
    {
        simulator.AddSimulationModifier(boost::make_shared<ContactRegionTaggingModifier<2>>());
    }

    // Add force laws
    auto p_boundary_force = boost::make_shared<VarAdhesionMorseMembraneForce<2>>();
    p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
    p_boundary_force->SetElementWellDepth(cor_spring_const);
    p_boundary_force->SetElementRestLength(cor_rest_length);
    p_boundary_force->SetLaminaWellDepth(2.0 * cor_spring_const);
    p_boundary_force->SetLaminaRestLength(cor_rest_length);
    p_boundary_force->SetApicalWellDepthMult(apical_lam_mult);
    p_boundary_force->SetStiffnessMult(stiffness_mult);
//    SetLaminaWellDepthMult(1.0 + 0.4 * (std::strtod(idString.c_str(), nullptr)));

    auto p_shear_force = boost::make_shared<ThreeRegionShearForce<2>>();
    p_main_modifier->AddImmersedBoundaryForce(p_shear_force);
    p_shear_force->SetSpringConst(kin_feedback_str);

//    auto p_cell_cell_force = boost::make_shared<ImmersedBoundaryMorseInteractionForce<2>>();
//    p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
//    p_cell_cell_force->SetWellDepth(tra_spring_const);
//    p_cell_cell_force->SetRestLength(interaction_dist * tra_rest_length);
//    p_cell_cell_force->SetLaminaRestLengthMult(2.0);
//    p_cell_cell_force->SetLaminaWellDepthMult(2.0);
//    p_cell_cell_force->SetAdditiveNormalNoise(true);
//    p_cell_cell_force->SetNormalNoiseMean(0.0);
//    p_cell_cell_force->SetNormalNoiseStdDev(normal_std);

    // Create and set an output directory that is different for each simulation
    std::stringstream output_directory;
    output_directory << "tooth_formation/Exe_BendingThreeRegion/sim/" << id_string;
    simulator.SetOutputDirectory(output_directory.str());

    // Calculate sampling multiple to have at least 5 frames per second on a 10 second video
    unsigned sampling_multiple = std::max(1u, static_cast<unsigned>(std::floor(num_time_steps / (10.0 * 2.0))));

    // Set simulation properties
    double dt = 0.01;
    simulator.SetDt(dt);
    simulator.SetSamplingTimestepMultiple(UINT_MAX);
    simulator.SetEndTime(num_time_steps * dt);
    p_svg_writer->SetSamplingMultiple(sampling_multiple);

    if (std::stoi(id_string) % std::thread::hardware_concurrency() == 0)
    {
        ProgressReporter &r_progress = simulator.rSetUpAndGetProgressReporter();
        r_progress.SetOutputToConsole(true);
    }

    double vol_start = p_mesh->GetVolumeOfElement(4u);

    try
    {
        simulator.Solve();
    }
    catch(const Exception& e)
    {
        std::cout << e.GetMessage() << std::endl;
    }

    // Calculate the maximal height variation along the lamina
    double max_y = 0.0;
    double min_y = 1.0;

    for (unsigned node_idx = 0; node_idx < p_mesh->GetLamina(0)->GetNumNodes(); node_idx++)
    {
        double this_y = p_mesh->GetLamina(0)->GetNode(node_idx)->rGetLocation()[1];

        if (this_y > max_y)
        {
            max_y = this_y;
        }
        if (this_y < min_y)
        {
            min_y = this_y;
        }
    }

    double vol_end = p_mesh->GetVolumeOfElement(4u);

    PRINT_VARIABLE(vol_end / vol_start);

    OutputFileHandler results_handler(output_directory.str(), false);
    out_stream results_file = results_handler.OpenOutputFile("results.csv");

    // Output summary statistics to results file
    (*results_file) << "id,"
                    << "max_y_var,"
                    << "asymmetry_measure" << std::endl;

    (*results_file) << id_string << ","
                    << boost::lexical_cast<std::string>(max_y - min_y) << ","
                    << p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(4, unit_vector<double>(2, 0));

    // Tidy up
    results_file->close();
}
