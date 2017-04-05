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

#include "AdditiveNormalLocationModifier.hpp"
#include "ApicalAndBasalTaggingModifier.hpp"
#include "CellId.hpp"
#include "CellRegionWriter.hpp"
#include "CellsGenerator.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ContactRegionTaggingModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ExecutableSupport.hpp"
#include "FluidSource.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryMorseInteractionForce.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "ThreeRegionInteractionForces.hpp"
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

/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(std::string idString, double corRestLength, double corSpringConst, double traRestLength,
                           double traSpringConst, double rhsAdhesionMod, double interactionDist, double stiffnessMult,
                           unsigned reMeshFreq, unsigned numTimeSteps, bool apicalLamina);
void OutputToConsole(std::string idString, std::string leading);

int main(int argc, char* argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

    // Define command line options
    boost::program_options::options_description general_options("This is a Chaste Immersed Boundary executable.\n");
    general_options.add_options()
                    ("help", "produce help message")
                    ("ID", boost::program_options::value<std::string>(),"ID string for the simulation")
                    ("CRL", boost::program_options::value<double>()->default_value(0.0),"Cortical rest length")
                    ("CSC", boost::program_options::value<double>()->default_value(0.0),"Cortical spring constant")
                    ("TRL", boost::program_options::value<double>()->default_value(0.0),"Transmembrane rest length")
                    ("TSC", boost::program_options::value<double>()->default_value(0.0),"Transmembrane spring constant")
                    ("AD", boost::program_options::value<double>()->default_value(0.0),"Adhesion modifier")
                    ("DI", boost::program_options::value<double>()->default_value(0.0),"Interaction distance for cell-cell forces")
                    ("SM", boost::program_options::value<double>()->default_value(0.0),"Stiffness multiplier for membrane forces")
                    ("RM", boost::program_options::value<unsigned>()->default_value(1u),"ReMesh frequency")
                    ("TS", boost::program_options::value<unsigned>()->default_value(1000u),"Number of time steps")
                    ("AL", boost::program_options::value<bool>()->default_value(false),"Whether to include apical lamina");

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
    double stiffness_mult = variables_map["SM"].as<double>();
    unsigned remesh_freq = variables_map["RM"].as<unsigned>();
    unsigned num_time_steps = variables_map["TS"].as<unsigned>();
    bool apical_lamina = variables_map["AL"].as<bool>();

    OutputToConsole(id_string, "Started");
    SetupSingletons();
    SetupAndRunSimulation(id_string, cor_rest_length, cor_spring_const, tra_rest_length, tra_spring_const,
                          rhs_adhesion_mod, interaction_dist, stiffness_mult, remesh_freq, num_time_steps,
                          apical_lamina);
    DestroySingletons();
    OutputToConsole(id_string, "Completed");
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

void OutputToConsole(std::string idString, std::string leading)
{
    // Compose the message
    std::stringstream message;
    message << leading << " simulation with ID string " << idString << std::endl;

    // Send it to the console
    std::cout << message.str() << std::flush;
}

void SetupAndRunSimulation(std::string idString, double corRestLength, double corSpringConst, double traRestLength,
                           double traSpringConst, double rhsAdhesionMod, double interactionDist, double stiffnessMult,
                           unsigned reMeshFreq, unsigned numTimeSteps, bool apicalLamina)
{
    /*
     * 1: Num cells
     * 2: Num nodes per cell
     * 3: Superellipse exponent
     * 4: Superellipse aspect ratio
     * 5: Random y-variation
     * 6: Include basal lamina
     * 7: Include apical lamina
     */
    ImmersedBoundaryPalisadeMeshGenerator gen(15, 128, 0.1, 2.0, 0.0, true, apicalLamina);
    ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

    p_mesh->SetNumGridPtsXAndY(256);

    std::cout << p_mesh->GetSpacingRatio() << std::endl;

    std::vector<CellPtr> cells;
    CellsGenerator<NoCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

    ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
    cell_population.SetIfPopulationHasActiveSources(false);
    cell_population.SetInteractionDistance(interactionDist);

    cell_population.SetReMeshFrequency(reMeshFreq);
    cell_population.SetOutputNodeRegionToVtk(true);

    OffLatticeSimulation<2> simulator(cell_population);
    simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2, 2> >());
    simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

    // Add normal location modifier first, so it happens before force calculation
    MAKE_PTR(AdditiveNormalLocationModifier<2>, p_noise);
    p_noise->SetStdDev(0.25 * p_mesh->GetCharacteristicNodeSpacing());
    simulator.AddSimulationModifier(p_noise);

    // Add main immersed boundary simulation modifier
    MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
    simulator.AddSimulationModifier(p_main_modifier);

    MAKE_PTR(ThreeRegionSvgWriter<2>, p_svg_writer);
    simulator.AddSimulationModifier(p_svg_writer);

    if (apicalLamina)
    {
        MAKE_PTR(ApicalAndBasalTaggingModifier<2>, p_tagger);
        simulator.AddSimulationModifier(p_tagger);
    }
    else
    {
        MAKE_PTR(ContactRegionTaggingModifier<2>, p_tagger);
        simulator.AddSimulationModifier(p_tagger);
    }

    // Add force laws
    MAKE_PTR(VarAdhesionMorseMembraneForce<2>, p_boundary_force);
    p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
    p_boundary_force->SetElementWellDepth(corSpringConst);
    p_boundary_force->SetElementRestLength(corRestLength);
    p_boundary_force->SetLaminaWellDepth(2.0 * corSpringConst);
    p_boundary_force->SetLaminaRestLength(corRestLength);
    p_boundary_force->SetStiffnessMult(stiffnessMult);

    MAKE_PTR(ImmersedBoundaryMorseInteractionForce<2>, p_cell_cell_force);
    p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
    p_cell_cell_force->SetWellDepth(traSpringConst);
    p_cell_cell_force->SetRestLength(0.25 * interactionDist * traRestLength);
    p_cell_cell_force->SetLaminaRestLengthMult(0.5);
    p_cell_cell_force->SetLaminaWellDepthMult(2.0);

    // Create and set an output directory that is different for each simulation
    std::stringstream output_directory;
    output_directory << "tooth_formation/Exe_BendingThreeRegion/sim/" << idString;
    simulator.SetOutputDirectory(output_directory.str());

    // Calculate sampling multiple to have at least 5 frames per second on a 15 second video
    unsigned sampling_multiple = std::max(1u, static_cast<unsigned>(std::floor(numTimeSteps / (15.0 * 2.0))));

    // Set simulation properties
    double dt = 0.01;
    simulator.SetDt(dt);
    simulator.SetSamplingTimestepMultiple(UINT_MAX);
    simulator.SetEndTime(numTimeSteps * dt);
    p_svg_writer->SetSamplingMultiple(sampling_multiple);

    simulator.Solve();

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

    OutputFileHandler results_handler(output_directory.str(), false);
    out_stream results_file = results_handler.OpenOutputFile("results.csv");

    // Output summary statistics to results file
    (*results_file) << "id,"
                    << "max_y_var,"
                    << "asymmetry_measure" << std::endl;

    (*results_file) << idString << ","
                    << boost::lexical_cast<std::string>(max_y - min_y) << ","
                    << p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(4, unit_vector<double>(2, 0));

    // Tidy up
    results_file->close();
}
