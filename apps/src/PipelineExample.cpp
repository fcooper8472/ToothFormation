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

#include <boost/make_shared.hpp>
#include <cxxtest/TestSuite.h>
#include <thread>

#include "CellId.hpp"
#include "CellsGenerator.hpp"
#include "ChasteMakeUnique.hpp"
#include "ExecutableSupport.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundarySvgWriter.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "OutputFileHandler.hpp"
#include "SuperellipseGenerator.hpp"

// Program option includes for handling command line arguments
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

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
        ("FS", bpo::value<double>()->default_value(1e2),"Interaction strength")
        ("RL", bpo::value<double>()->default_value(0.25),"Rest length")
        ("TS", bpo::value<unsigned>()->default_value(1000u),"Number of time steps");

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
    unsigned id_val = std::stoi(id_string);

    std::cout << "Starting simulation with ID string " << id_string << std::endl;
    SetupSingletons(id_val);
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
    double force_strength = gVariablesMap["FS"].as<double>();
    double rest_length = gVariablesMap["RL"].as<double>();
    unsigned num_time_steps = gVariablesMap["TS"].as<unsigned>();

    PRINT_4_VARIABLES(id_string, force_strength, rest_length, num_time_steps);

    /*
     * 1: num nodes
     * 2: superellipse exponent
     * 3: cell width
     * 4: cell height
     * 5: bottom left x
     * 6: bottom left y
     */
    SuperellipseGenerator gen(128, 1.0, 0.4, 0.6, 0.3, 0.2);
    const std::vector<c_vector<double, 2> > locations = gen.GetPointsAsVectors();

    std::vector<Node<2>* > nodes;
    std::vector<ImmersedBoundaryElement<2,2>* > elements;

    for (unsigned location = 0; location < locations.size(); location++)
    {
        nodes.push_back(new Node<2>(location, locations[location], true));
    }

    elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));

    auto p_mesh = our::make_unique<ImmersedBoundaryMesh<2,2>>(nodes, elements);
    p_mesh->SetNumGridPtsXAndY(64u);

    std::vector<CellPtr> cells;
    CellsGenerator<NoCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

    ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
    cell_population.SetIfPopulationHasActiveSources(false);

    OffLatticeSimulation<2> simulator(cell_population);
    simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2, 2>>());
    simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

    // Add svg writer
    auto p_svg_writer = boost::make_shared<ImmersedBoundarySvgWriter<2>>();
    simulator.AddSimulationModifier(p_svg_writer);
    p_svg_writer->SetSamplingMultiple(num_time_steps / 8);

    // Add main immersed boundary simulation modifier and random noise
    auto p_main_modifier = boost::make_shared<ImmersedBoundarySimulationModifier<2>>();
    simulator.AddSimulationModifier(p_main_modifier);

    // Add force law
    auto p_boundary_force = boost::make_shared<ImmersedBoundaryLinearMembraneForce<2>>();
    p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
    p_boundary_force->SetElementSpringConst(force_strength);
    p_boundary_force->SetElementRestLength(rest_length);

    // The output directory
    const std::string output_dir = "TestPipelineExample/sim/" + id_string;

    // Set simulation properties
    double dt = 0.01;
    simulator.SetOutputDirectory(output_dir);
    simulator.SetDt(dt);
    simulator.SetSamplingTimestepMultiple(UINT_MAX);
    simulator.SetEndTime(dt * num_time_steps);

    simulator.Solve();

    OutputFileHandler results_handler(output_dir, false);
    out_stream results_file = results_handler.OpenOutputFile("results.csv");

    // Output summary statistics to results file
    (*results_file) << "id,esf" << std::endl;
    (*results_file) << id_string << ',' << std::to_string(p_mesh->GetElongationShapeFactorOfElement(0u)) << std::endl;

    results_file->close();
}
