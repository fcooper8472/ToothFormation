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
#include "ImmersedBoundaryEnumerations.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryMorseInteractionForce.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "MathsCustomFunctions.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "ThreeRegionShearForce.hpp"
#include "ThreeRegionSvgWriter.hpp"
#include "ThreeRegionInteractionForces.hpp"
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
        ("CRL", bpo::value<double>()->default_value(0.25),"Cortical rest length")
        ("CSC", bpo::value<double>()->default_value(0.0),"Cortical spring constant")
        ("SUP", bpo::value<double>()->default_value(0.0),"Support strength")
        ("TRL", bpo::value<double>()->default_value(0.25),"Transmembrane rest length")
        ("TSC", bpo::value<double>()->default_value(0.0),"Transmembrane spring constant")
        ("KFS", bpo::value<double>()->default_value(0.0),"Kinematic Feedback strength")
        ("ALM", bpo::value<double>()->default_value(1.0),"Apical lamina stiffness multiplier")
        ("DI", bpo::value<double>()->default_value(0.0),"Interaction distance for cell-cell forces")
        ("SM", bpo::value<double>()->default_value(0.0),"Stiffness multiplier for membrane forces")
        ("AAM", bpo::value<double>()->default_value(0.0),"Apical-apical interaction multiplier")
        ("NS", bpo::value<double>()->default_value(0.0),"Standard deviation for normal noise")
        ("CYF", bpo::value<double>()->default_value(1.0),"Cyclic frequency for stiffness modulation")
        ("GOP", bpo::value<double>()->default_value(0.5),"Gradient on proportion")
        ("DF", bpo::value<double>()->default_value(0.0),"Diagonal fraction")
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
    double cor_rest_length = gVariablesMap["CRL"].as<double>();
    double cor_spring_const = gVariablesMap["CSC"].as<double>();
    double supp_strength = gVariablesMap["SUP"].as<double>();
    double tra_rest_length = gVariablesMap["TRL"].as<double>();
    double tra_spring_const = gVariablesMap["TSC"].as<double>();
    double kin_feedback_str = gVariablesMap["KFS"].as<double>();
    double apical_lam_mult = gVariablesMap["ALM"].as<double>();
    double interaction_dist = gVariablesMap["DI"].as<double>();
    double stiffness_mult = gVariablesMap["SM"].as<double>();
    double apical_apical_mult = gVariablesMap["AAM"].as<double>();
    double normal_std = gVariablesMap["NS"].as<double>();
    double cyclic_frequency = gVariablesMap["CYF"].as<double>();
    double grad_on_prop = gVariablesMap["GOP"].as<double>();
    double diagonal_fraction = gVariablesMap["DF"].as<double>();
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
     * 10: Absolute gap between cells
     */
    ImmersedBoundaryPalisadeMeshGenerator gen(15u, 128u, 0.05, 2.0, 0.0, true, apical_lamina, false, 192u, interaction_dist * tra_rest_length);
    ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

    std::array<unsigned, 3> region_sizes = {{6u, 3u, 6u}};

    std::cout << p_mesh->GetSpacingRatio() << std::endl;

    std::vector<CellPtr> cells;
    CellsGenerator<NoCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

    ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
    cell_population.SetIfPopulationHasActiveSources(true);
    cell_population.SetInteractionDistance(interaction_dist);
    cell_population.SetReMeshFrequency(remesh_freq);
    cell_population.SetOutputNodeRegionToVtk(true);

    OffLatticeSimulation<2> simulator(cell_population);
    simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2, 2>>());
    simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

    // Add main immersed boundary simulation modifier and random noise
    auto p_main_modifier = boost::make_shared<ImmersedBoundarySimulationModifier<2>>();
    p_main_modifier->SetNoiseLengthScale(0.03);
    p_main_modifier->SetNoiseSkip(4u);
    p_main_modifier->SetNoiseStrength(normal_std);
    p_main_modifier->SetAdditiveNormalNoise(true);
    p_main_modifier->SetZeroFieldSums(true);
    simulator.AddSimulationModifier(p_main_modifier);

    auto p_svg_writer = boost::make_shared<ThreeRegionSvgWriter<2>>();
    p_svg_writer->SetRegionSizes(region_sizes);
    simulator.AddSimulationModifier(p_svg_writer);

    if (apical_lamina)
    {
        simulator.AddSimulationModifier(boost::make_shared<ApicalAndBasalTaggingModifier<2>>());
    }
    else
    {
        simulator.AddSimulationModifier(boost::make_shared<ContactRegionTaggingModifier<2>>());
    }

    auto p_area_modifier = boost::make_shared<SimpleTargetAreaModifier<2>>();
    simulator.AddSimulationModifier(p_area_modifier);
    p_area_modifier->SetReferenceTargetArea(p_mesh->GetVolumeOfElement(4));
    p_area_modifier->SetGrowthDuration(0.0);

    // Add force laws
    auto p_boundary_force = boost::make_shared<VarAdhesionMorseMembraneForce<2>>();
    p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
    p_boundary_force->SetElementWellDepth(cor_spring_const);
    p_boundary_force->SetElementRestLength(cor_rest_length);
    p_boundary_force->SetLaminaWellDepth(2.0 * cor_spring_const);
    p_boundary_force->SetLaminaRestLength(cor_rest_length);
//    p_boundary_force->SetApicalWellDepthMult(apical_lam_mult);
    p_boundary_force->SetStiffnessMult(stiffness_mult);
    p_boundary_force->SetSupportStrength(supp_strength);
    p_boundary_force->SetDiagonalFraction(diagonal_fraction);
    p_boundary_force->SetRegionSizes(region_sizes);
    p_boundary_force->SetCyclicFrequency(cyclic_frequency);
    p_boundary_force->SetGradientOnProportion(grad_on_prop);

    auto p_cell_cell_force = boost::make_shared<ThreeRegionInteractionForces<2>>();
    p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
    p_cell_cell_force->SetBasicInteractionStrength(tra_spring_const);
    p_cell_cell_force->SetBasicInteractionDist(tra_rest_length);
    p_cell_cell_force->SetAdhesionMultiplier(apical_apical_mult);
    p_cell_cell_force->SetRegionSizes(region_sizes);


    // Create and set an output directory that is different for each simulation
    std::stringstream output_directory;
    output_directory << "tooth_formation/Exe_BendingThreeRegion/sim/" << id_string;
    simulator.SetOutputDirectory(output_directory.str());

    // Calculate sampling multiple to have at least 5 frames per second on a 15 second video
    const double num_secs = 15.0;
    const double num_fps = 5.0;
    unsigned sampling_multiple = std::max(1u, static_cast<unsigned>(std::floor(num_time_steps / (num_secs * num_fps))));

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
    std::vector<double> lamina_y_vals;
    for (unsigned node_idx = 0; node_idx < p_mesh->GetLamina(0)->GetNumNodes(); node_idx++)
    {
        lamina_y_vals.push_back(p_mesh->GetLamina(0)->GetNode(node_idx)->rGetLocation()[1]);
    }
    auto minmax_y_vals = std::minmax_element(lamina_y_vals.begin(), lamina_y_vals.end());

    double vol_end = p_mesh->GetVolumeOfElement(4u);

    // Calculate the average lean of cells
    std::vector<double> leans;
    for (unsigned elem_idx = 0; elem_idx < p_mesh->GetNumElements(); ++elem_idx)
    {
        // Pointer to element
        auto p_elem = p_mesh->GetElement(elem_idx);

        // Pointers to corner nodes
        auto p_lt_basal = p_elem->rGetCornerNodes()[LEFT_BASAL_CORNER];
        auto p_lt_apical = p_elem->rGetCornerNodes()[LEFT_APICAL_CORNER];
        auto p_rt_basal = p_elem->rGetCornerNodes()[RIGHT_BASAL_CORNER];
        auto p_rt_apical = p_elem->rGetCornerNodes()[RIGHT_APICAL_CORNER];

        if (p_lt_basal == nullptr || p_lt_apical == nullptr || p_rt_basal == nullptr || p_rt_apical == nullptr)
        {
            continue;
        }

        unsigned lb_idx = p_elem->GetNodeLocalIndex(p_lt_basal->GetIndex());
        unsigned la_idx = p_elem->GetNodeLocalIndex(p_lt_apical->GetIndex());
        unsigned rb_idx = p_elem->GetNodeLocalIndex(p_rt_basal->GetIndex());
        unsigned ra_idx = p_elem->GetNodeLocalIndex(p_rt_apical->GetIndex());

        unsigned lt_third = SmallDifferenceMod(lb_idx, la_idx, p_elem->GetNumNodes()) / 3;
        unsigned rt_third = SmallDifferenceMod(rb_idx, ra_idx, p_elem->GetNumNodes()) / 3;

        // Middle third, anticlockwise from left apical corner
        auto lt_vec = p_mesh->GetVectorFromAtoB(
                p_elem->GetNode(AdvanceMod(la_idx, lt_third + lt_third, p_elem->GetNumNodes()))->rGetLocation(),
                p_elem->GetNode(AdvanceMod(la_idx, lt_third, p_elem->GetNumNodes()))->rGetLocation()
        );

        // Middle third, clockwise from right apical corner
        auto rt_vec = p_mesh->GetVectorFromAtoB(
                p_elem->GetNode(AdvanceMod(ra_idx, -rt_third - rt_third, p_elem->GetNumNodes()))->rGetLocation(),
                p_elem->GetNode(AdvanceMod(ra_idx, -rt_third, p_elem->GetNumNodes()))->rGetLocation()
        );

        if (elem_idx < region_sizes[0])
        {
            leans.emplace_back(std::atan(lt_vec[0] / lt_vec[1]));
            leans.emplace_back(std::atan(rt_vec[0] / rt_vec[1]));
        }
        else if (elem_idx >= region_sizes[0] + region_sizes[1])
        {
            leans.emplace_back(std::atan(-lt_vec[0] / lt_vec[1]));
            leans.emplace_back(std::atan(-rt_vec[0] / rt_vec[1]));
        }
    }

    const double lean_mean = std::accumulate(leans.begin(), leans.end(), 0.0) / leans.size();
    const double lean_var = std::inner_product(leans.begin(), leans.end(), leans.begin(), 0.0) / leans.size() - lean_mean * lean_mean;

    // Calculate the average (absolute) skew of non-central elements about long axis
    std::vector<double> abs_skews;
    for (unsigned elem_idx = 0; elem_idx < p_mesh->GetNumElements(); ++elem_idx)
    {
        // Ignore the central cells
        if (elem_idx < region_sizes[0] || elem_idx >= region_sizes[0] + region_sizes[1])
        {
            // Get and orient the short axis
            const auto skew = p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(elem_idx, unit_vector<double>(2, 0));
            abs_skews.push_back(std::fabs(skew));
        }
    }

    const double skew_mean = std::accumulate(abs_skews.begin(), abs_skews.end(), 0.0) / abs_skews.size();
    const double skew_var = std::inner_product(abs_skews.begin(), abs_skews.end(), abs_skews.begin(), 0.0) / abs_skews.size() - skew_mean * skew_mean;

    // Lamina angle is angle of triangle with base 0.5 and height of max y variation
    const double lamina_angle = std::atan((*minmax_y_vals.second - *minmax_y_vals.first) / 0.5);
    const double lean_ratio = lean_mean / lamina_angle;

    PRINT_VECTOR(leans);
    PRINT_VARIABLE(vol_end / vol_start);

    OutputFileHandler results_handler(output_directory.str(), false);
    out_stream results_file = results_handler.OpenOutputFile("results.csv");

    // Generate mp4 from svg sequence, assuming ffmpeg has been built with librsvg enabled
    {
        const std::string svg_dir = results_handler.GetOutputDirectoryFullPath() + "results_from_time_0/";
        const std::string number = std::stoi(id_string) < 10 ? '0' + id_string : id_string;
        const std::string mp4_name = svg_dir + number + ".mp4";
        const std::string command =
                "ffmpeg -v 0 -r " +
                std::to_string(static_cast<unsigned>(num_fps)) +
                " -pattern_type glob -i \"" +
                svg_dir +
                "*.svg\" -c:v libx264 -pix_fmt yuv420p -crf 0 -preset slow -y " +
                mp4_name +
                "> /dev/null";

        std::cout << "C++: Generating mp4" << std::endl;
        std::system(command.c_str());
        std::cout << "C++: Finished generating mp4"<< std::endl;
    }

    // Output summary statistics to results file
    (*results_file) << "id,"
                    << "max_y_var,"
                    << "skew_mean,"
                    << "skew_std,"
                    << "lean_mean,"
                    << "lean_std,"
                    << "lean_ratio"
                    << std::endl;


    (*results_file) << id_string << ","
                    << std::to_string(*minmax_y_vals.second - *minmax_y_vals.first) << ","
                    << std::to_string(skew_mean) << ","
                    << std::to_string(std::sqrt(skew_var)) << ","
                    << std::to_string(lean_mean) << ","
                    << std::to_string(std::sqrt(lean_var)) << ","
                    << std::to_string(lean_ratio);
    // Tidy up
    results_file->close();


}
