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

#ifndef TESTTHREEPREGIONSIMULATION_HPP_
#define TESTTHREEPREGIONSIMULATION_HPP_

// Needed for the test environment
#include "AbstractCellBasedTestSuite.hpp"

#include "CellRegionWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FluidSource.hpp"
#include "ImmersedBoundaryCellCellInteractionForce.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "ThreeRegionInteractionForces.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"

#include "ThreeRegionMeshGenerator.hpp"

#include <boost/make_shared.hpp>
#include "ForwardEulerNumericalMethod.hpp"

#include "Debug.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestThreeRegionSimulation : public AbstractCellBasedTestSuite
{
public:
    void TestThreeRegionSim()
    {
        /*
         * 1: Num cells
         * 2: Num nodes per cell
         * 3: Superellipse exponent
         * 4: Superellipse aspect ratio
         * 5: Random y-variation
         * 6: Include membrane
         */
        ThreeRegionMeshGenerator gen(12, 128, 0.1, 2.0, 0.0, true);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

        p_mesh->SetNumGridPtsXAndY(256);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(false);
        cell_population.SetInteractionDistance(2.0 * cell_population.GetInteractionDistance());

        cell_population.SetReMeshFrequency(1u);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2, 2> >());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        MAKE_PTR(ImmersedBoundaryLinearMembraneForce<2>, p_boundary_force);
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetElementSpringConst(1.0 * 1e7);

        MAKE_PTR(ThreeRegionInteractionForces<2>, p_cell_cell_force);
        p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetBasicInteractionStrength(1e6);
        cell_population.AddCellWriter<CellRegionWriter>();

        // Set simulation properties
        double dt = 0.01;
        simulator.SetOutputDirectory("TestThreeRegionSim");
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(1000.0 * dt);

        simulator.Solve();
    }
};

#endif /*TESTTHREEPREGIONSIMULATION_HPP_*/