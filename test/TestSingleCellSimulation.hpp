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

#ifndef TESTSINGLECELLSIMULATION_HPP_
#define TESTSINGLECELLSIMULATION_HPP_

// Needed for the test environment
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "ImmersedBoundaryBiasedMembraneForce.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SuperellipseGenerator.hpp"

#include "ForwardEulerNumericalMethod.hpp"
#include <boost/make_shared.hpp>

#include "Debug.hpp"


// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestShortImmersedBoundarySimulations : public AbstractCellBasedTestSuite
{
public:

    void TestShortSingleCellSim()
    {
        /*
         * 1: num nodes
         * 2: superellipse exponent
         * 3: cell width
         * 4: cell height
         * 5: bottom left x
         * 6: bottom left y
         */
        std::unique_ptr<SuperellipseGenerator> p_gen(new SuperellipseGenerator(128, 0.1, 0.4, 0.6, 0.3, 0.2));
        const std::vector<c_vector<double, 2>> locations = p_gen->GetPointsAsVectors();

        std::vector<Node<2>* > nodes;
        std::vector<ImmersedBoundaryElement<2,2>* > elements;

        for (const auto& location : locations)
        {
            nodes.push_back(new Node<2>(nodes.size(), location, true));
        }

        elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));

        std::unique_ptr<ImmersedBoundaryMesh<2,2>> p_mesh(new ImmersedBoundaryMesh<2,2>(nodes, elements));
        p_mesh->SetNumGridPtsXAndY(64);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(false);
        cell_population.SetReMeshFrequency(20u);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2,2> >());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

        // Add main immersed boundary simulation modifier
        auto p_main_modifier = boost::make_shared<ImmersedBoundarySimulationModifier<2>>();
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force law
        auto p_boundary_force = boost::make_shared<ImmersedBoundaryBiasedMembraneForce<2>>();
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetElementSpringConst(1.0 * 1e7);

        // Set simulation properties
        double dt = 0.01;
        simulator.SetOutputDirectory("tooth_formation/TestSingleCellSimulation");
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(4u);
        simulator.SetEndTime(10000 * dt);

        simulator.Solve();
    }

};


#endif /*TESTSINGLECELLSIMULATION_HPP_*/
