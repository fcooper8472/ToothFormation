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

#include "SimulationStoppingModifier.hpp"

#include "Exception.hpp"

template <unsigned DIM>
SimulationStoppingModifier<DIM>::SimulationStoppingModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template <unsigned DIM>
void SimulationStoppingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    ImmersedBoundaryMesh<DIM, DIM>* p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>*>(&(rCellPopulation.rGetMesh()));

    if (CalculateLaminaAngle(p_mesh) >= mThreshold)
    {
        EXCEPTION("Lamina bend exceeded; finishing simulation.");
    }
}

template <unsigned DIM>
double SimulationStoppingModifier<DIM>::CalculateLaminaAngle(ImmersedBoundaryMesh<DIM, DIM>* const pIbMesh)
{
    // Calculate the maximal height variation along the lamina
    std::vector<double> lamina_y_vals;
    lamina_y_vals.reserve(pIbMesh->GetLamina(0)->GetNumNodes());

    for (unsigned node_idx = 0; node_idx < pIbMesh->GetLamina(0)->GetNumNodes(); node_idx++)
    {
        lamina_y_vals.push_back(pIbMesh->GetLamina(0)->GetNode(node_idx)->rGetLocation()[1]);
    }
    auto minmax_y_vals = std::minmax_element(lamina_y_vals.begin(), lamina_y_vals.end());

    return std::atan((*minmax_y_vals.second - *minmax_y_vals.first) / 0.5);
}

template <unsigned DIM>
void SimulationStoppingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
}

template <unsigned DIM>
void SimulationStoppingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned int DIM>
double SimulationStoppingModifier<DIM>::GetThreshold() const
{
    return mThreshold;
}

template<unsigned int DIM>
void SimulationStoppingModifier<DIM>::SetThreshold(double threshold)
{
    mThreshold = threshold;
}

// Explicit instantiation
template class SimulationStoppingModifier<1>;
template class SimulationStoppingModifier<2>;
template class SimulationStoppingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimulationStoppingModifier)
