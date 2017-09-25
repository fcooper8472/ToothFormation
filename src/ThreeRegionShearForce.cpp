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

#include "ThreeRegionShearForce.hpp"
#include "ImmersedBoundaryEnumerations.hpp"
#include "UblasCustomFunctions.hpp"

template <unsigned DIM>
ThreeRegionShearForce<DIM>::ThreeRegionShearForce()
        : AbstractImmersedBoundaryForce<DIM>(),
          mSpringConst(1e3),
          mExcludedRegions({LEFT_LATERAL_REGION, RIGHT_LATERAL_REGION, LEFT_BASAL_REGION, RIGHT_BASAL_REGION, LAMINA_REGION, NO_REGION})
{
}

template <unsigned DIM>
void ThreeRegionShearForce<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                                                                       ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    // Calculate which of the three regions the element is in
    auto middle_region = [&rCellPopulation](const unsigned& elem_idx) -> bool
    {
        EXCEPT_IF_NOT(rCellPopulation.GetNumElements() % 3 == 0);
        unsigned num_elems_per_region = rCellPopulation.GetNumElements() / 3;
        return elem_idx >= num_elems_per_region && elem_idx < 2 * num_elems_per_region;
    };

    // Calculate force only if neither node is in a lamina, and if nodes are in different elements
    auto condition_satisfied = [&](const std::pair<Node<DIM>*, Node<DIM>*>& pair) -> bool
    {
        // If any of the excluded regions match the region of either node, return false
        if (std::any_of(std::begin(mExcludedRegions), std::end(mExcludedRegions),
                        [&](unsigned i){return i == pair.first->GetRegion() || i == pair.second->GetRegion();}))
        {
            return false;
        }
        // This force only acts on nodes in different elements
        if (*(pair.first->ContainingElementsBegin()) == *(pair.second->ContainingElementsBegin()))
        {
            return false;
        }
        // If both nodes are in the middle region, return false
        if (middle_region(*(pair.first->ContainingElementsBegin())) &&
            middle_region(*(pair.second->ContainingElementsBegin())))
        {
            return false;
        }
        // Return true only if the nodes are within threshold distance
        auto vec_a2b = rCellPopulation.rGetMesh().GetVectorFromAtoB(pair.first->rGetLocation(), pair.second->rGetLocation());
        return norm_2(vec_a2b) < rCellPopulation.GetInteractionDistance();
    };

    for (auto&& pair : rNodePairs)
    {
        if (condition_satisfied(pair))
        {
            // Determine outer and inner nodes by finding element closest in index to the centre
            auto half_num_elems = static_cast<int>(rCellPopulation.GetNumElements() / 2);
            Node<DIM>* p_outer_node = std::abs(half_num_elems - static_cast<int>(*(pair.first->ContainingElementsBegin()))) >
                                      std::abs(half_num_elems - static_cast<int>(*(pair.second->ContainingElementsBegin()))) ?
                                      pair.first : pair.second;
            Node<DIM>* p_inner_node = p_outer_node == pair.first ? pair.second : pair.first;

            // Get a unit vector from the previous node to the next node, oriented in the positive y direction
            auto GetChord = [&rCellPopulation](const Node<DIM>& rNode) -> c_vector<double, DIM>
            {
                const auto& r_elem = *rCellPopulation.GetElement(*(rNode.ContainingElementsBegin()));
                const unsigned local_idx = r_elem.GetNodeLocalIndex(rNode.GetIndex());
                const unsigned next_idx = (local_idx + 1u) % r_elem.GetNumNodes();
                const unsigned prev_idx = ((local_idx + r_elem.GetNumNodes()) - 1) % r_elem.GetNumNodes();
                auto vec_a2b = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_elem.GetNode(prev_idx)->rGetLocation(),
                                                                            r_elem.GetNode(next_idx)->rGetLocation());
                if (vec_a2b[1] < 0.0)
                {
                    vec_a2b *= -1.0;
                }
                return vec_a2b;
            };

            const c_vector<double, DIM> average_chord = 0.5 * (GetChord(*p_inner_node) + GetChord(*p_outer_node));

            const unsigned outer_elem_idx = *(p_outer_node->ContainingElementsBegin());
            const unsigned inner_elem_idx = *(p_inner_node->ContainingElementsBegin());

            const double outer_elem_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(outer_elem_idx, false);
            const double inner_elem_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(inner_elem_idx, false);
            const double elem_spacing = 0.5 * (outer_elem_spacing + inner_elem_spacing);

            const double eff_spring_const = mSpringConst * elem_spacing / rCellPopulation.GetIntrinsicSpacing();

            /*
             * We must scale each applied force by a factor of elem_spacing / local spacing, so that forces
             * balance when spread to the grid later (where the multiplicative factor is the local spacing)
             */
            c_vector<double, DIM> outer_force = (eff_spring_const * elem_spacing / outer_elem_spacing) * average_chord;
            p_outer_node->AddAppliedForceContribution(outer_force);

            c_vector<double, DIM> inner_force = (-1.0 * eff_spring_const * elem_spacing / inner_elem_spacing) * average_chord;
            p_inner_node->AddAppliedForceContribution(inner_force);
        }
    }

    if (this->mAdditiveNormalNoise)
    {
        this->AddNormalNoiseToNodes(rCellPopulation);
    }
}

template<unsigned DIM>
void ThreeRegionShearForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SpringConstant>" << mSpringConst << "</SpringConstant>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template<unsigned DIM>
double ThreeRegionShearForce<DIM>::GetSpringConst() const
{
    return mSpringConst;
}

template<unsigned DIM>
void ThreeRegionShearForce<DIM>::SetSpringConst(double springConst)
{
    mSpringConst = springConst;
}


// Explicit instantiation
template class ThreeRegionShearForce<1>;
template class ThreeRegionShearForce<2>;
template class ThreeRegionShearForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ThreeRegionShearForce)
