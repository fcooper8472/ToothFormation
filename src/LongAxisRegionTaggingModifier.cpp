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

#include <ImmersedBoundaryEnumerations.hpp>
#include "LongAxisRegionTaggingModifier.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

template <unsigned DIM>
LongAxisRegionTaggingModifier<DIM>::LongAxisRegionTaggingModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template <unsigned DIM>
LongAxisRegionTaggingModifier<DIM>::~LongAxisRegionTaggingModifier()
{
}

template <unsigned DIM>
void LongAxisRegionTaggingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    ImmersedBoundaryMesh<DIM, DIM>* p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>*>(&(rCellPopulation.rGetMesh()));

    // Iterate over elements
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = p_mesh->GetElementIteratorBegin();
         elem_it != p_mesh->GetElementIteratorEnd();
         ++elem_it)
    {
        const c_vector<double, DIM> centroid = p_mesh->GetCentroidOfElement(elem_it->GetIndex());

        // Orient short axis to point in the positive x direction
        c_vector<double, DIM> short_axis = p_mesh->GetShortAxisOfElement(elem_it->GetIndex());
        short_axis = short_axis[0] > 0.0 ? short_axis : -short_axis;

        // Calculate the long axis, oriented upwards
        c_vector<double, DIM> long_axis;
        long_axis[0] = -short_axis[1];
        long_axis[1] = short_axis[0];

        std::vector<std::pair<unsigned, double> > node_dist_map;
        std::vector<bool> node_on_left(elem_it->GetNumNodes(), true);

        for (unsigned node_idx = 0; node_idx < elem_it->GetNumNodes(); ++node_idx)
        {
            const c_vector<double, DIM> vec_to_node = p_mesh->GetVectorFromAtoB(centroid, elem_it->GetNode(node_idx)->rGetLocation());
            const double dist = inner_prod(long_axis, vec_to_node);
            node_dist_map.push_back(std::make_pair(node_idx, dist));

            // Calculate whether the node is on the left or right of the element
            if (inner_prod(short_axis, vec_to_node) > 0.0)
            {
                node_on_left[node_idx] = false;
            }
        }

        std::sort(node_dist_map.begin(), node_dist_map.end(), ComparisonForDistanceMap);

        const double top = node_dist_map.back().second;
        const double bot = node_dist_map.front().second;

        const double threshold = 0.1 * (top - bot);

        for (unsigned pair_idx = 0; pair_idx < node_dist_map.size(); ++pair_idx)
        {
            const unsigned node_idx = node_dist_map[pair_idx].first;
            const double dist = node_dist_map[pair_idx].second;
            const bool on_left = node_on_left[node_idx];

            if (dist > top - threshold)
            {
                elem_it->GetNode(node_idx)->SetRegion(on_left ? LEFT_APICAL_REGION : RIGHT_APICAL_REGION);
            }
            else if (dist > top - 2.0 * threshold)
            {
                elem_it->GetNode(node_idx)->SetRegion(on_left ? LEFT_PERIAPICAL_REGION : RIGHT_PERIAPICAL_REGION);
            }
            else if (dist < bot + threshold)
            {
                elem_it->GetNode(node_idx)->SetRegion(on_left ? LEFT_BASAL_REGION : RIGHT_BASAL_REGION);
            }
            else
            {
                elem_it->GetNode(node_idx)->SetRegion(on_left ? LEFT_LATERAL_REGION : RIGHT_LATERAL_REGION);
            }
        }
    }

}

template <unsigned DIM>
void LongAxisRegionTaggingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
}

template <unsigned DIM>
bool LongAxisRegionTaggingModifier<DIM>::ComparisonForDistanceMap(std::pair<unsigned, double> idxDistPairA,
                                                                  std::pair<unsigned, double> idxDistPairB)
{
    return idxDistPairA.second < idxDistPairB.second;
}

template <unsigned DIM>
void LongAxisRegionTaggingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class LongAxisRegionTaggingModifier<1>;
template class LongAxisRegionTaggingModifier<2>;
template class LongAxisRegionTaggingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LongAxisRegionTaggingModifier)
