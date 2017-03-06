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

#include "ContactRegionTaggingModifier.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

template <unsigned DIM>
ContactRegionTaggingModifier<DIM>::ContactRegionTaggingModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template <unsigned DIM>
ContactRegionTaggingModifier<DIM>::~ContactRegionTaggingModifier()
{
}

template <unsigned DIM>
void ContactRegionTaggingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    ImmersedBoundaryMesh<DIM, DIM>* p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>*>(&(rCellPopulation.rGetMesh()));

    // Just get this once
    double nbr_dist = p_mesh->GetNeighbourDist();

    // Iterate over elements
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = p_mesh->GetElementIteratorBegin();
         elem_it != p_mesh->GetElementIteratorEnd();
         ++elem_it)
    {
        const c_vector<double, DIM> centroid = p_mesh->GetCentroidOfElement(elem_it->GetIndex());

        // Orient short axis to point in the positive x direction
        c_vector<double, DIM> short_axis = p_mesh->GetShortAxisOfElement(elem_it->GetIndex());
        short_axis = short_axis[0] > 0.0 ? short_axis : -short_axis;

        for (unsigned node_idx = 0; node_idx < elem_it->GetNumNodes(); ++node_idx)
        {
            // Get the current node, and determine whether it's on the left or right of the long axis
            Node<DIM>* p_this_node = p_mesh->GetNode(elem_it->GetNodeGlobalIndex(node_idx));
            bool node_on_left = inner_prod(short_axis, p_mesh->GetVectorFromAtoB(centroid, p_this_node->rGetLocation())) < 0.0;

            p_this_node->SetRegion(NO_REGION);

            // The vec of node neighbours includes all within neighbouring boxes: need to check against neighbour dist
            const std::vector<unsigned>& node_neighbours = elem_it->GetNode(node_idx)->rGetNeighbours();

            for(std::vector<unsigned>::const_iterator gbl_idx_it = node_neighbours.begin();
                gbl_idx_it != node_neighbours.end(); ++gbl_idx_it)
            {
                Node<DIM>* p_other_node = p_mesh->GetNode(*gbl_idx_it);

                if (p_mesh->GetDistanceBetweenNodes(p_this_node->GetIndex(), p_other_node->GetIndex()) < nbr_dist)
                {
                    if (p_other_node->GetRegion() == LAMINA_REGION)
                    {
                        // If it's below the centroid it will be the basal lamina
                        if (p_this_node->rGetLocation()[1] < centroid[1])
                        {
                            p_this_node->SetRegion(node_on_left ? LEFT_BASAL_REGION : RIGHT_BASAL_REGION);
                        }
                        else // it's the apical lamina (if present)
                        {
                            p_this_node->SetRegion(node_on_left ? LEFT_APICAL_REGION : RIGHT_APICAL_REGION);
                        }
                    }
                    else // it's a regular element
                    {
                        p_this_node->SetRegion(node_on_left ? LEFT_LATERAL_REGION : RIGHT_LATERAL_REGION);
                    }
                }
            }
        }
    }

}

template <unsigned DIM>
void ContactRegionTaggingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
}

template <unsigned DIM>
bool ContactRegionTaggingModifier<DIM>::ComparisonForDistanceMap(std::pair<unsigned, double> idxDistPairA,
                                                                  std::pair<unsigned, double> idxDistPairB)
{
    return idxDistPairA.second < idxDistPairB.second;
}

template <unsigned DIM>
void ContactRegionTaggingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ContactRegionTaggingModifier<1>;
template class ContactRegionTaggingModifier<2>;
template class ContactRegionTaggingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ContactRegionTaggingModifier)
