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
#include "UblasCustomFunctions.hpp"

template <unsigned DIM>
ContactRegionTaggingModifier<DIM>::ContactRegionTaggingModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}


template <unsigned DIM>
void ContactRegionTaggingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    auto p_cell_pop = static_cast<ImmersedBoundaryCellPopulation<DIM>*>(&rCellPopulation);
    auto p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>*>(&(rCellPopulation.rGetMesh()));

    /*
     * Annoyingly, we need a deep copy of the nodes for the private box collection because the neighbours are added
     * directly to the nodes, and we cannot overwrite the data used for mechanical interactions in the system.
     *
     * The neighbours of these copied nodes will be used when determining node regions, which requires a longer
     * length-scale than mechanical interactions.
     */
    std::vector<Node<DIM>*> copy_of_nodes;
    copy_of_nodes.reserve(p_mesh->GetNumNodes());
    for (const auto& p_node : p_mesh->rGetNodes())
    {
        copy_of_nodes.push_back(new Node<DIM>(p_node->GetIndex(), p_node->rGetLocation(), true));
    }

    // Verify (for now) that node indices are contiguous and ordered 0, ..., N-1
    std::sort(copy_of_nodes.begin(), copy_of_nodes.end(),
              [](Node<DIM>* a, Node<DIM>* b) -> bool {return a->GetIndex() < b->GetIndex();});

    EXCEPT_IF_NOT(copy_of_nodes.size() == copy_of_nodes.back()->GetIndex() + 1u);

    // Update the box collection
    std::vector<std::pair<Node<DIM>*, Node<DIM>*> > unused_node_pairs;
    mpBoxCollection->CalculateNodePairs(copy_of_nodes, unused_node_pairs);

    // Iterate over elements
    for (auto elem_it = p_mesh->GetElementIteratorBegin(); elem_it != p_mesh->GetElementIteratorEnd(); ++elem_it)
    {
        // Wipe the corners
        for (auto& p_node : elem_it->rGetCornerNodes())
        {
            p_node = nullptr;
        }

        const c_vector<double, DIM> centroid = p_mesh->GetCentroidOfElement(elem_it->GetIndex());
        const unsigned num_nodes_elem = elem_it->GetNumNodes();
        const unsigned max_basal_nodes = std::lround(num_nodes_elem / 6.0);

        // Orient short axis to point in the positive x direction
        c_vector<double, DIM> short_axis = p_mesh->GetShortAxisOfElement(elem_it->GetIndex());
        short_axis = short_axis[0] > 0.0 ? short_axis : -short_axis;
        const auto long_axis = Create_c_vector(-short_axis[1], short_axis[0]);

        // Lambda to compare two elements by distance along the long axis
        auto compare_long_axis = [&](const unsigned& a, const unsigned& b) ->bool
        {
            const auto& to_a = p_mesh->GetVectorFromAtoB(centroid, elem_it->GetNode(a)->rGetLocation());
            const auto& to_b = p_mesh->GetVectorFromAtoB(centroid, elem_it->GetNode(b)->rGetLocation());
            return inner_prod(to_a, long_axis) < inner_prod(to_b, long_axis);
        };

        // Lambda to compare two elements by distance along the short axis
        auto compare_short_axis = [&](const unsigned& a, const unsigned& b) ->bool
        {
            const auto& to_a = p_mesh->GetVectorFromAtoB(centroid, elem_it->GetNode(a)->rGetLocation());
            const auto& to_b = p_mesh->GetVectorFromAtoB(centroid, elem_it->GetNode(b)->rGetLocation());
            return inner_prod(to_a, short_axis) < inner_prod(to_b, short_axis);
        };

        // Create a vector filled with {0, ..., num_nodes_elem-1}, that can be used for partial sorting
        std::vector<unsigned> partial_sort(num_nodes_elem);
        std::iota(partial_sort.begin(), partial_sort.end(), 0u);

        // Partially sort so that first indices are those of nodes lowest down the long axis
        std::nth_element(partial_sort.begin(), partial_sort.begin() + max_basal_nodes,
                         partial_sort.end(), compare_long_axis);

        // Up to max_basal_nodes are all the nodes that are allowed to be basal (but they are not necessarily so)
        std::vector<unsigned> possible_basal_indices;
        for (unsigned i = 0; i < max_basal_nodes; ++i)
        {
            possible_basal_indices.emplace_back(partial_sort[i]);
        }

        // This will be filled in during the loop
        std::vector<unsigned> actual_basal_indices;

        // First loop over every node.  Set initially as all apical, then tag correct nodes as basal
        for (unsigned node_idx = 0; node_idx < num_nodes_elem; ++node_idx)
        {
            // Get the current node, and determine whether it's on the left or right of the long axis
            Node<DIM>* const p_this_node = p_mesh->GetNode(elem_it->GetNodeGlobalIndex(node_idx));
            const bool node_on_left = inner_prod(short_axis, p_mesh->GetVectorFromAtoB(centroid, p_this_node->rGetLocation())) < 0.0;

            p_this_node->SetRegion(node_on_left ? LEFT_APICAL_REGION : RIGHT_APICAL_REGION);

            // If this node might be basal, we find out
            if(std::any_of(possible_basal_indices.begin(), possible_basal_indices.end(), [&](const unsigned& i){return node_idx == i;}))
            {
                // The vec of node neighbours includes all within neighbouring boxes: need to check against neighbour dist
                // This is the mechanical interaction distance, so is the neighbours of the real nodes, not our copies
                const std::vector<unsigned>& node_neighbours = elem_it->GetNode(node_idx)->rGetNeighbours();

                for(const auto& global_idx : node_neighbours)
                {
                    Node<DIM> *const p_other_node = p_mesh->GetNode(global_idx);

                    // If the other node is in a lamina, and within the threshold distance, mark our node as basal
                    if (p_other_node->GetRegion() == LAMINA_REGION &&
                        p_mesh->GetDistanceBetweenNodes(p_this_node->GetIndex(), p_other_node->GetIndex()) < p_cell_pop->GetInteractionDistance())
                    {
                        p_this_node->SetRegion(node_on_left ? LEFT_BASAL_REGION : RIGHT_BASAL_REGION);
                        actual_basal_indices.emplace_back(node_idx);
                        break;
                    }
                }
            }
        }

        // If no basal nodes were tagged, put the lowest node along the long axis in the vector
        if (actual_basal_indices.empty())
        {
            std::nth_element(partial_sort.begin(), partial_sort.begin(), partial_sort.end(), compare_long_axis);
            actual_basal_indices.emplace_back(*partial_sort.begin());
        }

        auto left_right_basal = std::minmax_element(actual_basal_indices.begin(), actual_basal_indices.end(), compare_short_axis);
        const unsigned furthest_left_idx = *left_right_basal.first;
        const unsigned furthest_right_idx = *left_right_basal.second;

        if (furthest_left_idx != furthest_right_idx)
        {
            elem_it->rGetCornerNodes()[LEFT_BASAL_CORNER] = elem_it->GetNode(furthest_left_idx);
            elem_it->rGetCornerNodes()[RIGHT_BASAL_CORNER] = elem_it->GetNode(furthest_right_idx);
        }

        /*
         * We now go counterclockwise from the furthest right, to identify the right lateral / periapical nodes
         */

        // If we investigate too many nodes, we have gone too far; this must mean there are no 'free' apical nodes
        const unsigned gone_too_far = std::lround(0.75 * num_nodes_elem);

        // Define apical surface as a number of consecutive nodes that are not neighbours of nodes in other boundaries
        unsigned this_idx = furthest_right_idx;
        unsigned num_right_lat = 0;
        unsigned num_consecutive_misses = 0u;
        while (num_consecutive_misses < 5u)
        {
            this_idx = (this_idx + 1) % num_nodes_elem;
            num_consecutive_misses++;
            num_right_lat++;

            // Avoid problems near the basal corners by always accepting the fist few nodes
            if (num_right_lat < std::lround(0.25 * num_nodes_elem))
            {
                num_consecutive_misses = 0;
                continue;
            }

            // If we've gone too far, there's no apical surface, and we deal with the case separately
            if (num_right_lat == gone_too_far)
            {
                num_consecutive_misses = UNSIGNED_UNSET;
                break;
            }

            // If we get to here, we actually have to check the node neighbours!
            Node<DIM>* p_this_node = elem_it->GetNode(this_idx);

            // We need the larger box-collection information, as these are not mechanical interactions.
            const std::vector<unsigned>& node_neighbours = copy_of_nodes[p_this_node->GetIndex()]->rGetNeighbours();

            for(const unsigned& gbl_idx : node_neighbours)
            {
                Node<DIM>* p_other_node = p_mesh->GetNode(gbl_idx);

                if (p_mesh->NodesInDifferentElementOrLamina(p_this_node, p_other_node) &&
                    p_mesh->GetDistanceBetweenNodes(p_this_node->GetIndex(), p_other_node->GetIndex()) < mInteractionDist)
                {
                        num_consecutive_misses = 0;
                        break;
                }
            }
        }

        // If we didn't go too far, tag the nodes
        if (num_consecutive_misses != UNSIGNED_UNSET)
        {
            const unsigned num_lateral = std::lround(0.70 * num_right_lat);
            const unsigned num_pa = std::lround(0.8 * num_right_lat);
            for (unsigned i = 0; i < num_right_lat; ++i)
            {
                const unsigned local_idx = (furthest_right_idx + 1 + i) % num_nodes_elem;
                elem_it->GetNode(local_idx)->SetRegion(i <= num_lateral ? RIGHT_LATERAL_REGION : i <= num_pa ? RIGHT_PERIAPICAL_REGION : RIGHT_APICAL_REGION);
            }

            const unsigned right_apical_corner_idx = (furthest_right_idx + 1 + num_pa) % num_nodes_elem;
            elem_it->rGetCornerNodes()[RIGHT_APICAL_CORNER] = elem_it->GetNode(right_apical_corner_idx);
        }

        /*
         * We now go clockwise from the furthest left basal node, to identify the left lateral / periapical nodes
         */

        this_idx = furthest_left_idx;
        unsigned num_left_lat = 0;
        num_consecutive_misses = 0;
        while (num_consecutive_misses < 5u)
        {
            this_idx = (this_idx + num_nodes_elem - 1) % num_nodes_elem;
            num_consecutive_misses++;
            num_left_lat++;

            // Avoid problems near the basal corners by always accepting the fist few nodes
            if (num_left_lat < std::lround(0.25 * num_nodes_elem))
            {
                num_consecutive_misses = 0;
                continue;
            }

            // If we've gone too far, there's no apical surface, and we deal with the case separately
            if (num_left_lat == gone_too_far)
            {
                num_consecutive_misses = UNSIGNED_UNSET;
                break;
            }

            // If we get to here, we actually have to check the node neighbours!
            Node<DIM>* p_this_node = elem_it->GetNode(this_idx);

            // The vec of node neighbours includes all within neighbouring boxes: need to check against neighbour dist
            const std::vector<unsigned>& node_neighbours = copy_of_nodes[p_this_node->GetIndex()]->rGetNeighbours();

            for(const unsigned& gbl_idx : node_neighbours)
            {
                Node<DIM>* p_other_node = p_mesh->GetNode(gbl_idx);

                if (p_mesh->NodesInDifferentElementOrLamina(p_this_node, p_other_node) &&
                    p_mesh->GetDistanceBetweenNodes(p_this_node->GetIndex(), p_other_node->GetIndex()) < mInteractionDist)
                {
                    num_consecutive_misses = 0;
                    break;
                }
            }
        }

        // If we didn't go too far, tag the nodes
        if (num_consecutive_misses != UNSIGNED_UNSET)
        {
            const unsigned num_lateral = std::lround(0.7 * num_left_lat);
            const unsigned num_pa = std::lround(0.8 * num_right_lat);
            for (unsigned i = 0; i < num_left_lat; ++i)
            {
                unsigned local_idx = (furthest_left_idx + num_nodes_elem - i - 1) % num_nodes_elem;
                elem_it->GetNode(local_idx)->SetRegion(i <= num_lateral ? LEFT_LATERAL_REGION : i <= num_pa ? LEFT_PERIAPICAL_REGION : LEFT_APICAL_REGION);
            }

            const unsigned left_apical_corner_idx = (furthest_left_idx + num_nodes_elem - 1 - num_pa) % num_nodes_elem;
            elem_it->rGetCornerNodes()[LEFT_APICAL_CORNER] = elem_it->GetNode(left_apical_corner_idx);
        }
        // If we did go too far, tag all non-basal nodes as left or right based on their orientation about the long axis
        else
        {
            c_vector<double, DIM> axis;

            // If there were no basal nodes tagged, we just use the short axis.
            if (furthest_left_idx == furthest_right_idx)
            {
                axis = short_axis;
            }
            // Else, generate an axis based on the basal nodes (more robust to odd cell shapes).
            else
            {
                const c_vector<double, DIM> furthest_left = elem_it->GetNode(furthest_left_idx)->rGetLocation();
                const c_vector<double, DIM> furthest_right = elem_it->GetNode(furthest_right_idx)->rGetLocation();
                axis = p_mesh->GetVectorFromAtoB(furthest_left, furthest_right);
            }

            for (unsigned i = (furthest_right_idx + 1) % num_nodes_elem; i != furthest_left_idx; i = (i + 1) % num_nodes_elem)
            {
                Node<DIM>* p_this_node = elem_it->GetNode(i);
                bool node_on_left = inner_prod(axis, p_mesh->GetVectorFromAtoB(centroid, p_this_node->rGetLocation())) < 0.0;
                p_this_node->SetRegion(node_on_left ? LEFT_LATERAL_REGION : RIGHT_LATERAL_REGION);
            }
        }
    }

    // Tidy up out copied nodes
    for (const auto& p_node : copy_of_nodes)
    {
        delete p_node;
    }
}

template <unsigned DIM>
void ContactRegionTaggingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    // Fist, decide how big the interaction distance must be.  Make an educated guess based on number of cells.
    const double rough_cell_width = 1.0 / rCellPopulation.GetNumAllCells();
    mInteractionDist = 0.2 * rough_cell_width;  // this is a bit of a guess

    // Set up the box collection
    c_vector<double, 2 * 2> domain_size;
    domain_size(0) = 0.0;
    domain_size(1) = 1.0;
    domain_size(2) = 0.0;
    domain_size(3) = 1.0;

    mpBoxCollection = new ObsoleteBoxCollection<DIM>(mInteractionDist, domain_size, true, true);
    mpBoxCollection->SetupLocalBoxesHalfOnly();
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
