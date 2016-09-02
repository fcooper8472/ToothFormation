/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "ThreeRegionInteractionForces.hpp"
#include "ImmersedBoundaryElement.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ThreeRegionInteractionForces<DIM>::ThreeRegionInteractionForces()
        : AbstractImmersedBoundaryForce<DIM>(),
          mpMesh(NULL),
          mBasicInteractionStrength(1e3),
          mBasicInteractionDist(DOUBLE_UNSET)
{
}

template<unsigned DIM>
ThreeRegionInteractionForces<DIM>::~ThreeRegionInteractionForces()
{
}

template<unsigned DIM>
void ThreeRegionInteractionForces<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    /*
     * This force class calculates the force between pairs of nodes in different immersed boundaries.  Each node must
     * therefore store a dimensionless parameter representing the quantity of different transmembrane proteins at that
     * location.  We attach these quantities as node attributes, and keep track of where in the node attributes vector
     * each protein concentration is stored.
     */

    // This will be triggered only once - during simulation set up
    if (mpMesh == NULL)
    {
        mpMesh = &(rCellPopulation.rGetMesh());

        mBasicInteractionDist = 0.25 * rCellPopulation.GetInteractionDistance();

        // First verify that all nodes have the same number of attributes
        unsigned num_node_attributes = rCellPopulation.GetNode(0)->GetNumNodeAttributes();
        for (unsigned node_idx = 0; node_idx < rCellPopulation.GetNumNodes(); node_idx++ )
        {
            if (num_node_attributes != rCellPopulation.GetNode(node_idx)->GetNumNodeAttributes())
            {
                EXCEPTION("All nodes must have the same number of attributes to use this force class.");
            }
        }
    }

    MARK;

//    // Helper variables for loop
//    double normed_dist;
//    double protein_mult;
//
//    // The spring constant will be scaled by an amount determined by the intrinsic spacing
//    double intrinsic_spacing = rCellPopulation.GetIntrinsicSpacing();
//    double node_a_elem_spacing;
//    double node_b_elem_spacing;
//    double elem_spacing;
//
//    // The effective spring constant will be a scaled version of mBasicInteractionStrength
//    double effective_spring_const;
//
//    c_vector<double, DIM> vector_between_nodes;
//    c_vector<double, DIM> force_a_to_b;
//    c_vector<double, DIM> force_b_to_a;
//
//    Node<DIM>* p_node_a;
//    Node<DIM>* p_node_b;
//
//    // If using Morse potential, this can be pre-calculated
//    double well_width = 0.25 * rCellPopulation.GetInteractionDistance();
//
//    // Loop over all pairs of nodes that might be interacting
//    for (unsigned pair = 0; pair < rNodePairs.size(); pair++)
//    {
//        /*
//         * Interactions only occur between different elements.  Since each node is only ever in a single cell, we can
//         * test equality of the first ContainingElement.
//         */
//        if ( *(rNodePairs[pair].first->ContainingElementsBegin()) !=
//             *(rNodePairs[pair].second->ContainingElementsBegin()) )
//        {
//            p_node_a = rNodePairs[pair].first;
//            p_node_b = rNodePairs[pair].second;
//
//            std::vector<double>& r_a_attribs = p_node_a->rGetNodeAttributes();
//            std::vector<double>& r_b_attribs = p_node_b->rGetNodeAttributes();
//
//            vector_between_nodes = mpMesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
//            normed_dist = norm_2(vector_between_nodes);
//
//            if (normed_dist < rCellPopulation.GetInteractionDistance())
//            {
//                // Get the element spacing for each of the nodes concerned and calculate the effective spring constant
//                node_a_elem_spacing = mpMesh->GetAverageNodeSpacingOfElement(*(p_node_a->rGetContainingElementIndices().begin()), false);
//                node_b_elem_spacing = mpMesh->GetAverageNodeSpacingOfElement(*(p_node_b->rGetContainingElementIndices().begin()), false);
//                elem_spacing = 0.5 * (node_a_elem_spacing + node_b_elem_spacing);
//
//                effective_spring_const = mBasicInteractionStrength * elem_spacing / intrinsic_spacing;
//
//                // The protein multiplier is a function of the levels of each protein in the current and comparison nodes
//                protein_mult = std::min(r_a_attribs[e_cad_idx], r_b_attribs[e_cad_idx]) +
//                               std::min(r_a_attribs[p_cad_idx], r_b_attribs[p_cad_idx]) +
//                               std::max(r_a_attribs[integrin_idx], r_b_attribs[integrin_idx]);
//
//                /*
//                 * We must scale each applied force by a factor of elem_spacing / local spacing, so that forces
//                 * balance when spread to the grid later (where the multiplicative factor is the local spacing)
//                 */
//                double morse_exp = exp((mBasicInteractionDist - normed_dist) / well_width);
//                vector_between_nodes *= 2.0 * well_width * effective_spring_const * protein_mult * morse_exp *
//                                        (1.0 - morse_exp) / normed_dist;
//
//                force_a_to_b = vector_between_nodes * elem_spacing / node_a_elem_spacing;
//                p_node_a->AddAppliedForceContribution(force_a_to_b);
//
//                force_b_to_a = -1.0 * vector_between_nodes * elem_spacing / node_b_elem_spacing;
//                p_node_b->AddAppliedForceContribution(force_b_to_a);
//            }
//        }
//    }
}

template<unsigned DIM>
void ThreeRegionInteractionForces<DIM>::SetBasicInteractionStrength(double basicInteractionStrength)
{
    mBasicInteractionStrength = basicInteractionStrength;
}

template<unsigned DIM>
double ThreeRegionInteractionForces<DIM>::GetBasicInteractionStrength()
{
    return mBasicInteractionStrength;
}

template<unsigned DIM>
void ThreeRegionInteractionForces<DIM>::SetBasicInteractionDist(double basicInteractionDist)
{
    mBasicInteractionDist = basicInteractionDist;
}

template<unsigned DIM>
double ThreeRegionInteractionForces<DIM>::GetBasicInteractionDist()
{
    return mBasicInteractionDist;
}

template<unsigned DIM>
void ThreeRegionInteractionForces<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<BasicInteractionStrength>" << mBasicInteractionStrength << "</BasicInteractionStrength>\n";
    *rParamsFile << "\t\t\t<BasicInteractionDist>" << mBasicInteractionDist << "</BasicInteractionDist>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

// Explicit instantiation
template class ThreeRegionInteractionForces<1>;
template class ThreeRegionInteractionForces<2>;
template class ThreeRegionInteractionForces<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ThreeRegionInteractionForces)
