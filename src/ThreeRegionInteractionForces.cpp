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

#include "ThreeRegionInteractionForces.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

template <unsigned DIM>
ThreeRegionInteractionForces<DIM>::ThreeRegionInteractionForces()
        : AbstractImmersedBoundaryForce<DIM>(),
          mpMesh(NULL),
          mElemAttributeLocation(UNSIGNED_UNSET),
          mBasicInteractionStrength(1e3),
          mBasicInteractionDist(DOUBLE_UNSET),
          mAdhesionMultiplier(2.0),
          mRegionSizes()
{
}

template<unsigned DIM>
ThreeRegionInteractionForces<DIM>::~ThreeRegionInteractionForces()
{
}

template <unsigned DIM>
void ThreeRegionInteractionForces<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM> *, Node<DIM> *> > &rNodePairs,
                                                                             ImmersedBoundaryCellPopulation<DIM> &rCellPopulation)
{
    /*
     * This force class divides a palisade of cells into three regions: left, middle, and right.
     * In addition, one cell on the wrap-around behaves like those in the middle.
     *
     * The left and right cell properties are mirrored, and the middle cells are 'neutral'.
     *
     * It is assumed that there are (at least initially) 3n+1 cells.
     */

    // This will be triggered only once - during simulation set up
    if (mpMesh == NULL)
    {
        mpMesh = &(rCellPopulation.rGetMesh());

        if ((mpMesh->GetNumElements()) % 3 != 0 || mpMesh->GetNumElements() < 3)
        {
            EXCEPTION("This class assumes 3n cells arranged in a palisade.");
        }

        // First verify that all elements have the same number of attributes
        unsigned num_elem_attributes = rCellPopulation.GetElement(0)->GetNumElementAttributes();
        for (unsigned elem_idx = 0; elem_idx < rCellPopulation.GetNumElements(); elem_idx++)
        {
            if (num_elem_attributes != rCellPopulation.GetElement(elem_idx)->GetNumElementAttributes())
            {
                EXCEPTION("All elements must have the same number of attributes to use this force class.");
            }
        }

        // Define the location in the element attributes vector
        mElemAttributeLocation = num_elem_attributes;

        if (std::accumulate(mRegionSizes.begin(), mRegionSizes.end(), 0u) != mpMesh->GetNumElements())
        {
            EXCEPTION("mRegionSizes must be set to the correct number of elements");
        }

        std::array<unsigned, 3> region_size_sums;
        std::partial_sum(mRegionSizes.begin(), mRegionSizes.end(), region_size_sums.begin());

        for (unsigned elem_idx = 0; elem_idx < region_size_sums[0]; elem_idx++)
        {
            mpMesh->GetElement(elem_idx)->AddElementAttribute(THREE_REGION_LEFT);
        }
        for (unsigned elem_idx = region_size_sums[0]; elem_idx < region_size_sums[1]; elem_idx++)
        {
            mpMesh->GetElement(elem_idx)->AddElementAttribute(THREE_REGION_MID);
        }
        for (unsigned elem_idx = region_size_sums[1]; elem_idx < region_size_sums[2]; elem_idx++)
        {
            mpMesh->GetElement(elem_idx)->AddElementAttribute(THREE_REGION_RIGHT);
        }

        mBasicInteractionDist = 0.25 * rCellPopulation.GetInteractionDistance();
    }

    // Helper variables for loop
    double normed_dist;
    double elem_type_mult;

    // The spring constant will be scaled by an amount determined by the intrinsic spacing
    double intrinsic_spacing = rCellPopulation.GetIntrinsicSpacing();
    double node_a_elem_spacing;
    double node_b_elem_spacing;
    double elem_spacing;

    // The effective spring constant will be a scaled version of mBasicInteractionStrength
    double effective_spring_const;

    c_vector<double, DIM> vector_between_nodes;
    c_vector<double, DIM> force_a_to_b;
    c_vector<double, DIM> force_b_to_a;

    Node<DIM> *p_node_a;
    Node<DIM> *p_node_b;

    // If using Morse potential, this can be pre-calculated
    double well_width = 0.25 * rCellPopulation.GetInteractionDistance();

    // Loop over all pairs of nodes that might be interacting
    for (unsigned pair = 0; pair < rNodePairs.size(); pair++)
    {
        /*
        * Interactions only occur between different elements.  Since each node is only ever in a single cell, we can
        * test equality of the first ContainingElement.
        */
        if (mpMesh->NodesInDifferentElementOrLamina(rNodePairs[pair].first, rNodePairs[pair].second))
        {
            p_node_a = rNodePairs[pair].first;
            p_node_b = rNodePairs[pair].second;

            vector_between_nodes = mpMesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
            normed_dist = norm_2(vector_between_nodes);

            if (normed_dist < rCellPopulation.GetInteractionDistance())
            {
                // Get the element spacing for each of the nodes concerned and calculate the effective spring constant
                if (p_node_a->GetRegion() != LAMINA_REGION)
                {
                    node_a_elem_spacing = mpMesh->GetAverageNodeSpacingOfElement(*(p_node_a->rGetContainingElementIndices().begin()), false);
                }
                else
                {
                    node_a_elem_spacing = mpMesh->GetAverageNodeSpacingOfLamina(0, false);
                }
                if (p_node_b->GetRegion() != LAMINA_REGION)
                {
                    node_b_elem_spacing = mpMesh->GetAverageNodeSpacingOfElement(*(p_node_b->rGetContainingElementIndices().begin()), false);
                }
                else
                {
                    node_b_elem_spacing = mpMesh->GetAverageNodeSpacingOfLamina(0, false);
                }

                elem_spacing = 0.5 * (node_a_elem_spacing + node_b_elem_spacing);

                effective_spring_const = mBasicInteractionStrength * elem_spacing / intrinsic_spacing;

                // Here we calculate the adhesion variation due to node regions
                elem_type_mult = CalculateElementTypeMult(p_node_a, p_node_b);

                /*
                * We must scale each applied force by a factor of elem_spacing / local spacing, so that forces
                * balance when spread to the grid later (where the multiplicative factor is the local spacing)
                */
                double morse_exp = exp((mBasicInteractionDist - normed_dist) / well_width);
                vector_between_nodes *= 2.0 * well_width * effective_spring_const * elem_type_mult * morse_exp *
                                        (1.0 - morse_exp) / normed_dist;

                force_a_to_b = vector_between_nodes * elem_spacing / node_a_elem_spacing;
                p_node_a->AddAppliedForceContribution(force_a_to_b);

                force_b_to_a = -1.0 * vector_between_nodes * elem_spacing / node_b_elem_spacing;
                p_node_b->AddAppliedForceContribution(force_b_to_a);
            }
        }
    }

    if (this->mAdditiveNormalNoise)
    {
        this->AddNormalNoiseToNodes(rCellPopulation);
    }
}

template<unsigned DIM>
unsigned ThreeRegionInteractionForces<DIM>::GetElemType(unsigned elemIdx)
{
    unsigned elem_type = UINT_MAX;

    if (mElemAttributeLocation != UNSIGNED_UNSET)
    {
        elem_type = mpMesh->GetElement(elemIdx)->rGetElementAttributes()[mElemAttributeLocation];
    }

    return elem_type;
}

template<unsigned DIM>
unsigned ThreeRegionInteractionForces<DIM>::GetElemType(Node<DIM>* pNode)
{
    unsigned elem_type = UINT_MAX;

    if (mElemAttributeLocation != UNSIGNED_UNSET)
    {
        std::set<unsigned>& r_containing = pNode->rGetContainingElementIndices();

        if (r_containing.size() == 0)
        {
            elem_type = THREE_REGION_LAMINA;
        }
        else
        {
            unsigned elem_idx = *(r_containing.begin());
            elem_type = (unsigned)mpMesh->GetElement(elem_idx)->rGetElementAttributes()[mElemAttributeLocation];
        }
    }

    return elem_type;
}

template<unsigned DIM>
double ThreeRegionInteractionForces<DIM>::CalculateElementTypeMult(Node<DIM>* pNodeA, Node<DIM>* pNodeB)
{
    double type_mult = 1.0;

    unsigned elem_type_a = GetElemType(pNodeA);
    unsigned elem_type_b = GetElemType(pNodeB);

    unsigned region_a = pNodeA->GetRegion();
    unsigned region_b = pNodeB->GetRegion();

    if (region_a == RIGHT_APICAL_REGION || region_a == RIGHT_PERIAPICAL_REGION ||
        region_a == LEFT_APICAL_REGION || region_a == LEFT_PERIAPICAL_REGION)
    {
        type_mult = mAdhesionMultiplier;
    }

    if (region_b == RIGHT_APICAL_REGION || region_b == RIGHT_PERIAPICAL_REGION ||
        region_b == LEFT_APICAL_REGION || region_b == LEFT_PERIAPICAL_REGION)
    {
        type_mult = mAdhesionMultiplier;
    }
//    if (elem_type_a == THREE_REGION_MID)
//    {
//        if (region_a == LEFT_BASAL_REGION || region_a == RIGHT_BASAL_REGION)
//        {
//            type_mult = 2.0;
//        }
//    }
//
//    if (elem_type_b == THREE_REGION_MID)
//    {
//        if (region_b == LEFT_BASAL_REGION || region_b == RIGHT_BASAL_REGION)
//        {
//            type_mult = 2.0;
//        }
//    }

    return type_mult;
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
void ThreeRegionInteractionForces<DIM>::SetAdhesionMultiplier(double adhesionMultiplier)
{
    mAdhesionMultiplier = adhesionMultiplier;
}

template<unsigned DIM>
double ThreeRegionInteractionForces<DIM>::GetAdhesionMultiplier()
{
    return mAdhesionMultiplier;
}

template<unsigned DIM>
void ThreeRegionInteractionForces<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<BasicInteractionStrength>" << mBasicInteractionStrength << "</BasicInteractionStrength>\n";
    *rParamsFile << "\t\t\t<BasicInteractionDist>" << mBasicInteractionDist << "</BasicInteractionDist>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template<unsigned DIM>
const std::array<unsigned int, 3>& ThreeRegionInteractionForces<DIM>::GetRegionSizes() const
{
    return mRegionSizes;
}

template<unsigned DIM>
void ThreeRegionInteractionForces<DIM>::SetRegionSizes(const std::array<unsigned int, 3>& regionSizes)
{
    mRegionSizes = regionSizes;
}

// Explicit instantiation
template class ThreeRegionInteractionForces<1>;
template class ThreeRegionInteractionForces<2>;
template class ThreeRegionInteractionForces<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ThreeRegionInteractionForces)
