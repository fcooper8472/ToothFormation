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

#include "ThreeRegionInteractionForces.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

template <unsigned DIM>
ThreeRegionInteractionForces<DIM>::ThreeRegionInteractionForces()
        : AbstractImmersedBoundaryForce<DIM>(),
          mpMesh(NULL),
          mElemAttributeLocation(UNSIGNED_UNSET),
          mBasicInteractionStrength(1e3),
          mBasicInteractionDist(0.25),
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
    }

    // Pre-fetch some values associated with the cell populations
    const double intrinsic_spacing = rCellPopulation.GetIntrinsicSpacing();
    const double eff_rest_length = mBasicInteractionDist * rCellPopulation.GetInteractionDistance();
    const double eff_well_width = 0.25 * rCellPopulation.GetInteractionDistance();

    // Loop over all pairs of nodes that might be interacting
    for (const auto& node_pair : rNodePairs)
    {
        /*
        * Interactions only occur between different elements.  Since each node is only ever in a single cell, we can
        * test equality of the first ContainingElement.
        */
        if (mpMesh->NodesInDifferentElementOrLamina(node_pair.first, node_pair.second))
        {
            Node<DIM>* const p_node_a = node_pair.first;
            Node<DIM>* const p_node_b = node_pair.second;

            c_vector<double, DIM> vec_a2b = mpMesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
            const double normed_dist = norm_2(vec_a2b);

            if (normed_dist < rCellPopulation.GetInteractionDistance())
            {
                // Get the element spacing for each of the nodes concerned and calculate the effective spring constant
                const double node_a_elem_spacing = p_node_a->GetRegion() == LAMINA_REGION ?
                                                   mpMesh->GetAverageNodeSpacingOfLamina(0, false) :
                                                   mpMesh->GetAverageNodeSpacingOfElement(*(p_node_a->rGetContainingElementIndices().begin()), false);

                const double node_b_elem_spacing = p_node_a->GetRegion() == LAMINA_REGION ?
                                                   mpMesh->GetAverageNodeSpacingOfLamina(0, false) :
                                                   mpMesh->GetAverageNodeSpacingOfElement(*(p_node_a->rGetContainingElementIndices().begin()), false);

                const double elem_spacing = 0.5 * (node_a_elem_spacing + node_b_elem_spacing);
                const double eff_well_depth = mBasicInteractionStrength * elem_spacing / intrinsic_spacing;

                // Calculate the adhesion variation due to node regions
                const double elem_type_mult = CalculateElementTypeMult(p_node_a, p_node_b);


                const double morse_exp = std::exp((eff_rest_length - normed_dist) / eff_well_width);

                /*
                 * We must scale each applied force by a factor of elem_spacing / local spacing, so that forces
                 * balance when spread to the grid later (where the multiplicative factor is the local spacing)
                 */
                vec_a2b *= 2.0 * eff_well_width * eff_well_depth * morse_exp * (1.0 - morse_exp) / normed_dist;

                const c_vector<double, DIM> force_a2b = vec_a2b * (elem_spacing / node_a_elem_spacing);
                p_node_a->AddAppliedForceContribution(force_a2b);

                const c_vector<double, DIM> force_b2a = vec_a2b * (-1.0 * elem_spacing / node_b_elem_spacing);
                p_node_b->AddAppliedForceContribution(force_b2a);
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

    // Predicate to query whether a given node is right or left apical
    auto is_apical_node = [](Node<DIM>* const p_node)->bool
    {
        unsigned region = p_node->GetRegion();
        return region == RIGHT_APICAL_REGION || region == LEFT_APICAL_REGION;
    };

    // First, if either belong to the lamina, return
    for (const auto& p_node : {pNodeA, pNodeB})
    {
        unsigned elem_type = GetElemType(p_node);

        // If either of the nodes is lamina region, just set the relevant type multiplier and stop
        if (elem_type == THREE_REGION_LAMINA)
        {
            type_mult = 2.0;
            break;
        }
    }

    // If both are apical, the multiplier is the given member variable
    if (is_apical_node(pNodeA) && is_apical_node(pNodeB))
    {
        type_mult = mAdhesionMultiplier;
    }

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
