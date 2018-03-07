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

#include "VarAdhesionMorseMembraneForce.hpp"
#include "ImmersedBoundaryEnumerations.hpp"
#include "UblasVectorInclude.hpp"

template <unsigned DIM>
VarAdhesionMorseMembraneForce<DIM>::VarAdhesionMorseMembraneForce()
        : AbstractImmersedBoundaryForce<DIM>(),
          mElementWellDepth(1e6),
          mElementRestLength(0.5),
          mLaminaWellDepth(1e6),
          mLaminaRestLength(0.5),
          mApicalWellDepthMult(1.0),
          mWellWidth(0.25),
          mStiffnessMult(1.0)
{
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                                                              ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{

    auto p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM> *>(&(rCellPopulation.rGetMesh()));
    if (std::accumulate(mRegionSizes.begin(), mRegionSizes.end(), 0u) != p_mesh->GetNumElements())
    {
        EXCEPTION("mRegionSizes must be set to the correct number of elements");
    }

    // Data common across the entire cell population
    double intrinsicSpacingSquared = rCellPopulation.GetIntrinsicSpacing() * rCellPopulation.GetIntrinsicSpacing();

    // Loop over all elements ( <DIM, DIM> )
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = rCellPopulation.rGetMesh().GetElementIteratorBegin();
         elem_it != rCellPopulation.rGetMesh().GetElementIteratorEnd();
         ++elem_it)
    {
        CalculateForcesOnElement(*elem_it, rCellPopulation, intrinsicSpacingSquared);
    }

    // Loop over all laminas ( <DIM-1, DIM> )
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryLaminaIterator lam_it = rCellPopulation.rGetMesh().GetLaminaIteratorBegin();
         lam_it != rCellPopulation.rGetMesh().GetLaminaIteratorEnd();
         ++lam_it)
    {
        CalculateForcesOnElement(*lam_it, rCellPopulation, intrinsicSpacingSquared);
    }
}

#include "Debug.hpp"
template <unsigned DIM>
template <unsigned ELEMENT_DIM>
void VarAdhesionMorseMembraneForce<DIM>::CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, DIM>& rElement,
                                                                  ImmersedBoundaryCellPopulation<DIM>& rCellPopulation,
                                                                  double intrinsicSpacingSquared)
{
    auto& r_mesh = rCellPopulation.rGetMesh();

    // Get index and number of nodes of current element
    const unsigned elem_idx = rElement.GetIndex();
    const unsigned num_nodes = rElement.GetNumNodes();

    unsigned elem_region = 2u;
    if (rElement.GetIndex() < mRegionSizes[0])
    {
        elem_region = 0;
    }
    else if (rElement.GetIndex() < mRegionSizes[0] + mRegionSizes[1])
    {
        elem_region = 1;
    }

    if (ELEMENT_DIM != DIM)
    {
        elem_region = 3;
    }

    // Make a vector to store the force on node i+1 from node i
    std::vector<c_vector<double, DIM> > force_to_next(num_nodes);

    /*
     * Get the node spacing ratio for this element.  The rest length and spring constant are derived from this
     * characteristic length.
     *
     * The spring constant is derived with reference to the intrinsic spacing, so that with different node spacings
     * the user-defined parameters do not have to be updated.
     *
     * The correct factor to increase the spring constant by is (intrinsic spacing / node_spacing)^2.  One factor
     * takes into account the energy considerations of the elastic springs, and the other takes account of the
     * factor of node_spacing used in discretising the force relation.
     */

    double node_spacing = 0.0;
    double well_depth = 0.0;
    double rest_length = 0.0;
    double well_width = mWellWidth;

    // Determine if we're in a lamina or not
    if (ELEMENT_DIM < DIM)
    {
        node_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfLamina(elem_idx, false);

        well_depth = mLaminaWellDepth * intrinsicSpacingSquared / (node_spacing * node_spacing);
        rest_length = mLaminaRestLength * node_spacing;
        well_width *= node_spacing;

        // Alter the stiffness of the apical lamina, if required
        if (rElement.GetIndex() == 1)
        {
            well_depth *= mApicalWellDepthMult;
        }
    }
    else // regular element
    {
        node_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(elem_idx, false);

        well_depth = mElementWellDepth * intrinsicSpacingSquared / (node_spacing * node_spacing);
        rest_length = mElementRestLength * node_spacing;
        well_width *= node_spacing;
    }

    // Calculate quantity of gradient component
    auto elem_centroid = r_mesh.GetCentroidOfElement(elem_idx);
    auto mid_centroid = r_mesh.GetCentroidOfElement((rCellPopulation.rGetMesh().GetNumElements() - 1u) / 2);
    const double dist_from_centre = std::fabs(r_mesh.GetVectorFromAtoB(elem_centroid, mid_centroid)[0]);

    const double end_of_gradient = 1.0 * mStiffnessMult;
    const double mid_of_gradient = 0.5 * (1.0 + mStiffnessMult);

    const double gradient_weight = mid_of_gradient + 2.0 * (end_of_gradient - mid_of_gradient) * dist_from_centre;


    // Loop over nodes and calculate the force exerted on node i+1 by node i
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Index of the next node, calculated modulo number of nodes in this element
        unsigned next_idx = (node_idx + 1) % num_nodes;

        // If element rather than lamina, calculate the stiffness multiplier
        double stiffness_mult = CalculateStiffnessMult(elem_region, rElement.GetNode(node_idx)->GetRegion(), elem_idx, gradient_weight);

        // Morse force (derivative of Morse potential wrt distance between nodes
        force_to_next[node_idx] = rCellPopulation.rGetMesh().GetVectorFromAtoB(rElement.GetNodeLocation(node_idx),
                                                                               rElement.GetNodeLocation(next_idx));
        double normed_dist = norm_2(force_to_next[node_idx]);

        double morse_exp = exp((rest_length - normed_dist) / well_width);
        force_to_next[node_idx] *= 2.0 * well_width * well_depth * stiffness_mult * morse_exp * (1.0 - morse_exp) / normed_dist;
    }

    // Add the contributions of springs adjacent to each node
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Get index of previous node
        unsigned prev_idx = (node_idx + num_nodes - 1) % num_nodes;

        c_vector<double, DIM> aggregate_force = force_to_next[node_idx] - force_to_next[prev_idx];

        // Add the aggregate force contribution to the node
        rElement.GetNode(node_idx)->AddAppliedForceContribution(aggregate_force);
    }

    auto DistanceNodeToNode = [&](Node<DIM>* const a, Node<DIM>* const b) -> double
    {
        const unsigned a_idx = a->GetIndex();
        const unsigned b_idx = b->GetIndex();
        // Get the number of nodes in the shortest path around the element from a to b
        unsigned difference = a_idx > b_idx ? a_idx - b_idx : b_idx - a_idx;
        difference = std::min(difference, num_nodes - difference);

        return node_spacing * difference;
    };

    // Add the corner contributions
    if (DIM == ELEMENT_DIM)
    {
        // Square root of 5 is the diagonal of a 2:1 ratio rectangle, and 1/15 is rough width of a cell
        constexpr double width = 0.90 / 15.0;
        constexpr double diagonal = width * 2.23606797749978969640;
        constexpr double height = width * 1.92;

        const auto &rCorners = rElement.rGetCornerNodes();

        if (std::none_of(rCorners.begin(), rCorners.end(), [](Node<DIM> *a) { return a == nullptr; }))
        {

            Node<DIM>* const p_lt_ap = rCorners[LEFT_APICAL_CORNER];
            Node<DIM>* const p_rt_ap = rCorners[RIGHT_APICAL_CORNER];
            Node<DIM>* const p_lt_ba = rCorners[LEFT_BASAL_CORNER];
            Node<DIM>* const p_rt_ba = rCorners[RIGHT_BASAL_CORNER];

            const double left_height = DistanceNodeToNode(p_lt_ba, p_lt_ap);
            const double right_height = DistanceNodeToNode(p_rt_ba, p_rt_ap);
            const double top_height = 0.9 * DistanceNodeToNode(p_rt_ap, p_lt_ap);
            const double bot_height = 0.9 * DistanceNodeToNode(p_rt_ba, p_lt_ba);


//            // Left and mid region
//            if (elem_region == 0u || elem_region == 1u)
//            {
//                const double this_diagonal = elem_region == 1u ? diagonal : 1.0 * diagonal;
//                const double this_factor = factor == 1u ? factor : 1.0 * factor;
//
//                auto diag_vec = r_mesh.GetVectorFromAtoB(p_lt_ba->rGetLocation(), p_rt_ap->rGetLocation());
//                const double length = norm_2(diag_vec);
//                const auto diag_force = this_factor * (length - this_diagonal) * mElementWellDepth * diag_vec / length;
//                p_rt_ap->AddAppliedForceContribution(-diag_force);
//                p_lt_ba->AddAppliedForceContribution(diag_force);
//            }
//
//            // Mid and right region
//            if (elem_region == 1u || elem_region == 2u)
//            {
//                const double this_diagonal = elem_region == 1u ? diagonal : 1.0 * diagonal;
//                const double this_factor = factor == 1u ? factor : 1.0 * factor;
//
//                auto diag_vec = r_mesh.GetVectorFromAtoB(p_rt_ba->rGetLocation(), p_lt_ap->rGetLocation());
//                const double length = norm_2(diag_vec);
//                const auto diag_force = this_factor * (length - this_diagonal) * mElementWellDepth * diag_vec / length;
//                p_lt_ap->AddAppliedForceContribution(-diag_force);
//                p_rt_ba->AddAppliedForceContribution(diag_force);
//            }

            { // left
                const c_vector<double, DIM> left_vec = r_mesh.GetVectorFromAtoB(p_lt_ba->rGetLocation(), p_lt_ap->rGetLocation());
                const double len_1 = norm_2(left_vec);
                const c_vector<double, DIM> force_1 = (mSupportStrength * (len_1 - left_height) * well_depth / len_1) * left_vec;

                p_lt_ap->AddAppliedForceContribution(-force_1);
                p_lt_ba->AddAppliedForceContribution(force_1);
            }

            { // right
                const c_vector<double, DIM> right_vec = r_mesh.GetVectorFromAtoB(p_rt_ba->rGetLocation(), p_rt_ap->rGetLocation());
                const double len_2 = norm_2(right_vec);
                const c_vector<double, DIM> force_2 = (mSupportStrength * (len_2 - right_height) * well_depth / len_2) * right_vec;

                p_rt_ap->AddAppliedForceContribution(-force_2);
                p_rt_ba->AddAppliedForceContribution(force_2);
            }

            { // top
                const c_vector<double, DIM> top_vec = r_mesh.GetVectorFromAtoB(p_rt_ap->rGetLocation(), p_lt_ap->rGetLocation());
                const double len_3 = norm_2(top_vec);
                const c_vector<double, DIM> force_3 = (mSupportStrength * (len_3 - bot_height) * well_depth / len_3) * top_vec;

                p_lt_ap->AddAppliedForceContribution(-force_3);
                p_rt_ap->AddAppliedForceContribution(force_3);
            }

            { // bottom
                const c_vector<double, DIM> bot_vec = r_mesh.GetVectorFromAtoB(p_rt_ba->rGetLocation(), p_lt_ba->rGetLocation());
                const double len_4 = norm_2(bot_vec);
                const c_vector<double, DIM> force_4 = (mSupportStrength * (len_4 - bot_height) * well_depth / len_4) * bot_vec;

                p_lt_ba->AddAppliedForceContribution(-force_4);
                p_rt_ba->AddAppliedForceContribution(force_4);
            }
        }
    }


}

template<unsigned DIM>
double VarAdhesionMorseMembraneForce<DIM>::CalculateStiffnessMult(
        unsigned elem_region,
        unsigned node_region,
        unsigned elem_idx,
        double gradientWeight)
{
    double stiffness_mult = 1.0;

    // Left
    if (elem_region == 0)
    {
        if (node_region == RIGHT_APICAL_REGION)// || node_region == RIGHT_PERIAPICAL_REGION)
        {
            stiffness_mult = gradientWeight;
        }
    }
    // Right
    else if (elem_region == 2)
    {
        if (node_region == LEFT_APICAL_REGION)// || node_region == LEFT_PERIAPICAL_REGION)
        {
            stiffness_mult = gradientWeight;
        }
    }
    // Centre
    else if (elem_region == 1)
    {
        stiffness_mult = 1.0;
    }

    return stiffness_mult;
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ElementWellDepth>" << mElementWellDepth << "</ElementWellDepth>\n";
    *rParamsFile << "\t\t\t<ElementRestLength>" << mElementRestLength << "</ElementRestLength>\n";
    *rParamsFile << "\t\t\t<LaminaWellDepth>" << mLaminaWellDepth << "</LaminaWellDepth>\n";
    *rParamsFile << "\t\t\t<LaminaRestLength>" << mLaminaRestLength << "</LaminaRestLength>\n";
    *rParamsFile << "\t\t\t<WellWidth>" << mWellWidth << "</WellWidth>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template <unsigned DIM>
double VarAdhesionMorseMembraneForce<DIM>::GetElementWellDepth() const
{
    return mElementWellDepth;
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetElementWellDepth(double elementWellDepth)
{
    mElementWellDepth = elementWellDepth;
}

template <unsigned DIM>
double VarAdhesionMorseMembraneForce<DIM>::GetElementRestLength() const
{
    return mElementRestLength;
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetElementRestLength(double elementRestLength)
{
    mElementRestLength = elementRestLength;
}

template <unsigned DIM>
double VarAdhesionMorseMembraneForce<DIM>::GetLaminaWellDepth() const
{
    return mLaminaWellDepth;
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetLaminaWellDepth(double laminaWellDepth)
{
    mLaminaWellDepth = laminaWellDepth;
}

template <unsigned DIM>
double VarAdhesionMorseMembraneForce<DIM>::GetLaminaRestLength() const
{
    return mLaminaRestLength;
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetLaminaRestLength(double laminaRestLength)
{
    mLaminaRestLength = laminaRestLength;
}

template <unsigned DIM>
double VarAdhesionMorseMembraneForce<DIM>::GetApicalWellDepthMult() const
{
    return mLaminaRestLength;
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetApicalWellDepthMult(double apicalWellDepthMult)
{
    mApicalWellDepthMult = apicalWellDepthMult;
}

template <unsigned DIM>
double VarAdhesionMorseMembraneForce<DIM>::GetWellWidth() const
{
    return mWellWidth;
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetWellWidth(double wellWidth)
{
    mWellWidth = wellWidth;
}

template <unsigned DIM>
double VarAdhesionMorseMembraneForce<DIM>::GetStiffnessMult() const
{
    return mStiffnessMult;
}

template <unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetStiffnessMult(double stiffnessMult)
{
    mStiffnessMult = stiffnessMult;
}

template<unsigned DIM>
const std::array<unsigned int, 3>& VarAdhesionMorseMembraneForce<DIM>::GetRegionSizes() const
{
    return mRegionSizes;
}

template<unsigned DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetRegionSizes(const std::array<unsigned int, 3>& regionSizes)
{
    mRegionSizes = regionSizes;
}

template<unsigned int DIM>
double VarAdhesionMorseMembraneForce<DIM>::GetSupportStrength() const
{
    return mSupportStrength;
}

template<unsigned int DIM>
void VarAdhesionMorseMembraneForce<DIM>::SetSupportStrength(double supportStrength)
{
    mSupportStrength = supportStrength;
}


// Explicit instantiation
template class VarAdhesionMorseMembraneForce<1>;
template class VarAdhesionMorseMembraneForce<2>;
template class VarAdhesionMorseMembraneForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VarAdhesionMorseMembraneForce)
