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

#include "GradualChangeMorseMembraneForce.hpp"
#include "ImmersedBoundaryEnumerations.hpp"
#include "RandomNumberGenerator.hpp"

template <unsigned DIM>
GradualChangeMorseMembraneForce<DIM>::GradualChangeMorseMembraneForce()
        : AbstractImmersedBoundaryForce<DIM>(),
          mElementWellDepth(1e6),
          mElementRestLength(0.5),
          mLaminaWellDepth(1e6),
          mLaminaRestLength(0.5),
          mWellWidth(0.25),
          mStiffnessMult(1.0)
{
}

template <unsigned DIM>
GradualChangeMorseMembraneForce<DIM>::~GradualChangeMorseMembraneForce()
{
}

template <unsigned DIM>
void GradualChangeMorseMembraneForce<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                                                                ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    // Data common across the entire cell population
    double intrinsicSpacingSquared = rCellPopulation.GetIntrinsicSpacing() * rCellPopulation.GetIntrinsicSpacing();

    // Loop over all elements ( <DIM, DIM> )
    for (auto elem_it = rCellPopulation.rGetMesh().GetElementIteratorBegin();
         elem_it != rCellPopulation.rGetMesh().GetElementIteratorEnd();
         ++elem_it)
    {
        CalculateForcesOnElement(*elem_it, rCellPopulation, intrinsicSpacingSquared);
    }

    // Loop over all laminas ( <DIM-1, DIM> )
    for (auto lam_it = rCellPopulation.rGetMesh().GetLaminaIteratorBegin();
         lam_it != rCellPopulation.rGetMesh().GetLaminaIteratorEnd();
         ++lam_it)
    {
        CalculateForcesOnElement(*lam_it, rCellPopulation, intrinsicSpacingSquared);
    }

    if (this->mAdditiveNormalNoise)
    {
        this->AddNormalNoiseToNodes(rCellPopulation);
    }
}

template <unsigned DIM>
template <unsigned ELEMENT_DIM>
void GradualChangeMorseMembraneForce<DIM>::CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, DIM>& rElement,
                                                                    ImmersedBoundaryCellPopulation<DIM>& rCellPopulation,
                                                                    double intrinsicSpacingSquared)
{
    if (rCellPopulation.GetNumElements() % 2 != 1)
    {
        EXCEPTION("This force class assumes there are an odd number of elements");
    }

    // Get index and number of nodes of current element
    unsigned elem_idx = rElement.GetIndex();
    unsigned num_nodes = rElement.GetNumNodes();
    unsigned num_elems = rCellPopulation.rGetMesh().GetNumElements();
    unsigned middle_idx = (num_elems - 1) / 2;

    // If we need to generate random numbers, create a random number generator
    RandomNumberGenerator* p_gen = this->mAdditiveNormalNoise ? RandomNumberGenerator::Instance() : NULL;
    // Calculate the stiffness mult that applies to specific regions of elements
    double stiffness_mult = 1.0;
    if (ELEMENT_DIM == DIM)
    {
        double proportion_to_edge = fabs(elem_idx - static_cast<double>(middle_idx)) / static_cast<double>(middle_idx);
        stiffness_mult = 1.0 + proportion_to_edge * (mStiffnessMult - 1.0);
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

        // \todo Remove this hack for reducing apical lamina stiffness
        if (rElement.GetIndex() == 1)
        {
            well_depth *= 0.1;
        }
    }
    else // regular element
    {
        node_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(elem_idx, false);

        well_depth = mElementWellDepth * intrinsicSpacingSquared / (node_spacing * node_spacing);
        rest_length = mElementRestLength * node_spacing;
        well_width *= node_spacing;
    }

    // Loop over nodes and calculate the force exerted on node i+1 by node i
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Index of the next node, calculated modulo number of nodes in this element
        unsigned next_idx = (node_idx + 1) % num_nodes;

        // If element rather than lamina, calculate the stuffness multiplier
        double multiplier = UseStiffnessMult(elem_idx, middle_idx, rElement.GetNode(node_idx)->GetRegion()) ? stiffness_mult : 1.0;

        // Morse force (derivative of Morse potential wrt distance between nodes
        force_to_next[node_idx] = rCellPopulation.rGetMesh().GetVectorFromAtoB(rElement.GetNodeLocation(node_idx),
                                                                               rElement.GetNodeLocation(next_idx));
        double normed_dist = norm_2(force_to_next[node_idx]);

        double morse_exp = exp((rest_length - normed_dist) / well_width);
        force_to_next[node_idx] *= 2.0 * well_width * well_depth * multiplier * morse_exp * (1.0 - morse_exp) / normed_dist;
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
}

template <unsigned DIM>
bool GradualChangeMorseMembraneForce<DIM>::UseStiffnessMult(unsigned elemIdx, unsigned middleIdx, unsigned nodeRegion)
{
    // Left
    if (elemIdx < middleIdx)
    {
        if (nodeRegion == RIGHT_APICAL_REGION)// || nodeRegion == RIGHT_PERIAPICAL_REGION)
        {
            return true;
        }
    }
    // Right
    else if (elemIdx > middleIdx)
    {
        if (nodeRegion == LEFT_APICAL_REGION)// || nodeRegion == LEFT_PERIAPICAL_REGION)
        {
            return true;
        }
    }

    return false;
}

template <unsigned DIM>
void GradualChangeMorseMembraneForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
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
double GradualChangeMorseMembraneForce<DIM>::GetElementWellDepth() const
{
    return mElementWellDepth;
}

template <unsigned DIM>
void GradualChangeMorseMembraneForce<DIM>::SetElementWellDepth(double elementWellDepth)
{
    mElementWellDepth = elementWellDepth;
}

template <unsigned DIM>
double GradualChangeMorseMembraneForce<DIM>::GetElementRestLength() const
{
    return mElementRestLength;
}

template <unsigned DIM>
void GradualChangeMorseMembraneForce<DIM>::SetElementRestLength(double elementRestLength)
{
    mElementRestLength = elementRestLength;
}

template <unsigned DIM>
double GradualChangeMorseMembraneForce<DIM>::GetLaminaWellDepth() const
{
    return mLaminaWellDepth;
}

template <unsigned DIM>
void GradualChangeMorseMembraneForce<DIM>::SetLaminaWellDepth(double laminaWellDepth)
{
    mLaminaWellDepth = laminaWellDepth;
}

template <unsigned DIM>
double GradualChangeMorseMembraneForce<DIM>::GetLaminaRestLength() const
{
    return mLaminaRestLength;
}

template <unsigned DIM>
void GradualChangeMorseMembraneForce<DIM>::SetLaminaRestLength(double laminaRestLength)
{
    mLaminaRestLength = laminaRestLength;
}

template <unsigned DIM>
double GradualChangeMorseMembraneForce<DIM>::GetWellWidth() const
{
    return mWellWidth;
}

template <unsigned DIM>
void GradualChangeMorseMembraneForce<DIM>::SetWellWidth(double wellWidth)
{
    mWellWidth = wellWidth;
}

template <unsigned DIM>
double GradualChangeMorseMembraneForce<DIM>::GetStiffnessMult() const
{
    return mStiffnessMult;
}

template <unsigned DIM>
void GradualChangeMorseMembraneForce<DIM>::SetStiffnessMult(double stiffnessMult)
{
    mStiffnessMult = stiffnessMult;
}

// Explicit instantiation
template class GradualChangeMorseMembraneForce<1>;
template class GradualChangeMorseMembraneForce<2>;
template class GradualChangeMorseMembraneForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GradualChangeMorseMembraneForce)
