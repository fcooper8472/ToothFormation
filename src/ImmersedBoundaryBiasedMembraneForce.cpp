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

#include "ImmersedBoundaryBiasedMembraneForce.hpp"
#include "UblasVectorInclude.hpp"

template <unsigned DIM>
ImmersedBoundaryBiasedMembraneForce<DIM>::ImmersedBoundaryBiasedMembraneForce()
        : AbstractImmersedBoundaryForce<DIM>(),
          mElementSpringConst(1e6),
          mElementRestLength(0.5),
          mLaminaSpringConst(1e6),
          mLaminaRestLength(0.5)
{
}

template <unsigned DIM>
ImmersedBoundaryBiasedMembraneForce<DIM>::~ImmersedBoundaryBiasedMembraneForce()
{
}

template <unsigned DIM>
void ImmersedBoundaryBiasedMembraneForce<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                                                                    ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    // Data common across the entire cell population
    double intrinsicSpacingSquared = rCellPopulation.GetIntrinsicSpacing() * rCellPopulation.GetIntrinsicSpacing();

    // Loop over all elements ( <DIM, DIM> )
    for (auto elem_it = rCellPopulation.rGetMesh().GetElementIteratorBegin(); elem_it != rCellPopulation.rGetMesh().GetElementIteratorEnd(); ++elem_it)
    {
        CalculateForcesOnElement(*elem_it, rCellPopulation, intrinsicSpacingSquared);
    }

    // Loop over all laminas ( <DIM-1, DIM> )
    for (auto lam_it = rCellPopulation.rGetMesh().GetLaminaIteratorBegin(); lam_it != rCellPopulation.rGetMesh().GetLaminaIteratorEnd(); ++lam_it)
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
void ImmersedBoundaryBiasedMembraneForce<DIM>::CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, DIM>& rElement,
                                                                        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation,
                                                                        double intrinsicSpacingSquared)
{
    // Get index and number of nodes of current element
    unsigned elem_idx = rElement.GetIndex();
    unsigned num_nodes = rElement.GetNumNodes();

    // Make a vector to store the force on node i+1 from node i
    std::vector<c_vector<double, DIM> > elastic_Forceto_next_node(num_nodes);

    // Lambda to compare two nodes in this element by their x-position
    auto compare_x = [&rElement](const unsigned& a, const unsigned& b) -> bool
    {
        return rElement.GetNode(a)->rGetLocation()[0] < rElement.GetNode(b)->rGetLocation()[0];
    };

    // Lambda to compare two nodes in this element by their y-position
    auto compare_y = [&rElement](const unsigned& a, const unsigned& b) -> bool
    {
        return rElement.GetNode(a)->rGetLocation()[1] < rElement.GetNode(b)->rGetLocation()[1];
    };

    // Local indices of nodes in this element, sorted by x-pos
    std::vector<unsigned> sorted_indices(rElement.GetNumNodes());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0u);
    std::sort(sorted_indices.begin(), sorted_indices.end(), compare_x);

    // Create vectors containing indices of the left and right-most nodes in the element
    const unsigned num_special_nodes = std::lround(sorted_indices.size() * 0.3);

    std::vector<unsigned> lt_most;
    std::copy_n(sorted_indices.begin(), num_special_nodes, std::back_inserter(lt_most));

    std::vector<unsigned> rt_most;
    std::copy_n(sorted_indices.rbegin(), num_special_nodes, std::back_inserter(rt_most));

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
    double spring_constant = 0.0;
    double rest_length = 0.0;

    // Determine if we're in a lamina or not
    if (ELEMENT_DIM < DIM)
    {
        node_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfLamina(elem_idx, false);

        spring_constant = mLaminaSpringConst * intrinsicSpacingSquared / (node_spacing * node_spacing);
        rest_length = mLaminaRestLength * node_spacing;
    }
    else // regular element
    {
        node_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(elem_idx, false);

        spring_constant = mElementSpringConst * intrinsicSpacingSquared / (node_spacing * node_spacing);
        rest_length = mElementRestLength * node_spacing;
    }

    // Loop over nodes and calculate the force exerted on node i+1 by node i
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Index of the next node, calculated modulo number of nodes in this element
        unsigned next_idx = (node_idx + 1) % num_nodes;

        double modified_spring_constant = spring_constant;
        double modified_rest_length = rest_length;

        if (std::any_of(lt_most.begin(), lt_most.end(), [&](const unsigned& i){return node_idx == i;}))
        {
            modified_spring_constant *= 4.0;
        }
        else if (std::any_of(rt_most.begin(), rt_most.end(), [&](const unsigned& i){return node_idx == i;}))
        {
            modified_spring_constant *= 4.0;
        }

        // Hooke's law linear spring force
        elastic_Forceto_next_node[node_idx] = rCellPopulation.rGetMesh().GetVectorFromAtoB(rElement.GetNodeLocation(node_idx), rElement.GetNodeLocation(next_idx));
        double normed_dist = norm_2(elastic_Forceto_next_node[node_idx]);
        elastic_Forceto_next_node[node_idx] *= modified_spring_constant * (normed_dist - modified_rest_length) / normed_dist;
    }

    // Add the contributions of springs adjacent to each node
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Get index of previous node
        unsigned prev_idx = (node_idx + num_nodes - 1) % num_nodes;

        c_vector<double, DIM> aggregate_force = elastic_Forceto_next_node[node_idx] - elastic_Forceto_next_node[prev_idx];

        // Add the aggregate force contribution to the node
        rElement.GetNode(node_idx)->AddAppliedForceContribution(aggregate_force);
    }

    // Add on width/height springs
    auto lf_min_max = std::minmax_element(lt_most.begin(), lt_most.end(), compare_y);
    auto rt_min_max = std::minmax_element(rt_most.begin(), rt_most.end(), compare_y);

    const double width = 0.4;
    const double height = 0.6;
    const double factor = 0.25;

    Node<DIM>* p_top_lt = rElement.GetNode(*lf_min_max.second);
    Node<DIM>* p_top_rt = rElement.GetNode(*rt_min_max.second);
    Node<DIM>* p_bot_lt = rElement.GetNode(*lf_min_max.first);
    Node<DIM>* p_bot_rt = rElement.GetNode(*rt_min_max.first);

    auto left = p_top_lt->rGetLocation() - p_bot_lt->rGetLocation();
    const double len_1 = norm_2(left);
    const auto force_1 = factor * (len_1 - height) * mElementSpringConst * left / len_1;
    p_top_lt->AddAppliedForceContribution(-force_1);
    p_bot_lt->AddAppliedForceContribution(force_1);

    auto right = p_top_rt->rGetLocation() - p_bot_rt->rGetLocation();
    const double len_2 = norm_2(right);
    const auto force_2 = factor * (len_2 - height) * mElementSpringConst * right / len_2;
    p_top_rt->AddAppliedForceContribution(-force_2);
    p_bot_rt->AddAppliedForceContribution(force_2);

    auto top = p_top_lt->rGetLocation() - p_top_rt->rGetLocation();
    const double len_3 = norm_2(top);
    const auto force_3 = factor * (len_3 - width) * mElementSpringConst * top / len_3;
    p_top_lt->AddAppliedForceContribution(-force_3);
    p_top_rt->AddAppliedForceContribution(force_3);

    auto bot = p_bot_lt->rGetLocation() - p_bot_rt->rGetLocation();
    const double len_4 = norm_2(bot);
    const auto force_4 = factor * (len_4 - width) * mElementSpringConst * bot / len_4;
    p_bot_lt->AddAppliedForceContribution(-force_4);
    p_bot_rt->AddAppliedForceContribution(force_4);
}

template <unsigned DIM>
void ImmersedBoundaryBiasedMembraneForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ElementSpringConstant>" << mElementSpringConst << "</ElementSpringConstant>\n";
    *rParamsFile << "\t\t\t<ElementRestLength>" << mElementRestLength << "</ElementRestLength>\n";
    *rParamsFile << "\t\t\t<LaminaSpringConstant>" << mLaminaSpringConst << "</LaminaSpringConstant>\n";
    *rParamsFile << "\t\t\t<LaminaRestLength>" << mLaminaRestLength << "</LaminaRestLength>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template <unsigned DIM>
double ImmersedBoundaryBiasedMembraneForce<DIM>::GetElementSpringConst() const
{
    return mElementSpringConst;
}

template <unsigned DIM>
void ImmersedBoundaryBiasedMembraneForce<DIM>::SetElementSpringConst(double elementSpringConst)
{
    mElementSpringConst = elementSpringConst;
}

template <unsigned DIM>
double ImmersedBoundaryBiasedMembraneForce<DIM>::GetElementRestLength() const
{
    return mElementRestLength;
}

template <unsigned DIM>
void ImmersedBoundaryBiasedMembraneForce<DIM>::SetElementRestLength(double elementRestLength)
{
    mElementRestLength = elementRestLength;
}

template <unsigned DIM>
double ImmersedBoundaryBiasedMembraneForce<DIM>::GetLaminaSpringConst() const
{
    return mLaminaSpringConst;
}

template <unsigned DIM>
void ImmersedBoundaryBiasedMembraneForce<DIM>::SetLaminaSpringConst(double laminaSpringConst)
{
    mLaminaSpringConst = laminaSpringConst;
}

template <unsigned DIM>
double ImmersedBoundaryBiasedMembraneForce<DIM>::GetLaminaRestLength() const
{
    return mLaminaRestLength;
}

template <unsigned DIM>
void ImmersedBoundaryBiasedMembraneForce<DIM>::SetLaminaRestLength(double laminaRestLength)
{
    mLaminaRestLength = laminaRestLength;
}

// Explicit instantiation
template class ImmersedBoundaryBiasedMembraneForce<1>;
template class ImmersedBoundaryBiasedMembraneForce<2>;
template class ImmersedBoundaryBiasedMembraneForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryBiasedMembraneForce)
