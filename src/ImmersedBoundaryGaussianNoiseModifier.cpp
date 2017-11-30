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

#include "ImmersedBoundaryGaussianNoiseModifier.hpp"

#include "ChasteMakeUnique.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "UblasCustomFunctions.hpp"

#include "Exception.hpp"

template <unsigned DIM>
ImmersedBoundaryGaussianNoiseModifier<DIM>::ImmersedBoundaryGaussianNoiseModifier()
        : AbstractCellBasedSimulationModifier<DIM>(),
          mLowerCorner(),
          mUpperCorner(),
          mLengthScale(0.0),
          mpRandomField(nullptr)
{
}


template <unsigned DIM>
void ImmersedBoundaryGaussianNoiseModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
}

template <unsigned DIM>
void ImmersedBoundaryGaussianNoiseModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    // Verify members are set appropriately
    {
        for (unsigned dim = 0; dim < DIM; ++dim)
        {
            const double width_this_dim = mUpperCorner[dim] - mLowerCorner[dim];
            if (width_this_dim < DBL_EPSILON)
            {
                EXCEPTION("Domain width must be > 0. Check how mLowerCorner and mUpperCorner are set.");
            }

            if (mLowerCorner[dim] < 0.0 || mUpperCorner[dim] > 1.0)
            {
                EXCEPTION("Domain must fit within [0,1]x[0,1]. Check how mLowerCorner and mUpperCorner are set.");
            }
        }

        if (mLengthScale < DBL_EPSILON)
        {
            EXCEPTION("Correlation length scale must be > 0. Ensure mLengthScale is set correctly.");
        }
    }

    // Set other variables needed for the random field generator.
    auto p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>*>(&rCellPopulation.rGetMesh());
    const double grid_spacing = 1.0 / p_mesh->GetNumGridPtsX();

    std::array<unsigned, DIM> num_grid_pts;
    std::array<bool, DIM> periodicity;

    for (unsigned dim = 0; dim < DIM; ++dim)
    {
        const double width_this_dim = mUpperCorner[dim] - mLowerCorner[dim];
        num_grid_pts[dim] = std::lround(width_this_dim / grid_spacing);

        periodicity[dim] = std::fabs(1.0 - width_this_dim) < DBL_EPSILON;
    }

    // Calculate about a third of the eigenvalues \todo: remove this magic number
    const unsigned num_eigenvalues = 0.33 * std::accumulate(num_grid_pts.begin(), num_grid_pts.end(), 1, std::multiplies<unsigned>());


    // Set up the random field generator
    mpRandomField = our::make_unique<UniformGridRandomFieldGenerator<DIM>>(
            mLowerCorner, mUpperCorner, num_grid_pts, periodicity, num_eigenvalues, mLengthScale);
}

template <unsigned DIM>
void ImmersedBoundaryGaussianNoiseModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ImmersedBoundaryGaussianNoiseModifier<1>;
template class ImmersedBoundaryGaussianNoiseModifier<2>;
template class ImmersedBoundaryGaussianNoiseModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryGaussianNoiseModifier)
