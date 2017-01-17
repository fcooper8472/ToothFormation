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

#include "CellRegionWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellRegionWriter<ELEMENT_DIM, SPACE_DIM>::CellRegionWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cell_regions.dat")
{
    this->mVtkCellDataName = "Cell Regions";
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellRegionWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // We first need to test that the cell population is vertex-based
    ImmersedBoundaryCellPopulation<SPACE_DIM>* p_ibcp = dynamic_cast<ImmersedBoundaryCellPopulation<SPACE_DIM>*>(pCellPopulation);

    if (p_ibcp == NULL)
    {
        EXCEPTION("This class is intended for use only with ThreeRegion immersed boundary simulations");
    }

    // Bit of a fudge - we rely on the correct information being in location 2 of the vector.
    double cell_region = 0.0;
    if (p_ibcp->GetElementCorrespondingToCell(pCell)->GetNumElementAttributes() > 2)
    {
        cell_region = p_ibcp->GetElementCorrespondingToCell(pCell)->rGetElementAttributes()[2];
    }
    else
    {
        EXCEPTION("Expecting 3 or more element attributes...");
    }

    return cell_region;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellRegionWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double cell_region = this->GetCellDataForVtkOutput(pCell, pCellPopulation);

    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << " " << cell_region << " " << location_index;

    c_vector<double, SPACE_DIM> coords = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i = 0; i < SPACE_DIM; i++)
    {
        *this->mpOutStream << " " << coords[i];
    }
}

// Explicit instantiation
template class CellRegionWriter<1, 1>;
template class CellRegionWriter<1, 2>;
template class CellRegionWriter<2, 2>;
template class CellRegionWriter<1, 3>;
template class CellRegionWriter<2, 3>;
template class CellRegionWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellRegionWriter)
