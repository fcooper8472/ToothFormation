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

#ifndef SIMULATIONSTOPPINGMODIFIER_HPP_
#define SIMULATIONSTOPPINGMODIFIER_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"

/**
 * A modifier class which after each simulation time step exports the simulation geometry to a scalable vector graphic
 * (SVG) file.  Specialisation for "Three Region" tooth formation simulations.
 */
template <unsigned DIM>
class SimulationStoppingModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
    /** Threshold lamina angle, default value is 10 degrees (in radians) */
    double mThreshold = 0.174532925;


    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM> >(*this);
    }

    /**
     * Calculate the angle of the basal lamina, defined as:
     *
     *   arctan( (MAX - MIN) / 0.5 )
     *
     * where MAX is the maximum y-value of the lamina, and MIN is the minimum y-value of the lamina.
     *
     * @param pIbMesh the immersed boundary mesh
     * @return the lamina bend
     */
    double CalculateLaminaAngle(ImmersedBoundaryMesh<DIM, DIM>* const pIbMesh);

public:

    /**
     * Default constructor.
     */
    SimulationStoppingModifier();

    /**
     * Destructor.
     */
    virtual ~SimulationStoppingModifier() = default;

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /** @return mThreshold */
    double GetThreshold() const;

    /** @param threshold the new value for mThreshold */
    void SetThreshold(double threshold);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimulationStoppingModifier)

#endif /*SIMULATIONSTOPPINGMODIFIER_HPP_*/
