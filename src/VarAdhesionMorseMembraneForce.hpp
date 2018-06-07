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

#ifndef VARADHESIONMORSEMEMBRANEFORCE_HPP_
#define VARADHESIONMORSEMEMBRANEFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"

/**
 * A force class for use in immersed boundary simulations. This force implements Morse-potential-like links between
 * adjacent nodes in each immersed boundary. https://en.wikipedia.org/wiki/Morse_potential
 * The well width is a constant interaction strength, the rest length is an equilibrium bond distance, and the well
 * width is a parameter governing the profile of the curve.
 */
template <unsigned DIM>
class VarAdhesionMorseMembraneForce : public AbstractImmersedBoundaryForce<DIM>
{
private:
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
        archive& boost::serialization::base_object<AbstractImmersedBoundaryForce<DIM> >(*this);
        archive& mElementWellDepth;
        archive& mElementRestLength;
        archive& mLaminaWellDepth;
        archive& mLaminaRestLength;
        archive& mApicalWellDepthMult;
        archive& mWellWidth;
        archive& mStiffnessMult;
        archive& mSupportStrength;
    }

    /** The basic interaction strength */
    double mElementWellDepth;

    /** The rest length associated with each element as a fraction of the average node spacing */
    double mElementRestLength;

    /** The spring constant associated with each lamina */
    double mLaminaWellDepth;

    /** The rest length associated with each lamina as a fraction of the average node spacing */
    double mLaminaRestLength;

    /** Stiffness multiplier for the apical lamina */
    double mApicalWellDepthMult;

    /** The well width as a fraction of the average node spacing in either an element or lamina */
    double mWellWidth;

    /** The multiplicative factor to change the well depth by for particular nodes */
    double mStiffnessMult;

    /** Strength of the springs used to support the cell's shape */
    double mSupportStrength;

    /** Strength of the diagonal spring as a proportion of the support strength */
    double mDiagonalFraction;

    /** Frequency of oscillation for cyclic stiffness calculation */
    double mCyclicFrequency;

    /** The proportion of each cycle at which the gradient is applied */
    double mGradientOnProportion;

    /**
     * Array representing the number of elements in each of the three regions
     */
    std::array<unsigned, 3> mRegionSizes;

    /** Vector containing a random phase in [0, 2pi) for use in cyclic stiffness calculation */
    std::vector<double> mElementPhases;

    /**
     * Helper method for CalculateForcesOnElement().
     * Determine whether to modify stiffness, based on the elem region and node region.
     * @param elem_region the elem region (left, centre, right, membrane)
     * @param node_region the node region (RIGHT_APICAL, RIGHT_BASAL, etc)
     * @param elem_idx the index of the current element
     * @param gradientWeight the weighting due to the location of the element with index elem_idx
     * @return the stiffness multiplier altering interaction properties
     */
    double CalculateStiffnessMult(unsigned elem_region,
                                   unsigned node_region,
                                   unsigned elem_idx,
                                   double gradientWeight);

    /**
     * Helper method for AddImmersedBoundaryForceContribution.
     * Calculates forces, and can accept either an element or a lamina
     *
     * @tparam ELEMENT_DIM either DIM or DIM-1 depending on whether receiving an element or a lamina
     * @param rElement the element or lamina add forces to
     * @param rCellPopulation the immersed boundary cell population
     */
    template <unsigned ELEMENT_DIM>
    void CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, DIM>& rElement,
                                  ImmersedBoundaryCellPopulation<DIM>& rCellPopulation,
                                  double intrinsicSpacingSquared);

public:
    /** Constructor */
    VarAdhesionMorseMembraneForce();

    /** Destructor */
    virtual ~VarAdhesionMorseMembraneForce() = default;

    /**
     * Overridden AddImmersedBoundaryForceContribution() method.
     * Calculates basic elasticity in the membrane of each immersed boundary as a result of interactions.
     *
     * @param rNodePairs reference to a vector set of node pairs between which to contribute the force
     * @param rCellPopulation reference to the cell population
     */
    void AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                              ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputImmersedBoundaryForceParameters() method.
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile);

    /** @return mElementWellDepth */
    double GetElementWellDepth() const;

    /** @param elementWellDepth the new value of mElementWellDepth */
    void SetElementWellDepth(double elementWellDepth);

    /** @return mElementRestLength */
    double GetElementRestLength() const;

    /** @param elementRestLength the new value of mElementRestLength */
    void SetElementRestLength(double elementRestLength);

    /** @return mLaminaWellDepth */
    double GetLaminaWellDepth() const;

    /** @param laminaWellDepth the new value of mLaminaWellDepth */
    void SetLaminaWellDepth(double laminaWellDepth);

    /** @return mLaminaRestLength */
    double GetLaminaRestLength() const;

    /** @param laminaRestLength the new value of mLaminaRestLength */
    void SetLaminaRestLength(double laminaRestLength);

    /** @return mApicalWellDepthMult */
    double GetApicalWellDepthMult() const;

    /** @param apicalWellDepthMult the new value of mApicalWellDepthMult */
    void SetApicalWellDepthMult(double apicalWellDepthMult);

    /** @return mWellWidth */
    double GetWellWidth() const;

    /** @param wellWidth the new value of mWellWidth */
    void SetWellWidth(double wellWidth);

    /** @return mStiffnessMult */
    double GetStiffnessMult() const;

    /** @param stiffnessMult the new value of mStiffnessMult */
    void SetStiffnessMult(double stiffnessMult);

    /** @return mSupportStrength */
    double GetSupportStrength() const;

    /** @param supportStrength the new value of mSupportStrength */
    void SetSupportStrength(double supportStrength);

    /** @return mDiagonalFraction */
    double GetDiagonalFraction() const;

    /** @param diagonalFraction the new value of mDiagonalFraction */
    void SetDiagonalFraction(double diagonalFraction);

    /** @return mCyclicFrequency */
    double GetCyclicFrequency() const;

    /** @param cyclicFrequency the new value of mCyclicFrequency */
    void SetCyclicFrequency(double cyclicFrequency);

    /** @return mGradientOnProportion */
    double GetGradientOnProportion() const;

    /** @param gradientOnProportion the new value of mGradientOnProportion */
    void SetGradientOnProportion(double gradientOnProportion);

    /** @return mRegionSizes */
    const std::array<unsigned int, 3>& GetRegionSizes() const;

    /** @param regionSizes the new value of  mRegionSizes */
    void SetRegionSizes(const std::array<unsigned int, 3>& regionSizes);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VarAdhesionMorseMembraneForce)

#endif /*VARADHESIONMORSEMEMBRANEFORCE_HPP_*/
