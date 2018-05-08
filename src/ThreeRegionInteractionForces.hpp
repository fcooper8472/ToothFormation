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

#ifndef THREEREGIONINTERACTIONFORCES_HPP_
#define THREEREGIONINTERACTIONFORCES_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"

#include <iostream>

/**
 * A force class for use in Vertex-based simulations. This force is based on the
 * energy function proposed by Farhadifar et al in  Curr. Biol., 2007, 17, 2095-2104.
 */
template<unsigned DIM>
class ThreeRegionInteractionForces : public AbstractImmersedBoundaryForce<DIM>
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
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractImmersedBoundaryForce<DIM> >(*this);
        archive & mElemAttributeLocation;
        archive & mBasicInteractionStrength;
        archive & mBasicInteractionDist;
    }

protected:

    /**
     * The permitted types of element in these simulations
     */
    enum ThreeRegionElemType
    {
        THREE_REGION_LEFT,
        THREE_REGION_MID,
        THREE_REGION_RIGHT,
        THREE_REGION_LAMINA
    };

    /**
     * The immersed boundary mesh.
     */
    ImmersedBoundaryMesh<DIM,DIM>* mpMesh;

    /**
     * The location in the element attriutes vector of the three-region cell type
     */
    unsigned mElemAttributeLocation;

    /**
     * The cell-cell spring constant.
     *
     * Initialised to 1e3 during setup.
     */
    double mBasicInteractionStrength;

    /**
     * The cell-cell rest length as a proportion of the cell interaction distance.
     *
     * Initialised to 0.25 times the cell interaction distance during setup.
     */
    double mBasicInteractionDist;

    /**
     * The multiplyer by which we increase certain adhesions.
     */
    double mAdhesionMultiplier;

    /**
     * Array representing the number of elements in each of the three regions
     */
    std::array<unsigned, 3> mRegionSizes;

    /**
     * Get the type of requested element by interrogating its attribute vector
     *
     * @param elemIdx the index of the element
     * @return the element type
     */
    unsigned GetElemType(unsigned elemIdx);

    /**
     * Get the element type given a node, by interrogating the attribute vector
     * of its containing element
     *
     * @param pNode the node
     * @return the element type
     */
    unsigned GetElemType(Node<DIM>* pNode);

    /**
     * Calculate a multiplier based on the element type and node region of the
     * currently-interacting pair of nodes.
     *
     * @param pNodeA the first node in the pair
     * @param pNodeB the second node in the pair
     * @return the type multiplier
     */
    double CalculateElementTypeMult(Node<DIM>* pNodeA, Node<DIM>* pNodeB);

public:

    /**
     * Constructor.
     */
    ThreeRegionInteractionForces();

    /**
     * Destructor.
     */
    virtual ~ThreeRegionInteractionForces() = default;

    /**
     * Overridden AddImmersedBoundaryForceContribution() method.
     *
     * Calculates the force on each node in the immersed boundary cell population as a result of cell-cell interactions.
     *
     * @param rNodePairs reference to a vector set of node pairs between which to contribute the force
     * @param rCellPopulation reference to the cell population
     */
    void AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                              ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

    /**
     * @param basicInteractionDist the new value of mBasicInteractionStrength
     */
    void SetBasicInteractionStrength(double basicInteractionStrength);

    /**
     * @return mBasicInteractionStrength
     */
    double GetBasicInteractionStrength();

    /**
     * @param basicInteractionDist the new value of mBasicInteractionDist
     */
    void SetBasicInteractionDist(double basicInteractionDist);

    /**
     * @return mBasicInteractionDist
     */
    double GetBasicInteractionDist();

    /**
     * @param adhesionMultiplier the new value of mAdhesionMultiplier
     */
    void SetAdhesionMultiplier(double adhesionMultiplier);

    /**
     * @return mAdhesionMultiplier
     */
    double GetAdhesionMultiplier();

    /**
     * @return mRegionSizes
     */
    const std::array<unsigned int, 3>& GetRegionSizes() const;

    /**
     * @param regionSizes the new value of  mRegionSizes
     */
    void SetRegionSizes(const std::array<unsigned int, 3>& regionSizes);

    /**
     * Overridden OutputImmersedBoundaryForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ThreeRegionInteractionForces)

#endif /*THREEREGIONINTERACTIONFORCES_HPP_*/
