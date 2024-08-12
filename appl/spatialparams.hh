// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup FreeflowTests
 * \brief Definition of the spatial parameters for the 3D Channel Problems.
 */

#ifndef DUMUX_CHANNEL3D_SPATIAL_PARAMS_HH
#define DUMUX_CHANNEL3D_SPATIAL_PARAMS_HH

#include <dumux/freeflow/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the freeflow problems.
 */
template<class GridGeometry, class Scalar>
class Channel3DSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, Channel3DSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, Channel3DSpatialParams<GridGeometry, Scalar>>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename SubControlVolume::Traits::GlobalPosition;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr bool enablePseudoThreeDWallFriction = dim != 3;

public:
    Channel3DSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        maxHeight_ = getParam<Scalar>("Problem.Height");
        preExtrusionFactor_ = getParam<Scalar>("SpatialParams.PreExtrusionFactor", 0.6666667);

        if constexpr (enablePseudoThreeDWallFriction)
            extrusionFactor_ = preExtrusionFactor_ * maxHeight_;
        else
            extrusionFactor_ = 1.0;
        std::cout<< "extrusionFactor_ is "<< extrusionFactor_ << std::endl;
    }

    /*!
     * \brief Return the temperature in the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 283.15; }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return extrusionFactor_; }

    void setrelHeights(const std::shared_ptr<const std::vector<Scalar>>& relHeights)
    {
        relHeight_ = relHeights;
    }
    void setrelHeightGrads(const std::shared_ptr<const std::vector<Scalar>>& relHeightGradsX,
                           const std::shared_ptr<const std::vector<Scalar>>& relHeightGradsY)
    {
        relHeightGradX_ = relHeightGradsX;
        relHeightGradY_ = relHeightGradsY;
    }

private:
    Scalar extrusionFactor_;
    Scalar preExtrusionFactor_;
    Scalar maxHeight_;
    std::shared_ptr<const std::vector<Scalar>> relHeight_;
    std::shared_ptr<const std::vector<Scalar>> relHeightGradX_;
    std::shared_ptr<const std::vector<Scalar>> relHeightGradY_;
};

} // end namespace Dumux

#endif
