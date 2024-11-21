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
 * \ingroup Mapper
 * \brief defines a mapper to get for a velocity dof the connected elements.
 */
#ifndef DUMUX_FACE_DOF_TO_ELEMENT_INDEX_MAPPER_HH
#define DUMUX_FACE_DOF_TO_ELEMENT_INDEX_MAPPER_HH

#include <vector>
#include <array>
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup Mapper
 * \brief defines a mapper to get for a face dof the connected scvs.
 */
template<class GridGeometry>
class DofToEIdxMapper
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexType = typename SubControlVolume::Traits::GridIndexType;

public:
    DofToEIdxMapper()
    { }

    IndexType neighboringElementIdx(const SubControlVolume& scv) const
    {
        const auto& indices = dofToScvIdxMap_[scv.dofIndex()];
        return (indices[0] == scv.elementIndex()) ? indices[1] : indices[0];
    }

    void update(const GridGeometry& gg)
    { update_(gg); }

protected:
    void update_(const GridGeometry& gg)
    {
        std::cout << "gg.numDofs(): " << gg.numDofs() << std::endl;
        dofToScvIdxMap_.resize(gg.numDofs(), std::array<int,2>({-1,-1}));

        // assign scvs indices to each dof
        auto fvGeometry = localView(gg);
        for (const auto& element : elements(gg.gridView()))
        {
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto& indices = dofToScvIdxMap_[scv.dofIndex()];
                auto idx = (indices[0] == -1) ? 0 : 1;
                indices[idx] = scv.elementIndex();
            }
        }

    }
    std::vector<std::array<int,2>> dofToScvIdxMap_;
};

} // end namespace Dumux

#endif
