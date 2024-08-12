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

#ifndef ENABLECACHING
#define ENABLECACHING 1
#endif

#ifndef GRID_DIM
#define GRID_DIM 3
#endif

#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Pseudo3DStokesVariableHeight {};
struct Pseudo3DStokesVariableHeightMomentum { using InheritsFrom = std::tuple<Pseudo3DStokesVariableHeight, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct Pseudo3DStokesVariableHeightMass { using InheritsFrom = std::tuple<Pseudo3DStokesVariableHeight, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Pseudo3DStokesVariableHeight>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Pseudo3DStokesVariableHeight>
{
    static constexpr int dim = GRID_DIM;

    using HostGrid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, dim> >;

#if HAVE_DUNE_SUBGRID
    using type = Dune::SubGrid<HostGrid::dimension, HostGrid>;
#else
    using type = HostGrid;
#endif
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Pseudo3DStokesVariableHeight>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Channel3DSpatialParams<GridGeometry, Scalar>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Pseudo3DStokesVariableHeightMomentum>
{ using type = Pseudo3DStokesVariableHeightProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::Pseudo3DStokesVariableHeightMass>
{ using type = Pseudo3DStokesVariableHeightProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Pseudo3DStokesVariableHeight> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Pseudo3DStokesVariableHeight> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Pseudo3DStokesVariableHeight> { static constexpr bool value = ENABLECACHING; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Pseudo3DStokesVariableHeight>
{
    using Traits = MultiDomainTraits<TTag::Pseudo3DStokesVariableHeightMomentum, TTag::Pseudo3DStokesVariableHeightMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties
