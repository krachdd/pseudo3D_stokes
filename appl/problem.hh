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

#include <dune/common/float_cmp.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

#include "doftoelementmapper.hh"
namespace Dumux {

template <class TypeTag, class BaseProblem>
class Pseudo3DStokesVariableHeightProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr bool enablePseudoThreeDWallFriction = dim != 3;

public:
    Pseudo3DStokesVariableHeightProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        pRef_ = getParam<Scalar>("Problem.pRef");
        deltaP_ = getParam<Scalar>("Problem.DeltaP");
        maxHeight_ = getParam<Scalar>("Problem.Height");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        verticalFlow_ = getParam<bool>("Problem.VerticalFlow", false);
        std::cout << "Vertical flow is " << verticalFlow_;
        height_ = getParam<Scalar>("Problem.Height");
        if(dim == 3 && !Dune::FloatCmp::eq(height_, this->gridGeometry().bBoxMax()[2]))
            DUNE_THROW(Dune::InvalidStateException, "z-dimension must equal height");

        addPressureCorrection_ = getParam<bool>("Problem.AddPressureCorrection", true);

        if constexpr (ParentType::isMomentumProblem())
        {
            eIdxMap_ = DofToEIdxMapper<GridGeometry>();
            eIdxMap_.value().update(this->gridGeometry());
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume
     */
    template<class ElementVolumeVariables>
    Sources source(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume& scv) const
    {
        auto source = Sources(0.0);

        if constexpr (ParentType::isMomentumProblem() && enablePseudoThreeDWallFriction)
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            Scalar height = (*heights_)[eIdx];
            static const Scalar factor = getParam<Scalar>("Problem.PseudoWallFractionFactor", 12.0); //8.0 for vmax //12.0 for vmean

            // prefactor due to preFactorDrag
            Scalar preFactorDrag = 1;
            if(scv.dofAxis() == 0)
                preFactorDrag = (*preFactorDrag_x_)[eIdx];
            else
                preFactorDrag = (*preFactorDrag_y_)[eIdx];

            source = this->pseudo3DWallFriction(element, fvGeometry, elemVolVars, scv, height, factor*preFactorDrag);

            if(addPressureCorrection_)
            {
                const auto& scvf = (*scvfs(fvGeometry, scv).begin());
                Scalar correction(0.0);

                // For the frontal face we need to add the pressure correction
                const auto pressure = this->pressure(element, fvGeometry, scvf) - this->referencePressure(element, fvGeometry, scvf);
                correction = pressure*Extrusion::area(fvGeometry, scvf)*elemVolVars[scv].extrusionFactor();

                // Account for the orientation of the face.
                correction *= scvf.directionSign();

                // The actual factor which needs to be corrected due to porosity changes
                if(!scv.boundary())
                {
                    const auto eIdxI = scv.elementIndex();
                    const auto eIdxJ = scv.boundary() ? eIdxI : eIdxMap_.value().neighboringElementIdx(scv);
                    const auto h_avg = 0.5*((*heights_)[eIdxI] +(*heights_)[eIdxJ]);
                    correction *= (1.0 - h_avg/maxHeight_);
                }

                // Scale by volume
                correction /= Extrusion::volume(fvGeometry, scv)*elemVolVars[scv].extrusionFactor();

                source += correction;
            }

        }
        return source;
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values.setAllDirichlet();
            if (isOutlet_(globalPos) || isInlet_(globalPos))
                values.setAllNeumann();
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        // no-flow/no-slip
        return DirichletValues(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
        const auto& globalPos = scvf.ipGlobal();

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto p = isInlet_(globalPos) ? pRef_ + deltaP_ : pRef_;
            values = NavierStokesMomentumBoundaryFlux<typename GridGeometry::DiscretizationMethod>::fixedPressureMomentumFlux(
                *this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, p
            );
        }
        else
        {
            values = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>::scalarOutflowFlux(
                *this, element, fvGeometry, scvf, elemVolVars
            );
        }
        return values;
    }

    // \}

    //! Returns the analytical solution for the flux through the rectangular channel
    Scalar analyticalFlux() const
    {
        const Scalar h = height_;
        const Scalar w = this->gridGeometry().bBoxMax()[1];
        const Scalar L = this->gridGeometry().bBoxMax()[0];

        const Scalar mu = nu_*rho_;

        return h*h*h * w * deltaP_ / (12*mu*L) * (1.0 - 0.630 * h/w);
    }

    //Returns auxiliary calculation for darcypermeability
    Scalar darcyPermFactor() const
    {
        const Scalar mu = nu_*rho_;
        const Scalar h = maxHeight_;
        const Scalar w = this->gridGeometry().bBoxMax()[1];
        const Scalar L = this->gridGeometry().bBoxMax()[0];
        if (verticalFlow_)
            return mu/ ((deltaP_/w) * (L*h) );//TODO I think deltaP_ should be divided by w instead of L
        else
            return mu/ ((deltaP_/L) * (w*h) );
    }

    //Returns auxiliary calculation for darcypermeability if you use half of the domain
    // not ideal, should be replaced by some xmin, max stuff
    Scalar darcyPermFactorHalfDomain() const
    {
        return darcyPermFactor() * 2.0;
    }

    void setrelHeights(const std::shared_ptr<const std::vector<Scalar>>& relHeights)
    {
        heights_ = relHeights;
    }

    void setpreFactorsDrag(const std::shared_ptr<const std::vector<Scalar>>& preFactorDrag_x,
                           const std::shared_ptr<const std::vector<Scalar>>& preFactorDrag_y)
    {
        preFactorDrag_x_ = preFactorDrag_x;
        preFactorDrag_y_ = preFactorDrag_y;
    }

    bool isVerticalFlow() const
    {
        return verticalFlow_;
    }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    {
        const auto flowDir = flowDirection_();
        return globalPos[flowDir] < eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        const auto flowDir = flowDirection_();
        return globalPos[flowDir] > this->gridGeometry().bBoxMax()[flowDir] - eps_;
    }

    int flowDirection_() const
    { return (verticalFlow_ ? 1 : 0);}

    bool verticalFlow_;
    static constexpr Scalar eps_=1e-6;
    Scalar deltaP_;
    Scalar gradP_;
    Scalar pRef_;
    Scalar extrusionFactor_;
    Scalar height_;
    Scalar rho_;
    Scalar nu_;
    Scalar permeability_;
    std::shared_ptr<const std::vector<Scalar>> heights_;
    std::shared_ptr<const std::vector<Scalar>> preFactorDrag_x_;
    std::shared_ptr<const std::vector<Scalar>> preFactorDrag_y_;
    Scalar maxHeight_;
    Scalar frictionFactor_;
    bool addPressureCorrection_;
    std::optional<DofToEIdxMapper<GridGeometry>> eIdxMap_;
};

} // end namespace Dumux
