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

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/freeflow/navierstokes/fluxoveraxisalignedsurface.hh>

#include "properties.hh"


// Reads the image file into an array
void prepareGridWithPrefactor(
                const Dumux::Detail::RasterImageData::Result<std::size_t>& grayScaleImage,
                double maxHeight,
                 std::vector<unsigned char>& dataGrid, 
                 std::vector<double>& dataHeight,
                 std::vector<double>& preFactorDrag_x,
                 std::vector<double>& preFactorDrag_y)
{
    auto maxValue = max_element(std::begin(grayScaleImage), std::end(grayScaleImage));
    std::cout << "maximum value of GeometryFile is : " << *maxValue << std::endl;
    unsigned int j = 0;
    for(unsigned int i = 0; i < grayScaleImage.size(); i++) {
        if(grayScaleImage[i] != 0 )
        {
            dataHeight.push_back(( 1.*(grayScaleImage[i]) / *maxValue)*maxHeight);
            dataGrid.push_back(1);
        }
        else
        {
            preFactorDrag_x.erase(preFactorDrag_x.begin() + i-j);
            preFactorDrag_y.erase(preFactorDrag_y.begin() + i-j);
            dataGrid.push_back(0);
            j++;
        }
    }
}


template<class HostGridView, class SubGridView, class V>
std::shared_ptr<const std::vector<double>> getSubgridFields(const HostGridView& hostGridView,
                                                             const SubGridView& subGridView,
                                                             const V& hostGridHeights)
{
    std::vector<double> subgridHeights(subGridView.size(0));
    for (const auto& e : elements(subGridView))
    {
        const auto eIdx = subGridView.indexSet().index(e);
        subgridHeights[eIdx] = hostGridHeights[eIdx];   //
    }
    return std::make_shared<const std::vector<double>>(std::move(subgridHeights));
}


int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::Pseudo3DStokesVariableHeightMomentum;
    using MassTypeTag = Properties::TTag::Pseudo3DStokesVariableHeightMass;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid    
    using Grid = GetPropType<MassTypeTag, Properties::Grid>;
    Dumux::GridManager<Grid> gridManager;
    
    constexpr int dim = 2;
    using HostGrid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim> >;
    using SubGrid = Dune::SubGrid<dim, HostGrid>;
    using HostGridManager = Dumux::GridManager<HostGrid>;

    HostGridManager externalHostGridManager;
    externalHostGridManager.init();
    auto& hostGrid = externalHostGridManager.grid();
    const auto& gridView = externalHostGridManager.grid().leafGridView();
    
    std::string geomFile = getParam<std::string>("Grid.GeometryFile");
    std::string preFactorDrag_xFile = getParam<std::string>("Grid.PreFactorDragFileX", "noFile");
    std::string preFactorDrag_yFile = getParam<std::string>("Grid.PreFactorDragFileY", "noFile");

    const auto upperRight = getParam<std::vector<double>>("Grid.UpperRight");
    const auto cells = getParam<std::vector<double>>("Grid.Cells");
    const auto maxHeight = getParam<double>("Problem.Height");
    const auto writeVtu = getParam<bool>("Vtk.WriteVtuData");
    
    std::vector<double> pixelRes;
    pixelRes.push_back(upperRight[0]/cells[0]);
    pixelRes.push_back(upperRight[1]/cells[1]);
    
    std::cout << "Pixel Resolution is " << pixelRes[0] << ", " << pixelRes[1] << std::endl;
    std::cout << "maxHeight is " << maxHeight << std::endl;
    
    const auto grayScaleImage = NetPBMReader::readPGM<std::size_t>(geomFile);
    auto maxValueTest = *max_element(std::begin(grayScaleImage), std::end(grayScaleImage));
    std::cout << "Max value of gray scale image in main file is: " << maxValueTest << std::endl;
    std::cout << "grayScaleImage read"<< std::endl;

    std::vector<unsigned char> dataGrid;
    std::vector<double> dataHeight;
    std::vector<double> preFactorDrag_x (grayScaleImage.size(), 1);
    std::vector<double> preFactorDrag_y (grayScaleImage.size(), 1);
    
    std::cout << "before gird is prepared"<< std::endl;
    // Checks if there is a prefactor File given for x and y. If that is the case the values are read
    if(preFactorDrag_xFile.compare("noFile") != 0)
    {
        std::cout << "prefactor File is given for preFactorDrag_xFile" << std::endl;
        preFactorDrag_x = readFileToContainer<std::vector<double>>(preFactorDrag_xFile);
    }
    else
        std::cout << " NO prefactor File is given for preFactorDrag_xFile" << std::endl;
    if(preFactorDrag_yFile.compare("noFile") != 0)
    {
        std::cout << "prefactor File is given for preFactorDrag_yFile" << std::endl;
        preFactorDrag_y = readFileToContainer<std::vector<double>>(preFactorDrag_yFile);
    }
    else
        std::cout << " NO prefactor File is given for preFactorDrag_xFile" << std::endl;

    prepareGridWithPrefactor(grayScaleImage, maxHeight, dataGrid, dataHeight, preFactorDrag_x, preFactorDrag_y);
    std::cout << "gird is prepared"<< std::endl;


    using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;

    auto elementSelector = [&gridView, &dataGrid](const auto& element)
    {
        auto index = gridView.indexSet().index(element);
        bool isFluid = true;
        if(dataGrid[index]==0){ isFluid = false;}
        return isFluid;
    };
    Dumux::GridManager<SubGrid> subgridManager;
    subgridManager.init(hostGrid, elementSelector);

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = subgridManager.grid().leafGridView();


    const auto relHeights = getSubgridFields(gridView, leafGridView, dataHeight);
    const auto preFactorDrag_x_subGrid = getSubgridFields(gridView, leafGridView, preFactorDrag_x);
    const auto preFactorDrag_y_subGrid = getSubgridFields(gridView, leafGridView, preFactorDrag_y);

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    momentumProblem->setrelHeights(relHeights);
    momentumProblem->spatialParams().setrelHeights(relHeights);

    momentumProblem->setpreFactorsDrag(preFactorDrag_x_subGrid, preFactorDrag_y_subGrid);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);
    massProblem->setrelHeights(relHeights);
    massProblem->spatialParams().setrelHeights(relHeights);
    massProblem->setpreFactorsDrag(preFactorDrag_x_subGrid, preFactorDrag_y_subGrid);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // compute coupling stencil and afterwards initialize grid variables (need coupling information)
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.addField(dataHeight, "relheight");
    vtkWriter.addField(preFactorDrag_x, "preFactorDrag_x");
    vtkWriter.addField(preFactorDrag_y, "preFactorDrag_y");
    if(writeVtu)
        vtkWriter.write(0.0);

    // the linear solver
    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    //using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // set up 9 planes over which fluxes are calculated
    FluxOverAxisAlignedSurface flux(*massGridVariables, x[massIdx], assembler->localResidual(massIdx));

    using GridView = typename GetPropType<MassTypeTag, Properties::GridGeometry>::GridView;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    // Define relevant positions
    const Scalar xMin = massGridGeometry->bBoxMin()[0];
    const Scalar xMax = massGridGeometry->bBoxMax()[0];
    const Scalar xMax_2 = xMin + 0.5*(xMax - xMin);
    const Scalar xMax_4 = xMin + 0.25*(xMax - xMin);
    const Scalar yMin = massGridGeometry->bBoxMin()[1];
    const Scalar yMax = massGridGeometry->bBoxMax()[1];
    const Scalar yMax_2 = yMin + 0.5*(yMax - yMin);
    const Scalar yMax_4 = yMin + 0.25*(yMax - yMin);

//     std::size_t normalDirectionIndex = 0;
//     if(massProblem->isVerticalFlow())
//         normalDirectionIndex = 1;

#if GRID_DIM == 3
    const Scalar zMin = massGridGeometry->bBoxMin()[2];
    const Scalar zMax = massGridGeometry->bBoxMax()[2];
    const Scalar zMax_2 = zMin + 0.5*(zMax - zMin);
    const Scalar zMax_4 = zMin + 0.25*(zMax - zMin);
#endif

    // the first plane is at the inlet
#if GRID_DIM == 3
    const auto inletLowerLeft = GlobalPosition{xMin, yMin, zMin};
    const auto inletUpperRight = GlobalPosition{xMin, yMax, zMax};
    flux.addAxisAlignedSurface("inlet", inletLowerLeft, inletUpperRight);
#else
    auto inletLowerLeft = GlobalPosition{xMin, yMin};
    auto inletUpperRight = GlobalPosition{xMin, yMax};
    if(massProblem->isVerticalFlow())
    {
        inletUpperRight[0] = xMax;
        inletUpperRight[1] = yMin;
    }
    flux.addAxisAlignedSurface("inlet", inletLowerLeft, inletUpperRight);
#endif

    // the second plane is at the middle of the channel
    const Scalar planePosMiddleX = xMin + 0.5*(xMax - xMin);
    const Scalar planePosMiddleY = yMin + 0.5*(yMax - yMin);
#if GRID_DIM == 3
    const auto middleLowerLeft = GlobalPosition{planePosMiddleX, yMin, zMin};
    const auto middleUpperRight = GlobalPosition{planePosMiddleX, yMax, zMax};
    flux.addAxisAlignedSurface("middle", middleLowerLeft, middleUpperRight);
#else
    auto middleLowerLeft = GlobalPosition{planePosMiddleX, yMin};
    auto middleUpperRight = GlobalPosition{planePosMiddleX, yMax};
    if(massProblem->isVerticalFlow())
    {
        middleLowerLeft[0] = xMin;
        middleLowerLeft[1] = planePosMiddleY;
        middleUpperRight[0] = xMax;
        middleUpperRight[1] = planePosMiddleY;
    }
    flux.addAxisAlignedSurface("middle", middleLowerLeft, middleUpperRight);
#endif

    // The last plane is placed at the outlet of the channel.
#if GRID_DIM == 3
    const auto outletLowerLeft = GlobalPosition{xMax, yMin, zMin};
    const auto outletUpperRight = GlobalPosition{xMax, yMax, zMax};
    flux.addAxisAlignedSurface("outlet", outletLowerLeft, outletUpperRight);
#else
    auto outletLowerLeft = GlobalPosition{xMax, yMin};
    auto outletUpperRight = GlobalPosition{xMax, yMax};
    if(massProblem->isVerticalFlow())
    {
        outletLowerLeft[0] = xMin;
        outletLowerLeft[1] = yMax;
    }
    flux.addAxisAlignedSurface("outlet", outletLowerLeft, outletUpperRight);
#endif

    // Planes for aniso computation 4-9
    // 4-6 in main pressure gradient direction
    // 4 at in inlet, covering half of it
#if GRID_DIM == 3
    // inlet of the subdomain
    const auto aniso_main_inletLowerLeft     = GlobalPosition{ xMin,   yMin,   zMin   };
    const auto aniso_main_inletUpperRight    = GlobalPosition{ xMin,   yMax_2, zMax_2 };
    const auto aniso_perp_inletLowerLeft     = GlobalPosition{ xMin,   yMin,   zMin   };
    const auto aniso_perp_inletUpperRight    = GlobalPosition{ xMax_2, yMin,   zMax_2 };

    // center of the subdomain
    const auto aniso_main_middleLowerLeft    = GlobalPosition{ xMax_4, yMin,   zMin   };
    const auto aniso_main_middleUpperRight   = GlobalPosition{ xMax_4, yMax_2, zMax_2 };
    const auto aniso_perp_middleLowerLeft    = GlobalPosition{ xMin,   yMax_4, zMin   };
    const auto aniso_perp_middleUpperRight   = GlobalPosition{ xMax_2, yMax_4, zMax_2 };

    // outlet of the subdomain
    const auto aniso_main_outletLowerLeft    = GlobalPosition{ xMax_2, yMin,   zMin   };
    const auto aniso_main_outletUpperRight   = GlobalPosition{ xMax_2, yMax_2, zMax_2 };
    const auto aniso_perp_outletLowerLeft    = GlobalPosition{ xMin,   yMax_2, zMin   };
    const auto aniso_perp_outletUpperRight   = GlobalPosition{ xMax_2, yMax_2, zMax_2 };
#else
    // inlet of the subdomain
    const auto aniso_main_inletLowerLeft     = GlobalPosition{ xMin,   yMin   };
    const auto aniso_main_inletUpperRight    = GlobalPosition{ xMin,   yMax_2 };
    const auto aniso_perp_inletLowerLeft     = GlobalPosition{ xMin,   yMin   };
    const auto aniso_perp_inletUpperRight    = GlobalPosition{ xMax_2, yMin   };

    // center of the subdomain
    const auto aniso_main_middleLowerLeft    = GlobalPosition{ xMax_4, yMin  };
    const auto aniso_main_middleUpperRight   = GlobalPosition{ xMax_4, yMax_2 };
    const auto aniso_perp_middleLowerLeft    = GlobalPosition{ xMin,   yMax_4 };
    const auto aniso_perp_middleUpperRight   = GlobalPosition{ xMax_2, yMax_4 };

    // outlet of the subdomain
    const auto aniso_main_outletLowerLeft    = GlobalPosition{ xMax_2, yMin,  };
    const auto aniso_main_outletUpperRight   = GlobalPosition{ xMax_2, yMax_2 };
    const auto aniso_perp_outletLowerLeft    = GlobalPosition{ xMin,   yMax_2 };
    const auto aniso_perp_outletUpperRight   = GlobalPosition{ xMax_2, yMax_2 };

#endif 

if ( massProblem->isVerticalFlow())
    {
    // If the Flux is vertical perpendicular direction is main direction and vice versa
    flux.addAxisAlignedSurface( "aniso_main_inlet",  aniso_perp_inletLowerLeft,  aniso_perp_inletUpperRight  );
    flux.addAxisAlignedSurface( "aniso_main_middle", aniso_perp_middleLowerLeft, aniso_perp_middleUpperRight );
    flux.addAxisAlignedSurface( "aniso_main_outlet", aniso_perp_outletLowerLeft, aniso_perp_outletUpperRight );
    flux.addAxisAlignedSurface( "aniso_perp_inlet",  aniso_main_inletLowerLeft,  aniso_main_inletUpperRight  );
    flux.addAxisAlignedSurface( "aniso_perp_middle", aniso_main_middleLowerLeft, aniso_main_middleUpperRight );
    flux.addAxisAlignedSurface( "aniso_perp_outlet", aniso_main_outletLowerLeft, aniso_main_outletUpperRight );
    }
else 
    {
    flux.addAxisAlignedSurface( "aniso_main_inlet",  aniso_main_inletLowerLeft,  aniso_main_inletUpperRight  );
    flux.addAxisAlignedSurface( "aniso_main_middle", aniso_main_middleLowerLeft, aniso_main_middleUpperRight );
    flux.addAxisAlignedSurface( "aniso_main_outlet", aniso_main_outletLowerLeft, aniso_main_outletUpperRight );
    flux.addAxisAlignedSurface( "aniso_perp_inlet",  aniso_perp_inletLowerLeft,  aniso_perp_inletUpperRight  );
    flux.addAxisAlignedSurface( "aniso_perp_middle", aniso_perp_middleLowerLeft, aniso_perp_middleUpperRight );
    flux.addAxisAlignedSurface( "aniso_perp_outlet", aniso_perp_outletLowerLeft, aniso_perp_outletUpperRight );
    }


    // linearize & solve
    Dune::Timer timer;
    std::cout << "before solve"<< std::endl;
    nonLinearSolver.solve(x);
    std::cout << "after solve"<< std::endl;

    // write vtk output
    if(writeVtu)
        vtkWriter.write(1.0);

    // calculate and print mass fluxes over the planes
    const auto rho_ = getParam<Scalar>("Component.LiquidDensity");
    
    flux.calculateAllFluxes();
    std::cout << "mass flux at inlet is: " << flux.flux("inlet") << std::endl;
    std::cout << "mass flux at middle is: " << flux.flux("middle") << std::endl;
    std::cout << "mass flux at outlet is: " << flux.flux("outlet") << std::endl;
    std::cout << "volume flux at inlet is: " << flux.flux("inlet") / rho_<< std::endl;
    std::cout << "volume flux at middle is: " << flux.flux("middle") / rho_ << std::endl;
    std::cout << "volume flux at outlet is: " << flux.flux("outlet") / rho_<< std::endl;

    std::cout << "Fluxes over subdomain borders: " << std::endl;
    std::cout << "\tMass flux main direction inlet:     " << flux.flux( "aniso_main_inlet" )  << std::endl;
    std::cout << "\tMass flux main direction outlet:    " << flux.flux( "aniso_main_outlet" ) << std::endl;
    std::cout << "\tMass flux perp direction inlet:     " << flux.flux( "aniso_perp_inlet" )  << std::endl;
    std::cout << "\tMass flux perp direction outlet:    " << flux.flux( "aniso_perp_outlet" ) << std::endl;
    std::cout << "\tVolume flux main direction inlet:   " << flux.flux( "aniso_main_inlet" )  / rho_ << std::endl;
    std::cout << "\tVolume flux main direction outlet:  " << flux.flux( "aniso_main_outlet" ) / rho_ << std::endl;
    std::cout << "\tVolume flux perp direction inlet:   " << flux.flux( "aniso_perp_inlet" )  / rho_ << std::endl;
    std::cout << "\tVolume flux perp direction outlet:  " << flux.flux( "aniso_perp_outlet" ) / rho_ << std::endl;

    std::cout << "analyticalFlux: " << massProblem->analyticalFlux()*1e3 << std::endl;
    
    std::cout << "k11 is darcypermeability in main direction of the driving force!" << std::endl;
    std::cout << "darcypermeability:     " << massProblem->darcyPermFactor() * ( flux.flux( "outlet" ) / rho_ )            << std::endl; 
    std::cout << "k11_darcypermeability: " << massProblem->darcyPermFactorHalfDomain() * ( flux.flux( "aniso_main_outlet" ) / rho_ ) << std::endl; 
    std::cout << "k12_darcypermeability: " << massProblem->darcyPermFactorHalfDomain() * ( flux.flux( "aniso_perp_outlet" ) / rho_ ) << std::endl; 
    
    timer.stop();

    const auto& comm = Dune::MPIHelper::getCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
