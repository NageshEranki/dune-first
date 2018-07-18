// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <array>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid/grid.hh>

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/linearelasticity.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/function/callableadapter.hh>

#include <dune/first/first.hh>
#include "driver.hh"

int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    const double x1 = 10.0;
    const double y1 = 1.0;

    double mu = 100.0;
    double lambda = 10000.0;
    double g = 0.0;

    Dune::FieldVector<double,dim> L({x1,y1});

    typedef Dune::YaspGrid<dim> Grid;
    Grid grid(L,{10,1});

    grid.globalRefine(0);
    typedef Grid::LeafGridView GV;
    GV gv = grid.leafGridView();

    typedef DeformationFunction<GV> Defct;
    Defct defct(gv);

    typedef Dune::GeometryGrid<Grid,Defct> GeoGrid;
    GeoGrid geogrid(grid,defct);

    typedef GeoGrid::LeafGridView GeoGV;
    GeoGV geogv = geogrid.leafGridView();

    testp1<GeoGV,Defct>(geogv,mu,lambda,g,defct);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
