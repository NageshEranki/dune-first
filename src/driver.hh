template<class GV, class defunct>
void testp1 (const GV& gv, double mu, double lambda, double constG, defunct& defct )
{
  using Dune::PDELab::Backend::native;

  typedef typename GV::Grid::ctype DF;

  const int dim = GV::dimension;

  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
  FEM fem(gv);

  // make function space
  typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
  typedef Dune::PDELab::ISTL::VectorBackend<> ComponentVectorBackend;

  typedef Dune::PDELab::DefaultLeafOrderingTag Mapper;

  //Change to enti
  typedef Dune::PDELab::LexicographicOrderingTag OrderingTag;
  typedef Dune::PDELab::VectorGridFunctionSpace<
    GV,
    FEM,
    dim,
    Dune::PDELab::ISTL::VectorBackend<>,
    ComponentVectorBackend,
    Constraints,
    OrderingTag,
    Mapper
    > GFS;
  GFS gfs(gv,fem);
  gfs.name("displacement");

  // model description
  typedef ModelProblem<GV> Param;
  Dune::FieldVector<double, dim> G(0.0); G[dim-1] = -constG;
  Param param(G, mu, lambda);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<double>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(param,gfs,cg);

  std::cout << gfs.size() << " DOFs\n";
  std::cout << cg.size() << " constraint DOFs\n";

  // make local operator
  typedef Dune::PDELab::LinearElasticity<Param> LO;
  LO lo(param);

  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // 2D Q1

  // make grid operator
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,LO,
    MBE,
    double,double,double,
    C,C> GOS;
  GOS gos(gfs,cg,gfs,cg,lo,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GOS::Traits::Domain V;
  V x0(gfs);
  x0 = 0.0;
  typedef Dune::PDELab::LinearElasticityDirichletExtensionAdapter<Param> Displacement;
  Displacement u_fnkt(gv, param);
  Dune::PDELab::interpolate(u_fnkt,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  // represent operator as a matrix
  typedef typename GOS::Traits::Jacobian M;
  M m(gos);
  m = 0.0;
  gos.jacobian(x0,m);
  if (gfs.size() <= 16)
  {
    Dune::printmatrix(std::cout,native(m),"global stiffness matrix","row",9,3);
  }

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  gos.residual(x0,r);

  // make ISTL solver
  typedef typename M::Container ISTL_M;
  typedef typename V::Container ISTL_V;
  Dune::MatrixAdapter<ISTL_M,ISTL_V,ISTL_V> opa(native(m));
  Dune::SeqILU0<ISTL_M,ISTL_V,ISTL_V> ilu0(native(m),1e-2);

  Dune::CGSolver<ISTL_V> solver(opa,ilu0,1E-20,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  V x(gfs,0.0);
  Dune::PDELab::set_nonconstrained_dofs(cg,1.0,x);
  solver.apply(native(x),native(r),stat);
  x += x0;
  Dune::printvector(std::cout, native(x),"Solution Vector","row");

  typedef Dune::PDELab::VectorDiscreteGridFunction<GFS,V> XDGF;
  XDGF xdgf(gfs,x);

  defct.move(native(x));

  auto stationaryVTKWriter_u = std::make_shared<Dune::SubsamplingVTKWriter<GV>>(gv,0);

  std::string basename_u;
  basename_u =  "output";
  basename_u += "_u";

  Dune::VTKSequenceWriter<GV> vtkwriter_u(stationaryVTKWriter_u,basename_u,"","");
  vtkwriter_u.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<XDGF> >(xdgf,"x_h"));

  double time = 0.0;

  vtkwriter_u.write(time,Dune::VTK::appendedraw);

  // for(const auto& x: X)
  // {
  //   std::cout << x << '\n';
  // }
  // vtkwriter_u.write(time,Dune::VTK::appendedraw);
}
