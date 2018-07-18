#ifndef DUNE_FIRST_HH
#define DUNE_FIRST_HH

template<typename GV>
class ModelProblem
  : public Dune::PDELab::LinearElasticityParameterInterface<
  Dune::PDELab::LinearElasticityParameterTraits<GV, double>,
  ModelProblem<GV> >
{
public:

  typedef Dune::PDELab::LinearElasticityParameterTraits<GV, double> Traits;

  ModelProblem(typename Traits::RangeType G,
    typename Traits::RangeFieldType l,
    typename Traits::RangeFieldType m) :
    G_(G), lambda_(l), mu_(m)
  {}

  void
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
    typename Traits::RangeType & y) const
  {
    y = G_;
  }

  template<typename I>
  bool isDirichlet(const I & ig,
    const typename Traits::IntersectionDomainType & coord
    ) const
  {
    typename Traits::DomainType xg = ig.geometry().global( coord );
    return (xg[0] == 0.0 );//|| xg[0] == 10.0);  // Dirichlet b.c. on left & right boundary
  }

  void
  u (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
    typename Traits::RangeType & y) const
  {
    y = 0.0;
  }

  template<typename I>
  void
  g(const I& ig, const typename Traits::IntersectionDomainType& coord,
    typename Traits::RangeType& y, auto& dim) const
  {
    // y[dim-1] = 25.0;
    typename Traits::DomainType xg = ig.geometry().global(coord);
    if (xg[0] == 10.0)// && -2.0<xg[1]-5.0<2.0 && xg[2] == 2.5)
    {
      y[dim-1] = 3.0;
    } else {
      y = 0.0;
    }
  }

  typename Traits::RangeFieldType
  lambda (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return lambda_;
  }

  typename Traits::RangeFieldType
  mu (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return mu_;
  }

private:

  typename Traits::RangeType G_;
  typename Traits::RangeFieldType lambda_;
  typename Traits::RangeFieldType mu_;

};

template <class GV>
class DeformationFunction: public Dune::DiscreteCoordFunction< double, GV::dimension, DeformationFunction<GV> >
{
    static const int dim = GV::dimension;
    typedef DeformationFunction<GV> This;
    typedef Dune::DiscreteCoordFunction< double, GV::dimension, This > Base;

    typedef typename GV::template Codim<dim>::Entity HostVertex;
    typedef typename GV::template Codim<0>::Entity HostEntity;

public:

    DeformationFunction(const GV& gv_): gv(gv_)
    {
      positions.resize(gv.size(dim));
      const auto& index_set = gv.indexSet();
      for(const auto& vertex: vertices(gv))
      {
        auto geo = vertex.geometry();
        auto global = geo.corner(0);
        int vertex_id = index_set.index(vertex);
        positions[vertex_id] = global;
      }
    }

  // called with hostEntity, where hostEntity is a vertex
  void evaluate ( const HostVertex& hostVertex,unsigned int corner,
    Dune::FieldVector<double,dim>& y ) const
  {
    const auto& index_set = gv.indexSet();
    int corner_id = index_set.index(hostVertex);
    y = positions[corner_id];
  }


  // called with hostEntity, where hostEntity is a face
  void evaluate (const HostEntity& hostEntity, unsigned int corner,
    Dune::FieldVector<double,dim>& y) const
  {
    const auto& index_set = gv.indexSet();
    int corner_id = index_set.subIndex(hostEntity,corner,dim);
    y = positions[corner_id];
  }

  template <class V>
  void move(const V& x)
  {
    const auto& index_set = gv.indexSet();
    int num_corners = gv.size(dim);
    for(const auto& vertex: vertices(gv))
    {
      auto geo = vertex.geometry();
      int corner_id = index_set.index(vertex);
      auto& view = positions[corner_id];
      view[0] += x[corner_id];
      view[1] += x[corner_id + num_corners];
    }
  }

private:
  const GV& gv;
  std::vector< Dune::FieldVector< double,dim > > positions;
};


#endif // DUNE_FIRST_HH
