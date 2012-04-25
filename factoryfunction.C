

// libMesh includes
#include "dense_vector.h"
#include "factory.h"
#include "function_base.h"

// An example of a hand-coded function

class ExampleOneFunction : public FunctionBase<Number>
{
  virtual void init () { };

  /**
   * Clears the function.
   */
  virtual void clear () { };

  /**
   * Returns a new copy of the function.  The new copy should be as
   * ``deep'' as necessary to allow independent destruction and
   * simultaneous evaluations of the copies in different threads.
   */
  AutoPtr<FunctionBase<Number> > clone () { return (AutoPtr<FunctionBase<Number> >)NULL; };

  virtual Number operator() (const Point& p,
                             const Real time = 0)
    {
      return 1;
    }

  virtual void operator() (const Point& p,
                           const Real time,
                           DenseVector<Number>& output)
    {
      for (unsigned int i=0; i != output.size(); ++i)
        output(i) = 1;
    }

  virtual AutoPtr<FunctionBase<Number> > clone() const {
    return AutoPtr<FunctionBase<Number> >
      (new ExampleOneFunction());
  }
};

#ifdef LIBMESH_USE_SEPARATE_NAMESPACE
namespace libMesh {
#endif

//-------------------------------------------------
// Full specialization for the Factory<FunctionBase<Number> >
// So we can look up hand-coded functions by name string
template<>
std::map<std::string, Factory<FunctionBase<Number> >*>&
Factory<FunctionBase<Number> >::factory_map()
{
  static std::map<std::string, Factory<FunctionBase<Number> >*> _map;
  return _map;
}

FactoryImp<ExampleOneFunction, FunctionBase<Number> > example_one_factory ("example_one");

#ifdef LIBMESH_USE_SEPARATE_NAMESPACE
} // namespace libMesh
#endif

