#include <deal.II/base/function.h>

using namespace dealii;
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const; //boundary value
};
