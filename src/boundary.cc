
#include <deal.II/base/function.h>
#include "boundary.h"

namespace Geothermal
{

template <int dim>
double BoundaryValues<dim>::value(const dealii::Point<dim> &p,
                                  const unsigned int /*component*/) const
{
  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1));
  const double time = this->get_time();
  return 10. * sin(time * 3.1415926); // boundary value is set to zero in this case
}

}