#include "righthandside.h"

namespace Geothermal
{


template <int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
                                 const unsigned int component) const
{
    (void)component;
    Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
    Assert(dim == 3, ExcNotImplemented());

    const double time = this->get_time(); //get time
    return 0;
}

} // namespace Geothermal