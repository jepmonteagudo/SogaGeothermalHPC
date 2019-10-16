#pragma once
#include "righthandside.h"

template <int dim>
double RightHandSide<dim>::value(const dealii::Point<dim> &p,
                                 const unsigned int component) const
{
    (void)component;
    Assert(component == 0, dealii::ExcIndexRange(component, 0, 1)); // for debug
    Assert(dim == 3, dealii::ExcNotImplemented());

    const double time = this->get_time(); //get time
    return 0;
}
