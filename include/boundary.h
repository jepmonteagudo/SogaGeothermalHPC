#pragma once
#include <deal.II/base/function.h>
#include <iostream>
#include <fstream>

template <int dim>
class BoundaryValues : public dealii::Function<dim>
{
public:
    virtual double value(const dealii::Point<dim> &p,
                         const unsigned int component = 0) const; //boundary value
};

# include <boundary.tcc>