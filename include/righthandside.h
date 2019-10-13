#pragma once
#include <deal.II/base/function.h>

namespace Geothermal
{
template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide()
        : Function<dim>(),
          period(0.2)
    {
    }
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const; //assign rhs equal to zero at each point

private:
    const double period;
};

} // namespace Geothermal
