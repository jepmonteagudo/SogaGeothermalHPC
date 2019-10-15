#include <deal.II/base/function.h>

template <int dim>
class RightHandSide : public dealii::Function<dim>
{
public:
    RightHandSide()
        : dealii::Function<dim>(),
          period(0.2)
    {
    }
    virtual double value(const dealii::Point<dim> &p,
                         const unsigned int component = 0) const; //assign rhs equal to zero at each point

private:
    const double period;
};
