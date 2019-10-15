
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>
// #include <deal.II/dofs/dof_accessor.h>
// #include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

using namespace dealii;

template <int dim>
class HeatEquation
{
public:
    HeatEquation();
    void run();

private:
    void grid_input();
    void setup_system();
    void solve_time_step();
    void output_results() const;

    Triangulation<dim> triangulation; //grid
    FE_Q<dim> fe;                     //element
    DoFHandler<dim> dof_handler;      //grid<->eleemnt

    SparsityPattern sparsity_pattern;    // sparsity
    SparseMatrix<double> mass_matrix;    // M
    SparseMatrix<double> laplace_matrix; //A
    SparseMatrix<double> system_matrix;  //M + k*theta*A

    Vector<double> solution;     // solution at n
    Vector<double> old_solution; //solution at n-1
    Vector<double> system_rhs;   //rhs

    double time;
    double time_step;
    unsigned int timestep_number;

    const double theta;
};
