#pragma once
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <iostream>
#include <fstream>

class test
{
public:
    void func();
};

template <int dim>
class HeatEquation
{
public:
    HeatEquation();
    void run();

private:
    void grid_input();
    void setup_system();
    void assemble_system();
    void solve_time_step();
    void output_results() const;

    dealii::Triangulation<dim> triangulation; //grid
    dealii::FE_Q<dim> fe;                     //element
    dealii::DoFHandler<dim> dof_handler;      //grid<->eleemnt

    dealii::SparsityPattern sparsity_pattern;    // sparsity
    dealii::SparseMatrix<double> mass_matrix;    // M
    dealii::SparseMatrix<double> laplace_matrix; //A
    dealii::SparseMatrix<double> system_matrix;  //M + k*theta*A

    dealii::Vector<double> solution;     // solution at n
    dealii::Vector<double> old_solution; //solution at n-1
    dealii::Vector<double> system_rhs;   //rhs

    double time;
    double time_step;
    unsigned int timestep_number;

    const double theta;
};

#include <heateq.tcc>
