#include "heateq.h"
#include "boundary.h"
#include "righthandside.h"
#include <iostream>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/fe/fe_values.h>

// #include <deal.II/grid/tria_accessor.h>
// #include <deal.II/grid/tria_iterator.h>
// #include <deal.II/grid/grid_generator.h>
// #include <deal.II/grid/grid_tools.h>
// #include <deal.II/grid/grid_out.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/numerics/data_out.h>

void test::func()
{
    std::cout << "good" << std::endl;
}

template <int dim>
HeatEquation<dim>::HeatEquation() // initialization
    : fe(1),
      dof_handler(triangulation),
      time(0.0),
      time_step(1. / 20), //a time step constant at 1/500 (remember that one period of the source on the right hand side was set to 0.2 above,
                          //so we resolve each period with 100 time steps)
      timestep_number(0),
      theta(0.5)
{
}

template <int dim>
void HeatEquation<dim>::grid_input()
{
    dealii::GridIn<dim> gridin;
    gridin.attach_triangulation(triangulation);
    std::ifstream f("mesh.msh");
    gridin.read_msh(f);
}

template <int dim>
void HeatEquation<dim>::setup_system()
{
    dof_handler.distribute_dofs(fe); // distribute dofs to grid

    std::cout << std::endl
              << "==========================================="
              << std::endl
              << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl
              << std::endl;

    dealii::DynamicSparsityPattern dsp(dof_handler.n_dofs()); // sparsity
    dealii::DoFTools::make_sparsity_pattern(dof_handler,
                                            dsp);
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);    // initialize M using given sparsity parttern
    laplace_matrix.reinit(sparsity_pattern); //initialize A using given sparsity parttern
    system_matrix.reinit(sparsity_pattern);  // initialize M + k*theta*A using given sparsity parttern

    dealii::MatrixCreator::create_mass_matrix(dof_handler,
                                              dealii::QGauss<dim>(fe.degree + 1),
                                              mass_matrix); //Assemble the mass matrix and a right hand side vector.
                                                            //If no coefficient is given (i.e., if the pointer to a function object is zero as it is by default),
                                                            //the coefficient is taken as being constant and equal to one.
    dealii::MatrixCreator::create_laplace_matrix(dof_handler,
                                                 dealii::QGauss<dim>(fe.degree + 1),
                                                 laplace_matrix); //Assemble the Laplace matrix.
                                                                  //If no coefficient is given (i.e., if the pointer to a function object is zero as it is by default),
                                                                  //the coefficient is taken as being constant and equal to one.
                                                                  //In case you want to specify constraints and use the default argument for the coefficient you have to specify the (unused) coefficient argument as (const Function<spacedim> *const)nullptr.

    solution.reinit(dof_handler.n_dofs());
    old_solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim>
void HeatEquation<dim>::assemble_system()
{
    dealii::Vector<double> tmp;           // this vector is for
    dealii::Vector<double> forcing_terms; // this vector is for forcing_terms

    tmp.reinit(solution.size());           //initialize tmp
    forcing_terms.reinit(solution.size()); // initialize forcing_terms
    
    mass_matrix.vmult(system_rhs, old_solution); // matrix multiplication system_rhs = mass_matrix*old_solution

    laplace_matrix.vmult(tmp, old_solution);       // tmp = laplace_matrix*old_solution
    system_rhs.add(-(1 - theta) * time_step, tmp); // system_rhs = system_rhs -(1 - theta) * time_step*tmp  注意，这里system_rhs是一个vector，所以add是将两个元素相乘了

    RightHandSide<dim> rhs_function;
    rhs_function.set_time(time);
    dealii::VectorTools::create_right_hand_side(dof_handler,
                                                    dealii::QGauss<dim>(fe.degree + 1),
                                                    rhs_function,
                                                    tmp);
    forcing_terms = tmp;
    forcing_terms *= time_step * theta;

    rhs_function.set_time(time - time_step);
    dealii::VectorTools::create_right_hand_side(dof_handler,
                                                    dealii::QGauss<dim>(fe.degree + 1),
                                                    rhs_function,
                                                    tmp);

    forcing_terms.add(time_step * (1 - theta), tmp); // 形成forcing term = f(x,n-1)*(1-theta) + f(x,n)*theta

    system_rhs += forcing_terms; // system_rhs = system_rhs + forcing_term : sys_Old*U_Old + f(x,n-1)*(1-theta) + f(x,n)*theta

    system_matrix.copy_from(mass_matrix);
    system_matrix.add(theta * time_step, laplace_matrix); //sys = M + k*theta*A

    BoundaryValues<dim> boundary_values_function; // creat boundary value object
    boundary_values_function.set_time(time);      //set the proper time

    std::map<dealii::types::global_dof_index, double> boundary_values;
    dealii::VectorTools::interpolate_boundary_values(dof_handler, //evaluate value by interpolation
                                                         1,
                                                         boundary_values_function,
                                                         boundary_values);

    dealii::MatrixTools::apply_boundary_values(boundary_values,
                                                   system_matrix,
                                                   solution,
                                                   system_rhs);
}

template <int dim>
void HeatEquation<dim>::solve_time_step()
{
    dealii::SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm()); // setting for cg
    dealii::SolverCG<> cg(solver_control);                                   // config cg

    dealii::PreconditionSSOR<> preconditioner;     // precond
    preconditioner.initialize(system_matrix, 1.0); //initialize precond

    cg.solve(system_matrix, solution, system_rhs,
             preconditioner); // solve eq

    std::cout << "     " << solver_control.last_step()
              << " CG iterations." << std::endl;
}

// @sect4{<code>HeatEquation::output_results</code>}
//
// Neither is there anything new in generating graphical output:
template <int dim>
void HeatEquation<dim>::output_results() const
{
    dealii::DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "U");

    data_out.build_patches();

    const std::string filename = "solution-" + dealii::Utilities::int_to_string(timestep_number, 3) +
                                 ".vtk";
    std::ofstream output(filename.c_str());
    data_out.write_vtk(output);
}

template <int dim>
void HeatEquation<dim>::run()
{
    grid_input();
    setup_system();

    dealii::VectorTools::interpolate(dof_handler,
                                     dealii::ZeroFunction<dim>(),
                                     old_solution); // interpolate the old solution based on dof_handler, here using interpolation because we refine the global
    solution = old_solution;                        // updating the solutin with sinterpolated old solution

    output_results(); // output

    while (time <= 1.)
    {
        time += time_step; // get new time point
        ++timestep_number; // get the No. of the new time step

        std::cout << "Time step " << timestep_number << " at t=" << time
                  << std::endl;

        assemble_system();

        solve_time_step();

        output_results();

        old_solution = solution;
    }
}
