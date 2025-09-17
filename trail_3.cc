#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/numerics/data_postprocessor.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>

namespace Step8
{
  using namespace dealii;

  // ------------------------------------------------------------------
  // ElasticProblem class: sets up and solves plane stress elasticity
  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem();
    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void output_results() const;
    void print_node_displacements() const;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;
    FESystem<dim>      fe;

    AffineConstraints<double> constraints;
    SparsityPattern          sparsity_pattern;
    SparseMatrix<double>     system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    // Material properties
    double E, nu, lambda_ps, mu_ps;
  };
  // ------------------------------------------------------------------

  template <int dim>
  ElasticProblem<dim>::ElasticProblem()
    : dof_handler(triangulation)
    , fe(FE_Q<dim>(1), dim)
    , E(30e6)
    , nu(0.25)
  {
    lambda_ps = E * nu / (1.0 - nu * nu);
    mu_ps     = E / (2.0 * (1.0 + nu));
  }

  template <int dim>
  void right_hand_side(const std::vector<Point<dim>> &points,
                       std::vector<Tensor<1, dim>>    &values)
  {
    AssertDimension(values.size(), points.size());
    for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = Tensor<1, dim>();
  }

  template <int dim>
  void ElasticProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      dof_handler, 1, Functions::ZeroFunction<dim>(dim), constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  }

  template <int dim>
  void ElasticProblem<dim>::assemble_system()
  {
    const double thickness = 0.1;
    QGauss<dim>   quad(fe.degree+1);
    FEValues<dim> fe_values(fe, quad,
                            update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quad.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<double> lambda_vals(n_q_points), mu_vals(n_q_points);
    std::vector<Tensor<1, dim>> rhs_vals(n_q_points);

    for (auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      Functions::ConstantFunction<dim> lam(lambda_ps), mu(mu_ps);
      lam.value_list(fe_values.get_quadrature_points(), lambda_vals);
      mu.value_list(fe_values.get_quadrature_points(), mu_vals);
      right_hand_side(fe_values.get_quadrature_points(), rhs_vals);

      for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int ci = fe.system_to_component_index(i).first;
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const unsigned int cj = fe.system_to_component_index(j).first;
            cell_matrix(i,j) += (
              lambda_vals[q] * fe_values.shape_grad(i,q)[ci] * fe_values.shape_grad(j,q)[cj]
              + mu_vals[q] * fe_values.shape_grad(i,q)[cj] * fe_values.shape_grad(j,q)[ci]
              + (ci==cj ? mu_vals[q] * (fe_values.shape_grad(i,q) * fe_values.shape_grad(j,q)) : 0.0)
            ) * fe_values.JxW(q) * thickness;
          }
          cell_rhs(i) += fe_values.shape_value(i,q) * rhs_vals[q][ci] * fe_values.JxW(q) * thickness;
        }
      }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }

    // Neumann on boundary id=2
    QGauss<dim-1> face_quad(fe.degree+1);
    FEFaceValues<dim> fe_face(fe, face_quad,
                               update_values | update_normal_vectors | update_JxW_values);
    const double tx=0.0, ty=15.0;
    for (auto &cell : dof_handler.active_cell_iterators())
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id()==2)
        {
          fe_face.reinit(cell,f);
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int q=0; q<face_quad.size(); ++q)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int ci = fe.system_to_component_index(i).first;
              system_rhs(local_dof_indices[i]) += fe_face.shape_value(i,q)
                * (ci==0 ? tx : ty) * fe_face.JxW(q);
            }
        }
    constraints.distribute(system_rhs);
  }

  template <int dim>
  void ElasticProblem<dim>::solve()
  {
    SolverControl sc(100000, 1e-8 * system_rhs.l2_norm());
    SolverCG<Vector<double>> solver(sc);
    PreconditionSSOR<SparseMatrix<double>> pre;
    pre.initialize(system_matrix, 1.2);
    solver.solve(system_matrix, solution, system_rhs, pre);
    constraints.distribute(solution);
    std::cout<<"Converged in "<<sc.last_step()<<" steps, res="<<sc.last_value()<<std::endl;
  }

  // ------------------------------------------------------------------
  // StressPostprocessor class: computes stress and strain energy density
  template <int dim>
  class StressPostprocessor : public DataPostprocessor<dim>
  {
  public:
    StressPostprocessor(const double E, const double nu);

    virtual void evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &inputs,
      std::vector<Vector<double>> &computed_quantities) const override;

    virtual std::vector<std::string> get_names() const override
    {
      return {"sigma_xx","sigma_yy","sigma_xy","strain_energy_density"};
    }

    virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override
    {
      return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
        4, DataComponentInterpretation::component_is_scalar);
    }

    virtual UpdateFlags get_needed_update_flags() const override
    {
      return update_gradients;
    }

  private:
    const double E, nu, factor;
  };

  template <int dim>
  //constructor for the StressPostprocessor class
  StressPostprocessor<dim>::StressPostprocessor(const double E_, const double nu_)
    : DataPostprocessor<dim>(), E(E_), nu(nu_), factor(E_ / (1 - nu_*nu_)) {}

  template <int dim>
  void StressPostprocessor<dim>::evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const
  {
    for (unsigned int p=0; p<inputs.solution_gradients.size(); ++p)
    {
      Tensor<2, dim> grad;
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          grad[i][j] = inputs.solution_gradients[p][i][j];

      SymmetricTensor<2, dim> eps;
      for (unsigned int i=0;i<dim;++i)
        for (unsigned int j=0;j<dim;++j)
          eps[i][j] = 0.5*(grad[i][j]+grad[j][i]);
      SymmetricTensor<2, dim> sig;
      sig[0][0] = factor*(eps[0][0]+nu*eps[1][1]);
      sig[1][1] = factor*(eps[1][1]+nu*eps[0][0]);
      sig[0][1] = factor*(1-nu)*eps[0][1];
      double W = 0.5 * (sig * eps);
      computed_quantities[p][0] = sig[0][0];
      computed_quantities[p][1] = sig[1][1];
      computed_quantities[p][2] = sig[0][1];
      computed_quantities[p][3] = W;
    }
  }

  template <int dim>
  void ElasticProblem<dim>::output_results() const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    // displacement vector
    std::vector<std::string> disp_names(dim);
    for (unsigned int i=0;i<dim;++i)
      disp_names[i] = "disp_"+std::to_string(i);
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      disp_interp(dim, DataComponentInterpretation::component_is_part_of_vector);
    data_out.add_data_vector(solution, disp_names,
                             DataOut<dim>::type_dof_data,
                             disp_interp);

    // stress and energy
    StressPostprocessor<dim> stress_post(E, nu);
    data_out.add_data_vector(solution, stress_post);

    data_out.build_patches();
    std::ofstream out("solution_with_stress.vtk");
    data_out.write_vtk(out);
  }

  template <int dim>
  void ElasticProblem<dim>::print_node_displacements() const
  {
    const double eps=1e-8;
    std::vector<Point<dim>> corners = {Point<dim>(0,0),Point<dim>(12,0),Point<dim>(0,3),Point<dim>(12,3)};
    std::cout<<std::fixed<<std::setprecision(6)<<"Corner displacements:\n";
    for (auto &pt:corners)
    {
      bool found=false;
      for (auto &cell:dof_handler.active_cell_iterators())
        for (unsigned int v=0;v<GeometryInfo<dim>::vertices_per_cell;++v)
          if ((cell->vertex(v)-pt).norm()<eps)
          {
            std::cout<<pt<<": ";
            for (unsigned int d=0;d<dim;++d)
              std::cout<<solution[cell->vertex_dof_index(v,d)]<<" ";
            std::cout<<"\n";
            found=true; break;
          }
      if (!found) std::cout<<pt<<" not on mesh\n";
    }
  }

  template <int dim>
  void ElasticProblem<dim>::run()
  {
    GridGenerator::subdivided_hyper_rectangle(triangulation,{10,2},Point<dim>(0,0),Point<dim>(12,3));
    for (auto &cell:triangulation.active_cell_iterators())
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell;++f)
        if (cell->face(f)->at_boundary())
        {
          auto c=cell->face(f)->center();
          if (std::abs(c[0])<1e-12) cell->face(f)->set_boundary_id(2);
          else if (std::abs(c[0]-12)<1e-12) cell->face(f)->set_boundary_id(1);
          else cell->face(f)->set_boundary_id(0);
        }
    setup_system();
    assemble_system();
    solve();
    output_results();
    print_node_displacements();
  }
} // namespace Step8

int main()
{
  try
  {
    Step8::ElasticProblem<2> problem;
    problem.run();
  }
  catch (std::exception &exc)
  {
    std::cerr<<"Exception: "<<exc.what()<<std::endl;
    return 1;
  }   
  return 0;
}
