// Copyright (C) 2017 Chris Richardson and Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>

#include <ufc.h>

#include <dolfin/fem/fem_utils.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/assemble_local.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/fem/DiscreteOperators.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/LinearVariationalProblem.h>
#include <dolfin/fem/LinearVariationalSolver.h>
#include <dolfin/fem/LocalSolver.h>
#include <dolfin/fem/NonlinearVariationalProblem.h>
#include <dolfin/fem/NonlinearVariationalSolver.h>
#include <dolfin/fem/PointSource.h>
#include <dolfin/fem/SystemAssembler.h>
#include <dolfin/fem/PETScDMCollection.h>
#include <dolfin/fem/SparsityPatternBuilder.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/la/GenericTensor.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/SparsityPattern.h>

namespace py = pybind11;

namespace dolfin_wrappers
{
  void fem(py::module& m)
  {
    // Delcare UFC objects
    py::class_<ufc::finite_element, std::shared_ptr<ufc::finite_element>>
      (m, "ufc_finite_element", "UFC finite element object");
    py::class_<ufc::dofmap, std::shared_ptr<ufc::dofmap>>
      (m, "ufc_dofmap", "UFC dofmap object");
    py::class_<ufc::form, std::shared_ptr<ufc::form>>
      (m, "ufc_form", "UFC form object");

    // Function to convert pointers (from JIT usually) to UFC objects
    m.def("make_ufc_finite_element",
          [](std::uintptr_t e)
          {
            ufc::finite_element * p = reinterpret_cast<ufc::finite_element *>(e);
            return std::shared_ptr<const ufc::finite_element>(p);
          });

    m.def("make_ufc_dofmap",
          [](std::uintptr_t e)
          {
            ufc::dofmap * p = reinterpret_cast<ufc::dofmap *>(e);
            return std::shared_ptr<const ufc::dofmap>(p);
          });

    m.def("make_ufc_form",
          [](std::uintptr_t e)
          {
            ufc::form * p = reinterpret_cast<ufc::form *>(e);
            return std::shared_ptr<const ufc::form>(p);
          });

    // dolfin::FiniteElement class
    py::class_<dolfin::FiniteElement, std::shared_ptr<dolfin::FiniteElement>>
      (m, "FiniteElement", "DOLFIN FiniteElement object")
      .def(py::init<std::shared_ptr<const ufc::finite_element>>())
      .def("num_sub_elements", &dolfin::FiniteElement::num_sub_elements)
      .def("evaluate_dofs", [](const dolfin::FiniteElement& self, py::object f,
                               py::array_t<double> coordinate_dofs,
                               int cell_orientation, const dolfin::Cell& c)
           {
             const ufc::function* _f = nullptr;
             if (py::hasattr(f, "_cpp_object"))
               _f = f.attr("_cpp_object").cast<ufc::function*>();
             else
               _f = f.cast<ufc::function*>();

             ufc::cell ufc_cell;
             c.get_cell_data(ufc_cell);
             py::array_t<double, py::array::c_style> dofs(self.space_dimension());
             self.evaluate_dofs(dofs.mutable_data(), *_f,
                                coordinate_dofs.data(), cell_orientation, ufc_cell);
             return dofs;
           }, "Evaluate degrees of freedom on element for a given function")
      .def("tabulate_dof_coordinates", [](const dolfin::FiniteElement& self, const dolfin::Cell& cell)
           {
             // Get cell vertex coordinates
             std::vector<double> coordinate_dofs;
             cell.get_coordinate_dofs(coordinate_dofs);

             // Tabulate the coordinates
             boost::multi_array<double, 2> _dof_coords;
             self.tabulate_dof_coordinates(_dof_coords, coordinate_dofs, cell);

             // Copy data and return
             typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenArray;
             EigenArray dof_coords = Eigen::Map<EigenArray>(_dof_coords.data(),
                                                            _dof_coords.shape()[0],
                                                            _dof_coords.shape()[1]);
             return dof_coords;
           }, "Tabulate coordinates of dofs on cell")
      .def("space_dimension", &dolfin::FiniteElement::space_dimension)
      .def("value_dimension", &dolfin::FiniteElement::value_dimension)
      .def("signature", &dolfin::FiniteElement::signature);

    // dolfin::GenericDofMap class
    py::class_<dolfin::GenericDofMap, std::shared_ptr<dolfin::GenericDofMap>>
      (m, "GenericDofMap", "DOLFIN DofMap object")
      .def("index_map", &dolfin::GenericDofMap::index_map)
      .def("shared_nodes", &dolfin::GenericDofMap::shared_nodes)
      .def("cell_dofs", &dolfin::GenericDofMap::cell_dofs)
      .def("dofs", (std::vector<dolfin::la_index>(dolfin::GenericDofMap::*)() const)
           &dolfin::GenericDofMap::dofs)
      .def("dofs", (std::vector<dolfin::la_index>(dolfin::GenericDofMap::*)(const dolfin::Mesh&, std::size_t) const)
           &dolfin::GenericDofMap::dofs)
      .def("entity_dofs", (std::vector<dolfin::la_index>(dolfin::GenericDofMap::*)(const dolfin::Mesh&, std::size_t) const)
           &dolfin::GenericDofMap::entity_dofs)
      .def("entity_closure_dofs", (std::vector<dolfin::la_index>(dolfin::GenericDofMap::*)(const dolfin::Mesh&, std::size_t) const)
           &dolfin::GenericDofMap::entity_closure_dofs)
      .def("entity_dofs", (std::vector<dolfin::la_index>(dolfin::GenericDofMap::*)(const dolfin::Mesh&,
                                                                                   std::size_t,
                                                                                   const std::vector<std::size_t>&) const)
           &dolfin::GenericDofMap::entity_dofs)
      .def("entity_closure_dofs", (std::vector<dolfin::la_index>(dolfin::GenericDofMap::*)(const dolfin::Mesh&,
                                                                                           std::size_t,
                                                                                           const std::vector<std::size_t>&) const)
           &dolfin::GenericDofMap::entity_closure_dofs)
      .def("num_entity_dofs", &dolfin::GenericDofMap::num_entity_dofs)
      .def("tabulate_local_to_global_dofs", &dolfin::GenericDofMap::tabulate_local_to_global_dofs)
      .def("clear_sub_map_data", &dolfin::GenericDofMap::clear_sub_map_data)
      .def("tabulate_entity_dofs", [](const dolfin::GenericDofMap& instance, std::size_t entity_dim,
                                      std::size_t cell_entity_index)
           {
             std::vector<std::size_t> dofs(instance.num_entity_dofs(entity_dim));
             instance.tabulate_entity_dofs(dofs, entity_dim, cell_entity_index);
             return py::array_t<std::size_t>(dofs.size(), dofs.data());
           })
      .def("block_size", &dolfin::GenericDofMap::block_size)
      .def("tabulate_local_to_global_dofs", [](const dolfin::GenericDofMap& instance)
           {
             std::vector<std::size_t> dofs;
             instance.tabulate_local_to_global_dofs(dofs);
             return py::array_t<std::size_t>(dofs.size(), dofs.data());
           })
      .def("set", &dolfin::GenericDofMap::set);

    // dolfin::DofMap class
    py::class_<dolfin::DofMap, std::shared_ptr<dolfin::DofMap>, dolfin::GenericDofMap>
      (m, "DofMap", "DOLFIN DofMap object")
      .def(py::init<std::shared_ptr<const ufc::dofmap>, const dolfin::Mesh&>())
      .def(py::init<std::shared_ptr<const ufc::dofmap>, const dolfin::Mesh&, std::shared_ptr<const dolfin::SubDomain>>())
      .def("ownership_range", &dolfin::DofMap::ownership_range)
      .def("cell_dofs", &dolfin::DofMap::cell_dofs);

    // dolfin::SparsityPatternBuilder
    py::class_<dolfin::SparsityPatternBuilder>(m, "SparsityPatternBuilder")
      .def_static("build", &dolfin::SparsityPatternBuilder::build,
                  py::arg("sparsity_pattern"),py::arg("mesh"),
                  py::arg("dofmaps"), py::arg("cells"),
                  py::arg("interior_facets"), py::arg("exterior_facets"),
                  py::arg("vertices"), py::arg("diagonal"),
                  py::arg("init")=true, py::arg("finalize")=true);

    // dolfin::DirichletBC class
    py::class_<dolfin::DirichletBC, std::shared_ptr<dolfin::DirichletBC>>
      (m, "DirichletBC", "DOLFIN DirichletBC object")
      .def(py::init<std::shared_ptr<const dolfin::FunctionSpace>,
           std::shared_ptr<const dolfin::GenericFunction>,
           std::shared_ptr<const dolfin::SubDomain>>())
      .def("zero", &dolfin::DirichletBC::zero)
      .def("zero_columns", &dolfin::DirichletBC::zero_columns,
           py::arg("A"), py::arg("b"), py::arg("diagonal_value")=0.0)
      .def("get_boundary_values", [](const dolfin::DirichletBC& instance)
           {
             dolfin::DirichletBC::Map map;
             instance.get_boundary_values(map);
             return map;
           })
      .def("apply", (void (dolfin::DirichletBC::*)(dolfin::GenericVector&) const)
           &dolfin::DirichletBC::apply)
      .def("apply", (void (dolfin::DirichletBC::*)(dolfin::GenericMatrix&) const)
           &dolfin::DirichletBC::apply)
      .def("user_subdomain", &dolfin::DirichletBC::user_sub_domain)
      .def("set_value", &dolfin::DirichletBC::set_value)
      .def("set_value", [](dolfin::DirichletBC& self, py::object value)
           {
             auto _u = value.attr("_cpp_object").cast<std::shared_ptr<const dolfin::GenericFunction>>();
             self.set_value(_u);
           });

    // dolfin::AssemblerBase class
    py::class_<dolfin::AssemblerBase, std::shared_ptr<dolfin::AssemblerBase>>
      (m, "AssemblerBase")
      .def_readwrite("add_values", &dolfin::Assembler::add_values)
      .def_readwrite("keep_diagonal", &dolfin::Assembler::keep_diagonal)
      .def_readwrite("finalize_tensor", &dolfin::Assembler::finalize_tensor);

    // dolfin::Assembler class
    py::class_<dolfin::Assembler, std::shared_ptr<dolfin::Assembler>, dolfin::AssemblerBase>
      (m, "Assembler", "DOLFIN Assembler object")
      .def(py::init<>())
      .def("assemble", &dolfin::Assembler::assemble);

    // dolfin::SystemAssembler class
    py::class_<dolfin::SystemAssembler, std::shared_ptr<dolfin::SystemAssembler>, dolfin::AssemblerBase>
      (m, "SystemAssembler", "DOLFIN SystemAssembler object")
      .def(py::init<std::shared_ptr<const dolfin::Form>, std::shared_ptr<const dolfin::Form>,
           std::vector<std::shared_ptr<const dolfin::DirichletBC>>>())
      .def("assemble", (void (dolfin::SystemAssembler::*)(dolfin::GenericMatrix&, dolfin::GenericVector&))
           &dolfin::SystemAssembler::assemble)
      .def("assemble", (void (dolfin::SystemAssembler::*)(dolfin::GenericMatrix&, dolfin::GenericVector&,
                                                          const dolfin::GenericVector&))
           &dolfin::SystemAssembler::assemble)
      .def("assemble", (void (dolfin::SystemAssembler::*)(dolfin::GenericMatrix&))
           &dolfin::SystemAssembler::assemble)
      .def("assemble", (void (dolfin::SystemAssembler::*)(dolfin::GenericVector&))
           &dolfin::SystemAssembler::assemble);

    // dolfin::DiscreteOperators
    py::class_<dolfin::DiscreteOperators> (m, "DiscreteOperators")
      .def_static("build_gradient", &dolfin::DiscreteOperators::build_gradient);

    // dolfin::Form class
    py::class_<dolfin::Form, std::shared_ptr<dolfin::Form>>
      (m, "Form", "DOLFIN Form object")
      .def(py::init<std::shared_ptr<const ufc::form>,
                    std::vector<std::shared_ptr<const dolfin::FunctionSpace>>>())
      .def("num_coefficients", &dolfin::Form::num_coefficients, "Return number of coefficients in form")
      .def("original_coefficient_position", &dolfin::Form::original_coefficient_position)
      .def("set_coefficient", (void (dolfin::Form::*)(std::size_t, std::shared_ptr<const dolfin::GenericFunction>))
           &dolfin::Form::set_coefficient, "Doc")
      .def("set_coefficient", (void (dolfin::Form::*)(std::string, std::shared_ptr<const dolfin::GenericFunction>))
           &dolfin::Form::set_coefficient, "Doc")
      .def("set_mesh", &dolfin::Form::set_mesh)
      .def("rank", &dolfin::Form::rank)
      .def("mesh", &dolfin::Form::mesh);

    // dolfin::PointSource class
    py::class_<dolfin::PointSource, std::shared_ptr<dolfin::PointSource>>
      (m, "PointSource")
      .def(py::init<std::shared_ptr<const dolfin::FunctionSpace>, const dolfin::Point&, double>(),
           py::arg("V"), py::arg("p"), py::arg("magnitude")=1.0)
      .def(py::init<std::shared_ptr<const dolfin::FunctionSpace>, std::shared_ptr<const dolfin::FunctionSpace>, const dolfin::Point&, double>(),
           py::arg("V0"), py::arg("V1"), py::arg("p"), py::arg("magnitude")=1.0)
      .def(py::init<std::shared_ptr<const dolfin::FunctionSpace>, const std::vector<std::pair<const dolfin::Point*, double>>>())
      .def(py::init<std::shared_ptr<const dolfin::FunctionSpace>, std::shared_ptr<const dolfin::FunctionSpace>,
           const std::vector<std::pair<const dolfin::Point*, double>>>())
      .def("apply", (void (dolfin::PointSource::*)(dolfin::GenericVector&)) &dolfin::PointSource::apply)
      .def("apply", (void (dolfin::PointSource::*)(dolfin::GenericMatrix&)) &dolfin::PointSource::apply);

    // dolfin::LinearVariationalProblem class
    py::class_<dolfin::LinearVariationalProblem,
               std::shared_ptr<dolfin::LinearVariationalProblem>>
      (m, "LinearVariationalProblem")
      .def(py::init<std::shared_ptr<const dolfin::Form>,
           std::shared_ptr<const dolfin::Form>,
           std::shared_ptr<dolfin::Function>,
           std::vector<std::shared_ptr<const dolfin::DirichletBC>>>());

    // dolfin::LinearVariationalSolver class
    py::class_<dolfin::LinearVariationalSolver,
               std::shared_ptr<dolfin::LinearVariationalSolver>,
               dolfin::Variable>(m, "LinearVariationalSolver")
      .def(py::init<std::shared_ptr<dolfin::LinearVariationalProblem>>())
      .def("solve", &dolfin::LinearVariationalSolver::solve);

    // dolfin::NonlinearVariationalProblem class
    py::class_<dolfin::NonlinearVariationalProblem,
               std::shared_ptr<dolfin::NonlinearVariationalProblem>>
      (m, "NonlinearVariationalProblem")
      .def(py::init<std::shared_ptr<const dolfin::Form>,
           std::shared_ptr<dolfin::Function>,
           std::vector<std::shared_ptr<const dolfin::DirichletBC>>,
           std::shared_ptr<const dolfin::Form>>());

    // dolfin::NonlinearVariationalSolver class
    py::class_<dolfin::NonlinearVariationalSolver,
               std::shared_ptr<dolfin::NonlinearVariationalSolver>,
               dolfin::Variable>
      (m, "NonlinearVariationalSolver")
      .def(py::init<std::shared_ptr<dolfin::NonlinearVariationalProblem>>())
      .def("solve", &dolfin::NonlinearVariationalSolver::solve);

    // dolfin::LocalSolver class
    py::class_<dolfin::LocalSolver, std::shared_ptr<dolfin::LocalSolver>>
      local_solver(m, "LocalSolver");

    py::enum_<dolfin::LocalSolver::SolverType>(local_solver, "SolverType")
      .value("LU", dolfin::LocalSolver::SolverType::LU)
      .value("Cholesky", dolfin::LocalSolver::SolverType::Cholesky);

    local_solver.def(py::init<std::shared_ptr<const dolfin::Form>,
                     std::shared_ptr<const dolfin::Form>,
                     dolfin::LocalSolver::SolverType>())
      .def("solve_local_rhs", &dolfin::LocalSolver::solve_local_rhs)
      .def("solve_local_rhs", [](dolfin::LocalSolver& self, py::object u)
           {
             auto _u = u.attr("_cpp_object").cast<dolfin::Function*>();
             self.solve_local_rhs(*_u);
           });


#ifdef HAS_PETSC
    // dolfin::PETScDMCollection
    py::class_<dolfin::PETScDMCollection, std::shared_ptr<dolfin::PETScDMCollection>>
      (m, "PETScDMCollection")
      .def_static("create_transfer_matrix", &dolfin::PETScDMCollection::create_transfer_matrix);
#endif

    // Assemble functions

    m.def("assemble", (void (*)(dolfin::GenericTensor&, const dolfin::Form&)) &dolfin::assemble);
    m.def("assemble", (double (*)(const dolfin::Form&)) &dolfin::assemble);

    m.def("assemble_system", (void (*)(dolfin::GenericMatrix&, dolfin::GenericVector&,
                                       const dolfin::Form&, const dolfin::Form&,
                                       std::vector<std::shared_ptr<const dolfin::DirichletBC>>))
          &dolfin::assemble_system);

    m.def("assemble_system", (void (*)(dolfin::GenericMatrix&, dolfin::GenericVector&,
                                       const dolfin::Form&, const dolfin::Form&,
                                       std::vector<std::shared_ptr<const dolfin::DirichletBC>>,
                                       const dolfin::GenericVector&))
          &dolfin::assemble_system);

    m.def("assemble_local",
          (Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(*)(const dolfin::Form&, const dolfin::Cell&))
          &dolfin::assemble_local);

    // FEM utils functions
    m.def("set_coordinates", [](dolfin::MeshGeometry& geometry, const py::object u)
          {
            try
            {
              auto _u = u.attr("_cpp_object").cast<const dolfin::Function*>();
              dolfin::set_coordinates(geometry, *_u);
            }
            catch (const std::runtime_error& e)
            {
              // Do nothing, pybind11 will try next function
            }
          });
    m.def("set_coordinates", &dolfin::set_coordinates);

    m.def("get_coordinates", [](py::object u, const dolfin::MeshGeometry& geometry)
          {
            try
            {
              auto _u = u.attr("_cpp_object").cast<dolfin::Function*>();
              dolfin::get_coordinates(*_u, geometry);
            }
            catch (const std::runtime_error& e)
            {
              // Do nothing, pybind11 will try next function
            }
          });
    m.def("get_coordinates", &dolfin::get_coordinates);

    m.def("vertex_to_dof_map", &dolfin::vertex_to_dof_map);
    m.def("dof_to_vertex_map", &dolfin::dof_to_vertex_map);
  }

}
