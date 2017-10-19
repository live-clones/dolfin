DOLFIN pybind11 API differences
===============================

DOLFIN's Pyton bindings are being ported from SWIG to pybind11. The new Python wrappers are not complete, and some changes have been done to the API to clean it up and bring it more in line with the underlying DOLFIN C++ code.

Differences
-----------

.. |todo| replace:: TODO, not present in pybind11
.. |new| replace:: NEW, not present in SWIG
.. |wontfix| replace:: WONTFIX, will not be included in pybind11
.. |fixed| replace:: FIXED, pybind11 is updated to include this

.. list-table::
   :header-rows: 1

   * - Name
     - Status
   * - dolfin.AbstractCell
     - |todo|
   * - dolfin.AbstractCell.geometric_dimension
     - |todo|
   * - dolfin.AbstractCell.has_simplex_facets
     - |todo|
   * - dolfin.AbstractCell.is_simplex
     - |todo|
   * - dolfin.AbstractCell.topological_dimension
     - |todo|
   * - dolfin.AbstractDomain
     - |todo|
   * - dolfin.AbstractDomain.geometric_dimension
     - |todo|
   * - dolfin.AbstractDomain.topological_dimension
     - |todo|
   * - dolfin.AdaptiveLinearVariationalSolver.adapt_problem
     - |todo|
   * - dolfin.AdaptiveLinearVariationalSolver.adaptive_data
     - |todo|
   * - dolfin.AdaptiveLinearVariationalSolver.default_parameters
     - |todo|
   * - dolfin.AdaptiveLinearVariationalSolver.evaluate_goal
     - |todo|
   * - dolfin.AdaptiveLinearVariationalSolver.extract_bcs
     - |todo|
   * - dolfin.AdaptiveLinearVariationalSolver.solve_primal
     - |todo|
   * - dolfin.AdaptiveLinearVariationalSolver.str
     - |todo|
   * - dolfin.AdaptiveNonlinearVariationalSolver.adapt_problem
     - |todo|
   * - dolfin.AdaptiveNonlinearVariationalSolver.adaptive_data
     - |todo|
   * - dolfin.AdaptiveNonlinearVariationalSolver.default_parameters
     - |todo|
   * - dolfin.AdaptiveNonlinearVariationalSolver.evaluate_goal
     - |todo|
   * - dolfin.AdaptiveNonlinearVariationalSolver.extract_bcs
     - |todo|
   * - dolfin.AdaptiveNonlinearVariationalSolver.solve_primal
     - |todo|
   * - dolfin.AdaptiveNonlinearVariationalSolver.str
     - |todo|
   * - dolfin.And
     - |todo|
   * - dolfin.Argument
     - |todo|
   * - dolfin.Argument.T
     - |todo|
   * - dolfin.Argument.dx
     - |todo|
   * - dolfin.Argument.evaluate
     - |todo|
   * - dolfin.Argument.function_space
     - |todo|
   * - dolfin.Argument.geometric_dimension
     - |todo|
   * - dolfin.Argument.is_cellwise_constant
     - |todo|
   * - dolfin.Argument.number
     - |todo|
   * - dolfin.Argument.part
     - |todo|
   * - dolfin.Argument.ufl_disable_profiling
     - |todo|
   * - dolfin.Argument.ufl_domain
     - |todo|
   * - dolfin.Argument.ufl_domains
     - |todo|
   * - dolfin.Argument.ufl_element
     - |todo|
   * - dolfin.Argument.ufl_enable_profiling
     - |todo|
   * - dolfin.Argument.ufl_free_indices
     - |todo|
   * - dolfin.Argument.ufl_function_space
     - |todo|
   * - dolfin.Argument.ufl_index_dimensions
     - |todo|
   * - dolfin.Argument.ufl_operands
     - |todo|
   * - dolfin.Argument.ufl_shape
     - |todo|
   * - dolfin.Arguments
     - |todo|
   * - dolfin.Assembler.assemble_cells
     - |todo|
   * - dolfin.Assembler.assemble_exterior_facets
     - |todo|
   * - dolfin.Assembler.assemble_interior_facets
     - |todo|
   * - dolfin.Assembler.assemble_vertices
     - |todo|
   * - dolfin.Assembler.init_global_tensor
     - |todo|
   * - dolfin.AssemblerBase
     - |todo|
   * - dolfin.AssemblerBase.add_values
     - |todo|
   * - dolfin.AssemblerBase.finalize_tensor
     - |todo|
   * - dolfin.AssemblerBase.init_global_tensor
     - |todo|
   * - dolfin.AssemblerBase.keep_diagonal
     - |todo|
   * - dolfin.AutoSubDomain.geometric_dimension
     - |todo|
   * - dolfin.AutoSubDomain.map_tolerance
     - |todo|
   * - dolfin.AutoSubDomain.snap
     - |todo|
   * - dolfin.BDF1
     - |todo|
   * - dolfin.BDF1.bcs
     - |todo|
   * - dolfin.BDF1.dolfin_stage_forms
     - |todo|
   * - dolfin.BDF1.dt
     - |todo|
   * - dolfin.BDF1.dt_stage_offset
     - |todo|
   * - dolfin.BDF1.id
     - |todo|
   * - dolfin.BDF1.implicit
     - |todo|
   * - dolfin.BDF1.jacobian_index
     - |todo|
   * - dolfin.BDF1.label
     - |todo|
   * - dolfin.BDF1.last_stage
     - |todo|
   * - dolfin.BDF1.name
     - |todo|
   * - dolfin.BDF1.order
     - |todo|
   * - dolfin.BDF1.parameters
     - |todo|
   * - dolfin.BDF1.rename
     - |todo|
   * - dolfin.BDF1.rhs_form
     - |todo|
   * - dolfin.BDF1.solution
     - |todo|
   * - dolfin.BDF1.stage_solutions
     - |todo|
   * - dolfin.BDF1.str
     - |todo|
   * - dolfin.BDF1.t
     - |todo|
   * - dolfin.BDF1.to_adm
     - |todo|
   * - dolfin.BDF1.to_tlm
     - |todo|
   * - dolfin.BDF1.ufl_stage_forms
     - |todo|
   * - dolfin.BackwardEuler.bcs
     - |todo|
   * - dolfin.BackwardEuler.dt_stage_offset
     - |todo|
   * - dolfin.BackwardEuler.id
     - |todo|
   * - dolfin.BackwardEuler.implicit
     - |todo|
   * - dolfin.BackwardEuler.jacobian_index
     - |todo|
   * - dolfin.BackwardEuler.label
     - |todo|
   * - dolfin.BackwardEuler.name
     - |todo|
   * - dolfin.BackwardEuler.parameters
     - |todo|
   * - dolfin.BackwardEuler.rename
     - |todo|
   * - dolfin.BackwardEuler.str
     - |todo|
   * - dolfin.BasisFunction
     - |todo|
   * - dolfin.BasisFunction.eval
     - |todo|
   * - dolfin.BasisFunction.eval_derivatives
     - |todo|
   * - dolfin.BasisFunction.update_index
     - |todo|
   * - dolfin.BlockMatrix.apply
     - |todo|
   * - dolfin.BlockMatrix.get_block
     - |todo|
   * - dolfin.BlockMatrix.schur_approximation
     - |todo|
   * - dolfin.BlockMatrix.set_block
     - |todo|
   * - dolfin.BlockMatrix.size
     - |todo|
   * - dolfin.BlockMatrix.str
     - |todo|
   * - dolfin.BlockMatrix.zero
     - |todo|
   * - dolfin.BlockVector.axpy
     - |todo|
   * - dolfin.BlockVector.copy
     - |todo|
   * - dolfin.BlockVector.empty
     - |todo|
   * - dolfin.BlockVector.get_block
     - |todo|
   * - dolfin.BlockVector.inner
     - |todo|
   * - dolfin.BlockVector.max
     - |todo|
   * - dolfin.BlockVector.min
     - |todo|
   * - dolfin.BlockVector.norm
     - |todo|
   * - dolfin.BlockVector.set_block
     - |todo|
   * - dolfin.BlockVector.size
     - |todo|
   * - dolfin.BlockVector.str
     - |todo|
   * - dolfin.BoostGraphOrdering
     - |todo|
   * - dolfin.BoostGraphOrdering.compute_cuthill_mckee
     - |todo|
   * - dolfin.BoundaryMesh.cell_name
     - |new|
   * - dolfin.BoundaryMesh.child
     - |todo|
   * - dolfin.BoundaryMesh.clean
     - |todo|
   * - dolfin.BoundaryMesh.clear_child
     - |todo|
   * - dolfin.BoundaryMesh.depth
     - |todo|
   * - dolfin.BoundaryMesh.entity_map
     - |todo|
   * - dolfin.BoundaryMesh.geometric_dimension
     - |new|
   * - dolfin.BoundaryMesh.ghost_mode
     - |todo|
   * - dolfin.BoundaryMesh.has_child
     - |todo|
   * - dolfin.BoundaryMesh.has_parent
     - |todo|
   * - dolfin.BoundaryMesh.label
     - |todo|
   * - dolfin.BoundaryMesh.leaf_node
     - |todo|
   * - dolfin.BoundaryMesh.name
     - |todo|
   * - dolfin.BoundaryMesh.order
     - |todo|
   * - dolfin.BoundaryMesh.parameters
     - |todo|
   * - dolfin.BoundaryMesh.parent
     - |todo|
   * - dolfin.BoundaryMesh.rename
     - |todo|
   * - dolfin.BoundaryMesh.renumber_by_color
     - |todo|
   * - dolfin.BoundaryMesh.root_node
     - |todo|
   * - dolfin.BoundaryMesh.scale
     - |todo|
   * - dolfin.BoundaryMesh.set_child
     - |todo|
   * - dolfin.BoundaryMesh.set_parent
     - |todo|
   * - dolfin.BoundaryMesh.str
     - |todo|
   * - dolfin.BoundingBoxTree.collides
     - |todo|
   * - dolfin.BoundingBoxTree.collides_entity
     - |todo|
   * - dolfin.BoundingBoxTree.compute_closest_point
     - |todo|
   * - dolfin.BoundingBoxTree.compute_process_collisions
     - |todo|
   * - dolfin.BoundingBoxTree3D
     - |todo|
   * - dolfin.BoundingBoxTree3D.build
     - |todo|
   * - dolfin.BoundingBoxTree3D.compute_closest_entity
     - |todo|
   * - dolfin.BoundingBoxTree3D.compute_closest_point
     - |todo|
   * - dolfin.BoundingBoxTree3D.compute_collisions
     - |todo|
   * - dolfin.BoundingBoxTree3D.compute_entity_collisions
     - |todo|
   * - dolfin.BoundingBoxTree3D.compute_first_collision
     - |todo|
   * - dolfin.BoundingBoxTree3D.compute_first_entity_collision
     - |todo|
   * - dolfin.BoundingBoxTree3D.compute_process_collisions
     - |todo|
   * - dolfin.BoundingBoxTree3D.create
     - |todo|
   * - dolfin.BoundingBoxTree3D.str
     - |todo|
   * - dolfin.BoxMesh.cell_name
     - |new|
   * - dolfin.BoxMesh.child
     - |todo|
   * - dolfin.BoxMesh.clean
     - |todo|
   * - dolfin.BoxMesh.clear_child
     - |todo|
   * - dolfin.BoxMesh.create
     - |todo|
   * - dolfin.BoxMesh.depth
     - |todo|
   * - dolfin.BoxMesh.geometric_dimension
     - |new|
   * - dolfin.BoxMesh.ghost_mode
     - |todo|
   * - dolfin.BoxMesh.has_child
     - |todo|
   * - dolfin.BoxMesh.has_parent
     - |todo|
   * - dolfin.BoxMesh.label
     - |todo|
   * - dolfin.BoxMesh.leaf_node
     - |todo|
   * - dolfin.BoxMesh.name
     - |todo|
   * - dolfin.BoxMesh.order
     - |todo|
   * - dolfin.BoxMesh.parameters
     - |todo|
   * - dolfin.BoxMesh.parent
     - |todo|
   * - dolfin.BoxMesh.rename
     - |todo|
   * - dolfin.BoxMesh.renumber_by_color
     - |todo|
   * - dolfin.BoxMesh.root_node
     - |todo|
   * - dolfin.BoxMesh.scale
     - |todo|
   * - dolfin.BoxMesh.set_child
     - |todo|
   * - dolfin.BoxMesh.set_parent
     - |todo|
   * - dolfin.BoxMesh.str
     - |todo|
   * - dolfin.BrokenElement
     - |todo|
   * - dolfin.BrokenElement.cell
     - |todo|
   * - dolfin.BrokenElement.degree
     - |todo|
   * - dolfin.BrokenElement.extract_component
     - |todo|
   * - dolfin.BrokenElement.extract_reference_component
     - |todo|
   * - dolfin.BrokenElement.extract_subelement_component
     - |todo|
   * - dolfin.BrokenElement.extract_subelement_reference_component
     - |todo|
   * - dolfin.BrokenElement.family
     - |todo|
   * - dolfin.BrokenElement.is_cellwise_constant
     - |todo|
   * - dolfin.BrokenElement.mapping
     - |todo|
   * - dolfin.BrokenElement.num_sub_elements
     - |todo|
   * - dolfin.BrokenElement.quadrature_scheme
     - |todo|
   * - dolfin.BrokenElement.reconstruct
     - |todo|
   * - dolfin.BrokenElement.reference_value_shape
     - |todo|
   * - dolfin.BrokenElement.reference_value_size
     - |todo|
   * - dolfin.BrokenElement.shortstr
     - |todo|
   * - dolfin.BrokenElement.sub_elements
     - |todo|
   * - dolfin.BrokenElement.symmetry
     - |todo|
   * - dolfin.BrokenElement.value_shape
     - |todo|
   * - dolfin.BrokenElement.value_size
     - |todo|
   * - dolfin.ButcherMultiStageScheme
     - |todo|
   * - dolfin.ButcherMultiStageScheme.bcs
     - |todo|
   * - dolfin.ButcherMultiStageScheme.dolfin_stage_forms
     - |todo|
   * - dolfin.ButcherMultiStageScheme.dt
     - |todo|
   * - dolfin.ButcherMultiStageScheme.dt_stage_offset
     - |todo|
   * - dolfin.ButcherMultiStageScheme.id
     - |todo|
   * - dolfin.ButcherMultiStageScheme.implicit
     - |todo|
   * - dolfin.ButcherMultiStageScheme.jacobian_index
     - |todo|
   * - dolfin.ButcherMultiStageScheme.label
     - |todo|
   * - dolfin.ButcherMultiStageScheme.last_stage
     - |todo|
   * - dolfin.ButcherMultiStageScheme.name
     - |todo|
   * - dolfin.ButcherMultiStageScheme.order
     - |todo|
   * - dolfin.ButcherMultiStageScheme.parameters
     - |todo|
   * - dolfin.ButcherMultiStageScheme.rename
     - |todo|
   * - dolfin.ButcherMultiStageScheme.rhs_form
     - |todo|
   * - dolfin.ButcherMultiStageScheme.solution
     - |todo|
   * - dolfin.ButcherMultiStageScheme.stage_solutions
     - |todo|
   * - dolfin.ButcherMultiStageScheme.str
     - |todo|
   * - dolfin.ButcherMultiStageScheme.t
     - |todo|
   * - dolfin.ButcherMultiStageScheme.to_adm
     - |todo|
   * - dolfin.ButcherMultiStageScheme.to_tlm
     - |todo|
   * - dolfin.ButcherMultiStageScheme.ufl_stage_forms
     - |todo|
   * - dolfin.CN2.bcs
     - |todo|
   * - dolfin.CN2.dt_stage_offset
     - |todo|
   * - dolfin.CN2.id
     - |todo|
   * - dolfin.CN2.implicit
     - |todo|
   * - dolfin.CN2.jacobian_index
     - |todo|
   * - dolfin.CN2.label
     - |todo|
   * - dolfin.CN2.name
     - |todo|
   * - dolfin.CN2.parameters
     - |todo|
   * - dolfin.CN2.rename
     - |todo|
   * - dolfin.CN2.str
     - |todo|
   * - dolfin.CRITICAL
     - |todo|
   * - dolfin.Cell.cell_normal
     - |todo|
   * - dolfin.Cell.get_cell_data
     - |todo|
   * - dolfin.Cell.get_cell_topology
     - |todo|
   * - dolfin.Cell.get_coordinate_dofs
     - |todo|
   * - dolfin.Cell.incident
     - |todo|
   * - dolfin.Cell.init
     - |todo|
   * - dolfin.Cell.intersection
     - |todo|
   * - dolfin.Cell.mesh_id
     - |todo|
   * - dolfin.Cell.num_vertices
     - |todo|
   * - dolfin.Cell.order
     - |todo|
   * - dolfin.Cell.ordered
     - |todo|
   * - dolfin.Cell.owner
     - |todo|
   * - dolfin.Cell.squared_distance
     - |todo|
   * - dolfin.Cell.str
     - |todo|
   * - dolfin.Cell.type
     - |todo|
   * - dolfin.CellFunctionBool
     - |todo|
   * - dolfin.CellFunctionBool.array
     - |todo|
   * - dolfin.CellFunctionBool.child
     - |todo|
   * - dolfin.CellFunctionBool.clear_child
     - |todo|
   * - dolfin.CellFunctionBool.cpp_value_type
     - |todo|
   * - dolfin.CellFunctionBool.depth
     - |todo|
   * - dolfin.CellFunctionBool.dim
     - |todo|
   * - dolfin.CellFunctionBool.empty
     - |todo|
   * - dolfin.CellFunctionBool.has_child
     - |todo|
   * - dolfin.CellFunctionBool.has_parent
     - |todo|
   * - dolfin.CellFunctionBool.id
     - |todo|
   * - dolfin.CellFunctionBool.init
     - |todo|
   * - dolfin.CellFunctionBool.label
     - |todo|
   * - dolfin.CellFunctionBool.leaf_node
     - |todo|
   * - dolfin.CellFunctionBool.mesh
     - |todo|
   * - dolfin.CellFunctionBool.name
     - |todo|
   * - dolfin.CellFunctionBool.parameters
     - |todo|
   * - dolfin.CellFunctionBool.parent
     - |todo|
   * - dolfin.CellFunctionBool.rename
     - |todo|
   * - dolfin.CellFunctionBool.root_node
     - |todo|
   * - dolfin.CellFunctionBool.set_all
     - |todo|
   * - dolfin.CellFunctionBool.set_child
     - |todo|
   * - dolfin.CellFunctionBool.set_parent
     - |todo|
   * - dolfin.CellFunctionBool.set_value
     - |todo|
   * - dolfin.CellFunctionBool.set_values
     - |todo|
   * - dolfin.CellFunctionBool.size
     - |todo|
   * - dolfin.CellFunctionBool.str
     - |todo|
   * - dolfin.CellFunctionBool.ufl_id
     - |todo|
   * - dolfin.CellFunctionBool.where_equal
     - |todo|
   * - dolfin.CellFunctionDouble
     - |todo|
   * - dolfin.CellFunctionDouble.array
     - |todo|
   * - dolfin.CellFunctionDouble.child
     - |todo|
   * - dolfin.CellFunctionDouble.clear_child
     - |todo|
   * - dolfin.CellFunctionDouble.cpp_value_type
     - |todo|
   * - dolfin.CellFunctionDouble.depth
     - |todo|
   * - dolfin.CellFunctionDouble.dim
     - |todo|
   * - dolfin.CellFunctionDouble.empty
     - |todo|
   * - dolfin.CellFunctionDouble.has_child
     - |todo|
   * - dolfin.CellFunctionDouble.has_parent
     - |todo|
   * - dolfin.CellFunctionDouble.id
     - |todo|
   * - dolfin.CellFunctionDouble.init
     - |todo|
   * - dolfin.CellFunctionDouble.label
     - |todo|
   * - dolfin.CellFunctionDouble.leaf_node
     - |todo|
   * - dolfin.CellFunctionDouble.mesh
     - |todo|
   * - dolfin.CellFunctionDouble.name
     - |todo|
   * - dolfin.CellFunctionDouble.parameters
     - |todo|
   * - dolfin.CellFunctionDouble.parent
     - |todo|
   * - dolfin.CellFunctionDouble.rename
     - |todo|
   * - dolfin.CellFunctionDouble.root_node
     - |todo|
   * - dolfin.CellFunctionDouble.set_all
     - |todo|
   * - dolfin.CellFunctionDouble.set_child
     - |todo|
   * - dolfin.CellFunctionDouble.set_parent
     - |todo|
   * - dolfin.CellFunctionDouble.set_value
     - |todo|
   * - dolfin.CellFunctionDouble.set_values
     - |todo|
   * - dolfin.CellFunctionDouble.size
     - |todo|
   * - dolfin.CellFunctionDouble.str
     - |todo|
   * - dolfin.CellFunctionDouble.ufl_id
     - |todo|
   * - dolfin.CellFunctionDouble.where_equal
     - |todo|
   * - dolfin.CellFunctionInt
     - |todo|
   * - dolfin.CellFunctionInt.array
     - |todo|
   * - dolfin.CellFunctionInt.child
     - |todo|
   * - dolfin.CellFunctionInt.clear_child
     - |todo|
   * - dolfin.CellFunctionInt.cpp_value_type
     - |todo|
   * - dolfin.CellFunctionInt.depth
     - |todo|
   * - dolfin.CellFunctionInt.dim
     - |todo|
   * - dolfin.CellFunctionInt.empty
     - |todo|
   * - dolfin.CellFunctionInt.has_child
     - |todo|
   * - dolfin.CellFunctionInt.has_parent
     - |todo|
   * - dolfin.CellFunctionInt.id
     - |todo|
   * - dolfin.CellFunctionInt.init
     - |todo|
   * - dolfin.CellFunctionInt.label
     - |todo|
   * - dolfin.CellFunctionInt.leaf_node
     - |todo|
   * - dolfin.CellFunctionInt.mesh
     - |todo|
   * - dolfin.CellFunctionInt.name
     - |todo|
   * - dolfin.CellFunctionInt.parameters
     - |todo|
   * - dolfin.CellFunctionInt.parent
     - |todo|
   * - dolfin.CellFunctionInt.rename
     - |todo|
   * - dolfin.CellFunctionInt.root_node
     - |todo|
   * - dolfin.CellFunctionInt.set_all
     - |todo|
   * - dolfin.CellFunctionInt.set_child
     - |todo|
   * - dolfin.CellFunctionInt.set_parent
     - |todo|
   * - dolfin.CellFunctionInt.set_value
     - |todo|
   * - dolfin.CellFunctionInt.set_values
     - |todo|
   * - dolfin.CellFunctionInt.size
     - |todo|
   * - dolfin.CellFunctionInt.str
     - |todo|
   * - dolfin.CellFunctionInt.ufl_id
     - |todo|
   * - dolfin.CellFunctionInt.where_equal
     - |todo|
   * - dolfin.CellFunctionSizet
     - |todo|
   * - dolfin.CellFunctionSizet.array
     - |todo|
   * - dolfin.CellFunctionSizet.child
     - |todo|
   * - dolfin.CellFunctionSizet.clear_child
     - |todo|
   * - dolfin.CellFunctionSizet.cpp_value_type
     - |todo|
   * - dolfin.CellFunctionSizet.depth
     - |todo|
   * - dolfin.CellFunctionSizet.dim
     - |todo|
   * - dolfin.CellFunctionSizet.empty
     - |todo|
   * - dolfin.CellFunctionSizet.has_child
     - |todo|
   * - dolfin.CellFunctionSizet.has_parent
     - |todo|
   * - dolfin.CellFunctionSizet.id
     - |todo|
   * - dolfin.CellFunctionSizet.init
     - |todo|
   * - dolfin.CellFunctionSizet.label
     - |todo|
   * - dolfin.CellFunctionSizet.leaf_node
     - |todo|
   * - dolfin.CellFunctionSizet.mesh
     - |todo|
   * - dolfin.CellFunctionSizet.name
     - |todo|
   * - dolfin.CellFunctionSizet.parameters
     - |todo|
   * - dolfin.CellFunctionSizet.parent
     - |todo|
   * - dolfin.CellFunctionSizet.rename
     - |todo|
   * - dolfin.CellFunctionSizet.root_node
     - |todo|
   * - dolfin.CellFunctionSizet.set_all
     - |todo|
   * - dolfin.CellFunctionSizet.set_child
     - |todo|
   * - dolfin.CellFunctionSizet.set_parent
     - |todo|
   * - dolfin.CellFunctionSizet.set_value
     - |todo|
   * - dolfin.CellFunctionSizet.set_values
     - |todo|
   * - dolfin.CellFunctionSizet.size
     - |todo|
   * - dolfin.CellFunctionSizet.str
     - |todo|
   * - dolfin.CellFunctionSizet.ufl_id
     - |todo|
   * - dolfin.CellFunctionSizet.where_equal
     - |todo|
   * - dolfin.CellSize
     - |todo|
   * - dolfin.CellType.Type
     - |new|
   * - dolfin.CellType.cell_normal
     - |todo|
   * - dolfin.CellType.circumradius
     - |todo|
   * - dolfin.CellType.collides
     - |todo|
   * - dolfin.CellType.create
     - |todo|
   * - dolfin.CellType.create_entities
     - |todo|
   * - dolfin.CellType.dim
     - |todo|
   * - dolfin.CellType.entity_type
     - |todo|
   * - dolfin.CellType.facet_area
     - |todo|
   * - dolfin.CellType.facet_type
     - |todo|
   * - dolfin.CellType.h
     - |todo|
   * - dolfin.CellType.hexahedron
     - |todo|
   * - dolfin.CellType.inradius
     - |todo|
   * - dolfin.CellType.interval
     - |todo|
   * - dolfin.CellType.is_simplex
     - |todo|
   * - dolfin.CellType.normal
     - |todo|
   * - dolfin.CellType.num_entities
     - |todo|
   * - dolfin.CellType.num_vertices
     - |todo|
   * - dolfin.CellType.order
     - |todo|
   * - dolfin.CellType.ordered
     - |todo|
   * - dolfin.CellType.orientation
     - |todo|
   * - dolfin.CellType.point
     - |todo|
   * - dolfin.CellType.quadrilateral
     - |todo|
   * - dolfin.CellType.radius_ratio
     - |todo|
   * - dolfin.CellType.squared_distance
     - |todo|
   * - dolfin.CellType.string2type
     - |todo|
   * - dolfin.CellType.tetrahedron
     - |todo|
   * - dolfin.CellType.triangle
     - |todo|
   * - dolfin.CellType.volume
     - |todo|
   * - dolfin.CellType.vtk_mapping
     - |todo|
   * - dolfin.Coefficient
     - |todo|
   * - dolfin.Coefficient.T
     - |todo|
   * - dolfin.Coefficient.count
     - |todo|
   * - dolfin.Coefficient.dx
     - |todo|
   * - dolfin.Coefficient.evaluate
     - |todo|
   * - dolfin.Coefficient.geometric_dimension
     - |todo|
   * - dolfin.Coefficient.is_cellwise_constant
     - |todo|
   * - dolfin.Coefficient.ufl_disable_profiling
     - |todo|
   * - dolfin.Coefficient.ufl_domain
     - |todo|
   * - dolfin.Coefficient.ufl_domains
     - |todo|
   * - dolfin.Coefficient.ufl_element
     - |todo|
   * - dolfin.Coefficient.ufl_enable_profiling
     - |todo|
   * - dolfin.Coefficient.ufl_free_indices
     - |todo|
   * - dolfin.Coefficient.ufl_function_space
     - |todo|
   * - dolfin.Coefficient.ufl_index_dimensions
     - |todo|
   * - dolfin.Coefficient.ufl_operands
     - |todo|
   * - dolfin.Coefficient.ufl_shape
     - |todo|
   * - dolfin.Coefficients
     - |todo|
   * - dolfin.CollisionPredicates
     - |todo|
   * - dolfin.CollisionPredicates.collides
     - |todo|
   * - dolfin.CollisionPredicates.collides_segment_point
     - |todo|
   * - dolfin.CollisionPredicates.collides_segment_point_1d
     - |todo|
   * - dolfin.CollisionPredicates.collides_segment_point_2d
     - |todo|
   * - dolfin.CollisionPredicates.collides_segment_point_3d
     - |todo|
   * - dolfin.CollisionPredicates.collides_segment_segment
     - |todo|
   * - dolfin.CollisionPredicates.collides_segment_segment_1d
     - |todo|
   * - dolfin.CollisionPredicates.collides_segment_segment_2d
     - |todo|
   * - dolfin.CollisionPredicates.collides_segment_segment_3d
     - |todo|
   * - dolfin.CollisionPredicates.collides_tetrahedron_point_3d
     - |todo|
   * - dolfin.CollisionPredicates.collides_tetrahedron_segment_3d
     - |todo|
   * - dolfin.CollisionPredicates.collides_tetrahedron_tetrahedron_3d
     - |todo|
   * - dolfin.CollisionPredicates.collides_tetrahedron_triangle_3d
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_point
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_point_2d
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_point_3d
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_segment
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_segment_2d
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_segment_3d
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_triangle
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_triangle_2d
     - |todo|
   * - dolfin.CollisionPredicates.collides_triangle_triangle_3d
     - |todo|
   * - dolfin.CompiledExpression
     - |new|
   * - dolfin.CompiledExpression.T
     - |new|
   * - dolfin.CompiledExpression.compute_vertex_values
     - |new|
   * - dolfin.CompiledExpression.count
     - |new|
   * - dolfin.CompiledExpression.cpp_object
     - |new|
   * - dolfin.CompiledExpression.dx
     - |new|
   * - dolfin.CompiledExpression.evaluate
     - |new|
   * - dolfin.CompiledExpression.geometric_dimension
     - |new|
   * - dolfin.CompiledExpression.id
     - |new|
   * - dolfin.CompiledExpression.is_cellwise_constant
     - |new|
   * - dolfin.CompiledExpression.label
     - |new|
   * - dolfin.CompiledExpression.name
     - |new|
   * - dolfin.CompiledExpression.ufl_disable_profiling
     - |new|
   * - dolfin.CompiledExpression.ufl_domain
     - |new|
   * - dolfin.CompiledExpression.ufl_domains
     - |new|
   * - dolfin.CompiledExpression.ufl_element
     - |new|
   * - dolfin.CompiledExpression.ufl_enable_profiling
     - |new|
   * - dolfin.CompiledExpression.ufl_evaluate
     - |new|
   * - dolfin.CompiledExpression.ufl_free_indices
     - |new|
   * - dolfin.CompiledExpression.ufl_function_space
     - |new|
   * - dolfin.CompiledExpression.ufl_index_dimensions
     - |new|
   * - dolfin.CompiledExpression.ufl_operands
     - |new|
   * - dolfin.CompiledExpression.ufl_shape
     - |new|
   * - dolfin.CompiledExpression.value_dimension
     - |new|
   * - dolfin.CompiledExpression.value_rank
     - |new|
   * - dolfin.CompiledSubDomain.get_property
     - |new|
   * - dolfin.CompiledSubDomain.inside
     - |new|
   * - dolfin.CompiledSubDomain.map
     - |new|
   * - dolfin.CompiledSubDomain.mark
     - |new|
   * - dolfin.CompiledSubDomain.mark_cells
     - |new|
   * - dolfin.CompiledSubDomain.mark_facets
     - |new|
   * - dolfin.CompiledSubDomain.set_property
     - |new|
   * - dolfin.Constant.cpp_object
     - |new|
   * - dolfin.Constant.eval
     - |todo|
   * - dolfin.Constant.eval_cell
     - |todo|
   * - dolfin.Constant.get_generic_function
     - |todo|
   * - dolfin.Constant.get_property
     - |todo|
   * - dolfin.Constant.label
     - |todo|
   * - dolfin.Constant.parameters
     - |todo|
   * - dolfin.Constant.restrict
     - |todo|
   * - dolfin.Constant.set_generic_function
     - |todo|
   * - dolfin.Constant.set_property
     - |todo|
   * - dolfin.Constant.update
     - |todo|
   * - dolfin.Constant.value_dimension
     - |todo|
   * - dolfin.Constant.value_rank
     - |todo|
   * - dolfin.Constant.value_shape
     - |todo|
   * - dolfin.CoordinateMatrix
     - |todo|
   * - dolfin.CoordinateMatrix.base_one
     - |todo|
   * - dolfin.CoordinateMatrix.columns
     - |todo|
   * - dolfin.CoordinateMatrix.mpi_comm
     - |todo|
   * - dolfin.CoordinateMatrix.norm
     - |todo|
   * - dolfin.CoordinateMatrix.rows
     - |todo|
   * - dolfin.CoordinateMatrix.size
     - |todo|
   * - dolfin.CoordinateMatrix.values
     - |todo|
   * - dolfin.CrankNicolson
     - |todo|
   * - dolfin.CrankNicolson.bcs
     - |todo|
   * - dolfin.CrankNicolson.dolfin_stage_forms
     - |todo|
   * - dolfin.CrankNicolson.dt
     - |todo|
   * - dolfin.CrankNicolson.dt_stage_offset
     - |todo|
   * - dolfin.CrankNicolson.id
     - |todo|
   * - dolfin.CrankNicolson.implicit
     - |todo|
   * - dolfin.CrankNicolson.jacobian_index
     - |todo|
   * - dolfin.CrankNicolson.label
     - |todo|
   * - dolfin.CrankNicolson.last_stage
     - |todo|
   * - dolfin.CrankNicolson.name
     - |todo|
   * - dolfin.CrankNicolson.order
     - |todo|
   * - dolfin.CrankNicolson.parameters
     - |todo|
   * - dolfin.CrankNicolson.rename
     - |todo|
   * - dolfin.CrankNicolson.rhs_form
     - |todo|
   * - dolfin.CrankNicolson.solution
     - |todo|
   * - dolfin.CrankNicolson.stage_solutions
     - |todo|
   * - dolfin.CrankNicolson.str
     - |todo|
   * - dolfin.CrankNicolson.t
     - |todo|
   * - dolfin.CrankNicolson.to_adm
     - |todo|
   * - dolfin.CrankNicolson.to_tlm
     - |todo|
   * - dolfin.CrankNicolson.ufl_stage_forms
     - |todo|
   * - dolfin.DBG
     - |todo|
   * - dolfin.DEBUG
     - |todo|
   * - dolfin.DOLFIN_EPS_LARGE
     - |todo|
   * - dolfin.DOLFIN_LINELENGTH
     - |todo|
   * - dolfin.DOLFIN_SQRT_EPS
     - |todo|
   * - dolfin.DOLFIN_TERM_WIDTH
     - |todo|
   * - dolfin.DOLFIN_VERSION_GIT
     - |todo|
   * - dolfin.DOLFIN_VERSION_MAJOR
     - |todo|
   * - dolfin.DOLFIN_VERSION_MICRO
     - |todo|
   * - dolfin.DOLFIN_VERSION_MINOR
     - |todo|
   * - dolfin.DOLFIN_VERSION_RELEASE
     - |todo|
   * - dolfin.DOLFIN_VERSION_STRING
     - |todo|
   * - dolfin.DefaultFactory.create_krylov_solver
     - |todo|
   * - dolfin.DefaultFactory.create_layout
     - |todo|
   * - dolfin.DefaultFactory.create_linear_operator
     - |todo|
   * - dolfin.DefaultFactory.create_lu_solver
     - |todo|
   * - dolfin.DefaultFactory.krylov_solver_methods
     - |todo|
   * - dolfin.DefaultFactory.krylov_solver_preconditioners
     - |todo|
   * - dolfin.DefaultFactory.lu_solver_methods
     - |todo|
   * - dolfin.DirichletBC.child
     - |todo|
   * - dolfin.DirichletBC.clear_child
     - |todo|
   * - dolfin.DirichletBC.default_parameters
     - |todo|
   * - dolfin.DirichletBC.depth
     - |todo|
   * - dolfin.DirichletBC.has_child
     - |todo|
   * - dolfin.DirichletBC.has_parent
     - |todo|
   * - dolfin.DirichletBC.id
     - |todo|
   * - dolfin.DirichletBC.label
     - |todo|
   * - dolfin.DirichletBC.leaf_node
     - |todo|
   * - dolfin.DirichletBC.markers
     - |todo|
   * - dolfin.DirichletBC.name
     - |todo|
   * - dolfin.DirichletBC.parameters
     - |todo|
   * - dolfin.DirichletBC.parent
     - |todo|
   * - dolfin.DirichletBC.rename
     - |todo|
   * - dolfin.DirichletBC.root_node
     - |todo|
   * - dolfin.DirichletBC.set_child
     - |todo|
   * - dolfin.DirichletBC.set_parent
     - |todo|
   * - dolfin.DirichletBC.str
     - |todo|
   * - dolfin.DirichletBC.user_sub_domain
     - |todo|
   * - dolfin.DirichletBC.user_subdomain
     - |new|
   * - dolfin.DirichletBC.value
     - |todo|
   * - dolfin.Dn
     - |todo|
   * - dolfin.DofMap.cell_dimension
     - |todo|
   * - dolfin.DofMap.collapse
     - |todo|
   * - dolfin.DofMap.copy
     - |todo|
   * - dolfin.DofMap.create
     - |todo|
   * - dolfin.DofMap.extract_sub_dofmap
     - |todo|
   * - dolfin.DofMap.global_dimension
     - |todo|
   * - dolfin.DofMap.id
     - |todo|
   * - dolfin.DofMap.is_view
     - |todo|
   * - dolfin.DofMap.label
     - |todo|
   * - dolfin.DofMap.local_to_global_unowned
     - |todo|
   * - dolfin.DofMap.max_cell_dimension
     - |todo|
   * - dolfin.DofMap.max_element_dofs
     - |todo|
   * - dolfin.DofMap.name
     - |todo|
   * - dolfin.DofMap.num_element_dofs
     - |todo|
   * - dolfin.DofMap.num_entity_closure_dofs
     - |todo|
   * - dolfin.DofMap.num_facet_dofs
     - |todo|
   * - dolfin.DofMap.parameters
     - |todo|
   * - dolfin.DofMap.rename
     - |todo|
   * - dolfin.DofMap.str
     - |todo|
   * - dolfin.DofMap.tabulate_entity_closure_dofs
     - |todo|
   * - dolfin.DofMap.tabulate_facet_dofs
     - |todo|
   * - dolfin.DofMap.tabulate_global_dofs
     - |todo|
   * - dolfin.DomainBoundary.geometric_dimension
     - |todo|
   * - dolfin.DomainBoundary.map_tolerance
     - |todo|
   * - dolfin.DomainBoundary.snap
     - |todo|
   * - dolfin.DoubleArray
     - |todo|
   * - dolfin.DoubleArray.array
     - |todo|
   * - dolfin.DoubleArray.data
     - |todo|
   * - dolfin.DoubleArray.size
     - |todo|
   * - dolfin.DoubleArray.str
     - |todo|
   * - dolfin.DynamicMeshEditor
     - |todo|
   * - dolfin.DynamicMeshEditor.add_cell
     - |todo|
   * - dolfin.DynamicMeshEditor.add_vertex
     - |todo|
   * - dolfin.DynamicMeshEditor.close
     - |todo|
   * - dolfin.DynamicMeshEditor.open
     - |todo|
   * - dolfin.ERK
     - |todo|
   * - dolfin.ERK.bcs
     - |todo|
   * - dolfin.ERK.dolfin_stage_forms
     - |todo|
   * - dolfin.ERK.dt
     - |todo|
   * - dolfin.ERK.dt_stage_offset
     - |todo|
   * - dolfin.ERK.id
     - |todo|
   * - dolfin.ERK.implicit
     - |todo|
   * - dolfin.ERK.jacobian_index
     - |todo|
   * - dolfin.ERK.label
     - |todo|
   * - dolfin.ERK.last_stage
     - |todo|
   * - dolfin.ERK.name
     - |todo|
   * - dolfin.ERK.order
     - |todo|
   * - dolfin.ERK.parameters
     - |todo|
   * - dolfin.ERK.rename
     - |todo|
   * - dolfin.ERK.rhs_form
     - |todo|
   * - dolfin.ERK.solution
     - |todo|
   * - dolfin.ERK.stage_solutions
     - |todo|
   * - dolfin.ERK.str
     - |todo|
   * - dolfin.ERK.t
     - |todo|
   * - dolfin.ERK.to_adm
     - |todo|
   * - dolfin.ERK.to_tlm
     - |todo|
   * - dolfin.ERK.ufl_stage_forms
     - |todo|
   * - dolfin.ERK1
     - |todo|
   * - dolfin.ERK1.bcs
     - |todo|
   * - dolfin.ERK1.dolfin_stage_forms
     - |todo|
   * - dolfin.ERK1.dt
     - |todo|
   * - dolfin.ERK1.dt_stage_offset
     - |todo|
   * - dolfin.ERK1.id
     - |todo|
   * - dolfin.ERK1.implicit
     - |todo|
   * - dolfin.ERK1.jacobian_index
     - |todo|
   * - dolfin.ERK1.label
     - |todo|
   * - dolfin.ERK1.last_stage
     - |todo|
   * - dolfin.ERK1.name
     - |todo|
   * - dolfin.ERK1.order
     - |todo|
   * - dolfin.ERK1.parameters
     - |todo|
   * - dolfin.ERK1.rename
     - |todo|
   * - dolfin.ERK1.rhs_form
     - |todo|
   * - dolfin.ERK1.solution
     - |todo|
   * - dolfin.ERK1.stage_solutions
     - |todo|
   * - dolfin.ERK1.str
     - |todo|
   * - dolfin.ERK1.t
     - |todo|
   * - dolfin.ERK1.to_adm
     - |todo|
   * - dolfin.ERK1.to_tlm
     - |todo|
   * - dolfin.ERK1.ufl_stage_forms
     - |todo|
   * - dolfin.ERK4
     - |todo|
   * - dolfin.ERK4.bcs
     - |todo|
   * - dolfin.ERK4.dolfin_stage_forms
     - |todo|
   * - dolfin.ERK4.dt
     - |todo|
   * - dolfin.ERK4.dt_stage_offset
     - |todo|
   * - dolfin.ERK4.id
     - |todo|
   * - dolfin.ERK4.implicit
     - |todo|
   * - dolfin.ERK4.jacobian_index
     - |todo|
   * - dolfin.ERK4.label
     - |todo|
   * - dolfin.ERK4.last_stage
     - |todo|
   * - dolfin.ERK4.name
     - |todo|
   * - dolfin.ERK4.order
     - |todo|
   * - dolfin.ERK4.parameters
     - |todo|
   * - dolfin.ERK4.rename
     - |todo|
   * - dolfin.ERK4.rhs_form
     - |todo|
   * - dolfin.ERK4.solution
     - |todo|
   * - dolfin.ERK4.stage_solutions
     - |todo|
   * - dolfin.ERK4.str
     - |todo|
   * - dolfin.ERK4.t
     - |todo|
   * - dolfin.ERK4.to_adm
     - |todo|
   * - dolfin.ERK4.to_tlm
     - |todo|
   * - dolfin.ERK4.ufl_stage_forms
     - |todo|
   * - dolfin.ERROR
     - |todo|
   * - dolfin.ESDIRK3.bcs
     - |todo|
   * - dolfin.ESDIRK3.dt_stage_offset
     - |todo|
   * - dolfin.ESDIRK3.id
     - |todo|
   * - dolfin.ESDIRK3.implicit
     - |todo|
   * - dolfin.ESDIRK3.jacobian_index
     - |todo|
   * - dolfin.ESDIRK3.label
     - |todo|
   * - dolfin.ESDIRK3.name
     - |todo|
   * - dolfin.ESDIRK3.parameters
     - |todo|
   * - dolfin.ESDIRK3.rename
     - |todo|
   * - dolfin.ESDIRK3.str
     - |todo|
   * - dolfin.ESDIRK4.bcs
     - |todo|
   * - dolfin.ESDIRK4.dt_stage_offset
     - |todo|
   * - dolfin.ESDIRK4.id
     - |todo|
   * - dolfin.ESDIRK4.implicit
     - |todo|
   * - dolfin.ESDIRK4.jacobian_index
     - |todo|
   * - dolfin.ESDIRK4.label
     - |todo|
   * - dolfin.ESDIRK4.name
     - |todo|
   * - dolfin.ESDIRK4.parameters
     - |todo|
   * - dolfin.ESDIRK4.rename
     - |todo|
   * - dolfin.ESDIRK4.str
     - |todo|
   * - dolfin.Edge.incident
     - |todo|
   * - dolfin.Edge.init
     - |todo|
   * - dolfin.Edge.mesh_id
     - |todo|
   * - dolfin.Edge.owner
     - |todo|
   * - dolfin.Edge.str
     - |todo|
   * - dolfin.EdgeFunctionBool
     - |todo|
   * - dolfin.EdgeFunctionBool.array
     - |todo|
   * - dolfin.EdgeFunctionBool.child
     - |todo|
   * - dolfin.EdgeFunctionBool.clear_child
     - |todo|
   * - dolfin.EdgeFunctionBool.cpp_value_type
     - |todo|
   * - dolfin.EdgeFunctionBool.depth
     - |todo|
   * - dolfin.EdgeFunctionBool.dim
     - |todo|
   * - dolfin.EdgeFunctionBool.empty
     - |todo|
   * - dolfin.EdgeFunctionBool.has_child
     - |todo|
   * - dolfin.EdgeFunctionBool.has_parent
     - |todo|
   * - dolfin.EdgeFunctionBool.id
     - |todo|
   * - dolfin.EdgeFunctionBool.init
     - |todo|
   * - dolfin.EdgeFunctionBool.label
     - |todo|
   * - dolfin.EdgeFunctionBool.leaf_node
     - |todo|
   * - dolfin.EdgeFunctionBool.mesh
     - |todo|
   * - dolfin.EdgeFunctionBool.name
     - |todo|
   * - dolfin.EdgeFunctionBool.parameters
     - |todo|
   * - dolfin.EdgeFunctionBool.parent
     - |todo|
   * - dolfin.EdgeFunctionBool.rename
     - |todo|
   * - dolfin.EdgeFunctionBool.root_node
     - |todo|
   * - dolfin.EdgeFunctionBool.set_all
     - |todo|
   * - dolfin.EdgeFunctionBool.set_child
     - |todo|
   * - dolfin.EdgeFunctionBool.set_parent
     - |todo|
   * - dolfin.EdgeFunctionBool.set_value
     - |todo|
   * - dolfin.EdgeFunctionBool.set_values
     - |todo|
   * - dolfin.EdgeFunctionBool.size
     - |todo|
   * - dolfin.EdgeFunctionBool.str
     - |todo|
   * - dolfin.EdgeFunctionBool.ufl_id
     - |todo|
   * - dolfin.EdgeFunctionBool.where_equal
     - |todo|
   * - dolfin.EdgeFunctionDouble
     - |todo|
   * - dolfin.EdgeFunctionDouble.array
     - |todo|
   * - dolfin.EdgeFunctionDouble.child
     - |todo|
   * - dolfin.EdgeFunctionDouble.clear_child
     - |todo|
   * - dolfin.EdgeFunctionDouble.cpp_value_type
     - |todo|
   * - dolfin.EdgeFunctionDouble.depth
     - |todo|
   * - dolfin.EdgeFunctionDouble.dim
     - |todo|
   * - dolfin.EdgeFunctionDouble.empty
     - |todo|
   * - dolfin.EdgeFunctionDouble.has_child
     - |todo|
   * - dolfin.EdgeFunctionDouble.has_parent
     - |todo|
   * - dolfin.EdgeFunctionDouble.id
     - |todo|
   * - dolfin.EdgeFunctionDouble.init
     - |todo|
   * - dolfin.EdgeFunctionDouble.label
     - |todo|
   * - dolfin.EdgeFunctionDouble.leaf_node
     - |todo|
   * - dolfin.EdgeFunctionDouble.mesh
     - |todo|
   * - dolfin.EdgeFunctionDouble.name
     - |todo|
   * - dolfin.EdgeFunctionDouble.parameters
     - |todo|
   * - dolfin.EdgeFunctionDouble.parent
     - |todo|
   * - dolfin.EdgeFunctionDouble.rename
     - |todo|
   * - dolfin.EdgeFunctionDouble.root_node
     - |todo|
   * - dolfin.EdgeFunctionDouble.set_all
     - |todo|
   * - dolfin.EdgeFunctionDouble.set_child
     - |todo|
   * - dolfin.EdgeFunctionDouble.set_parent
     - |todo|
   * - dolfin.EdgeFunctionDouble.set_value
     - |todo|
   * - dolfin.EdgeFunctionDouble.set_values
     - |todo|
   * - dolfin.EdgeFunctionDouble.size
     - |todo|
   * - dolfin.EdgeFunctionDouble.str
     - |todo|
   * - dolfin.EdgeFunctionDouble.ufl_id
     - |todo|
   * - dolfin.EdgeFunctionDouble.where_equal
     - |todo|
   * - dolfin.EdgeFunctionInt
     - |todo|
   * - dolfin.EdgeFunctionInt.array
     - |todo|
   * - dolfin.EdgeFunctionInt.child
     - |todo|
   * - dolfin.EdgeFunctionInt.clear_child
     - |todo|
   * - dolfin.EdgeFunctionInt.cpp_value_type
     - |todo|
   * - dolfin.EdgeFunctionInt.depth
     - |todo|
   * - dolfin.EdgeFunctionInt.dim
     - |todo|
   * - dolfin.EdgeFunctionInt.empty
     - |todo|
   * - dolfin.EdgeFunctionInt.has_child
     - |todo|
   * - dolfin.EdgeFunctionInt.has_parent
     - |todo|
   * - dolfin.EdgeFunctionInt.id
     - |todo|
   * - dolfin.EdgeFunctionInt.init
     - |todo|
   * - dolfin.EdgeFunctionInt.label
     - |todo|
   * - dolfin.EdgeFunctionInt.leaf_node
     - |todo|
   * - dolfin.EdgeFunctionInt.mesh
     - |todo|
   * - dolfin.EdgeFunctionInt.name
     - |todo|
   * - dolfin.EdgeFunctionInt.parameters
     - |todo|
   * - dolfin.EdgeFunctionInt.parent
     - |todo|
   * - dolfin.EdgeFunctionInt.rename
     - |todo|
   * - dolfin.EdgeFunctionInt.root_node
     - |todo|
   * - dolfin.EdgeFunctionInt.set_all
     - |todo|
   * - dolfin.EdgeFunctionInt.set_child
     - |todo|
   * - dolfin.EdgeFunctionInt.set_parent
     - |todo|
   * - dolfin.EdgeFunctionInt.set_value
     - |todo|
   * - dolfin.EdgeFunctionInt.set_values
     - |todo|
   * - dolfin.EdgeFunctionInt.size
     - |todo|
   * - dolfin.EdgeFunctionInt.str
     - |todo|
   * - dolfin.EdgeFunctionInt.ufl_id
     - |todo|
   * - dolfin.EdgeFunctionInt.where_equal
     - |todo|
   * - dolfin.EdgeFunctionSizet
     - |todo|
   * - dolfin.EdgeFunctionSizet.array
     - |todo|
   * - dolfin.EdgeFunctionSizet.child
     - |todo|
   * - dolfin.EdgeFunctionSizet.clear_child
     - |todo|
   * - dolfin.EdgeFunctionSizet.cpp_value_type
     - |todo|
   * - dolfin.EdgeFunctionSizet.depth
     - |todo|
   * - dolfin.EdgeFunctionSizet.dim
     - |todo|
   * - dolfin.EdgeFunctionSizet.empty
     - |todo|
   * - dolfin.EdgeFunctionSizet.has_child
     - |todo|
   * - dolfin.EdgeFunctionSizet.has_parent
     - |todo|
   * - dolfin.EdgeFunctionSizet.id
     - |todo|
   * - dolfin.EdgeFunctionSizet.init
     - |todo|
   * - dolfin.EdgeFunctionSizet.label
     - |todo|
   * - dolfin.EdgeFunctionSizet.leaf_node
     - |todo|
   * - dolfin.EdgeFunctionSizet.mesh
     - |todo|
   * - dolfin.EdgeFunctionSizet.name
     - |todo|
   * - dolfin.EdgeFunctionSizet.parameters
     - |todo|
   * - dolfin.EdgeFunctionSizet.parent
     - |todo|
   * - dolfin.EdgeFunctionSizet.rename
     - |todo|
   * - dolfin.EdgeFunctionSizet.root_node
     - |todo|
   * - dolfin.EdgeFunctionSizet.set_all
     - |todo|
   * - dolfin.EdgeFunctionSizet.set_child
     - |todo|
   * - dolfin.EdgeFunctionSizet.set_parent
     - |todo|
   * - dolfin.EdgeFunctionSizet.set_value
     - |todo|
   * - dolfin.EdgeFunctionSizet.set_values
     - |todo|
   * - dolfin.EdgeFunctionSizet.size
     - |todo|
   * - dolfin.EdgeFunctionSizet.str
     - |todo|
   * - dolfin.EdgeFunctionSizet.ufl_id
     - |todo|
   * - dolfin.EdgeFunctionSizet.where_equal
     - |todo|
   * - dolfin.EigenFactory.create_krylov_solver
     - |todo|
   * - dolfin.EigenFactory.create_layout
     - |todo|
   * - dolfin.EigenFactory.create_linear_operator
     - |todo|
   * - dolfin.EigenFactory.create_lu_solver
     - |todo|
   * - dolfin.EigenFactory.krylov_solver_methods
     - |todo|
   * - dolfin.EigenFactory.krylov_solver_preconditioners
     - |todo|
   * - dolfin.EigenFactory.lu_solver_methods
     - |todo|
   * - dolfin.EigenKrylovSolver
     - |todo|
   * - dolfin.EigenKrylovSolver.default_parameters
     - |todo|
   * - dolfin.EigenKrylovSolver.get_operator
     - |todo|
   * - dolfin.EigenKrylovSolver.id
     - |todo|
   * - dolfin.EigenKrylovSolver.label
     - |todo|
   * - dolfin.EigenKrylovSolver.methods
     - |todo|
   * - dolfin.EigenKrylovSolver.name
     - |todo|
   * - dolfin.EigenKrylovSolver.parameter_type
     - |todo|
   * - dolfin.EigenKrylovSolver.parameters
     - |todo|
   * - dolfin.EigenKrylovSolver.preconditioners
     - |todo|
   * - dolfin.EigenKrylovSolver.rename
     - |todo|
   * - dolfin.EigenKrylovSolver.set_operator
     - |todo|
   * - dolfin.EigenKrylovSolver.set_operators
     - |todo|
   * - dolfin.EigenKrylovSolver.solve
     - |todo|
   * - dolfin.EigenKrylovSolver.str
     - |todo|
   * - dolfin.EigenKrylovSolver.update_parameters
     - |todo|
   * - dolfin.EigenLUSolver
     - |todo|
   * - dolfin.EigenLUSolver.default_parameters
     - |todo|
   * - dolfin.EigenLUSolver.get_operator
     - |todo|
   * - dolfin.EigenLUSolver.id
     - |todo|
   * - dolfin.EigenLUSolver.label
     - |todo|
   * - dolfin.EigenLUSolver.methods
     - |todo|
   * - dolfin.EigenLUSolver.name
     - |todo|
   * - dolfin.EigenLUSolver.parameter_type
     - |todo|
   * - dolfin.EigenLUSolver.parameters
     - |todo|
   * - dolfin.EigenLUSolver.rename
     - |todo|
   * - dolfin.EigenLUSolver.set_operator
     - |todo|
   * - dolfin.EigenLUSolver.set_operators
     - |todo|
   * - dolfin.EigenLUSolver.solve
     - |todo|
   * - dolfin.EigenLUSolver.str
     - |todo|
   * - dolfin.EigenLUSolver.update_parameters
     - |todo|
   * - dolfin.EigenMatrix.add
     - |todo|
   * - dolfin.EigenMatrix.add_local
     - |todo|
   * - dolfin.EigenMatrix.assign
     - |todo|
   * - dolfin.EigenMatrix.compress
     - |todo|
   * - dolfin.EigenMatrix.data_view
     - |new|
   * - dolfin.EigenMatrix.ident_local
     - |todo|
   * - dolfin.EigenMatrix.is_symmetric
     - |todo|
   * - dolfin.EigenMatrix.mat
     - |todo|
   * - dolfin.EigenMatrix.resize
     - |todo|
   * - dolfin.EigenMatrix.set_local
     - |todo|
   * - dolfin.EigenMatrix.setrow
     - |todo|
   * - dolfin.EigenMatrix.shared_instance
     - |todo|
   * - dolfin.EigenMatrix.zero_local
     - |todo|
   * - dolfin.EigenVector.abs
     - |todo|
   * - dolfin.EigenVector.add
     - |todo|
   * - dolfin.EigenVector.data
     - |todo|
   * - dolfin.EigenVector.resize
     - |todo|
   * - dolfin.EigenVector.shared_instance
     - |todo|
   * - dolfin.EigenVector.vec
     - |todo|
   * - dolfin.EnrichedElement
     - |todo|
   * - dolfin.EnrichedElement.cell
     - |todo|
   * - dolfin.EnrichedElement.degree
     - |todo|
   * - dolfin.EnrichedElement.extract_component
     - |todo|
   * - dolfin.EnrichedElement.extract_reference_component
     - |todo|
   * - dolfin.EnrichedElement.extract_subelement_component
     - |todo|
   * - dolfin.EnrichedElement.extract_subelement_reference_component
     - |todo|
   * - dolfin.EnrichedElement.family
     - |todo|
   * - dolfin.EnrichedElement.is_cellwise_constant
     - |todo|
   * - dolfin.EnrichedElement.mapping
     - |todo|
   * - dolfin.EnrichedElement.num_sub_elements
     - |todo|
   * - dolfin.EnrichedElement.quadrature_scheme
     - |todo|
   * - dolfin.EnrichedElement.reconstruct
     - |todo|
   * - dolfin.EnrichedElement.reference_value_shape
     - |todo|
   * - dolfin.EnrichedElement.reference_value_size
     - |todo|
   * - dolfin.EnrichedElement.shortstr
     - |todo|
   * - dolfin.EnrichedElement.sobolev_space
     - |todo|
   * - dolfin.EnrichedElement.sub_elements
     - |todo|
   * - dolfin.EnrichedElement.symmetry
     - |todo|
   * - dolfin.EnrichedElement.value_shape
     - |todo|
   * - dolfin.EnrichedElement.value_size
     - |todo|
   * - dolfin.Equation
     - |todo|
   * - dolfin.Equation.is_linear
     - |todo|
   * - dolfin.Equation.lhs
     - |todo|
   * - dolfin.Equation.rhs
     - |todo|
   * - dolfin.Equation.rhs_int
     - |todo|
   * - dolfin.ErrorControl
     - |todo|
   * - dolfin.ErrorControl.child
     - |todo|
   * - dolfin.ErrorControl.clear_child
     - |todo|
   * - dolfin.ErrorControl.compute_cell_residual
     - |todo|
   * - dolfin.ErrorControl.compute_dual
     - |todo|
   * - dolfin.ErrorControl.compute_extrapolation
     - |todo|
   * - dolfin.ErrorControl.compute_facet_residual
     - |todo|
   * - dolfin.ErrorControl.compute_indicators
     - |todo|
   * - dolfin.ErrorControl.default_parameters
     - |todo|
   * - dolfin.ErrorControl.depth
     - |todo|
   * - dolfin.ErrorControl.estimate_error
     - |todo|
   * - dolfin.ErrorControl.has_child
     - |todo|
   * - dolfin.ErrorControl.has_parent
     - |todo|
   * - dolfin.ErrorControl.id
     - |todo|
   * - dolfin.ErrorControl.label
     - |todo|
   * - dolfin.ErrorControl.leaf_node
     - |todo|
   * - dolfin.ErrorControl.name
     - |todo|
   * - dolfin.ErrorControl.parameters
     - |todo|
   * - dolfin.ErrorControl.parent
     - |todo|
   * - dolfin.ErrorControl.rename
     - |todo|
   * - dolfin.ErrorControl.residual_representation
     - |todo|
   * - dolfin.ErrorControl.root_node
     - |todo|
   * - dolfin.ErrorControl.set_child
     - |todo|
   * - dolfin.ErrorControl.set_parent
     - |todo|
   * - dolfin.ErrorControl.str
     - |todo|
   * - dolfin.Event
     - |todo|
   * - dolfin.Event.count
     - |todo|
   * - dolfin.Event.maxcount
     - |todo|
   * - dolfin.ExplicitEuler
     - |todo|
   * - dolfin.ExplicitEuler.bcs
     - |todo|
   * - dolfin.ExplicitEuler.dolfin_stage_forms
     - |todo|
   * - dolfin.ExplicitEuler.dt
     - |todo|
   * - dolfin.ExplicitEuler.dt_stage_offset
     - |todo|
   * - dolfin.ExplicitEuler.id
     - |todo|
   * - dolfin.ExplicitEuler.implicit
     - |todo|
   * - dolfin.ExplicitEuler.jacobian_index
     - |todo|
   * - dolfin.ExplicitEuler.label
     - |todo|
   * - dolfin.ExplicitEuler.last_stage
     - |todo|
   * - dolfin.ExplicitEuler.name
     - |todo|
   * - dolfin.ExplicitEuler.order
     - |todo|
   * - dolfin.ExplicitEuler.parameters
     - |todo|
   * - dolfin.ExplicitEuler.rename
     - |todo|
   * - dolfin.ExplicitEuler.rhs_form
     - |todo|
   * - dolfin.ExplicitEuler.solution
     - |todo|
   * - dolfin.ExplicitEuler.stage_solutions
     - |todo|
   * - dolfin.ExplicitEuler.str
     - |todo|
   * - dolfin.ExplicitEuler.t
     - |todo|
   * - dolfin.ExplicitEuler.to_adm
     - |todo|
   * - dolfin.ExplicitEuler.to_tlm
     - |todo|
   * - dolfin.ExplicitEuler.ufl_stage_forms
     - |todo|
   * - dolfin.ExplicitMidPoint.bcs
     - |todo|
   * - dolfin.ExplicitMidPoint.dt_stage_offset
     - |todo|
   * - dolfin.ExplicitMidPoint.id
     - |todo|
   * - dolfin.ExplicitMidPoint.implicit
     - |todo|
   * - dolfin.ExplicitMidPoint.jacobian_index
     - |todo|
   * - dolfin.ExplicitMidPoint.label
     - |todo|
   * - dolfin.ExplicitMidPoint.name
     - |todo|
   * - dolfin.ExplicitMidPoint.parameters
     - |todo|
   * - dolfin.ExplicitMidPoint.rename
     - |todo|
   * - dolfin.ExplicitMidPoint.str
     - |todo|
   * - dolfin.Expression.T
     - |new|
   * - dolfin.Expression.compute_vertex_values
     - |new|
   * - dolfin.Expression.count
     - |new|
   * - dolfin.Expression.cpp_object
     - |new|
   * - dolfin.Expression.dx
     - |new|
   * - dolfin.Expression.evaluate
     - |new|
   * - dolfin.Expression.geometric_dimension
     - |new|
   * - dolfin.Expression.id
     - |new|
   * - dolfin.Expression.is_cellwise_constant
     - |new|
   * - dolfin.Expression.label
     - |new|
   * - dolfin.Expression.name
     - |new|
   * - dolfin.Expression.str
     - |todo|
   * - dolfin.Expression.ufl_disable_profiling
     - |new|
   * - dolfin.Expression.ufl_domain
     - |new|
   * - dolfin.Expression.ufl_domains
     - |new|
   * - dolfin.Expression.ufl_element
     - |new|
   * - dolfin.Expression.ufl_enable_profiling
     - |new|
   * - dolfin.Expression.ufl_free_indices
     - |new|
   * - dolfin.Expression.ufl_function_space
     - |new|
   * - dolfin.Expression.ufl_index_dimensions
     - |new|
   * - dolfin.Expression.ufl_operands
     - |new|
   * - dolfin.Expression.ufl_shape
     - |new|
   * - dolfin.Expression.value_dimension
     - |new|
   * - dolfin.Expression.value_rank
     - |new|
   * - dolfin.Expression.value_shape
     - |todo|
   * - dolfin.Extrapolation
     - |todo|
   * - dolfin.Extrapolation.extrapolate
     - |todo|
   * - dolfin.FENICS_EPS
     - |todo|
   * - dolfin.FENICS_EPS_LARGE
     - |todo|
   * - dolfin.FENICS_PI
     - |todo|
   * - dolfin.FENICS_SQRT_EPS
     - |todo|
   * - dolfin.Face.incident
     - |todo|
   * - dolfin.Face.init
     - |todo|
   * - dolfin.Face.mesh_id
     - |todo|
   * - dolfin.Face.owner
     - |todo|
   * - dolfin.Face.str
     - |todo|
   * - dolfin.FaceFunctionBool
     - |todo|
   * - dolfin.FaceFunctionBool.array
     - |todo|
   * - dolfin.FaceFunctionBool.child
     - |todo|
   * - dolfin.FaceFunctionBool.clear_child
     - |todo|
   * - dolfin.FaceFunctionBool.cpp_value_type
     - |todo|
   * - dolfin.FaceFunctionBool.depth
     - |todo|
   * - dolfin.FaceFunctionBool.dim
     - |todo|
   * - dolfin.FaceFunctionBool.empty
     - |todo|
   * - dolfin.FaceFunctionBool.has_child
     - |todo|
   * - dolfin.FaceFunctionBool.has_parent
     - |todo|
   * - dolfin.FaceFunctionBool.id
     - |todo|
   * - dolfin.FaceFunctionBool.init
     - |todo|
   * - dolfin.FaceFunctionBool.label
     - |todo|
   * - dolfin.FaceFunctionBool.leaf_node
     - |todo|
   * - dolfin.FaceFunctionBool.mesh
     - |todo|
   * - dolfin.FaceFunctionBool.name
     - |todo|
   * - dolfin.FaceFunctionBool.parameters
     - |todo|
   * - dolfin.FaceFunctionBool.parent
     - |todo|
   * - dolfin.FaceFunctionBool.rename
     - |todo|
   * - dolfin.FaceFunctionBool.root_node
     - |todo|
   * - dolfin.FaceFunctionBool.set_all
     - |todo|
   * - dolfin.FaceFunctionBool.set_child
     - |todo|
   * - dolfin.FaceFunctionBool.set_parent
     - |todo|
   * - dolfin.FaceFunctionBool.set_value
     - |todo|
   * - dolfin.FaceFunctionBool.set_values
     - |todo|
   * - dolfin.FaceFunctionBool.size
     - |todo|
   * - dolfin.FaceFunctionBool.str
     - |todo|
   * - dolfin.FaceFunctionBool.ufl_id
     - |todo|
   * - dolfin.FaceFunctionBool.where_equal
     - |todo|
   * - dolfin.FaceFunctionDouble
     - |todo|
   * - dolfin.FaceFunctionDouble.array
     - |todo|
   * - dolfin.FaceFunctionDouble.child
     - |todo|
   * - dolfin.FaceFunctionDouble.clear_child
     - |todo|
   * - dolfin.FaceFunctionDouble.cpp_value_type
     - |todo|
   * - dolfin.FaceFunctionDouble.depth
     - |todo|
   * - dolfin.FaceFunctionDouble.dim
     - |todo|
   * - dolfin.FaceFunctionDouble.empty
     - |todo|
   * - dolfin.FaceFunctionDouble.has_child
     - |todo|
   * - dolfin.FaceFunctionDouble.has_parent
     - |todo|
   * - dolfin.FaceFunctionDouble.id
     - |todo|
   * - dolfin.FaceFunctionDouble.init
     - |todo|
   * - dolfin.FaceFunctionDouble.label
     - |todo|
   * - dolfin.FaceFunctionDouble.leaf_node
     - |todo|
   * - dolfin.FaceFunctionDouble.mesh
     - |todo|
   * - dolfin.FaceFunctionDouble.name
     - |todo|
   * - dolfin.FaceFunctionDouble.parameters
     - |todo|
   * - dolfin.FaceFunctionDouble.parent
     - |todo|
   * - dolfin.FaceFunctionDouble.rename
     - |todo|
   * - dolfin.FaceFunctionDouble.root_node
     - |todo|
   * - dolfin.FaceFunctionDouble.set_all
     - |todo|
   * - dolfin.FaceFunctionDouble.set_child
     - |todo|
   * - dolfin.FaceFunctionDouble.set_parent
     - |todo|
   * - dolfin.FaceFunctionDouble.set_value
     - |todo|
   * - dolfin.FaceFunctionDouble.set_values
     - |todo|
   * - dolfin.FaceFunctionDouble.size
     - |todo|
   * - dolfin.FaceFunctionDouble.str
     - |todo|
   * - dolfin.FaceFunctionDouble.ufl_id
     - |todo|
   * - dolfin.FaceFunctionDouble.where_equal
     - |todo|
   * - dolfin.FaceFunctionInt
     - |todo|
   * - dolfin.FaceFunctionInt.array
     - |todo|
   * - dolfin.FaceFunctionInt.child
     - |todo|
   * - dolfin.FaceFunctionInt.clear_child
     - |todo|
   * - dolfin.FaceFunctionInt.cpp_value_type
     - |todo|
   * - dolfin.FaceFunctionInt.depth
     - |todo|
   * - dolfin.FaceFunctionInt.dim
     - |todo|
   * - dolfin.FaceFunctionInt.empty
     - |todo|
   * - dolfin.FaceFunctionInt.has_child
     - |todo|
   * - dolfin.FaceFunctionInt.has_parent
     - |todo|
   * - dolfin.FaceFunctionInt.id
     - |todo|
   * - dolfin.FaceFunctionInt.init
     - |todo|
   * - dolfin.FaceFunctionInt.label
     - |todo|
   * - dolfin.FaceFunctionInt.leaf_node
     - |todo|
   * - dolfin.FaceFunctionInt.mesh
     - |todo|
   * - dolfin.FaceFunctionInt.name
     - |todo|
   * - dolfin.FaceFunctionInt.parameters
     - |todo|
   * - dolfin.FaceFunctionInt.parent
     - |todo|
   * - dolfin.FaceFunctionInt.rename
     - |todo|
   * - dolfin.FaceFunctionInt.root_node
     - |todo|
   * - dolfin.FaceFunctionInt.set_all
     - |todo|
   * - dolfin.FaceFunctionInt.set_child
     - |todo|
   * - dolfin.FaceFunctionInt.set_parent
     - |todo|
   * - dolfin.FaceFunctionInt.set_value
     - |todo|
   * - dolfin.FaceFunctionInt.set_values
     - |todo|
   * - dolfin.FaceFunctionInt.size
     - |todo|
   * - dolfin.FaceFunctionInt.str
     - |todo|
   * - dolfin.FaceFunctionInt.ufl_id
     - |todo|
   * - dolfin.FaceFunctionInt.where_equal
     - |todo|
   * - dolfin.FaceFunctionSizet
     - |todo|
   * - dolfin.FaceFunctionSizet.array
     - |todo|
   * - dolfin.FaceFunctionSizet.child
     - |todo|
   * - dolfin.FaceFunctionSizet.clear_child
     - |todo|
   * - dolfin.FaceFunctionSizet.cpp_value_type
     - |todo|
   * - dolfin.FaceFunctionSizet.depth
     - |todo|
   * - dolfin.FaceFunctionSizet.dim
     - |todo|
   * - dolfin.FaceFunctionSizet.empty
     - |todo|
   * - dolfin.FaceFunctionSizet.has_child
     - |todo|
   * - dolfin.FaceFunctionSizet.has_parent
     - |todo|
   * - dolfin.FaceFunctionSizet.id
     - |todo|
   * - dolfin.FaceFunctionSizet.init
     - |todo|
   * - dolfin.FaceFunctionSizet.label
     - |todo|
   * - dolfin.FaceFunctionSizet.leaf_node
     - |todo|
   * - dolfin.FaceFunctionSizet.mesh
     - |todo|
   * - dolfin.FaceFunctionSizet.name
     - |todo|
   * - dolfin.FaceFunctionSizet.parameters
     - |todo|
   * - dolfin.FaceFunctionSizet.parent
     - |todo|
   * - dolfin.FaceFunctionSizet.rename
     - |todo|
   * - dolfin.FaceFunctionSizet.root_node
     - |todo|
   * - dolfin.FaceFunctionSizet.set_all
     - |todo|
   * - dolfin.FaceFunctionSizet.set_child
     - |todo|
   * - dolfin.FaceFunctionSizet.set_parent
     - |todo|
   * - dolfin.FaceFunctionSizet.set_value
     - |todo|
   * - dolfin.FaceFunctionSizet.set_values
     - |todo|
   * - dolfin.FaceFunctionSizet.size
     - |todo|
   * - dolfin.FaceFunctionSizet.str
     - |todo|
   * - dolfin.FaceFunctionSizet.ufl_id
     - |todo|
   * - dolfin.FaceFunctionSizet.where_equal
     - |todo|
   * - dolfin.Facet.distance
     - |todo|
   * - dolfin.Facet.incident
     - |todo|
   * - dolfin.Facet.init
     - |todo|
   * - dolfin.Facet.mesh_id
     - |todo|
   * - dolfin.Facet.owner
     - |todo|
   * - dolfin.Facet.squared_distance
     - |todo|
   * - dolfin.Facet.str
     - |todo|
   * - dolfin.FacetArea.cpp_object
     - |new|
   * - dolfin.FacetArea.eval
     - |todo|
   * - dolfin.FacetArea.eval_cell
     - |todo|
   * - dolfin.FacetArea.get_generic_function
     - |todo|
   * - dolfin.FacetArea.get_property
     - |todo|
   * - dolfin.FacetArea.parameters
     - |todo|
   * - dolfin.FacetArea.rename
     - |todo|
   * - dolfin.FacetArea.restrict
     - |todo|
   * - dolfin.FacetArea.set_generic_function
     - |todo|
   * - dolfin.FacetArea.set_property
     - |todo|
   * - dolfin.FacetArea.str
     - |todo|
   * - dolfin.FacetArea.update
     - |todo|
   * - dolfin.FacetArea.value_shape
     - |todo|
   * - dolfin.FacetArea.value_size
     - |todo|
   * - dolfin.FacetCell
     - |todo|
   * - dolfin.FacetCell.cell_normal
     - |todo|
   * - dolfin.FacetCell.circumradius
     - |todo|
   * - dolfin.FacetCell.collides
     - |todo|
   * - dolfin.FacetCell.contains
     - |todo|
   * - dolfin.FacetCell.dim
     - |todo|
   * - dolfin.FacetCell.distance
     - |todo|
   * - dolfin.FacetCell.entities
     - |todo|
   * - dolfin.FacetCell.facet_area
     - |todo|
   * - dolfin.FacetCell.facet_index
     - |todo|
   * - dolfin.FacetCell.get_cell_data
     - |todo|
   * - dolfin.FacetCell.get_cell_topology
     - |todo|
   * - dolfin.FacetCell.get_coordinate_dofs
     - |todo|
   * - dolfin.FacetCell.get_vertex_coordinates
     - |todo|
   * - dolfin.FacetCell.global_index
     - |todo|
   * - dolfin.FacetCell.h
     - |todo|
   * - dolfin.FacetCell.incident
     - |todo|
   * - dolfin.FacetCell.index
     - |todo|
   * - dolfin.FacetCell.init
     - |todo|
   * - dolfin.FacetCell.inradius
     - |todo|
   * - dolfin.FacetCell.intersection
     - |todo|
   * - dolfin.FacetCell.is_ghost
     - |todo|
   * - dolfin.FacetCell.is_shared
     - |todo|
   * - dolfin.FacetCell.mesh
     - |todo|
   * - dolfin.FacetCell.mesh_id
     - |todo|
   * - dolfin.FacetCell.midpoint
     - |todo|
   * - dolfin.FacetCell.normal
     - |todo|
   * - dolfin.FacetCell.num_entities
     - |todo|
   * - dolfin.FacetCell.num_global_entities
     - |todo|
   * - dolfin.FacetCell.num_vertices
     - |todo|
   * - dolfin.FacetCell.order
     - |todo|
   * - dolfin.FacetCell.ordered
     - |todo|
   * - dolfin.FacetCell.orientation
     - |todo|
   * - dolfin.FacetCell.owner
     - |todo|
   * - dolfin.FacetCell.radius_ratio
     - |todo|
   * - dolfin.FacetCell.sharing_processes
     - |todo|
   * - dolfin.FacetCell.squared_distance
     - |todo|
   * - dolfin.FacetCell.str
     - |todo|
   * - dolfin.FacetCell.type
     - |todo|
   * - dolfin.FacetCell.volume
     - |todo|
   * - dolfin.FacetElement
     - |todo|
   * - dolfin.FacetFunctionBool
     - |todo|
   * - dolfin.FacetFunctionBool.array
     - |todo|
   * - dolfin.FacetFunctionBool.child
     - |todo|
   * - dolfin.FacetFunctionBool.clear_child
     - |todo|
   * - dolfin.FacetFunctionBool.cpp_value_type
     - |todo|
   * - dolfin.FacetFunctionBool.depth
     - |todo|
   * - dolfin.FacetFunctionBool.dim
     - |todo|
   * - dolfin.FacetFunctionBool.empty
     - |todo|
   * - dolfin.FacetFunctionBool.has_child
     - |todo|
   * - dolfin.FacetFunctionBool.has_parent
     - |todo|
   * - dolfin.FacetFunctionBool.id
     - |todo|
   * - dolfin.FacetFunctionBool.init
     - |todo|
   * - dolfin.FacetFunctionBool.label
     - |todo|
   * - dolfin.FacetFunctionBool.leaf_node
     - |todo|
   * - dolfin.FacetFunctionBool.mesh
     - |todo|
   * - dolfin.FacetFunctionBool.name
     - |todo|
   * - dolfin.FacetFunctionBool.parameters
     - |todo|
   * - dolfin.FacetFunctionBool.parent
     - |todo|
   * - dolfin.FacetFunctionBool.rename
     - |todo|
   * - dolfin.FacetFunctionBool.root_node
     - |todo|
   * - dolfin.FacetFunctionBool.set_all
     - |todo|
   * - dolfin.FacetFunctionBool.set_child
     - |todo|
   * - dolfin.FacetFunctionBool.set_parent
     - |todo|
   * - dolfin.FacetFunctionBool.set_value
     - |todo|
   * - dolfin.FacetFunctionBool.set_values
     - |todo|
   * - dolfin.FacetFunctionBool.size
     - |todo|
   * - dolfin.FacetFunctionBool.str
     - |todo|
   * - dolfin.FacetFunctionBool.ufl_id
     - |todo|
   * - dolfin.FacetFunctionBool.where_equal
     - |todo|
   * - dolfin.FacetFunctionDouble
     - |todo|
   * - dolfin.FacetFunctionDouble.array
     - |todo|
   * - dolfin.FacetFunctionDouble.child
     - |todo|
   * - dolfin.FacetFunctionDouble.clear_child
     - |todo|
   * - dolfin.FacetFunctionDouble.cpp_value_type
     - |todo|
   * - dolfin.FacetFunctionDouble.depth
     - |todo|
   * - dolfin.FacetFunctionDouble.dim
     - |todo|
   * - dolfin.FacetFunctionDouble.empty
     - |todo|
   * - dolfin.FacetFunctionDouble.has_child
     - |todo|
   * - dolfin.FacetFunctionDouble.has_parent
     - |todo|
   * - dolfin.FacetFunctionDouble.id
     - |todo|
   * - dolfin.FacetFunctionDouble.init
     - |todo|
   * - dolfin.FacetFunctionDouble.label
     - |todo|
   * - dolfin.FacetFunctionDouble.leaf_node
     - |todo|
   * - dolfin.FacetFunctionDouble.mesh
     - |todo|
   * - dolfin.FacetFunctionDouble.name
     - |todo|
   * - dolfin.FacetFunctionDouble.parameters
     - |todo|
   * - dolfin.FacetFunctionDouble.parent
     - |todo|
   * - dolfin.FacetFunctionDouble.rename
     - |todo|
   * - dolfin.FacetFunctionDouble.root_node
     - |todo|
   * - dolfin.FacetFunctionDouble.set_all
     - |todo|
   * - dolfin.FacetFunctionDouble.set_child
     - |todo|
   * - dolfin.FacetFunctionDouble.set_parent
     - |todo|
   * - dolfin.FacetFunctionDouble.set_value
     - |todo|
   * - dolfin.FacetFunctionDouble.set_values
     - |todo|
   * - dolfin.FacetFunctionDouble.size
     - |todo|
   * - dolfin.FacetFunctionDouble.str
     - |todo|
   * - dolfin.FacetFunctionDouble.ufl_id
     - |todo|
   * - dolfin.FacetFunctionDouble.where_equal
     - |todo|
   * - dolfin.FacetFunctionInt
     - |todo|
   * - dolfin.FacetFunctionInt.array
     - |todo|
   * - dolfin.FacetFunctionInt.child
     - |todo|
   * - dolfin.FacetFunctionInt.clear_child
     - |todo|
   * - dolfin.FacetFunctionInt.cpp_value_type
     - |todo|
   * - dolfin.FacetFunctionInt.depth
     - |todo|
   * - dolfin.FacetFunctionInt.dim
     - |todo|
   * - dolfin.FacetFunctionInt.empty
     - |todo|
   * - dolfin.FacetFunctionInt.has_child
     - |todo|
   * - dolfin.FacetFunctionInt.has_parent
     - |todo|
   * - dolfin.FacetFunctionInt.id
     - |todo|
   * - dolfin.FacetFunctionInt.init
     - |todo|
   * - dolfin.FacetFunctionInt.label
     - |todo|
   * - dolfin.FacetFunctionInt.leaf_node
     - |todo|
   * - dolfin.FacetFunctionInt.mesh
     - |todo|
   * - dolfin.FacetFunctionInt.name
     - |todo|
   * - dolfin.FacetFunctionInt.parameters
     - |todo|
   * - dolfin.FacetFunctionInt.parent
     - |todo|
   * - dolfin.FacetFunctionInt.rename
     - |todo|
   * - dolfin.FacetFunctionInt.root_node
     - |todo|
   * - dolfin.FacetFunctionInt.set_all
     - |todo|
   * - dolfin.FacetFunctionInt.set_child
     - |todo|
   * - dolfin.FacetFunctionInt.set_parent
     - |todo|
   * - dolfin.FacetFunctionInt.set_value
     - |todo|
   * - dolfin.FacetFunctionInt.set_values
     - |todo|
   * - dolfin.FacetFunctionInt.size
     - |todo|
   * - dolfin.FacetFunctionInt.str
     - |todo|
   * - dolfin.FacetFunctionInt.ufl_id
     - |todo|
   * - dolfin.FacetFunctionInt.where_equal
     - |todo|
   * - dolfin.FacetFunctionSizet
     - |todo|
   * - dolfin.FacetFunctionSizet.array
     - |todo|
   * - dolfin.FacetFunctionSizet.child
     - |todo|
   * - dolfin.FacetFunctionSizet.clear_child
     - |todo|
   * - dolfin.FacetFunctionSizet.cpp_value_type
     - |todo|
   * - dolfin.FacetFunctionSizet.depth
     - |todo|
   * - dolfin.FacetFunctionSizet.dim
     - |todo|
   * - dolfin.FacetFunctionSizet.empty
     - |todo|
   * - dolfin.FacetFunctionSizet.has_child
     - |todo|
   * - dolfin.FacetFunctionSizet.has_parent
     - |todo|
   * - dolfin.FacetFunctionSizet.id
     - |todo|
   * - dolfin.FacetFunctionSizet.init
     - |todo|
   * - dolfin.FacetFunctionSizet.label
     - |todo|
   * - dolfin.FacetFunctionSizet.leaf_node
     - |todo|
   * - dolfin.FacetFunctionSizet.mesh
     - |todo|
   * - dolfin.FacetFunctionSizet.name
     - |todo|
   * - dolfin.FacetFunctionSizet.parameters
     - |todo|
   * - dolfin.FacetFunctionSizet.parent
     - |todo|
   * - dolfin.FacetFunctionSizet.rename
     - |todo|
   * - dolfin.FacetFunctionSizet.root_node
     - |todo|
   * - dolfin.FacetFunctionSizet.set_all
     - |todo|
   * - dolfin.FacetFunctionSizet.set_child
     - |todo|
   * - dolfin.FacetFunctionSizet.set_parent
     - |todo|
   * - dolfin.FacetFunctionSizet.set_value
     - |todo|
   * - dolfin.FacetFunctionSizet.set_values
     - |todo|
   * - dolfin.FacetFunctionSizet.size
     - |todo|
   * - dolfin.FacetFunctionSizet.str
     - |todo|
   * - dolfin.FacetFunctionSizet.ufl_id
     - |todo|
   * - dolfin.FacetFunctionSizet.where_equal
     - |todo|
   * - dolfin.File.Type_binary
     - |todo|
   * - dolfin.File.Type_raw
     - |todo|
   * - dolfin.File.Type_svg
     - |todo|
   * - dolfin.File.Type_vtk
     - |todo|
   * - dolfin.File.Type_x3d
     - |todo|
   * - dolfin.File.Type_xml
     - |todo|
   * - dolfin.File.Type_xyz
     - |todo|
   * - dolfin.File.create_parent_path
     - |todo|
   * - dolfin.File.exists
     - |todo|
   * - dolfin.File.read
     - |new|
   * - dolfin.File.write
     - |new|
   * - dolfin.FiniteElementBase
     - |todo|
   * - dolfin.FiniteElementBase.cell
     - |todo|
   * - dolfin.FiniteElementBase.degree
     - |todo|
   * - dolfin.FiniteElementBase.extract_component
     - |todo|
   * - dolfin.FiniteElementBase.extract_reference_component
     - |todo|
   * - dolfin.FiniteElementBase.extract_subelement_component
     - |todo|
   * - dolfin.FiniteElementBase.extract_subelement_reference_component
     - |todo|
   * - dolfin.FiniteElementBase.family
     - |todo|
   * - dolfin.FiniteElementBase.is_cellwise_constant
     - |todo|
   * - dolfin.FiniteElementBase.mapping
     - |todo|
   * - dolfin.FiniteElementBase.num_sub_elements
     - |todo|
   * - dolfin.FiniteElementBase.quadrature_scheme
     - |todo|
   * - dolfin.FiniteElementBase.reference_value_shape
     - |todo|
   * - dolfin.FiniteElementBase.reference_value_size
     - |todo|
   * - dolfin.FiniteElementBase.sub_elements
     - |todo|
   * - dolfin.FiniteElementBase.symmetry
     - |todo|
   * - dolfin.FiniteElementBase.value_shape
     - |todo|
   * - dolfin.FiniteElementBase.value_size
     - |todo|
   * - dolfin.Form.cell_domains
     - |todo|
   * - dolfin.Form.check
     - |todo|
   * - dolfin.Form.child
     - |todo|
   * - dolfin.Form.clear_child
     - |todo|
   * - dolfin.Form.coefficient
     - |todo|
   * - dolfin.Form.coefficient_name
     - |todo|
   * - dolfin.Form.coefficient_number
     - |todo|
   * - dolfin.Form.coefficients
     - |todo|
   * - dolfin.Form.coloring
     - |todo|
   * - dolfin.Form.depth
     - |todo|
   * - dolfin.Form.exterior_facet_domains
     - |todo|
   * - dolfin.Form.function_space
     - |todo|
   * - dolfin.Form.function_spaces
     - |todo|
   * - dolfin.Form.has_child
     - |todo|
   * - dolfin.Form.has_parent
     - |todo|
   * - dolfin.Form.interior_facet_domains
     - |todo|
   * - dolfin.Form.leaf_node
     - |todo|
   * - dolfin.Form.parent
     - |todo|
   * - dolfin.Form.root_node
     - |todo|
   * - dolfin.Form.set_child
     - |todo|
   * - dolfin.Form.set_coefficients
     - |todo|
   * - dolfin.Form.set_parent
     - |todo|
   * - dolfin.Form.set_some_coefficients
     - |todo|
   * - dolfin.Form.ufc_form
     - |todo|
   * - dolfin.Form.vertex_domains
     - |todo|
   * - dolfin.ForwardEuler.bcs
     - |todo|
   * - dolfin.ForwardEuler.dt_stage_offset
     - |todo|
   * - dolfin.ForwardEuler.id
     - |todo|
   * - dolfin.ForwardEuler.implicit
     - |todo|
   * - dolfin.ForwardEuler.jacobian_index
     - |todo|
   * - dolfin.ForwardEuler.label
     - |todo|
   * - dolfin.ForwardEuler.name
     - |todo|
   * - dolfin.ForwardEuler.parameters
     - |todo|
   * - dolfin.ForwardEuler.rename
     - |todo|
   * - dolfin.ForwardEuler.str
     - |todo|
   * - dolfin.Function.child
     - |todo|
   * - dolfin.Function.clear_child
     - |todo|
   * - dolfin.Function.cpp_object
     - |new|
   * - dolfin.Function.depth
     - |todo|
   * - dolfin.Function.has_child
     - |todo|
   * - dolfin.Function.has_parent
     - |todo|
   * - dolfin.Function.label
     - |todo|
   * - dolfin.Function.parameters
     - |todo|
   * - dolfin.Function.parent
     - |todo|
   * - dolfin.Function.restrict
     - |todo|
   * - dolfin.Function.set_child
     - |todo|
   * - dolfin.Function.set_parent
     - |todo|
   * - dolfin.Function.str
     - |todo|
   * - dolfin.Function.update
     - |todo|
   * - dolfin.Function.value_size
     - |todo|
   * - dolfin.FunctionAXPY.Direction
     - |new|
   * - dolfin.FunctionAXPY.Direction_ADD_ADD
     - |todo|
   * - dolfin.FunctionAXPY.Direction_ADD_SUB
     - |todo|
   * - dolfin.FunctionAXPY.Direction_SUB_ADD
     - |todo|
   * - dolfin.FunctionAXPY.Direction_SUB_SUB
     - |todo|
   * - dolfin.FunctionAssigner.num_assigning_functions
     - |todo|
   * - dolfin.FunctionAssigner.num_receiving_functions
     - |todo|
   * - dolfin.FunctionSpace.assign
     - |todo|
   * - dolfin.FunctionSpace.child
     - |todo|
   * - dolfin.FunctionSpace.clear_child
     - |todo|
   * - dolfin.FunctionSpace.depth
     - |todo|
   * - dolfin.FunctionSpace.has_cell
     - |todo|
   * - dolfin.FunctionSpace.has_child
     - |todo|
   * - dolfin.FunctionSpace.has_element
     - |todo|
   * - dolfin.FunctionSpace.has_parent
     - |todo|
   * - dolfin.FunctionSpace.interpolate
     - |todo|
   * - dolfin.FunctionSpace.label
     - |todo|
   * - dolfin.FunctionSpace.leaf_node
     - |todo|
   * - dolfin.FunctionSpace.name
     - |todo|
   * - dolfin.FunctionSpace.parameters
     - |todo|
   * - dolfin.FunctionSpace.parent
     - |todo|
   * - dolfin.FunctionSpace.print_dofmap
     - |todo|
   * - dolfin.FunctionSpace.rename
     - |todo|
   * - dolfin.FunctionSpace.root_node
     - |todo|
   * - dolfin.FunctionSpace.set_child
     - |todo|
   * - dolfin.FunctionSpace.set_parent
     - |todo|
   * - dolfin.FunctionSpace.split
     - |todo|
   * - dolfin.FunctionSpace.str
     - |todo|
   * - dolfin.GRL1.bcs
     - |todo|
   * - dolfin.GRL1.dt_stage_offset
     - |todo|
   * - dolfin.GRL1.id
     - |todo|
   * - dolfin.GRL1.implicit
     - |todo|
   * - dolfin.GRL1.jacobian_index
     - |todo|
   * - dolfin.GRL1.label
     - |todo|
   * - dolfin.GRL1.name
     - |todo|
   * - dolfin.GRL1.parameters
     - |todo|
   * - dolfin.GRL1.rename
     - |todo|
   * - dolfin.GRL1.str
     - |todo|
   * - dolfin.GRL2.bcs
     - |todo|
   * - dolfin.GRL2.dt_stage_offset
     - |todo|
   * - dolfin.GRL2.id
     - |todo|
   * - dolfin.GRL2.implicit
     - |todo|
   * - dolfin.GRL2.jacobian_index
     - |todo|
   * - dolfin.GRL2.label
     - |todo|
   * - dolfin.GRL2.name
     - |todo|
   * - dolfin.GRL2.parameters
     - |todo|
   * - dolfin.GRL2.rename
     - |todo|
   * - dolfin.GRL2.str
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.adapt_problem
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.adaptive_data
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.default_parameters
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.evaluate_goal
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.extract_bcs
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.id
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.label
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.name
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.parameters
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.rename
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.solve
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.solve_primal
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.str
     - |todo|
   * - dolfin.GenericAdaptiveVariationalSolver.summary
     - |todo|
   * - dolfin.GenericBoundingBoxTree
     - |todo|
   * - dolfin.GenericBoundingBoxTree.build
     - |todo|
   * - dolfin.GenericBoundingBoxTree.compute_closest_entity
     - |todo|
   * - dolfin.GenericBoundingBoxTree.compute_closest_point
     - |todo|
   * - dolfin.GenericBoundingBoxTree.compute_collisions
     - |todo|
   * - dolfin.GenericBoundingBoxTree.compute_entity_collisions
     - |todo|
   * - dolfin.GenericBoundingBoxTree.compute_first_collision
     - |todo|
   * - dolfin.GenericBoundingBoxTree.compute_first_entity_collision
     - |todo|
   * - dolfin.GenericBoundingBoxTree.compute_process_collisions
     - |todo|
   * - dolfin.GenericBoundingBoxTree.create
     - |todo|
   * - dolfin.GenericBoundingBoxTree.str
     - |todo|
   * - dolfin.GenericDofMap
     - |todo|
   * - dolfin.GenericDofMap.block_size
     - |todo|
   * - dolfin.GenericDofMap.cell_dimension
     - |todo|
   * - dolfin.GenericDofMap.cell_dofs
     - |todo|
   * - dolfin.GenericDofMap.clear_sub_map_data
     - |todo|
   * - dolfin.GenericDofMap.collapse
     - |todo|
   * - dolfin.GenericDofMap.constrained_domain
     - |todo|
   * - dolfin.GenericDofMap.copy
     - |todo|
   * - dolfin.GenericDofMap.create
     - |todo|
   * - dolfin.GenericDofMap.dofs
     - |todo|
   * - dolfin.GenericDofMap.entity_closure_dofs
     - |todo|
   * - dolfin.GenericDofMap.entity_dofs
     - |todo|
   * - dolfin.GenericDofMap.extract_sub_dofmap
     - |todo|
   * - dolfin.GenericDofMap.global_dimension
     - |todo|
   * - dolfin.GenericDofMap.id
     - |todo|
   * - dolfin.GenericDofMap.index_map
     - |todo|
   * - dolfin.GenericDofMap.is_view
     - |todo|
   * - dolfin.GenericDofMap.label
     - |todo|
   * - dolfin.GenericDofMap.local_to_global_index
     - |todo|
   * - dolfin.GenericDofMap.local_to_global_unowned
     - |todo|
   * - dolfin.GenericDofMap.max_cell_dimension
     - |todo|
   * - dolfin.GenericDofMap.max_element_dofs
     - |todo|
   * - dolfin.GenericDofMap.name
     - |todo|
   * - dolfin.GenericDofMap.neighbours
     - |todo|
   * - dolfin.GenericDofMap.num_element_dofs
     - |todo|
   * - dolfin.GenericDofMap.num_entity_closure_dofs
     - |todo|
   * - dolfin.GenericDofMap.num_entity_dofs
     - |todo|
   * - dolfin.GenericDofMap.num_facet_dofs
     - |todo|
   * - dolfin.GenericDofMap.off_process_owner
     - |todo|
   * - dolfin.GenericDofMap.ownership_range
     - |todo|
   * - dolfin.GenericDofMap.parameters
     - |todo|
   * - dolfin.GenericDofMap.rename
     - |todo|
   * - dolfin.GenericDofMap.set
     - |todo|
   * - dolfin.GenericDofMap.shared_nodes
     - |todo|
   * - dolfin.GenericDofMap.str
     - |todo|
   * - dolfin.GenericDofMap.tabulate_entity_closure_dofs
     - |todo|
   * - dolfin.GenericDofMap.tabulate_entity_dofs
     - |todo|
   * - dolfin.GenericDofMap.tabulate_facet_dofs
     - |todo|
   * - dolfin.GenericDofMap.tabulate_global_dofs
     - |todo|
   * - dolfin.GenericDofMap.tabulate_local_to_global_dofs
     - |todo|
   * - dolfin.GenericFile
     - |todo|
   * - dolfin.GenericFile.name
     - |todo|
   * - dolfin.GenericFile.read
     - |todo|
   * - dolfin.GenericFile.write
     - |todo|
   * - dolfin.GenericFunction
     - |todo|
   * - dolfin.GenericFunction.compute_vertex_values
     - |todo|
   * - dolfin.GenericFunction.eval
     - |todo|
   * - dolfin.GenericFunction.eval_cell
     - |todo|
   * - dolfin.GenericFunction.id
     - |todo|
   * - dolfin.GenericFunction.label
     - |todo|
   * - dolfin.GenericFunction.name
     - |todo|
   * - dolfin.GenericFunction.parameters
     - |todo|
   * - dolfin.GenericFunction.rename
     - |todo|
   * - dolfin.GenericFunction.restrict
     - |todo|
   * - dolfin.GenericFunction.str
     - |todo|
   * - dolfin.GenericFunction.update
     - |todo|
   * - dolfin.GenericFunction.value_dimension
     - |todo|
   * - dolfin.GenericFunction.value_rank
     - |todo|
   * - dolfin.GenericFunction.value_shape
     - |todo|
   * - dolfin.GenericFunction.value_size
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.create_krylov_solver
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.create_layout
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.create_linear_operator
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.create_lu_solver
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.create_matrix
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.create_vector
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.krylov_solver_methods
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.krylov_solver_preconditioners
     - |todo|
   * - dolfin.GenericLinearAlgebraFactory.lu_solver_methods
     - |todo|
   * - dolfin.GenericLinearOperator
     - |todo|
   * - dolfin.GenericLinearOperator.id
     - |todo|
   * - dolfin.GenericLinearOperator.label
     - |todo|
   * - dolfin.GenericLinearOperator.mpi_comm
     - |todo|
   * - dolfin.GenericLinearOperator.mult
     - |todo|
   * - dolfin.GenericLinearOperator.name
     - |todo|
   * - dolfin.GenericLinearOperator.parameters
     - |todo|
   * - dolfin.GenericLinearOperator.rename
     - |todo|
   * - dolfin.GenericLinearOperator.shared_instance
     - |todo|
   * - dolfin.GenericLinearOperator.size
     - |todo|
   * - dolfin.GenericLinearOperator.str
     - |todo|
   * - dolfin.GenericLinearSolver
     - |todo|
   * - dolfin.GenericLinearSolver.id
     - |todo|
   * - dolfin.GenericLinearSolver.label
     - |todo|
   * - dolfin.GenericLinearSolver.name
     - |todo|
   * - dolfin.GenericLinearSolver.parameter_type
     - |todo|
   * - dolfin.GenericLinearSolver.parameters
     - |todo|
   * - dolfin.GenericLinearSolver.rename
     - |todo|
   * - dolfin.GenericLinearSolver.set_operator
     - |todo|
   * - dolfin.GenericLinearSolver.set_operators
     - |todo|
   * - dolfin.GenericLinearSolver.solve
     - |todo|
   * - dolfin.GenericLinearSolver.str
     - |todo|
   * - dolfin.GenericLinearSolver.update_parameters
     - |todo|
   * - dolfin.GenericMatrix
     - |todo|
   * - dolfin.GenericMatrix.add
     - |todo|
   * - dolfin.GenericMatrix.add_local
     - |todo|
   * - dolfin.GenericMatrix.apply
     - |todo|
   * - dolfin.GenericMatrix.array
     - |todo|
   * - dolfin.GenericMatrix.assign
     - |todo|
   * - dolfin.GenericMatrix.axpy
     - |todo|
   * - dolfin.GenericMatrix.copy
     - |todo|
   * - dolfin.GenericMatrix.empty
     - |todo|
   * - dolfin.GenericMatrix.factory
     - |todo|
   * - dolfin.GenericMatrix.get
     - |todo|
   * - dolfin.GenericMatrix.get_diagonal
     - |todo|
   * - dolfin.GenericMatrix.getrow
     - |todo|
   * - dolfin.GenericMatrix.id
     - |todo|
   * - dolfin.GenericMatrix.ident
     - |todo|
   * - dolfin.GenericMatrix.ident_local
     - |todo|
   * - dolfin.GenericMatrix.ident_zeros
     - |todo|
   * - dolfin.GenericMatrix.init
     - |todo|
   * - dolfin.GenericMatrix.init_vector
     - |todo|
   * - dolfin.GenericMatrix.is_symmetric
     - |todo|
   * - dolfin.GenericMatrix.label
     - |todo|
   * - dolfin.GenericMatrix.local_range
     - |todo|
   * - dolfin.GenericMatrix.mpi_comm
     - |todo|
   * - dolfin.GenericMatrix.mult
     - |todo|
   * - dolfin.GenericMatrix.name
     - |todo|
   * - dolfin.GenericMatrix.nnz
     - |todo|
   * - dolfin.GenericMatrix.norm
     - |todo|
   * - dolfin.GenericMatrix.parameters
     - |todo|
   * - dolfin.GenericMatrix.rank
     - |todo|
   * - dolfin.GenericMatrix.rename
     - |todo|
   * - dolfin.GenericMatrix.set
     - |todo|
   * - dolfin.GenericMatrix.set_diagonal
     - |todo|
   * - dolfin.GenericMatrix.set_local
     - |todo|
   * - dolfin.GenericMatrix.setrow
     - |todo|
   * - dolfin.GenericMatrix.shared_instance
     - |todo|
   * - dolfin.GenericMatrix.size
     - |todo|
   * - dolfin.GenericMatrix.str
     - |todo|
   * - dolfin.GenericMatrix.transpmult
     - |todo|
   * - dolfin.GenericMatrix.zero
     - |todo|
   * - dolfin.GenericMatrix.zero_local
     - |todo|
   * - dolfin.GenericTensor
     - |todo|
   * - dolfin.GenericTensor.add
     - |todo|
   * - dolfin.GenericTensor.add_local
     - |todo|
   * - dolfin.GenericTensor.apply
     - |todo|
   * - dolfin.GenericTensor.empty
     - |todo|
   * - dolfin.GenericTensor.factory
     - |todo|
   * - dolfin.GenericTensor.id
     - |todo|
   * - dolfin.GenericTensor.init
     - |todo|
   * - dolfin.GenericTensor.label
     - |todo|
   * - dolfin.GenericTensor.local_range
     - |todo|
   * - dolfin.GenericTensor.mpi_comm
     - |todo|
   * - dolfin.GenericTensor.name
     - |todo|
   * - dolfin.GenericTensor.parameters
     - |todo|
   * - dolfin.GenericTensor.rank
     - |todo|
   * - dolfin.GenericTensor.rename
     - |todo|
   * - dolfin.GenericTensor.set_local
     - |todo|
   * - dolfin.GenericTensor.shared_instance
     - |todo|
   * - dolfin.GenericTensor.size
     - |todo|
   * - dolfin.GenericTensor.str
     - |todo|
   * - dolfin.GenericTensor.zero
     - |todo|
   * - dolfin.GenericVector.abs
     - |todo|
   * - dolfin.GenericVector.add
     - |todo|
   * - dolfin.GenericVector.shared_instance
     - |todo|
   * - dolfin.GlobalParameters
     - |todo|
   * - dolfin.GlobalParameters.add
     - |todo|
   * - dolfin.GlobalParameters.add_unset
     - |todo|
   * - dolfin.GlobalParameters.assign
     - |todo|
   * - dolfin.GlobalParameters.clear
     - |todo|
   * - dolfin.GlobalParameters.copy
     - |todo|
   * - dolfin.GlobalParameters.default_parameters
     - |todo|
   * - dolfin.GlobalParameters.find_parameter
     - |todo|
   * - dolfin.GlobalParameters.find_parameter_set
     - |todo|
   * - dolfin.GlobalParameters.get
     - |todo|
   * - dolfin.GlobalParameters.get_range
     - |todo|
   * - dolfin.GlobalParameters.has_key
     - |todo|
   * - dolfin.GlobalParameters.has_parameter
     - |todo|
   * - dolfin.GlobalParameters.has_parameter_set
     - |todo|
   * - dolfin.GlobalParameters.items
     - |todo|
   * - dolfin.GlobalParameters.iterdata
     - |todo|
   * - dolfin.GlobalParameters.iteritems
     - |todo|
   * - dolfin.GlobalParameters.iterkeys
     - |todo|
   * - dolfin.GlobalParameters.itervalues
     - |todo|
   * - dolfin.GlobalParameters.keys
     - |todo|
   * - dolfin.GlobalParameters.name
     - |todo|
   * - dolfin.GlobalParameters.option_string
     - |todo|
   * - dolfin.GlobalParameters.parse
     - |todo|
   * - dolfin.GlobalParameters.remove
     - |todo|
   * - dolfin.GlobalParameters.rename
     - |todo|
   * - dolfin.GlobalParameters.set_range
     - |todo|
   * - dolfin.GlobalParameters.size
     - |todo|
   * - dolfin.GlobalParameters.str
     - |todo|
   * - dolfin.GlobalParameters.to_dict
     - |todo|
   * - dolfin.GlobalParameters.update
     - |todo|
   * - dolfin.GlobalParameters.values
     - |todo|
   * - dolfin.Graph
     - |todo|
   * - dolfin.H1
     - |todo|
   * - dolfin.H2
     - |todo|
   * - dolfin.HCurl
     - |todo|
   * - dolfin.HCurlElement
     - |todo|
   * - dolfin.HCurlElement.cell
     - |todo|
   * - dolfin.HCurlElement.degree
     - |todo|
   * - dolfin.HCurlElement.extract_component
     - |todo|
   * - dolfin.HCurlElement.extract_reference_component
     - |todo|
   * - dolfin.HCurlElement.extract_subelement_component
     - |todo|
   * - dolfin.HCurlElement.extract_subelement_reference_component
     - |todo|
   * - dolfin.HCurlElement.family
     - |todo|
   * - dolfin.HCurlElement.is_cellwise_constant
     - |todo|
   * - dolfin.HCurlElement.mapping
     - |todo|
   * - dolfin.HCurlElement.num_sub_elements
     - |todo|
   * - dolfin.HCurlElement.quadrature_scheme
     - |todo|
   * - dolfin.HCurlElement.reconstruct
     - |todo|
   * - dolfin.HCurlElement.reference_value_shape
     - |todo|
   * - dolfin.HCurlElement.reference_value_size
     - |todo|
   * - dolfin.HCurlElement.shortstr
     - |todo|
   * - dolfin.HCurlElement.sobolev_space
     - |todo|
   * - dolfin.HCurlElement.sub_elements
     - |todo|
   * - dolfin.HCurlElement.symmetry
     - |todo|
   * - dolfin.HCurlElement.value_shape
     - |todo|
   * - dolfin.HCurlElement.value_size
     - |todo|
   * - dolfin.HDF5Attribute
     - |todo|
   * - dolfin.HDF5Attribute.exists
     - |todo|
   * - dolfin.HDF5Attribute.items
     - |todo|
   * - dolfin.HDF5Attribute.keys
     - |todo|
   * - dolfin.HDF5Attribute.list_attributes
     - |todo|
   * - dolfin.HDF5Attribute.str
     - |todo|
   * - dolfin.HDF5Attribute.to_dict
     - |todo|
   * - dolfin.HDF5Attribute.type_str
     - |todo|
   * - dolfin.HDF5Attribute.values
     - |todo|
   * - dolfin.HDF5File.flush
     - |todo|
   * - dolfin.HDF5File.h5_id
     - |todo|
   * - dolfin.HDF5File.str
     - |todo|
   * - dolfin.HDiv
     - |todo|
   * - dolfin.HDivElement
     - |todo|
   * - dolfin.HDivElement.cell
     - |todo|
   * - dolfin.HDivElement.degree
     - |todo|
   * - dolfin.HDivElement.extract_component
     - |todo|
   * - dolfin.HDivElement.extract_reference_component
     - |todo|
   * - dolfin.HDivElement.extract_subelement_component
     - |todo|
   * - dolfin.HDivElement.extract_subelement_reference_component
     - |todo|
   * - dolfin.HDivElement.family
     - |todo|
   * - dolfin.HDivElement.is_cellwise_constant
     - |todo|
   * - dolfin.HDivElement.mapping
     - |todo|
   * - dolfin.HDivElement.num_sub_elements
     - |todo|
   * - dolfin.HDivElement.quadrature_scheme
     - |todo|
   * - dolfin.HDivElement.reconstruct
     - |todo|
   * - dolfin.HDivElement.reference_value_shape
     - |todo|
   * - dolfin.HDivElement.reference_value_size
     - |todo|
   * - dolfin.HDivElement.shortstr
     - |todo|
   * - dolfin.HDivElement.sobolev_space
     - |todo|
   * - dolfin.HDivElement.sub_elements
     - |todo|
   * - dolfin.HDivElement.symmetry
     - |todo|
   * - dolfin.HDivElement.value_shape
     - |todo|
   * - dolfin.HDivElement.value_size
     - |todo|
   * - dolfin.HierarchicalDirichletBC
     - |todo|
   * - dolfin.HierarchicalDirichletBC.child
     - |todo|
   * - dolfin.HierarchicalDirichletBC.clear_child
     - |todo|
   * - dolfin.HierarchicalDirichletBC.depth
     - |todo|
   * - dolfin.HierarchicalDirichletBC.has_child
     - |todo|
   * - dolfin.HierarchicalDirichletBC.has_parent
     - |todo|
   * - dolfin.HierarchicalDirichletBC.leaf_node
     - |todo|
   * - dolfin.HierarchicalDirichletBC.parent
     - |todo|
   * - dolfin.HierarchicalDirichletBC.root_node
     - |todo|
   * - dolfin.HierarchicalDirichletBC.set_child
     - |todo|
   * - dolfin.HierarchicalDirichletBC.set_parent
     - |todo|
   * - dolfin.HierarchicalErrorControl
     - |todo|
   * - dolfin.HierarchicalErrorControl.child
     - |todo|
   * - dolfin.HierarchicalErrorControl.clear_child
     - |todo|
   * - dolfin.HierarchicalErrorControl.depth
     - |todo|
   * - dolfin.HierarchicalErrorControl.has_child
     - |todo|
   * - dolfin.HierarchicalErrorControl.has_parent
     - |todo|
   * - dolfin.HierarchicalErrorControl.leaf_node
     - |todo|
   * - dolfin.HierarchicalErrorControl.parent
     - |todo|
   * - dolfin.HierarchicalErrorControl.root_node
     - |todo|
   * - dolfin.HierarchicalErrorControl.set_child
     - |todo|
   * - dolfin.HierarchicalErrorControl.set_parent
     - |todo|
   * - dolfin.HierarchicalForm
     - |todo|
   * - dolfin.HierarchicalForm.child
     - |todo|
   * - dolfin.HierarchicalForm.clear_child
     - |todo|
   * - dolfin.HierarchicalForm.depth
     - |todo|
   * - dolfin.HierarchicalForm.has_child
     - |todo|
   * - dolfin.HierarchicalForm.has_parent
     - |todo|
   * - dolfin.HierarchicalForm.leaf_node
     - |todo|
   * - dolfin.HierarchicalForm.parent
     - |todo|
   * - dolfin.HierarchicalForm.root_node
     - |todo|
   * - dolfin.HierarchicalForm.set_child
     - |todo|
   * - dolfin.HierarchicalForm.set_parent
     - |todo|
   * - dolfin.HierarchicalFunction
     - |todo|
   * - dolfin.HierarchicalFunction.clear_child
     - |todo|
   * - dolfin.HierarchicalFunction.depth
     - |todo|
   * - dolfin.HierarchicalFunction.has_child
     - |todo|
   * - dolfin.HierarchicalFunction.has_parent
     - |todo|
   * - dolfin.HierarchicalFunction.set_child
     - |todo|
   * - dolfin.HierarchicalFunction.set_parent
     - |todo|
   * - dolfin.HierarchicalFunctionSpace
     - |todo|
   * - dolfin.HierarchicalFunctionSpace.clear_child
     - |todo|
   * - dolfin.HierarchicalFunctionSpace.depth
     - |todo|
   * - dolfin.HierarchicalFunctionSpace.has_child
     - |todo|
   * - dolfin.HierarchicalFunctionSpace.has_parent
     - |todo|
   * - dolfin.HierarchicalFunctionSpace.set_child
     - |todo|
   * - dolfin.HierarchicalFunctionSpace.set_parent
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.child
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.clear_child
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.depth
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.has_child
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.has_parent
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.leaf_node
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.parent
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.root_node
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.set_child
     - |todo|
   * - dolfin.HierarchicalLinearVariationalProblem.set_parent
     - |todo|
   * - dolfin.HierarchicalMesh
     - |todo|
   * - dolfin.HierarchicalMesh.clear_child
     - |todo|
   * - dolfin.HierarchicalMesh.depth
     - |todo|
   * - dolfin.HierarchicalMesh.has_child
     - |todo|
   * - dolfin.HierarchicalMesh.has_parent
     - |todo|
   * - dolfin.HierarchicalMesh.set_child
     - |todo|
   * - dolfin.HierarchicalMesh.set_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.clear_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.depth
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.has_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.has_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.leaf_node
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.root_node
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.set_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionBool.set_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.clear_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.depth
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.has_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.has_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.leaf_node
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.root_node
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.set_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionDouble.set_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.clear_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.depth
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.has_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.has_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.leaf_node
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.root_node
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.set_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionInt.set_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.clear_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.depth
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.has_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.has_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.leaf_node
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.root_node
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.set_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionSizet.set_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionUInt
     - |todo|
   * - dolfin.HierarchicalMeshFunctionUInt.clear_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionUInt.depth
     - |todo|
   * - dolfin.HierarchicalMeshFunctionUInt.has_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionUInt.has_parent
     - |todo|
   * - dolfin.HierarchicalMeshFunctionUInt.set_child
     - |todo|
   * - dolfin.HierarchicalMeshFunctionUInt.set_parent
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.child
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.clear_child
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.depth
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.has_child
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.has_parent
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.leaf_node
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.parent
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.root_node
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.set_child
     - |todo|
   * - dolfin.HierarchicalNonlinearVariationalProblem.set_parent
     - |todo|
   * - dolfin.INFO
     - |todo|
   * - dolfin.ImplicitEuler
     - |todo|
   * - dolfin.ImplicitEuler.bcs
     - |todo|
   * - dolfin.ImplicitEuler.dolfin_stage_forms
     - |todo|
   * - dolfin.ImplicitEuler.dt
     - |todo|
   * - dolfin.ImplicitEuler.dt_stage_offset
     - |todo|
   * - dolfin.ImplicitEuler.id
     - |todo|
   * - dolfin.ImplicitEuler.implicit
     - |todo|
   * - dolfin.ImplicitEuler.jacobian_index
     - |todo|
   * - dolfin.ImplicitEuler.label
     - |todo|
   * - dolfin.ImplicitEuler.last_stage
     - |todo|
   * - dolfin.ImplicitEuler.name
     - |todo|
   * - dolfin.ImplicitEuler.order
     - |todo|
   * - dolfin.ImplicitEuler.parameters
     - |todo|
   * - dolfin.ImplicitEuler.rename
     - |todo|
   * - dolfin.ImplicitEuler.rhs_form
     - |todo|
   * - dolfin.ImplicitEuler.solution
     - |todo|
   * - dolfin.ImplicitEuler.stage_solutions
     - |todo|
   * - dolfin.ImplicitEuler.str
     - |todo|
   * - dolfin.ImplicitEuler.t
     - |todo|
   * - dolfin.ImplicitEuler.to_adm
     - |todo|
   * - dolfin.ImplicitEuler.to_tlm
     - |todo|
   * - dolfin.ImplicitEuler.ufl_stage_forms
     - |todo|
   * - dolfin.Index
     - |todo|
   * - dolfin.Index.count
     - |todo|
   * - dolfin.IndexMap.MapSize
     - |new|
   * - dolfin.IndexMap.MapSize_GLOBAL
     - |todo|
   * - dolfin.IndexMap.block_size
     - |todo|
   * - dolfin.IndexMap.global_index_owner
     - |todo|
   * - dolfin.IndexMap.init
     - |todo|
   * - dolfin.IndexMap.local_to_global
     - |todo|
   * - dolfin.IndexMap.local_to_global_unowned
     - |todo|
   * - dolfin.IndexMap.mpi_comm
     - |todo|
   * - dolfin.IndexMap.off_process_owner
     - |todo|
   * - dolfin.IndexMap.set_local_to_global
     - |todo|
   * - dolfin.IndexSet
     - |todo|
   * - dolfin.IndexSet.clear
     - |todo|
   * - dolfin.IndexSet.empty
     - |todo|
   * - dolfin.IndexSet.fill
     - |todo|
   * - dolfin.IndexSet.find
     - |todo|
   * - dolfin.IndexSet.has_index
     - |todo|
   * - dolfin.IndexSet.insert
     - |todo|
   * - dolfin.IndexSet.size
     - |todo|
   * - dolfin.IntArray
     - |todo|
   * - dolfin.IntArray.array
     - |todo|
   * - dolfin.IntArray.data
     - |todo|
   * - dolfin.IntArray.size
     - |todo|
   * - dolfin.IntArray.str
     - |todo|
   * - dolfin.Integral
     - |todo|
   * - dolfin.Integral.integral_type
     - |todo|
   * - dolfin.Integral.integrand
     - |todo|
   * - dolfin.Integral.metadata
     - |todo|
   * - dolfin.Integral.reconstruct
     - |todo|
   * - dolfin.Integral.subdomain_data
     - |todo|
   * - dolfin.Integral.subdomain_id
     - |todo|
   * - dolfin.Integral.ufl_domain
     - |todo|
   * - dolfin.InteriorElement
     - |todo|
   * - dolfin.IntervalMesh.cell_name
     - |new|
   * - dolfin.IntervalMesh.child
     - |todo|
   * - dolfin.IntervalMesh.clean
     - |todo|
   * - dolfin.IntervalMesh.clear_child
     - |todo|
   * - dolfin.IntervalMesh.create
     - |todo|
   * - dolfin.IntervalMesh.depth
     - |todo|
   * - dolfin.IntervalMesh.geometric_dimension
     - |new|
   * - dolfin.IntervalMesh.ghost_mode
     - |todo|
   * - dolfin.IntervalMesh.has_child
     - |todo|
   * - dolfin.IntervalMesh.has_parent
     - |todo|
   * - dolfin.IntervalMesh.label
     - |todo|
   * - dolfin.IntervalMesh.leaf_node
     - |todo|
   * - dolfin.IntervalMesh.name
     - |todo|
   * - dolfin.IntervalMesh.order
     - |todo|
   * - dolfin.IntervalMesh.parameters
     - |todo|
   * - dolfin.IntervalMesh.parent
     - |todo|
   * - dolfin.IntervalMesh.rename
     - |todo|
   * - dolfin.IntervalMesh.renumber_by_color
     - |todo|
   * - dolfin.IntervalMesh.root_node
     - |todo|
   * - dolfin.IntervalMesh.scale
     - |todo|
   * - dolfin.IntervalMesh.set_child
     - |todo|
   * - dolfin.IntervalMesh.set_parent
     - |todo|
   * - dolfin.IntervalMesh.str
     - |todo|
   * - dolfin.Jacobian
     - |todo|
   * - dolfin.Jacobian.T
     - |todo|
   * - dolfin.Jacobian.dx
     - |todo|
   * - dolfin.Jacobian.evaluate
     - |todo|
   * - dolfin.Jacobian.geometric_dimension
     - |todo|
   * - dolfin.Jacobian.is_cellwise_constant
     - |todo|
   * - dolfin.Jacobian.name
     - |todo|
   * - dolfin.Jacobian.ufl_disable_profiling
     - |todo|
   * - dolfin.Jacobian.ufl_domain
     - |todo|
   * - dolfin.Jacobian.ufl_domains
     - |todo|
   * - dolfin.Jacobian.ufl_enable_profiling
     - |todo|
   * - dolfin.Jacobian.ufl_free_indices
     - |todo|
   * - dolfin.Jacobian.ufl_index_dimensions
     - |todo|
   * - dolfin.Jacobian.ufl_operands
     - |todo|
   * - dolfin.Jacobian.ufl_shape
     - |todo|
   * - dolfin.JacobianDeterminant
     - |todo|
   * - dolfin.JacobianDeterminant.T
     - |todo|
   * - dolfin.JacobianDeterminant.dx
     - |todo|
   * - dolfin.JacobianDeterminant.evaluate
     - |todo|
   * - dolfin.JacobianDeterminant.geometric_dimension
     - |todo|
   * - dolfin.JacobianDeterminant.is_cellwise_constant
     - |todo|
   * - dolfin.JacobianDeterminant.name
     - |todo|
   * - dolfin.JacobianDeterminant.ufl_disable_profiling
     - |todo|
   * - dolfin.JacobianDeterminant.ufl_domain
     - |todo|
   * - dolfin.JacobianDeterminant.ufl_domains
     - |todo|
   * - dolfin.JacobianDeterminant.ufl_enable_profiling
     - |todo|
   * - dolfin.JacobianDeterminant.ufl_free_indices
     - |todo|
   * - dolfin.JacobianDeterminant.ufl_index_dimensions
     - |todo|
   * - dolfin.JacobianDeterminant.ufl_operands
     - |todo|
   * - dolfin.JacobianDeterminant.ufl_shape
     - |todo|
   * - dolfin.JacobianInverse
     - |todo|
   * - dolfin.JacobianInverse.T
     - |todo|
   * - dolfin.JacobianInverse.dx
     - |todo|
   * - dolfin.JacobianInverse.evaluate
     - |todo|
   * - dolfin.JacobianInverse.geometric_dimension
     - |todo|
   * - dolfin.JacobianInverse.is_cellwise_constant
     - |todo|
   * - dolfin.JacobianInverse.name
     - |todo|
   * - dolfin.JacobianInverse.ufl_disable_profiling
     - |todo|
   * - dolfin.JacobianInverse.ufl_domain
     - |todo|
   * - dolfin.JacobianInverse.ufl_domains
     - |todo|
   * - dolfin.JacobianInverse.ufl_enable_profiling
     - |todo|
   * - dolfin.JacobianInverse.ufl_free_indices
     - |todo|
   * - dolfin.JacobianInverse.ufl_index_dimensions
     - |todo|
   * - dolfin.JacobianInverse.ufl_operands
     - |todo|
   * - dolfin.JacobianInverse.ufl_shape
     - |todo|
   * - dolfin.KrylovSolver.default_parameters
     - |todo|
   * - dolfin.KrylovSolver.parameter_type
     - |todo|
   * - dolfin.KrylovSolver.str
     - |todo|
   * - dolfin.KrylovSolver.update_parameters
     - |todo|
   * - dolfin.L2
     - |todo|
   * - dolfin.LUSolver.default_parameters
     - |todo|
   * - dolfin.LUSolver.parameter_type
     - |todo|
   * - dolfin.LUSolver.set_operators
     - |todo|
   * - dolfin.LUSolver.str
     - |todo|
   * - dolfin.LUSolver.update_parameters
     - |todo|
   * - dolfin.Lagrange
     - |todo|
   * - dolfin.Lagrange.ddx
     - |todo|
   * - dolfin.Lagrange.degree
     - |todo|
   * - dolfin.Lagrange.dqdx
     - |todo|
   * - dolfin.Lagrange.eval
     - |todo|
   * - dolfin.Lagrange.id
     - |todo|
   * - dolfin.Lagrange.label
     - |todo|
   * - dolfin.Lagrange.name
     - |todo|
   * - dolfin.Lagrange.parameters
     - |todo|
   * - dolfin.Lagrange.point
     - |todo|
   * - dolfin.Lagrange.rename
     - |todo|
   * - dolfin.Lagrange.set
     - |todo|
   * - dolfin.Lagrange.size
     - |todo|
   * - dolfin.Lagrange.str
     - |todo|
   * - dolfin.Legendre
     - |todo|
   * - dolfin.Legendre.d2dx
     - |todo|
   * - dolfin.Legendre.ddx
     - |todo|
   * - dolfin.Legendre.eval
     - |todo|
   * - dolfin.LinearAlgebraObject
     - |todo|
   * - dolfin.LinearAlgebraObject.id
     - |todo|
   * - dolfin.LinearAlgebraObject.label
     - |todo|
   * - dolfin.LinearAlgebraObject.mpi_comm
     - |todo|
   * - dolfin.LinearAlgebraObject.name
     - |todo|
   * - dolfin.LinearAlgebraObject.parameters
     - |todo|
   * - dolfin.LinearAlgebraObject.rename
     - |todo|
   * - dolfin.LinearAlgebraObject.shared_instance
     - |todo|
   * - dolfin.LinearAlgebraObject.str
     - |todo|
   * - dolfin.LinearOperator.init_layout
     - |todo|
   * - dolfin.LinearOperator.instance
     - |new|
   * - dolfin.LinearOperator.mult
     - |todo|
   * - dolfin.LinearOperator.shared_instance
     - |todo|
   * - dolfin.LinearOperator.size
     - |todo|
   * - dolfin.LinearOperator.str
     - |todo|
   * - dolfin.LinearSolver
     - |todo|
   * - dolfin.LinearSolver.default_parameters
     - |todo|
   * - dolfin.LinearSolver.id
     - |todo|
   * - dolfin.LinearSolver.label
     - |todo|
   * - dolfin.LinearSolver.name
     - |todo|
   * - dolfin.LinearSolver.parameter_type
     - |todo|
   * - dolfin.LinearSolver.parameters
     - |todo|
   * - dolfin.LinearSolver.rename
     - |todo|
   * - dolfin.LinearSolver.set_operator
     - |todo|
   * - dolfin.LinearSolver.set_operators
     - |todo|
   * - dolfin.LinearSolver.solve
     - |todo|
   * - dolfin.LinearSolver.str
     - |todo|
   * - dolfin.LinearSolver.update_parameters
     - |todo|
   * - dolfin.LinearVariationalProblem.bilinear_form
     - |todo|
   * - dolfin.LinearVariationalProblem.child
     - |todo|
   * - dolfin.LinearVariationalProblem.clear_child
     - |todo|
   * - dolfin.LinearVariationalProblem.depth
     - |todo|
   * - dolfin.LinearVariationalProblem.has_child
     - |todo|
   * - dolfin.LinearVariationalProblem.has_parent
     - |todo|
   * - dolfin.LinearVariationalProblem.leaf_node
     - |todo|
   * - dolfin.LinearVariationalProblem.linear_form
     - |todo|
   * - dolfin.LinearVariationalProblem.parent
     - |todo|
   * - dolfin.LinearVariationalProblem.root_node
     - |todo|
   * - dolfin.LinearVariationalProblem.set_child
     - |todo|
   * - dolfin.LinearVariationalProblem.set_parent
     - |todo|
   * - dolfin.LinearVariationalProblem.solution
     - |todo|
   * - dolfin.LinearVariationalProblem.test_space
     - |todo|
   * - dolfin.LinearVariationalProblem.trial_space
     - |todo|
   * - dolfin.LinearVariationalSolver.default_parameters
     - |todo|
   * - dolfin.LinearVariationalSolver.str
     - |todo|
   * - dolfin.LocalAssembler
     - |todo|
   * - dolfin.LocalAssembler.assemble
     - |todo|
   * - dolfin.LocalAssembler.assemble_cell
     - |todo|
   * - dolfin.LocalAssembler.assemble_exterior_facet
     - |todo|
   * - dolfin.LocalAssembler.assemble_interior_facet
     - |todo|
   * - dolfin.LocalMeshData
     - |todo|
   * - dolfin.LocalMeshData.broadcast_mesh_data
     - |todo|
   * - dolfin.LocalMeshData.check
     - |todo|
   * - dolfin.LocalMeshData.clear
     - |todo|
   * - dolfin.LocalMeshData.domain_data
     - |todo|
   * - dolfin.LocalMeshData.extract_mesh_data
     - |todo|
   * - dolfin.LocalMeshData.geometry
     - |todo|
   * - dolfin.LocalMeshData.id
     - |todo|
   * - dolfin.LocalMeshData.label
     - |todo|
   * - dolfin.LocalMeshData.mpi_comm
     - |todo|
   * - dolfin.LocalMeshData.name
     - |todo|
   * - dolfin.LocalMeshData.parameters
     - |todo|
   * - dolfin.LocalMeshData.receive_mesh_data
     - |todo|
   * - dolfin.LocalMeshData.rename
     - |todo|
   * - dolfin.LocalMeshData.reorder
     - |todo|
   * - dolfin.LocalMeshData.str
     - |todo|
   * - dolfin.LocalMeshData.topology
     - |todo|
   * - dolfin.LocalSolver.SolverType
     - |new|
   * - dolfin.LogLevel
     - |new|
   * - dolfin.LogLevel.CRITICAL
     - |new|
   * - dolfin.LogLevel.DEBUG
     - |new|
   * - dolfin.LogLevel.ERROR
     - |new|
   * - dolfin.LogLevel.INFO
     - |new|
   * - dolfin.LogLevel.PROGRESS
     - |new|
   * - dolfin.LogLevel.TRACE
     - |new|
   * - dolfin.LogLevel.WARNING
     - |new|
   * - dolfin.MPI.MPI_AVG
     - |todo|
   * - dolfin.MPI.all_gather
     - |todo|
   * - dolfin.MPI.comm_null
     - |new|
   * - dolfin.MPI.comm_self
     - |new|
   * - dolfin.MPI.comm_world
     - |new|
   * - dolfin.MPI.compute_local_range
     - |todo|
   * - dolfin.MPI.finalized
     - |new|
   * - dolfin.MPI.gather
     - |todo|
   * - dolfin.MPI.global_offset
     - |todo|
   * - dolfin.MPI.index_owner
     - |todo|
   * - dolfin.MPI.init
     - |new|
   * - dolfin.MPI.initialized
     - |new|
   * - dolfin.MPI.is_broadcaster
     - |todo|
   * - dolfin.MPI.is_receiver
     - |todo|
   * - dolfin.MPI.responsible
     - |new|
   * - dolfin.MPICH_IGNORE_CXX_SEEK
     - |todo|
   * - dolfin.MPIInfo
     - |todo|
   * - dolfin.MPI_Comm
     - |todo|
   * - dolfin.Matrix.add
     - |todo|
   * - dolfin.Matrix.add_local
     - |todo|
   * - dolfin.Matrix.assign
     - |todo|
   * - dolfin.Matrix.ident_local
     - |todo|
   * - dolfin.Matrix.instance
     - |new|
   * - dolfin.Matrix.is_symmetric
     - |todo|
   * - dolfin.Matrix.set_local
     - |todo|
   * - dolfin.Matrix.setrow
     - |todo|
   * - dolfin.Matrix.shared_instance
     - |todo|
   * - dolfin.Matrix.zero_local
     - |todo|
   * - dolfin.Max
     - |todo|
   * - dolfin.Mesh.cell_name
     - |new|
   * - dolfin.Mesh.child
     - |todo|
   * - dolfin.Mesh.clean
     - |todo|
   * - dolfin.Mesh.clear_child
     - |todo|
   * - dolfin.Mesh.depth
     - |todo|
   * - dolfin.Mesh.geometric_dimension
     - |new|
   * - dolfin.Mesh.ghost_mode
     - |todo|
   * - dolfin.Mesh.has_child
     - |todo|
   * - dolfin.Mesh.has_parent
     - |todo|
   * - dolfin.Mesh.label
     - |todo|
   * - dolfin.Mesh.leaf_node
     - |todo|
   * - dolfin.Mesh.name
     - |todo|
   * - dolfin.Mesh.order
     - |todo|
   * - dolfin.Mesh.parameters
     - |todo|
   * - dolfin.Mesh.parent
     - |todo|
   * - dolfin.Mesh.rename
     - |todo|
   * - dolfin.Mesh.renumber_by_color
     - |todo|
   * - dolfin.Mesh.root_node
     - |todo|
   * - dolfin.Mesh.scale
     - |todo|
   * - dolfin.Mesh.set_child
     - |todo|
   * - dolfin.Mesh.set_parent
     - |todo|
   * - dolfin.Mesh.str
     - |todo|
   * - dolfin.MeshColoring.color
     - |todo|
   * - dolfin.MeshColoring.compute_colors
     - |todo|
   * - dolfin.MeshColoring.type_to_dim
     - |todo|
   * - dolfin.MeshConnectivity
     - |todo|
   * - dolfin.MeshConnectivity.clear
     - |todo|
   * - dolfin.MeshConnectivity.empty
     - |todo|
   * - dolfin.MeshConnectivity.hash
     - |todo|
   * - dolfin.MeshConnectivity.init
     - |todo|
   * - dolfin.MeshConnectivity.set_global_size
     - |todo|
   * - dolfin.MeshConnectivity.size
     - |todo|
   * - dolfin.MeshConnectivity.size_global
     - |todo|
   * - dolfin.MeshConnectivity.str
     - |todo|
   * - dolfin.MeshCoordinates.cpp_object
     - |new|
   * - dolfin.MeshCoordinates.eval
     - |todo|
   * - dolfin.MeshCoordinates.eval_cell
     - |todo|
   * - dolfin.MeshCoordinates.get_generic_function
     - |todo|
   * - dolfin.MeshCoordinates.get_property
     - |todo|
   * - dolfin.MeshCoordinates.parameters
     - |todo|
   * - dolfin.MeshCoordinates.rename
     - |todo|
   * - dolfin.MeshCoordinates.restrict
     - |todo|
   * - dolfin.MeshCoordinates.set_generic_function
     - |todo|
   * - dolfin.MeshCoordinates.set_property
     - |todo|
   * - dolfin.MeshCoordinates.str
     - |todo|
   * - dolfin.MeshCoordinates.update
     - |todo|
   * - dolfin.MeshCoordinates.value_shape
     - |todo|
   * - dolfin.MeshCoordinates.value_size
     - |todo|
   * - dolfin.MeshData
     - |todo|
   * - dolfin.MeshData.array
     - |todo|
   * - dolfin.MeshData.clear
     - |todo|
   * - dolfin.MeshData.create_array
     - |todo|
   * - dolfin.MeshData.erase_array
     - |todo|
   * - dolfin.MeshData.exists
     - |todo|
   * - dolfin.MeshData.id
     - |todo|
   * - dolfin.MeshData.label
     - |todo|
   * - dolfin.MeshData.name
     - |todo|
   * - dolfin.MeshData.parameters
     - |todo|
   * - dolfin.MeshData.rename
     - |todo|
   * - dolfin.MeshData.str
     - |todo|
   * - dolfin.MeshDisplacement
     - |todo|
   * - dolfin.MeshDisplacement.compute_vertex_values
     - |todo|
   * - dolfin.MeshDisplacement.eval
     - |todo|
   * - dolfin.MeshDisplacement.eval_cell
     - |todo|
   * - dolfin.MeshDisplacement.get_generic_function
     - |todo|
   * - dolfin.MeshDisplacement.get_property
     - |todo|
   * - dolfin.MeshDisplacement.id
     - |todo|
   * - dolfin.MeshDisplacement.label
     - |todo|
   * - dolfin.MeshDisplacement.name
     - |todo|
   * - dolfin.MeshDisplacement.parameters
     - |todo|
   * - dolfin.MeshDisplacement.rename
     - |todo|
   * - dolfin.MeshDisplacement.restrict
     - |todo|
   * - dolfin.MeshDisplacement.set_generic_function
     - |todo|
   * - dolfin.MeshDisplacement.set_property
     - |todo|
   * - dolfin.MeshDisplacement.str
     - |todo|
   * - dolfin.MeshDisplacement.sub
     - |todo|
   * - dolfin.MeshDisplacement.update
     - |todo|
   * - dolfin.MeshDisplacement.value_dimension
     - |todo|
   * - dolfin.MeshDisplacement.value_rank
     - |todo|
   * - dolfin.MeshDisplacement.value_shape
     - |todo|
   * - dolfin.MeshDisplacement.value_size
     - |todo|
   * - dolfin.MeshDomains
     - |todo|
   * - dolfin.MeshDomains.clear
     - |todo|
   * - dolfin.MeshDomains.get_marker
     - |todo|
   * - dolfin.MeshDomains.init
     - |todo|
   * - dolfin.MeshDomains.is_empty
     - |todo|
   * - dolfin.MeshDomains.markers
     - |todo|
   * - dolfin.MeshDomains.max_dim
     - |todo|
   * - dolfin.MeshDomains.num_marked
     - |todo|
   * - dolfin.MeshDomains.set_marker
     - |todo|
   * - dolfin.MeshEditor.add_entity_point
     - |todo|
   * - dolfin.MeshEditor.add_vertex_global
     - |todo|
   * - dolfin.MeshEditor.init_entities
     - |todo|
   * - dolfin.MeshEntity.incident
     - |todo|
   * - dolfin.MeshEntity.init
     - |todo|
   * - dolfin.MeshEntity.mesh_id
     - |todo|
   * - dolfin.MeshEntity.owner
     - |todo|
   * - dolfin.MeshEntity.str
     - |todo|
   * - dolfin.MeshFunctionBool
     - |todo|
   * - dolfin.MeshFunctionBool.array
     - |todo|
   * - dolfin.MeshFunctionBool.child
     - |todo|
   * - dolfin.MeshFunctionBool.clear_child
     - |todo|
   * - dolfin.MeshFunctionBool.cpp_value_type
     - |todo|
   * - dolfin.MeshFunctionBool.depth
     - |todo|
   * - dolfin.MeshFunctionBool.dim
     - |todo|
   * - dolfin.MeshFunctionBool.empty
     - |todo|
   * - dolfin.MeshFunctionBool.has_child
     - |todo|
   * - dolfin.MeshFunctionBool.has_parent
     - |todo|
   * - dolfin.MeshFunctionBool.id
     - |todo|
   * - dolfin.MeshFunctionBool.init
     - |todo|
   * - dolfin.MeshFunctionBool.label
     - |todo|
   * - dolfin.MeshFunctionBool.leaf_node
     - |todo|
   * - dolfin.MeshFunctionBool.mesh
     - |todo|
   * - dolfin.MeshFunctionBool.name
     - |todo|
   * - dolfin.MeshFunctionBool.parameters
     - |todo|
   * - dolfin.MeshFunctionBool.parent
     - |todo|
   * - dolfin.MeshFunctionBool.rename
     - |todo|
   * - dolfin.MeshFunctionBool.root_node
     - |todo|
   * - dolfin.MeshFunctionBool.set_all
     - |todo|
   * - dolfin.MeshFunctionBool.set_child
     - |todo|
   * - dolfin.MeshFunctionBool.set_parent
     - |todo|
   * - dolfin.MeshFunctionBool.set_value
     - |todo|
   * - dolfin.MeshFunctionBool.set_values
     - |todo|
   * - dolfin.MeshFunctionBool.size
     - |todo|
   * - dolfin.MeshFunctionBool.str
     - |todo|
   * - dolfin.MeshFunctionBool.ufl_id
     - |todo|
   * - dolfin.MeshFunctionBool.where_equal
     - |todo|
   * - dolfin.MeshFunctionDouble
     - |todo|
   * - dolfin.MeshFunctionDouble.array
     - |todo|
   * - dolfin.MeshFunctionDouble.child
     - |todo|
   * - dolfin.MeshFunctionDouble.clear_child
     - |todo|
   * - dolfin.MeshFunctionDouble.cpp_value_type
     - |todo|
   * - dolfin.MeshFunctionDouble.depth
     - |todo|
   * - dolfin.MeshFunctionDouble.dim
     - |todo|
   * - dolfin.MeshFunctionDouble.empty
     - |todo|
   * - dolfin.MeshFunctionDouble.has_child
     - |todo|
   * - dolfin.MeshFunctionDouble.has_parent
     - |todo|
   * - dolfin.MeshFunctionDouble.id
     - |todo|
   * - dolfin.MeshFunctionDouble.init
     - |todo|
   * - dolfin.MeshFunctionDouble.label
     - |todo|
   * - dolfin.MeshFunctionDouble.leaf_node
     - |todo|
   * - dolfin.MeshFunctionDouble.mesh
     - |todo|
   * - dolfin.MeshFunctionDouble.name
     - |todo|
   * - dolfin.MeshFunctionDouble.parameters
     - |todo|
   * - dolfin.MeshFunctionDouble.parent
     - |todo|
   * - dolfin.MeshFunctionDouble.rename
     - |todo|
   * - dolfin.MeshFunctionDouble.root_node
     - |todo|
   * - dolfin.MeshFunctionDouble.set_all
     - |todo|
   * - dolfin.MeshFunctionDouble.set_child
     - |todo|
   * - dolfin.MeshFunctionDouble.set_parent
     - |todo|
   * - dolfin.MeshFunctionDouble.set_value
     - |todo|
   * - dolfin.MeshFunctionDouble.set_values
     - |todo|
   * - dolfin.MeshFunctionDouble.size
     - |todo|
   * - dolfin.MeshFunctionDouble.str
     - |todo|
   * - dolfin.MeshFunctionDouble.ufl_id
     - |todo|
   * - dolfin.MeshFunctionDouble.where_equal
     - |todo|
   * - dolfin.MeshFunctionInt
     - |todo|
   * - dolfin.MeshFunctionInt.array
     - |todo|
   * - dolfin.MeshFunctionInt.child
     - |todo|
   * - dolfin.MeshFunctionInt.clear_child
     - |todo|
   * - dolfin.MeshFunctionInt.cpp_value_type
     - |todo|
   * - dolfin.MeshFunctionInt.depth
     - |todo|
   * - dolfin.MeshFunctionInt.dim
     - |todo|
   * - dolfin.MeshFunctionInt.empty
     - |todo|
   * - dolfin.MeshFunctionInt.has_child
     - |todo|
   * - dolfin.MeshFunctionInt.has_parent
     - |todo|
   * - dolfin.MeshFunctionInt.id
     - |todo|
   * - dolfin.MeshFunctionInt.init
     - |todo|
   * - dolfin.MeshFunctionInt.label
     - |todo|
   * - dolfin.MeshFunctionInt.leaf_node
     - |todo|
   * - dolfin.MeshFunctionInt.mesh
     - |todo|
   * - dolfin.MeshFunctionInt.name
     - |todo|
   * - dolfin.MeshFunctionInt.parameters
     - |todo|
   * - dolfin.MeshFunctionInt.parent
     - |todo|
   * - dolfin.MeshFunctionInt.rename
     - |todo|
   * - dolfin.MeshFunctionInt.root_node
     - |todo|
   * - dolfin.MeshFunctionInt.set_all
     - |todo|
   * - dolfin.MeshFunctionInt.set_child
     - |todo|
   * - dolfin.MeshFunctionInt.set_parent
     - |todo|
   * - dolfin.MeshFunctionInt.set_value
     - |todo|
   * - dolfin.MeshFunctionInt.set_values
     - |todo|
   * - dolfin.MeshFunctionInt.size
     - |todo|
   * - dolfin.MeshFunctionInt.str
     - |todo|
   * - dolfin.MeshFunctionInt.ufl_id
     - |todo|
   * - dolfin.MeshFunctionInt.where_equal
     - |todo|
   * - dolfin.MeshFunctionSizet
     - |todo|
   * - dolfin.MeshFunctionSizet.array
     - |todo|
   * - dolfin.MeshFunctionSizet.child
     - |todo|
   * - dolfin.MeshFunctionSizet.clear_child
     - |todo|
   * - dolfin.MeshFunctionSizet.cpp_value_type
     - |todo|
   * - dolfin.MeshFunctionSizet.depth
     - |todo|
   * - dolfin.MeshFunctionSizet.dim
     - |todo|
   * - dolfin.MeshFunctionSizet.empty
     - |todo|
   * - dolfin.MeshFunctionSizet.has_child
     - |todo|
   * - dolfin.MeshFunctionSizet.has_parent
     - |todo|
   * - dolfin.MeshFunctionSizet.id
     - |todo|
   * - dolfin.MeshFunctionSizet.init
     - |todo|
   * - dolfin.MeshFunctionSizet.label
     - |todo|
   * - dolfin.MeshFunctionSizet.leaf_node
     - |todo|
   * - dolfin.MeshFunctionSizet.mesh
     - |todo|
   * - dolfin.MeshFunctionSizet.name
     - |todo|
   * - dolfin.MeshFunctionSizet.parameters
     - |todo|
   * - dolfin.MeshFunctionSizet.parent
     - |todo|
   * - dolfin.MeshFunctionSizet.rename
     - |todo|
   * - dolfin.MeshFunctionSizet.root_node
     - |todo|
   * - dolfin.MeshFunctionSizet.set_all
     - |todo|
   * - dolfin.MeshFunctionSizet.set_child
     - |todo|
   * - dolfin.MeshFunctionSizet.set_parent
     - |todo|
   * - dolfin.MeshFunctionSizet.set_value
     - |todo|
   * - dolfin.MeshFunctionSizet.set_values
     - |todo|
   * - dolfin.MeshFunctionSizet.size
     - |todo|
   * - dolfin.MeshFunctionSizet.str
     - |todo|
   * - dolfin.MeshFunctionSizet.ufl_id
     - |todo|
   * - dolfin.MeshFunctionSizet.where_equal
     - |todo|
   * - dolfin.MeshGeometry.get_entity_index
     - |todo|
   * - dolfin.MeshGeometry.hash
     - |todo|
   * - dolfin.MeshGeometry.init
     - |todo|
   * - dolfin.MeshGeometry.init_entities
     - |todo|
   * - dolfin.MeshGeometry.num_entity_coordinates
     - |todo|
   * - dolfin.MeshGeometry.num_points
     - |todo|
   * - dolfin.MeshGeometry.num_vertices
     - |todo|
   * - dolfin.MeshGeometry.point
     - |todo|
   * - dolfin.MeshGeometry.point_coordinates
     - |todo|
   * - dolfin.MeshGeometry.set
     - |todo|
   * - dolfin.MeshGeometry.str
     - |todo|
   * - dolfin.MeshGeometry.vertex_coordinates
     - |todo|
   * - dolfin.MeshGeometry.x
     - |todo|
   * - dolfin.MeshHierarchy
     - |todo|
   * - dolfin.MeshHierarchy.coarsen
     - |todo|
   * - dolfin.MeshHierarchy.coarsest
     - |todo|
   * - dolfin.MeshHierarchy.finest
     - |todo|
   * - dolfin.MeshHierarchy.rebalance
     - |todo|
   * - dolfin.MeshHierarchy.refine
     - |todo|
   * - dolfin.MeshHierarchy.size
     - |todo|
   * - dolfin.MeshHierarchy.unrefine
     - |todo|
   * - dolfin.MeshHierarchy.weight
     - |todo|
   * - dolfin.MeshPartitioning
     - |todo|
   * - dolfin.MeshPartitioning.build_distributed_mesh
     - |todo|
   * - dolfin.MeshQuality.dihedral_angles
     - |todo|
   * - dolfin.MeshQuality.dihedral_angles_histogram_data
     - |todo|
   * - dolfin.MeshRenumbering
     - |todo|
   * - dolfin.MeshRenumbering.renumber_by_color
     - |todo|
   * - dolfin.MeshTopology.clear
     - |todo|
   * - dolfin.MeshTopology.coloring
     - |todo|
   * - dolfin.MeshTopology.init
     - |todo|
   * - dolfin.MeshTopology.init_ghost
     - |todo|
   * - dolfin.MeshTopology.init_global_indices
     - |todo|
   * - dolfin.MeshTopology.set_global_index
     - |todo|
   * - dolfin.MeshTopology.size_global
     - |todo|
   * - dolfin.MeshTopology.str
     - |todo|
   * - dolfin.MeshTransformation.scale
     - |todo|
   * - dolfin.MeshValueCollectionBool
     - |todo|
   * - dolfin.MeshValueCollectionBool.assign
     - |todo|
   * - dolfin.MeshValueCollectionBool.clear
     - |todo|
   * - dolfin.MeshValueCollectionBool.dim
     - |todo|
   * - dolfin.MeshValueCollectionBool.empty
     - |todo|
   * - dolfin.MeshValueCollectionBool.get_value
     - |todo|
   * - dolfin.MeshValueCollectionBool.id
     - |todo|
   * - dolfin.MeshValueCollectionBool.init
     - |todo|
   * - dolfin.MeshValueCollectionBool.label
     - |todo|
   * - dolfin.MeshValueCollectionBool.mesh
     - |todo|
   * - dolfin.MeshValueCollectionBool.name
     - |todo|
   * - dolfin.MeshValueCollectionBool.parameters
     - |todo|
   * - dolfin.MeshValueCollectionBool.rename
     - |todo|
   * - dolfin.MeshValueCollectionBool.set_value
     - |todo|
   * - dolfin.MeshValueCollectionBool.size
     - |todo|
   * - dolfin.MeshValueCollectionBool.str
     - |todo|
   * - dolfin.MeshValueCollectionBool.values
     - |todo|
   * - dolfin.MeshValueCollectionDouble
     - |todo|
   * - dolfin.MeshValueCollectionDouble.assign
     - |todo|
   * - dolfin.MeshValueCollectionDouble.clear
     - |todo|
   * - dolfin.MeshValueCollectionDouble.dim
     - |todo|
   * - dolfin.MeshValueCollectionDouble.empty
     - |todo|
   * - dolfin.MeshValueCollectionDouble.get_value
     - |todo|
   * - dolfin.MeshValueCollectionDouble.id
     - |todo|
   * - dolfin.MeshValueCollectionDouble.init
     - |todo|
   * - dolfin.MeshValueCollectionDouble.label
     - |todo|
   * - dolfin.MeshValueCollectionDouble.mesh
     - |todo|
   * - dolfin.MeshValueCollectionDouble.name
     - |todo|
   * - dolfin.MeshValueCollectionDouble.parameters
     - |todo|
   * - dolfin.MeshValueCollectionDouble.rename
     - |todo|
   * - dolfin.MeshValueCollectionDouble.set_value
     - |todo|
   * - dolfin.MeshValueCollectionDouble.size
     - |todo|
   * - dolfin.MeshValueCollectionDouble.str
     - |todo|
   * - dolfin.MeshValueCollectionDouble.values
     - |todo|
   * - dolfin.MeshValueCollectionInt
     - |todo|
   * - dolfin.MeshValueCollectionInt.assign
     - |todo|
   * - dolfin.MeshValueCollectionInt.clear
     - |todo|
   * - dolfin.MeshValueCollectionInt.dim
     - |todo|
   * - dolfin.MeshValueCollectionInt.empty
     - |todo|
   * - dolfin.MeshValueCollectionInt.get_value
     - |todo|
   * - dolfin.MeshValueCollectionInt.id
     - |todo|
   * - dolfin.MeshValueCollectionInt.init
     - |todo|
   * - dolfin.MeshValueCollectionInt.label
     - |todo|
   * - dolfin.MeshValueCollectionInt.mesh
     - |todo|
   * - dolfin.MeshValueCollectionInt.name
     - |todo|
   * - dolfin.MeshValueCollectionInt.parameters
     - |todo|
   * - dolfin.MeshValueCollectionInt.rename
     - |todo|
   * - dolfin.MeshValueCollectionInt.set_value
     - |todo|
   * - dolfin.MeshValueCollectionInt.size
     - |todo|
   * - dolfin.MeshValueCollectionInt.str
     - |todo|
   * - dolfin.MeshValueCollectionInt.values
     - |todo|
   * - dolfin.MeshValueCollectionSizet
     - |todo|
   * - dolfin.MeshValueCollectionSizet.assign
     - |todo|
   * - dolfin.MeshValueCollectionSizet.clear
     - |todo|
   * - dolfin.MeshValueCollectionSizet.dim
     - |todo|
   * - dolfin.MeshValueCollectionSizet.empty
     - |todo|
   * - dolfin.MeshValueCollectionSizet.get_value
     - |todo|
   * - dolfin.MeshValueCollectionSizet.id
     - |todo|
   * - dolfin.MeshValueCollectionSizet.init
     - |todo|
   * - dolfin.MeshValueCollectionSizet.label
     - |todo|
   * - dolfin.MeshValueCollectionSizet.mesh
     - |todo|
   * - dolfin.MeshValueCollectionSizet.name
     - |todo|
   * - dolfin.MeshValueCollectionSizet.parameters
     - |todo|
   * - dolfin.MeshValueCollectionSizet.rename
     - |todo|
   * - dolfin.MeshValueCollectionSizet.set_value
     - |todo|
   * - dolfin.MeshValueCollectionSizet.size
     - |todo|
   * - dolfin.MeshValueCollectionSizet.str
     - |todo|
   * - dolfin.MeshValueCollectionSizet.values
     - |todo|
   * - dolfin.MeshView
     - |todo|
   * - dolfin.MeshView.geometric_dimension
     - |todo|
   * - dolfin.MeshView.is_piecewise_linear_simplex_domain
     - |todo|
   * - dolfin.MeshView.topological_dimension
     - |todo|
   * - dolfin.MeshView.ufl_cell
     - |todo|
   * - dolfin.MeshView.ufl_id
     - |todo|
   * - dolfin.MeshView.ufl_mesh
     - |todo|
   * - dolfin.Min
     - |todo|
   * - dolfin.MultiMesh.bounding_box_tree
     - |todo|
   * - dolfin.MultiMesh.bounding_box_tree_boundary
     - |todo|
   * - dolfin.MultiMesh.clear
     - |todo|
   * - dolfin.MultiMesh.collision_map_cut_cells
     - |todo|
   * - dolfin.MultiMesh.compute_area
     - |todo|
   * - dolfin.MultiMesh.default_parameters
     - |todo|
   * - dolfin.MultiMesh.facet_normals
     - |todo|
   * - dolfin.MultiMesh.is_built
     - |todo|
   * - dolfin.MultiMesh.mpi_comm
     - |todo|
   * - dolfin.MultiMesh.plot_matplotlib
     - |todo|
   * - dolfin.MultiMesh.quadrature_rules_interface
     - |todo|
   * - dolfin.MultiMesh.quadrature_rules_overlap
     - |todo|
   * - dolfin.MultiMesh.str
     - |todo|
   * - dolfin.MultiMesh.type
     - |todo|
   * - dolfin.MultiMesh.ufl_cell
     - |todo|
   * - dolfin.MultiMesh.ufl_coordinate_element
     - |todo|
   * - dolfin.MultiMesh.ufl_domain
     - |todo|
   * - dolfin.MultiMesh.ufl_id
     - |todo|
   * - dolfin.MultiMeshAssembler
     - |todo|
   * - dolfin.MultiMeshAssembler.add_values
     - |todo|
   * - dolfin.MultiMeshAssembler.assemble
     - |todo|
   * - dolfin.MultiMeshAssembler.extend_cut_cell_integration
     - |todo|
   * - dolfin.MultiMeshAssembler.finalize_tensor
     - |todo|
   * - dolfin.MultiMeshAssembler.init_global_tensor
     - |todo|
   * - dolfin.MultiMeshAssembler.keep_diagonal
     - |todo|
   * - dolfin.MultiMeshDirichletBC
     - |todo|
   * - dolfin.MultiMeshDirichletBC.apply
     - |todo|
   * - dolfin.MultiMeshDirichletBC.homogenize
     - |todo|
   * - dolfin.MultiMeshDirichletBC.zero
     - |todo|
   * - dolfin.MultiMeshDofMap
     - |todo|
   * - dolfin.MultiMeshDofMap.add
     - |todo|
   * - dolfin.MultiMeshDofMap.build
     - |todo|
   * - dolfin.MultiMeshDofMap.clear
     - |todo|
   * - dolfin.MultiMeshDofMap.global_dimension
     - |todo|
   * - dolfin.MultiMeshDofMap.inactive_dofs
     - |todo|
   * - dolfin.MultiMeshDofMap.index_map
     - |todo|
   * - dolfin.MultiMeshDofMap.num_parts
     - |todo|
   * - dolfin.MultiMeshDofMap.off_process_owner
     - |todo|
   * - dolfin.MultiMeshDofMap.ownership_range
     - |todo|
   * - dolfin.MultiMeshDofMap.part
     - |todo|
   * - dolfin.MultiMeshDofMap.str
     - |todo|
   * - dolfin.MultiMeshForm
     - |todo|
   * - dolfin.MultiMeshForm.add
     - |todo|
   * - dolfin.MultiMeshForm.build
     - |todo|
   * - dolfin.MultiMeshForm.clear
     - |todo|
   * - dolfin.MultiMeshForm.function_space
     - |todo|
   * - dolfin.MultiMeshForm.multimesh
     - |todo|
   * - dolfin.MultiMeshForm.multimesh_coefficient
     - |todo|
   * - dolfin.MultiMeshForm.multimesh_coefficient_keys
     - |todo|
   * - dolfin.MultiMeshForm.multimesh_coefficients
     - |todo|
   * - dolfin.MultiMeshForm.num_parts
     - |todo|
   * - dolfin.MultiMeshForm.part
     - |todo|
   * - dolfin.MultiMeshForm.rank
     - |todo|
   * - dolfin.MultiMeshForm.set_multimesh_coefficient
     - |todo|
   * - dolfin.MultiMeshFunction.T
     - |todo|
   * - dolfin.MultiMeshFunction.assign
     - |todo|
   * - dolfin.MultiMeshFunction.assign_part
     - |todo|
   * - dolfin.MultiMeshFunction.copy
     - |todo|
   * - dolfin.MultiMeshFunction.count
     - |todo|
   * - dolfin.MultiMeshFunction.dx
     - |todo|
   * - dolfin.MultiMeshFunction.eval
     - |todo|
   * - dolfin.MultiMeshFunction.evaluate
     - |todo|
   * - dolfin.MultiMeshFunction.function_space
     - |todo|
   * - dolfin.MultiMeshFunction.geometric_dimension
     - |todo|
   * - dolfin.MultiMeshFunction.id
     - |todo|
   * - dolfin.MultiMeshFunction.interpolate
     - |todo|
   * - dolfin.MultiMeshFunction.is_cellwise_constant
     - |todo|
   * - dolfin.MultiMeshFunction.label
     - |todo|
   * - dolfin.MultiMeshFunction.name
     - |todo|
   * - dolfin.MultiMeshFunction.parameters
     - |todo|
   * - dolfin.MultiMeshFunction.part
     - |todo|
   * - dolfin.MultiMeshFunction.parts
     - |todo|
   * - dolfin.MultiMeshFunction.rename
     - |todo|
   * - dolfin.MultiMeshFunction.restrict
     - |todo|
   * - dolfin.MultiMeshFunction.restrict_as_ufc_function
     - |todo|
   * - dolfin.MultiMeshFunction.str
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_disable_profiling
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_domain
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_domains
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_element
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_enable_profiling
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_free_indices
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_function_space
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_index_dimensions
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_operands
     - |todo|
   * - dolfin.MultiMeshFunction.ufl_shape
     - |todo|
   * - dolfin.MultiMeshFunctionSpace.add
     - |new|
   * - dolfin.MultiMeshFunctionSpace.build
     - |new|
   * - dolfin.MultiMeshSubSpace
     - |todo|
   * - dolfin.MultiMeshSubSpace.add
     - |todo|
   * - dolfin.MultiMeshSubSpace.build
     - |todo|
   * - dolfin.MultiMeshSubSpace.dim
     - |todo|
   * - dolfin.MultiMeshSubSpace.dofmap
     - |todo|
   * - dolfin.MultiMeshSubSpace.id
     - |todo|
   * - dolfin.MultiMeshSubSpace.label
     - |todo|
   * - dolfin.MultiMeshSubSpace.lock_inactive_dofs
     - |todo|
   * - dolfin.MultiMeshSubSpace.multimesh
     - |todo|
   * - dolfin.MultiMeshSubSpace.name
     - |todo|
   * - dolfin.MultiMeshSubSpace.num_parts
     - |todo|
   * - dolfin.MultiMeshSubSpace.parameters
     - |todo|
   * - dolfin.MultiMeshSubSpace.part
     - |todo|
   * - dolfin.MultiMeshSubSpace.rename
     - |todo|
   * - dolfin.MultiMeshSubSpace.str
     - |todo|
   * - dolfin.MultiMeshSubSpace.view
     - |todo|
   * - dolfin.MultiStageScheme
     - |todo|
   * - dolfin.MultiStageScheme.bcs
     - |todo|
   * - dolfin.MultiStageScheme.dolfin_stage_forms
     - |todo|
   * - dolfin.MultiStageScheme.dt
     - |todo|
   * - dolfin.MultiStageScheme.dt_stage_offset
     - |todo|
   * - dolfin.MultiStageScheme.id
     - |todo|
   * - dolfin.MultiStageScheme.implicit
     - |todo|
   * - dolfin.MultiStageScheme.jacobian_index
     - |todo|
   * - dolfin.MultiStageScheme.label
     - |todo|
   * - dolfin.MultiStageScheme.last_stage
     - |todo|
   * - dolfin.MultiStageScheme.name
     - |todo|
   * - dolfin.MultiStageScheme.order
     - |todo|
   * - dolfin.MultiStageScheme.parameters
     - |todo|
   * - dolfin.MultiStageScheme.rename
     - |todo|
   * - dolfin.MultiStageScheme.rhs_form
     - |todo|
   * - dolfin.MultiStageScheme.solution
     - |todo|
   * - dolfin.MultiStageScheme.stage_solutions
     - |todo|
   * - dolfin.MultiStageScheme.str
     - |todo|
   * - dolfin.MultiStageScheme.t
     - |todo|
   * - dolfin.MultiStageScheme.to_adm
     - |todo|
   * - dolfin.MultiStageScheme.to_tlm
     - |todo|
   * - dolfin.MultiStageScheme.ufl_stage_forms
     - |todo|
   * - dolfin.NewtonSolver.default_parameters
     - |todo|
   * - dolfin.NewtonSolver.get_relaxation_parameter
     - |todo|
   * - dolfin.NewtonSolver.iteration
     - |todo|
   * - dolfin.NewtonSolver.krylov_iterations
     - |todo|
   * - dolfin.NewtonSolver.linear_solver
     - |todo|
   * - dolfin.NewtonSolver.relative_residual
     - |todo|
   * - dolfin.NewtonSolver.residual
     - |todo|
   * - dolfin.NewtonSolver.residual0
     - |todo|
   * - dolfin.NewtonSolver.set_relaxation_parameter
     - |todo|
   * - dolfin.NewtonSolver.str
     - |todo|
   * - dolfin.NodalEnrichedElement
     - |todo|
   * - dolfin.NodalEnrichedElement.cell
     - |todo|
   * - dolfin.NodalEnrichedElement.degree
     - |todo|
   * - dolfin.NodalEnrichedElement.extract_component
     - |todo|
   * - dolfin.NodalEnrichedElement.extract_reference_component
     - |todo|
   * - dolfin.NodalEnrichedElement.extract_subelement_component
     - |todo|
   * - dolfin.NodalEnrichedElement.extract_subelement_reference_component
     - |todo|
   * - dolfin.NodalEnrichedElement.family
     - |todo|
   * - dolfin.NodalEnrichedElement.is_cellwise_constant
     - |todo|
   * - dolfin.NodalEnrichedElement.mapping
     - |todo|
   * - dolfin.NodalEnrichedElement.num_sub_elements
     - |todo|
   * - dolfin.NodalEnrichedElement.quadrature_scheme
     - |todo|
   * - dolfin.NodalEnrichedElement.reconstruct
     - |todo|
   * - dolfin.NodalEnrichedElement.reference_value_shape
     - |todo|
   * - dolfin.NodalEnrichedElement.reference_value_size
     - |todo|
   * - dolfin.NodalEnrichedElement.shortstr
     - |todo|
   * - dolfin.NodalEnrichedElement.sobolev_space
     - |todo|
   * - dolfin.NodalEnrichedElement.sub_elements
     - |todo|
   * - dolfin.NodalEnrichedElement.symmetry
     - |todo|
   * - dolfin.NodalEnrichedElement.value_shape
     - |todo|
   * - dolfin.NodalEnrichedElement.value_size
     - |todo|
   * - dolfin.NonlinearProblem.J_pc
     - |todo|
   * - dolfin.NonlinearVariationalProblem.bcs
     - |todo|
   * - dolfin.NonlinearVariationalProblem.child
     - |todo|
   * - dolfin.NonlinearVariationalProblem.clear_child
     - |todo|
   * - dolfin.NonlinearVariationalProblem.depth
     - |todo|
   * - dolfin.NonlinearVariationalProblem.has_child
     - |todo|
   * - dolfin.NonlinearVariationalProblem.has_jacobian
     - |todo|
   * - dolfin.NonlinearVariationalProblem.has_lower_bound
     - |todo|
   * - dolfin.NonlinearVariationalProblem.has_parent
     - |todo|
   * - dolfin.NonlinearVariationalProblem.has_upper_bound
     - |todo|
   * - dolfin.NonlinearVariationalProblem.jacobian_form
     - |todo|
   * - dolfin.NonlinearVariationalProblem.leaf_node
     - |todo|
   * - dolfin.NonlinearVariationalProblem.lower_bound
     - |todo|
   * - dolfin.NonlinearVariationalProblem.parent
     - |todo|
   * - dolfin.NonlinearVariationalProblem.residual_form
     - |todo|
   * - dolfin.NonlinearVariationalProblem.root_node
     - |todo|
   * - dolfin.NonlinearVariationalProblem.set_child
     - |todo|
   * - dolfin.NonlinearVariationalProblem.set_parent
     - |todo|
   * - dolfin.NonlinearVariationalProblem.solution
     - |todo|
   * - dolfin.NonlinearVariationalProblem.test_space
     - |todo|
   * - dolfin.NonlinearVariationalProblem.trial_space
     - |todo|
   * - dolfin.NonlinearVariationalProblem.upper_bound
     - |todo|
   * - dolfin.NonlinearVariationalSolver.default_parameters
     - |todo|
   * - dolfin.NonlinearVariationalSolver.str
     - |todo|
   * - dolfin.Not
     - |todo|
   * - dolfin.OptimisationProblem.J_pc
     - |todo|
   * - dolfin.OptimisationProblem.form
     - |todo|
   * - dolfin.Or
     - |todo|
   * - dolfin.PETScBaseMatrix
     - |todo|
   * - dolfin.PETScBaseMatrix.id
     - |todo|
   * - dolfin.PETScBaseMatrix.init_vector
     - |todo|
   * - dolfin.PETScBaseMatrix.label
     - |todo|
   * - dolfin.PETScBaseMatrix.local_range
     - |todo|
   * - dolfin.PETScBaseMatrix.mat
     - |todo|
   * - dolfin.PETScBaseMatrix.mpi_comm
     - |todo|
   * - dolfin.PETScBaseMatrix.name
     - |todo|
   * - dolfin.PETScBaseMatrix.parameters
     - |todo|
   * - dolfin.PETScBaseMatrix.petsc_error
     - |todo|
   * - dolfin.PETScBaseMatrix.rename
     - |todo|
   * - dolfin.PETScBaseMatrix.size
     - |todo|
   * - dolfin.PETScBaseMatrix.str
     - |todo|
   * - dolfin.PETScDMCollection.petsc_error
     - |todo|
   * - dolfin.PETScDMCollection.reset
     - |todo|
   * - dolfin.PETScFactory.create_krylov_solver
     - |todo|
   * - dolfin.PETScFactory.create_layout
     - |todo|
   * - dolfin.PETScFactory.create_linear_operator
     - |todo|
   * - dolfin.PETScFactory.create_lu_solver
     - |todo|
   * - dolfin.PETScFactory.krylov_solver_methods
     - |todo|
   * - dolfin.PETScFactory.krylov_solver_preconditioners
     - |todo|
   * - dolfin.PETScFactory.lu_solver_methods
     - |todo|
   * - dolfin.PETScKrylovSolver.default_parameters
     - |todo|
   * - dolfin.PETScKrylovSolver.methods
     - |todo|
   * - dolfin.PETScKrylovSolver.monitor
     - |todo|
   * - dolfin.PETScKrylovSolver.mpi_comm
     - |todo|
   * - dolfin.PETScKrylovSolver.norm_type
     - |new|
   * - dolfin.PETScKrylovSolver.parameter_type
     - |todo|
   * - dolfin.PETScKrylovSolver.petsc_error
     - |todo|
   * - dolfin.PETScKrylovSolver.preconditioners
     - |todo|
   * - dolfin.PETScKrylovSolver.set_nonzero_guess
     - |todo|
   * - dolfin.PETScKrylovSolver.set_tolerances
     - |todo|
   * - dolfin.PETScKrylovSolver.str
     - |todo|
   * - dolfin.PETScKrylovSolver.update_parameters
     - |todo|
   * - dolfin.PETScLUSolver.default_parameters
     - |todo|
   * - dolfin.PETScLUSolver.methods
     - |todo|
   * - dolfin.PETScLUSolver.mpi_comm
     - |todo|
   * - dolfin.PETScLUSolver.parameter_type
     - |todo|
   * - dolfin.PETScLUSolver.set_from_options
     - |todo|
   * - dolfin.PETScLUSolver.set_operator
     - |todo|
   * - dolfin.PETScLUSolver.set_operators
     - |todo|
   * - dolfin.PETScLUSolver.str
     - |todo|
   * - dolfin.PETScLUSolver.update_parameters
     - |todo|
   * - dolfin.PETScLinearOperator
     - |todo|
   * - dolfin.PETScLinearOperator.id
     - |todo|
   * - dolfin.PETScLinearOperator.init_layout
     - |todo|
   * - dolfin.PETScLinearOperator.init_vector
     - |todo|
   * - dolfin.PETScLinearOperator.label
     - |todo|
   * - dolfin.PETScLinearOperator.local_range
     - |todo|
   * - dolfin.PETScLinearOperator.mat
     - |todo|
   * - dolfin.PETScLinearOperator.mpi_comm
     - |todo|
   * - dolfin.PETScLinearOperator.mult
     - |todo|
   * - dolfin.PETScLinearOperator.name
     - |todo|
   * - dolfin.PETScLinearOperator.parameters
     - |todo|
   * - dolfin.PETScLinearOperator.petsc_error
     - |todo|
   * - dolfin.PETScLinearOperator.rename
     - |todo|
   * - dolfin.PETScLinearOperator.shared_instance
     - |todo|
   * - dolfin.PETScLinearOperator.size
     - |todo|
   * - dolfin.PETScLinearOperator.str
     - |todo|
   * - dolfin.PETScMatrix.add
     - |todo|
   * - dolfin.PETScMatrix.add_local
     - |todo|
   * - dolfin.PETScMatrix.assign
     - |todo|
   * - dolfin.PETScMatrix.binary_dump
     - |todo|
   * - dolfin.PETScMatrix.ident_local
     - |todo|
   * - dolfin.PETScMatrix.is_symmetric
     - |todo|
   * - dolfin.PETScMatrix.petsc_error
     - |todo|
   * - dolfin.PETScMatrix.set_from_options
     - |todo|
   * - dolfin.PETScMatrix.set_local
     - |todo|
   * - dolfin.PETScMatrix.setrow
     - |todo|
   * - dolfin.PETScMatrix.shared_instance
     - |todo|
   * - dolfin.PETScMatrix.zero_local
     - |todo|
   * - dolfin.PETScObject
     - |todo|
   * - dolfin.PETScObject.petsc_error
     - |todo|
   * - dolfin.PETScPreconditioner.id
     - |todo|
   * - dolfin.PETScPreconditioner.label
     - |todo|
   * - dolfin.PETScPreconditioner.name
     - |todo|
   * - dolfin.PETScPreconditioner.parameters
     - |todo|
   * - dolfin.PETScPreconditioner.petsc_error
     - |todo|
   * - dolfin.PETScPreconditioner.rename
     - |todo|
   * - dolfin.PETScPreconditioner.set
     - |todo|
   * - dolfin.PETScPreconditioner.set_coordinates
     - |todo|
   * - dolfin.PETScPreconditioner.set_fieldsplit
     - |todo|
   * - dolfin.PETScPreconditioner.set_type
     - |todo|
   * - dolfin.PETScPreconditioner.str
     - |todo|
   * - dolfin.PETScSNESSolver.default_parameters
     - |todo|
   * - dolfin.PETScSNESSolver.get_options_prefix
     - |todo|
   * - dolfin.PETScSNESSolver.init
     - |todo|
   * - dolfin.PETScSNESSolver.methods
     - |todo|
   * - dolfin.PETScSNESSolver.mpi_comm
     - |todo|
   * - dolfin.PETScSNESSolver.petsc_error
     - |todo|
   * - dolfin.PETScSNESSolver.set_from_options
     - |todo|
   * - dolfin.PETScSNESSolver.set_options_prefix
     - |todo|
   * - dolfin.PETScSNESSolver.snes
     - |todo|
   * - dolfin.PETScTAOSolver.default_parameters
     - |todo|
   * - dolfin.PETScTAOSolver.init
     - |todo|
   * - dolfin.PETScTAOSolver.methods
     - |todo|
   * - dolfin.PETScTAOSolver.mpi_comm
     - |todo|
   * - dolfin.PETScTAOSolver.petsc_error
     - |todo|
   * - dolfin.PETScTAOSolver.tao
     - |todo|
   * - dolfin.PETScVector.abs
     - |todo|
   * - dolfin.PETScVector.add
     - |todo|
   * - dolfin.PETScVector.petsc_error
     - |todo|
   * - dolfin.PETScVector.reset
     - |todo|
   * - dolfin.PETScVector.set_from_options
     - |todo|
   * - dolfin.PETScVector.shared_instance
     - |todo|
   * - dolfin.PROGRESS
     - |todo|
   * - dolfin.ParameterValue
     - |todo|
   * - dolfin.ParameterValue.Type_Bool
     - |todo|
   * - dolfin.ParameterValue.Type_Float
     - |todo|
   * - dolfin.ParameterValue.Type_Int
     - |todo|
   * - dolfin.ParameterValue.Type_String
     - |todo|
   * - dolfin.ParameterValue.access_count
     - |todo|
   * - dolfin.ParameterValue.change_count
     - |todo|
   * - dolfin.ParameterValue.check_key
     - |todo|
   * - dolfin.ParameterValue.data
     - |todo|
   * - dolfin.ParameterValue.description
     - |todo|
   * - dolfin.ParameterValue.get_range
     - |todo|
   * - dolfin.ParameterValue.is_set
     - |todo|
   * - dolfin.ParameterValue.key
     - |todo|
   * - dolfin.ParameterValue.range_str
     - |todo|
   * - dolfin.ParameterValue.reset
     - |todo|
   * - dolfin.ParameterValue.set_range
     - |todo|
   * - dolfin.ParameterValue.str
     - |todo|
   * - dolfin.ParameterValue.type_str
     - |todo|
   * - dolfin.ParameterValue.value
     - |todo|
   * - dolfin.ParameterValue.value_str
     - |todo|
   * - dolfin.ParameterValue.warn_once
     - |todo|
   * - dolfin.Parameters.add_unset
     - |todo|
   * - dolfin.Parameters.clear
     - |todo|
   * - dolfin.Parameters.find_parameter
     - |todo|
   * - dolfin.Parameters.find_parameter_set
     - |todo|
   * - dolfin.Parameters.get
     - |todo|
   * - dolfin.Parameters.has_key
     - |todo|
   * - dolfin.Parameters.iterdata
     - |todo|
   * - dolfin.Parameters.iteritems
     - |todo|
   * - dolfin.Parameters.iterkeys
     - |todo|
   * - dolfin.Parameters.itervalues
     - |todo|
   * - dolfin.Parameters.option_string
     - |todo|
   * - dolfin.Parameters.remove
     - |todo|
   * - dolfin.Parameters.size
     - |todo|
   * - dolfin.Parameters.to_dict
     - |todo|
   * - dolfin.Parameters.values
     - |todo|
   * - dolfin.PermutationSymbol
     - |todo|
   * - dolfin.PermutationSymbol.T
     - |todo|
   * - dolfin.PermutationSymbol.dx
     - |todo|
   * - dolfin.PermutationSymbol.evaluate
     - |todo|
   * - dolfin.PermutationSymbol.geometric_dimension
     - |todo|
   * - dolfin.PermutationSymbol.is_cellwise_constant
     - |todo|
   * - dolfin.PermutationSymbol.ufl_disable_profiling
     - |todo|
   * - dolfin.PermutationSymbol.ufl_domain
     - |todo|
   * - dolfin.PermutationSymbol.ufl_domains
     - |todo|
   * - dolfin.PermutationSymbol.ufl_enable_profiling
     - |todo|
   * - dolfin.PermutationSymbol.ufl_free_indices
     - |todo|
   * - dolfin.PermutationSymbol.ufl_index_dimensions
     - |todo|
   * - dolfin.PermutationSymbol.ufl_operands
     - |todo|
   * - dolfin.PermutationSymbol.ufl_shape
     - |todo|
   * - dolfin.Point.cross
     - |todo|
   * - dolfin.Point.dot
     - |todo|
   * - dolfin.Point.rotate
     - |todo|
   * - dolfin.Point.squared_distance
     - |todo|
   * - dolfin.Point.squared_norm
     - |todo|
   * - dolfin.Point.str
     - |todo|
   * - dolfin.PointIntegralSolver.default_parameters
     - |todo|
   * - dolfin.PointIntegralSolver.id
     - |todo|
   * - dolfin.PointIntegralSolver.label
     - |todo|
   * - dolfin.PointIntegralSolver.name
     - |todo|
   * - dolfin.PointIntegralSolver.num_jacobian_computations
     - |todo|
   * - dolfin.PointIntegralSolver.parameters
     - |todo|
   * - dolfin.PointIntegralSolver.rename
     - |todo|
   * - dolfin.PointIntegralSolver.str
     - |todo|
   * - dolfin.Progress
     - |todo|
   * - dolfin.Progress.update
     - |todo|
   * - dolfin.RK4.bcs
     - |todo|
   * - dolfin.RK4.dt_stage_offset
     - |todo|
   * - dolfin.RK4.id
     - |todo|
   * - dolfin.RK4.implicit
     - |todo|
   * - dolfin.RK4.jacobian_index
     - |todo|
   * - dolfin.RK4.label
     - |todo|
   * - dolfin.RK4.name
     - |todo|
   * - dolfin.RK4.parameters
     - |todo|
   * - dolfin.RK4.rename
     - |todo|
   * - dolfin.RK4.str
     - |todo|
   * - dolfin.RKSolver.step
     - |todo|
   * - dolfin.RL1.bcs
     - |todo|
   * - dolfin.RL1.dt_stage_offset
     - |todo|
   * - dolfin.RL1.id
     - |todo|
   * - dolfin.RL1.implicit
     - |todo|
   * - dolfin.RL1.jacobian_index
     - |todo|
   * - dolfin.RL1.label
     - |todo|
   * - dolfin.RL1.name
     - |todo|
   * - dolfin.RL1.parameters
     - |todo|
   * - dolfin.RL1.rename
     - |todo|
   * - dolfin.RL1.str
     - |todo|
   * - dolfin.RL2.bcs
     - |todo|
   * - dolfin.RL2.dt_stage_offset
     - |todo|
   * - dolfin.RL2.id
     - |todo|
   * - dolfin.RL2.implicit
     - |todo|
   * - dolfin.RL2.jacobian_index
     - |todo|
   * - dolfin.RL2.label
     - |todo|
   * - dolfin.RL2.name
     - |todo|
   * - dolfin.RL2.parameters
     - |todo|
   * - dolfin.RL2.rename
     - |todo|
   * - dolfin.RL2.str
     - |todo|
   * - dolfin.RTLD_GLOBAL
     - |new|
   * - dolfin.RTLD_NOW
     - |new|
   * - dolfin.RectangleMesh.cell_name
     - |new|
   * - dolfin.RectangleMesh.child
     - |todo|
   * - dolfin.RectangleMesh.clean
     - |todo|
   * - dolfin.RectangleMesh.clear_child
     - |todo|
   * - dolfin.RectangleMesh.create
     - |todo|
   * - dolfin.RectangleMesh.depth
     - |todo|
   * - dolfin.RectangleMesh.geometric_dimension
     - |new|
   * - dolfin.RectangleMesh.ghost_mode
     - |todo|
   * - dolfin.RectangleMesh.has_child
     - |todo|
   * - dolfin.RectangleMesh.has_parent
     - |todo|
   * - dolfin.RectangleMesh.label
     - |todo|
   * - dolfin.RectangleMesh.leaf_node
     - |todo|
   * - dolfin.RectangleMesh.name
     - |todo|
   * - dolfin.RectangleMesh.order
     - |todo|
   * - dolfin.RectangleMesh.parameters
     - |todo|
   * - dolfin.RectangleMesh.parent
     - |todo|
   * - dolfin.RectangleMesh.rename
     - |todo|
   * - dolfin.RectangleMesh.renumber_by_color
     - |todo|
   * - dolfin.RectangleMesh.root_node
     - |todo|
   * - dolfin.RectangleMesh.scale
     - |todo|
   * - dolfin.RectangleMesh.set_child
     - |todo|
   * - dolfin.RectangleMesh.set_parent
     - |todo|
   * - dolfin.RectangleMesh.str
     - |todo|
   * - dolfin.RestrictedElement
     - |todo|
   * - dolfin.RestrictedElement.cell
     - |todo|
   * - dolfin.RestrictedElement.degree
     - |todo|
   * - dolfin.RestrictedElement.extract_component
     - |todo|
   * - dolfin.RestrictedElement.extract_reference_component
     - |todo|
   * - dolfin.RestrictedElement.extract_subelement_component
     - |todo|
   * - dolfin.RestrictedElement.extract_subelement_reference_component
     - |todo|
   * - dolfin.RestrictedElement.family
     - |todo|
   * - dolfin.RestrictedElement.is_cellwise_constant
     - |todo|
   * - dolfin.RestrictedElement.mapping
     - |todo|
   * - dolfin.RestrictedElement.num_restricted_sub_elements
     - |todo|
   * - dolfin.RestrictedElement.num_sub_elements
     - |todo|
   * - dolfin.RestrictedElement.quadrature_scheme
     - |todo|
   * - dolfin.RestrictedElement.reconstruct
     - |todo|
   * - dolfin.RestrictedElement.reference_value_shape
     - |todo|
   * - dolfin.RestrictedElement.reference_value_size
     - |todo|
   * - dolfin.RestrictedElement.restricted_sub_elements
     - |todo|
   * - dolfin.RestrictedElement.restriction_domain
     - |todo|
   * - dolfin.RestrictedElement.shortstr
     - |todo|
   * - dolfin.RestrictedElement.sub_element
     - |todo|
   * - dolfin.RestrictedElement.sub_elements
     - |todo|
   * - dolfin.RestrictedElement.symmetry
     - |todo|
   * - dolfin.RestrictedElement.value_shape
     - |todo|
   * - dolfin.RestrictedElement.value_size
     - |todo|
   * - dolfin.RushLarsenScheme
     - |todo|
   * - dolfin.RushLarsenScheme.bcs
     - |todo|
   * - dolfin.RushLarsenScheme.dolfin_stage_forms
     - |todo|
   * - dolfin.RushLarsenScheme.dt
     - |todo|
   * - dolfin.RushLarsenScheme.dt_stage_offset
     - |todo|
   * - dolfin.RushLarsenScheme.id
     - |todo|
   * - dolfin.RushLarsenScheme.implicit
     - |todo|
   * - dolfin.RushLarsenScheme.jacobian_index
     - |todo|
   * - dolfin.RushLarsenScheme.label
     - |todo|
   * - dolfin.RushLarsenScheme.last_stage
     - |todo|
   * - dolfin.RushLarsenScheme.name
     - |todo|
   * - dolfin.RushLarsenScheme.order
     - |todo|
   * - dolfin.RushLarsenScheme.parameters
     - |todo|
   * - dolfin.RushLarsenScheme.rename
     - |todo|
   * - dolfin.RushLarsenScheme.rhs_form
     - |todo|
   * - dolfin.RushLarsenScheme.solution
     - |todo|
   * - dolfin.RushLarsenScheme.stage_solutions
     - |todo|
   * - dolfin.RushLarsenScheme.str
     - |todo|
   * - dolfin.RushLarsenScheme.t
     - |todo|
   * - dolfin.RushLarsenScheme.to_adm
     - |todo|
   * - dolfin.RushLarsenScheme.to_tlm
     - |todo|
   * - dolfin.RushLarsenScheme.ufl_stage_forms
     - |todo|
   * - dolfin.SCOTCH
     - |todo|
   * - dolfin.SCOTCH.compute_gps
     - |todo|
   * - dolfin.SCOTCH.compute_reordering
     - |todo|
   * - dolfin.SLEPcEigenSolver.default_parameters
     - |todo|
   * - dolfin.SLEPcEigenSolver.eps
     - |todo|
   * - dolfin.SLEPcEigenSolver.get_iteration_number
     - |todo|
   * - dolfin.SLEPcEigenSolver.petsc_error
     - |todo|
   * - dolfin.SLEPcEigenSolver.str
     - |todo|
   * - dolfin.Scalar.add
     - |todo|
   * - dolfin.Scalar.add_local
     - |todo|
   * - dolfin.Scalar.copy
     - |todo|
   * - dolfin.Scalar.set_local
     - |todo|
   * - dolfin.Scalar.shared_instance
     - |todo|
   * - dolfin.SparsityPattern
     - |todo|
   * - dolfin.SparsityPattern.Type_sorted
     - |todo|
   * - dolfin.SparsityPattern.Type_unsorted
     - |todo|
   * - dolfin.SparsityPattern.apply
     - |todo|
   * - dolfin.SparsityPattern.init
     - |todo|
   * - dolfin.SparsityPattern.insert_full_rows_local
     - |todo|
   * - dolfin.SparsityPattern.insert_global
     - |todo|
   * - dolfin.SparsityPattern.insert_local
     - |todo|
   * - dolfin.SparsityPattern.insert_local_global
     - |todo|
   * - dolfin.SparsityPattern.local_range
     - |todo|
   * - dolfin.SparsityPattern.mpi_comm
     - |todo|
   * - dolfin.SparsityPattern.num_local_nonzeros
     - |todo|
   * - dolfin.SparsityPattern.num_nonzeros
     - |todo|
   * - dolfin.SparsityPattern.num_nonzeros_diagonal
     - |todo|
   * - dolfin.SparsityPattern.num_nonzeros_off_diagonal
     - |todo|
   * - dolfin.SparsityPattern.primary_dim
     - |todo|
   * - dolfin.SparsityPattern.rank
     - |todo|
   * - dolfin.SparsityPattern.str
     - |todo|
   * - dolfin.SparsityPatternBuilder.build_multimesh_sparsity_pattern
     - |todo|
   * - dolfin.SpecialFacetFunction
     - |todo|
   * - dolfin.SpecialFacetFunction.compute_vertex_values
     - |todo|
   * - dolfin.SpecialFacetFunction.eval
     - |todo|
   * - dolfin.SpecialFacetFunction.eval_cell
     - |todo|
   * - dolfin.SpecialFacetFunction.get_generic_function
     - |todo|
   * - dolfin.SpecialFacetFunction.get_property
     - |todo|
   * - dolfin.SpecialFacetFunction.id
     - |todo|
   * - dolfin.SpecialFacetFunction.label
     - |todo|
   * - dolfin.SpecialFacetFunction.name
     - |todo|
   * - dolfin.SpecialFacetFunction.parameters
     - |todo|
   * - dolfin.SpecialFacetFunction.rename
     - |todo|
   * - dolfin.SpecialFacetFunction.restrict
     - |todo|
   * - dolfin.SpecialFacetFunction.set_generic_function
     - |todo|
   * - dolfin.SpecialFacetFunction.set_property
     - |todo|
   * - dolfin.SpecialFacetFunction.str
     - |todo|
   * - dolfin.SpecialFacetFunction.update
     - |todo|
   * - dolfin.SpecialFacetFunction.value_dimension
     - |todo|
   * - dolfin.SpecialFacetFunction.value_rank
     - |todo|
   * - dolfin.SpecialFacetFunction.value_shape
     - |todo|
   * - dolfin.SpecialFacetFunction.value_size
     - |todo|
   * - dolfin.SubDomain.geometric_dimension
     - |todo|
   * - dolfin.SubDomain.map_tolerance
     - |todo|
   * - dolfin.SubDomain.snap
     - |todo|
   * - dolfin.SubMesh.cell_name
     - |new|
   * - dolfin.SubMesh.child
     - |todo|
   * - dolfin.SubMesh.clean
     - |todo|
   * - dolfin.SubMesh.clear_child
     - |todo|
   * - dolfin.SubMesh.depth
     - |todo|
   * - dolfin.SubMesh.geometric_dimension
     - |new|
   * - dolfin.SubMesh.ghost_mode
     - |todo|
   * - dolfin.SubMesh.has_child
     - |todo|
   * - dolfin.SubMesh.has_parent
     - |todo|
   * - dolfin.SubMesh.label
     - |todo|
   * - dolfin.SubMesh.leaf_node
     - |todo|
   * - dolfin.SubMesh.name
     - |todo|
   * - dolfin.SubMesh.order
     - |todo|
   * - dolfin.SubMesh.parameters
     - |todo|
   * - dolfin.SubMesh.parent
     - |todo|
   * - dolfin.SubMesh.rename
     - |todo|
   * - dolfin.SubMesh.renumber_by_color
     - |todo|
   * - dolfin.SubMesh.root_node
     - |todo|
   * - dolfin.SubMesh.scale
     - |todo|
   * - dolfin.SubMesh.set_child
     - |todo|
   * - dolfin.SubMesh.set_parent
     - |todo|
   * - dolfin.SubMesh.str
     - |todo|
   * - dolfin.SubSystemsManager.init_mpi
     - |todo|
   * - dolfin.SubSystemsManager.petsc_err_msg
     - |todo|
   * - dolfin.SubSystemsManager.singleton
     - |todo|
   * - dolfin.SubsetIterator.end
     - |todo|
   * - dolfin.SubsetIterator.end_iterator
     - |todo|
   * - dolfin.SubsetIterator.next
     - |todo|
   * - dolfin.SystemAssembler.init_global_tensor
     - |todo|
   * - dolfin.TAOLinearBoundSolver.default_parameters
     - |todo|
   * - dolfin.TAOLinearBoundSolver.get_matrix
     - |todo|
   * - dolfin.TAOLinearBoundSolver.get_vector
     - |todo|
   * - dolfin.TAOLinearBoundSolver.id
     - |todo|
   * - dolfin.TAOLinearBoundSolver.krylov_solvers
     - |todo|
   * - dolfin.TAOLinearBoundSolver.label
     - |todo|
   * - dolfin.TAOLinearBoundSolver.methods
     - |todo|
   * - dolfin.TAOLinearBoundSolver.name
     - |todo|
   * - dolfin.TAOLinearBoundSolver.parameters
     - |todo|
   * - dolfin.TAOLinearBoundSolver.petsc_error
     - |todo|
   * - dolfin.TAOLinearBoundSolver.preconditioners
     - |todo|
   * - dolfin.TAOLinearBoundSolver.rename
     - |todo|
   * - dolfin.TAOLinearBoundSolver.set_ksp
     - |todo|
   * - dolfin.TAOLinearBoundSolver.set_solver
     - |todo|
   * - dolfin.TAOLinearBoundSolver.str
     - |todo|
   * - dolfin.TAOLinearBoundSolver.tao
     - |todo|
   * - dolfin.TRACE
     - |todo|
   * - dolfin.Table.get
     - |todo|
   * - dolfin.Table.get_value
     - |todo|
   * - dolfin.Table.id
     - |todo|
   * - dolfin.Table.label
     - |todo|
   * - dolfin.Table.name
     - |todo|
   * - dolfin.Table.parameters
     - |todo|
   * - dolfin.Table.rename
     - |todo|
   * - dolfin.Table.set
     - |todo|
   * - dolfin.Table.str_latex
     - |todo|
   * - dolfin.TableEntry
     - |todo|
   * - dolfin.TensorConstant
     - |todo|
   * - dolfin.TensorLayout.Ghosts
     - |new|
   * - dolfin.TensorLayout.Sparsity
     - |new|
   * - dolfin.TensorLayout.id
     - |todo|
   * - dolfin.TensorLayout.index_map
     - |todo|
   * - dolfin.TensorLayout.is_ghosted
     - |todo|
   * - dolfin.TensorLayout.label
     - |todo|
   * - dolfin.TensorLayout.local_range
     - |todo|
   * - dolfin.TensorLayout.mpi_comm
     - |todo|
   * - dolfin.TensorLayout.name
     - |todo|
   * - dolfin.TensorLayout.parameters
     - |todo|
   * - dolfin.TensorLayout.primary_dim
     - |todo|
   * - dolfin.TensorLayout.rank
     - |todo|
   * - dolfin.TensorLayout.rename
     - |todo|
   * - dolfin.TensorLayout.size
     - |todo|
   * - dolfin.TensorLayout.str
     - |todo|
   * - dolfin.TensorProductCell
     - |todo|
   * - dolfin.TensorProductCell.cellname
     - |todo|
   * - dolfin.TensorProductCell.geometric_dimension
     - |todo|
   * - dolfin.TensorProductCell.has_simplex_facets
     - |todo|
   * - dolfin.TensorProductCell.is_simplex
     - |todo|
   * - dolfin.TensorProductCell.num_edges
     - |todo|
   * - dolfin.TensorProductCell.num_facets
     - |todo|
   * - dolfin.TensorProductCell.num_vertices
     - |todo|
   * - dolfin.TensorProductCell.reconstruct
     - |todo|
   * - dolfin.TensorProductCell.sub_cells
     - |todo|
   * - dolfin.TensorProductCell.topological_dimension
     - |todo|
   * - dolfin.TensorProductElement
     - |todo|
   * - dolfin.TensorProductElement.cell
     - |todo|
   * - dolfin.TensorProductElement.degree
     - |todo|
   * - dolfin.TensorProductElement.extract_component
     - |todo|
   * - dolfin.TensorProductElement.extract_reference_component
     - |todo|
   * - dolfin.TensorProductElement.extract_subelement_component
     - |todo|
   * - dolfin.TensorProductElement.extract_subelement_reference_component
     - |todo|
   * - dolfin.TensorProductElement.family
     - |todo|
   * - dolfin.TensorProductElement.is_cellwise_constant
     - |todo|
   * - dolfin.TensorProductElement.mapping
     - |todo|
   * - dolfin.TensorProductElement.num_sub_elements
     - |todo|
   * - dolfin.TensorProductElement.quadrature_scheme
     - |todo|
   * - dolfin.TensorProductElement.reconstruct
     - |todo|
   * - dolfin.TensorProductElement.reference_value_shape
     - |todo|
   * - dolfin.TensorProductElement.reference_value_size
     - |todo|
   * - dolfin.TensorProductElement.shortstr
     - |todo|
   * - dolfin.TensorProductElement.sobolev_space
     - |todo|
   * - dolfin.TensorProductElement.sub_elements
     - |todo|
   * - dolfin.TensorProductElement.symmetry
     - |todo|
   * - dolfin.TensorProductElement.value_shape
     - |todo|
   * - dolfin.TensorProductElement.value_size
     - |todo|
   * - dolfin.TensorProductMesh
     - |todo|
   * - dolfin.TensorProductMesh.geometric_dimension
     - |todo|
   * - dolfin.TensorProductMesh.is_piecewise_linear_simplex_domain
     - |todo|
   * - dolfin.TensorProductMesh.topological_dimension
     - |todo|
   * - dolfin.TensorProductMesh.ufl_cell
     - |todo|
   * - dolfin.TensorProductMesh.ufl_coordinate_element
     - |todo|
   * - dolfin.TensorProductMesh.ufl_id
     - |todo|
   * - dolfin.Time
     - |todo|
   * - dolfin.TimeSeries.clear
     - |todo|
   * - dolfin.TimeSeries.default_parameters
     - |todo|
   * - dolfin.TimeSeries.id
     - |todo|
   * - dolfin.TimeSeries.label
     - |todo|
   * - dolfin.TimeSeries.name
     - |todo|
   * - dolfin.TimeSeries.parameters
     - |todo|
   * - dolfin.TimeSeries.rename
     - |todo|
   * - dolfin.TimeSeries.str
     - |todo|
   * - dolfin.TimingClear
     - |new|
   * - dolfin.TimingClear.clear
     - |new|
   * - dolfin.TimingClear.keep
     - |new|
   * - dolfin.TimingType
     - |new|
   * - dolfin.TimingType.system
     - |new|
   * - dolfin.TimingType.user
     - |new|
   * - dolfin.TimingType.wall
     - |new|
   * - dolfin.UFC_SIGNATURE
     - |todo|
   * - dolfin.UFC_VERSION
     - |todo|
   * - dolfin.UFC_VERSION_MAINTENANCE
     - |todo|
   * - dolfin.UFC_VERSION_MAJOR
     - |todo|
   * - dolfin.UFC_VERSION_MINOR
     - |todo|
   * - dolfin.UFC_VERSION_RELEASE
     - |todo|
   * - dolfin.UFLException
     - |todo|
   * - dolfin.UFLException.args
     - |todo|
   * - dolfin.UFLException.with_traceback
     - |todo|
   * - dolfin.UIntArray
     - |todo|
   * - dolfin.UIntArray.array
     - |todo|
   * - dolfin.UIntArray.data
     - |todo|
   * - dolfin.UIntArray.size
     - |todo|
   * - dolfin.UIntArray.str
     - |todo|
   * - dolfin.UnitCubeMesh.cell_name
     - |new|
   * - dolfin.UnitCubeMesh.child
     - |todo|
   * - dolfin.UnitCubeMesh.clean
     - |todo|
   * - dolfin.UnitCubeMesh.clear_child
     - |todo|
   * - dolfin.UnitCubeMesh.create
     - |todo|
   * - dolfin.UnitCubeMesh.depth
     - |todo|
   * - dolfin.UnitCubeMesh.geometric_dimension
     - |new|
   * - dolfin.UnitCubeMesh.ghost_mode
     - |todo|
   * - dolfin.UnitCubeMesh.has_child
     - |todo|
   * - dolfin.UnitCubeMesh.has_parent
     - |todo|
   * - dolfin.UnitCubeMesh.label
     - |todo|
   * - dolfin.UnitCubeMesh.leaf_node
     - |todo|
   * - dolfin.UnitCubeMesh.name
     - |todo|
   * - dolfin.UnitCubeMesh.order
     - |todo|
   * - dolfin.UnitCubeMesh.parameters
     - |todo|
   * - dolfin.UnitCubeMesh.parent
     - |todo|
   * - dolfin.UnitCubeMesh.rename
     - |todo|
   * - dolfin.UnitCubeMesh.renumber_by_color
     - |todo|
   * - dolfin.UnitCubeMesh.root_node
     - |todo|
   * - dolfin.UnitCubeMesh.scale
     - |todo|
   * - dolfin.UnitCubeMesh.set_child
     - |todo|
   * - dolfin.UnitCubeMesh.set_parent
     - |todo|
   * - dolfin.UnitCubeMesh.str
     - |todo|
   * - dolfin.UnitIntervalMesh.cell_name
     - |new|
   * - dolfin.UnitIntervalMesh.child
     - |todo|
   * - dolfin.UnitIntervalMesh.clean
     - |todo|
   * - dolfin.UnitIntervalMesh.clear_child
     - |todo|
   * - dolfin.UnitIntervalMesh.depth
     - |todo|
   * - dolfin.UnitIntervalMesh.geometric_dimension
     - |new|
   * - dolfin.UnitIntervalMesh.ghost_mode
     - |todo|
   * - dolfin.UnitIntervalMesh.has_child
     - |todo|
   * - dolfin.UnitIntervalMesh.has_parent
     - |todo|
   * - dolfin.UnitIntervalMesh.label
     - |todo|
   * - dolfin.UnitIntervalMesh.leaf_node
     - |todo|
   * - dolfin.UnitIntervalMesh.name
     - |todo|
   * - dolfin.UnitIntervalMesh.order
     - |todo|
   * - dolfin.UnitIntervalMesh.parameters
     - |todo|
   * - dolfin.UnitIntervalMesh.parent
     - |todo|
   * - dolfin.UnitIntervalMesh.rename
     - |todo|
   * - dolfin.UnitIntervalMesh.renumber_by_color
     - |todo|
   * - dolfin.UnitIntervalMesh.root_node
     - |todo|
   * - dolfin.UnitIntervalMesh.scale
     - |todo|
   * - dolfin.UnitIntervalMesh.set_child
     - |todo|
   * - dolfin.UnitIntervalMesh.set_parent
     - |todo|
   * - dolfin.UnitIntervalMesh.str
     - |todo|
   * - dolfin.UnitQuadMesh.bounding_box_tree
     - |todo|
   * - dolfin.UnitQuadMesh.cell_orientations
     - |todo|
   * - dolfin.UnitQuadMesh.cells
     - |todo|
   * - dolfin.UnitQuadMesh.child
     - |todo|
   * - dolfin.UnitQuadMesh.clean
     - |todo|
   * - dolfin.UnitQuadMesh.clear_child
     - |todo|
   * - dolfin.UnitQuadMesh.color
     - |todo|
   * - dolfin.UnitQuadMesh.coordinates
     - |todo|
   * - dolfin.UnitQuadMesh.data
     - |todo|
   * - dolfin.UnitQuadMesh.depth
     - |todo|
   * - dolfin.UnitQuadMesh.domains
     - |todo|
   * - dolfin.UnitQuadMesh.geometry
     - |todo|
   * - dolfin.UnitQuadMesh.ghost_mode
     - |todo|
   * - dolfin.UnitQuadMesh.has_child
     - |todo|
   * - dolfin.UnitQuadMesh.has_parent
     - |todo|
   * - dolfin.UnitQuadMesh.hash
     - |todo|
   * - dolfin.UnitQuadMesh.hmax
     - |todo|
   * - dolfin.UnitQuadMesh.hmin
     - |todo|
   * - dolfin.UnitQuadMesh.id
     - |todo|
   * - dolfin.UnitQuadMesh.init
     - |todo|
   * - dolfin.UnitQuadMesh.init_cell_orientations
     - |todo|
   * - dolfin.UnitQuadMesh.init_global
     - |todo|
   * - dolfin.UnitQuadMesh.label
     - |todo|
   * - dolfin.UnitQuadMesh.leaf_node
     - |todo|
   * - dolfin.UnitQuadMesh.mpi_comm
     - |todo|
   * - dolfin.UnitQuadMesh.name
     - |todo|
   * - dolfin.UnitQuadMesh.num_cells
     - |todo|
   * - dolfin.UnitQuadMesh.num_edges
     - |todo|
   * - dolfin.UnitQuadMesh.num_entities
     - |todo|
   * - dolfin.UnitQuadMesh.num_faces
     - |todo|
   * - dolfin.UnitQuadMesh.num_facets
     - |todo|
   * - dolfin.UnitQuadMesh.num_vertices
     - |todo|
   * - dolfin.UnitQuadMesh.order
     - |todo|
   * - dolfin.UnitQuadMesh.ordered
     - |todo|
   * - dolfin.UnitQuadMesh.parameters
     - |todo|
   * - dolfin.UnitQuadMesh.parent
     - |todo|
   * - dolfin.UnitQuadMesh.rename
     - |todo|
   * - dolfin.UnitQuadMesh.renumber_by_color
     - |todo|
   * - dolfin.UnitQuadMesh.rmax
     - |todo|
   * - dolfin.UnitQuadMesh.rmin
     - |todo|
   * - dolfin.UnitQuadMesh.root_node
     - |todo|
   * - dolfin.UnitQuadMesh.rotate
     - |todo|
   * - dolfin.UnitQuadMesh.scale
     - |todo|
   * - dolfin.UnitQuadMesh.set_child
     - |todo|
   * - dolfin.UnitQuadMesh.set_parent
     - |todo|
   * - dolfin.UnitQuadMesh.size
     - |todo|
   * - dolfin.UnitQuadMesh.size_global
     - |todo|
   * - dolfin.UnitQuadMesh.smooth
     - |todo|
   * - dolfin.UnitQuadMesh.smooth_boundary
     - |todo|
   * - dolfin.UnitQuadMesh.snap_boundary
     - |todo|
   * - dolfin.UnitQuadMesh.str
     - |todo|
   * - dolfin.UnitQuadMesh.topology
     - |todo|
   * - dolfin.UnitQuadMesh.translate
     - |todo|
   * - dolfin.UnitQuadMesh.type
     - |todo|
   * - dolfin.UnitQuadMesh.ufl_cell
     - |todo|
   * - dolfin.UnitQuadMesh.ufl_coordinate_element
     - |todo|
   * - dolfin.UnitQuadMesh.ufl_domain
     - |todo|
   * - dolfin.UnitQuadMesh.ufl_id
     - |todo|
   * - dolfin.UnitSquareMesh.cell_name
     - |new|
   * - dolfin.UnitSquareMesh.child
     - |todo|
   * - dolfin.UnitSquareMesh.clean
     - |todo|
   * - dolfin.UnitSquareMesh.clear_child
     - |todo|
   * - dolfin.UnitSquareMesh.create
     - |todo|
   * - dolfin.UnitSquareMesh.depth
     - |todo|
   * - dolfin.UnitSquareMesh.geometric_dimension
     - |new|
   * - dolfin.UnitSquareMesh.ghost_mode
     - |todo|
   * - dolfin.UnitSquareMesh.has_child
     - |todo|
   * - dolfin.UnitSquareMesh.has_parent
     - |todo|
   * - dolfin.UnitSquareMesh.label
     - |todo|
   * - dolfin.UnitSquareMesh.leaf_node
     - |todo|
   * - dolfin.UnitSquareMesh.name
     - |todo|
   * - dolfin.UnitSquareMesh.order
     - |todo|
   * - dolfin.UnitSquareMesh.parameters
     - |todo|
   * - dolfin.UnitSquareMesh.parent
     - |todo|
   * - dolfin.UnitSquareMesh.rename
     - |todo|
   * - dolfin.UnitSquareMesh.renumber_by_color
     - |todo|
   * - dolfin.UnitSquareMesh.root_node
     - |todo|
   * - dolfin.UnitSquareMesh.scale
     - |todo|
   * - dolfin.UnitSquareMesh.set_child
     - |todo|
   * - dolfin.UnitSquareMesh.set_parent
     - |todo|
   * - dolfin.UnitSquareMesh.str
     - |todo|
   * - dolfin.UnitTetrahedronMesh
     - |todo|
   * - dolfin.UnitTetrahedronMesh.create
     - |todo|
   * - dolfin.UserExpression.T
     - |new|
   * - dolfin.UserExpression.compute_vertex_values
     - |new|
   * - dolfin.UserExpression.count
     - |new|
   * - dolfin.UserExpression.cpp_object
     - |new|
   * - dolfin.UserExpression.dx
     - |new|
   * - dolfin.UserExpression.evaluate
     - |new|
   * - dolfin.UserExpression.geometric_dimension
     - |new|
   * - dolfin.UserExpression.id
     - |new|
   * - dolfin.UserExpression.is_cellwise_constant
     - |new|
   * - dolfin.UserExpression.label
     - |new|
   * - dolfin.UserExpression.name
     - |new|
   * - dolfin.UserExpression.str
     - |todo|
   * - dolfin.UserExpression.ufl_disable_profiling
     - |new|
   * - dolfin.UserExpression.ufl_domain
     - |new|
   * - dolfin.UserExpression.ufl_domains
     - |new|
   * - dolfin.UserExpression.ufl_element
     - |new|
   * - dolfin.UserExpression.ufl_enable_profiling
     - |new|
   * - dolfin.UserExpression.ufl_free_indices
     - |new|
   * - dolfin.UserExpression.ufl_function_space
     - |new|
   * - dolfin.UserExpression.ufl_index_dimensions
     - |new|
   * - dolfin.UserExpression.ufl_operands
     - |new|
   * - dolfin.UserExpression.ufl_shape
     - |new|
   * - dolfin.UserExpression.value_dimension
     - |new|
   * - dolfin.UserExpression.value_rank
     - |new|
   * - dolfin.UserExpression.value_shape
     - |todo|
   * - dolfin.VTKFile
     - |new|
   * - dolfin.VTKFile.write
     - |new|
   * - dolfin.Variable.str
     - |todo|
   * - dolfin.Vector.add
     - |todo|
   * - dolfin.Vector.instance
     - |new|
   * - dolfin.Vector.shared_instance
     - |todo|
   * - dolfin.VectorConstant
     - |todo|
   * - dolfin.Vertex.incident
     - |todo|
   * - dolfin.Vertex.init
     - |todo|
   * - dolfin.Vertex.mesh_id
     - |todo|
   * - dolfin.Vertex.owner
     - |todo|
   * - dolfin.Vertex.str
     - |todo|
   * - dolfin.Vertex.x
     - |todo|
   * - dolfin.VertexFunctionBool
     - |todo|
   * - dolfin.VertexFunctionBool.array
     - |todo|
   * - dolfin.VertexFunctionBool.child
     - |todo|
   * - dolfin.VertexFunctionBool.clear_child
     - |todo|
   * - dolfin.VertexFunctionBool.cpp_value_type
     - |todo|
   * - dolfin.VertexFunctionBool.depth
     - |todo|
   * - dolfin.VertexFunctionBool.dim
     - |todo|
   * - dolfin.VertexFunctionBool.empty
     - |todo|
   * - dolfin.VertexFunctionBool.has_child
     - |todo|
   * - dolfin.VertexFunctionBool.has_parent
     - |todo|
   * - dolfin.VertexFunctionBool.id
     - |todo|
   * - dolfin.VertexFunctionBool.init
     - |todo|
   * - dolfin.VertexFunctionBool.label
     - |todo|
   * - dolfin.VertexFunctionBool.leaf_node
     - |todo|
   * - dolfin.VertexFunctionBool.mesh
     - |todo|
   * - dolfin.VertexFunctionBool.name
     - |todo|
   * - dolfin.VertexFunctionBool.parameters
     - |todo|
   * - dolfin.VertexFunctionBool.parent
     - |todo|
   * - dolfin.VertexFunctionBool.rename
     - |todo|
   * - dolfin.VertexFunctionBool.root_node
     - |todo|
   * - dolfin.VertexFunctionBool.set_all
     - |todo|
   * - dolfin.VertexFunctionBool.set_child
     - |todo|
   * - dolfin.VertexFunctionBool.set_parent
     - |todo|
   * - dolfin.VertexFunctionBool.set_value
     - |todo|
   * - dolfin.VertexFunctionBool.set_values
     - |todo|
   * - dolfin.VertexFunctionBool.size
     - |todo|
   * - dolfin.VertexFunctionBool.str
     - |todo|
   * - dolfin.VertexFunctionBool.ufl_id
     - |todo|
   * - dolfin.VertexFunctionBool.where_equal
     - |todo|
   * - dolfin.VertexFunctionDouble
     - |todo|
   * - dolfin.VertexFunctionDouble.array
     - |todo|
   * - dolfin.VertexFunctionDouble.child
     - |todo|
   * - dolfin.VertexFunctionDouble.clear_child
     - |todo|
   * - dolfin.VertexFunctionDouble.cpp_value_type
     - |todo|
   * - dolfin.VertexFunctionDouble.depth
     - |todo|
   * - dolfin.VertexFunctionDouble.dim
     - |todo|
   * - dolfin.VertexFunctionDouble.empty
     - |todo|
   * - dolfin.VertexFunctionDouble.has_child
     - |todo|
   * - dolfin.VertexFunctionDouble.has_parent
     - |todo|
   * - dolfin.VertexFunctionDouble.id
     - |todo|
   * - dolfin.VertexFunctionDouble.init
     - |todo|
   * - dolfin.VertexFunctionDouble.label
     - |todo|
   * - dolfin.VertexFunctionDouble.leaf_node
     - |todo|
   * - dolfin.VertexFunctionDouble.mesh
     - |todo|
   * - dolfin.VertexFunctionDouble.name
     - |todo|
   * - dolfin.VertexFunctionDouble.parameters
     - |todo|
   * - dolfin.VertexFunctionDouble.parent
     - |todo|
   * - dolfin.VertexFunctionDouble.rename
     - |todo|
   * - dolfin.VertexFunctionDouble.root_node
     - |todo|
   * - dolfin.VertexFunctionDouble.set_all
     - |todo|
   * - dolfin.VertexFunctionDouble.set_child
     - |todo|
   * - dolfin.VertexFunctionDouble.set_parent
     - |todo|
   * - dolfin.VertexFunctionDouble.set_value
     - |todo|
   * - dolfin.VertexFunctionDouble.set_values
     - |todo|
   * - dolfin.VertexFunctionDouble.size
     - |todo|
   * - dolfin.VertexFunctionDouble.str
     - |todo|
   * - dolfin.VertexFunctionDouble.ufl_id
     - |todo|
   * - dolfin.VertexFunctionDouble.where_equal
     - |todo|
   * - dolfin.VertexFunctionInt
     - |todo|
   * - dolfin.VertexFunctionInt.array
     - |todo|
   * - dolfin.VertexFunctionInt.child
     - |todo|
   * - dolfin.VertexFunctionInt.clear_child
     - |todo|
   * - dolfin.VertexFunctionInt.cpp_value_type
     - |todo|
   * - dolfin.VertexFunctionInt.depth
     - |todo|
   * - dolfin.VertexFunctionInt.dim
     - |todo|
   * - dolfin.VertexFunctionInt.empty
     - |todo|
   * - dolfin.VertexFunctionInt.has_child
     - |todo|
   * - dolfin.VertexFunctionInt.has_parent
     - |todo|
   * - dolfin.VertexFunctionInt.id
     - |todo|
   * - dolfin.VertexFunctionInt.init
     - |todo|
   * - dolfin.VertexFunctionInt.label
     - |todo|
   * - dolfin.VertexFunctionInt.leaf_node
     - |todo|
   * - dolfin.VertexFunctionInt.mesh
     - |todo|
   * - dolfin.VertexFunctionInt.name
     - |todo|
   * - dolfin.VertexFunctionInt.parameters
     - |todo|
   * - dolfin.VertexFunctionInt.parent
     - |todo|
   * - dolfin.VertexFunctionInt.rename
     - |todo|
   * - dolfin.VertexFunctionInt.root_node
     - |todo|
   * - dolfin.VertexFunctionInt.set_all
     - |todo|
   * - dolfin.VertexFunctionInt.set_child
     - |todo|
   * - dolfin.VertexFunctionInt.set_parent
     - |todo|
   * - dolfin.VertexFunctionInt.set_value
     - |todo|
   * - dolfin.VertexFunctionInt.set_values
     - |todo|
   * - dolfin.VertexFunctionInt.size
     - |todo|
   * - dolfin.VertexFunctionInt.str
     - |todo|
   * - dolfin.VertexFunctionInt.ufl_id
     - |todo|
   * - dolfin.VertexFunctionInt.where_equal
     - |todo|
   * - dolfin.VertexFunctionSizet
     - |todo|
   * - dolfin.VertexFunctionSizet.array
     - |todo|
   * - dolfin.VertexFunctionSizet.child
     - |todo|
   * - dolfin.VertexFunctionSizet.clear_child
     - |todo|
   * - dolfin.VertexFunctionSizet.cpp_value_type
     - |todo|
   * - dolfin.VertexFunctionSizet.depth
     - |todo|
   * - dolfin.VertexFunctionSizet.dim
     - |todo|
   * - dolfin.VertexFunctionSizet.empty
     - |todo|
   * - dolfin.VertexFunctionSizet.has_child
     - |todo|
   * - dolfin.VertexFunctionSizet.has_parent
     - |todo|
   * - dolfin.VertexFunctionSizet.id
     - |todo|
   * - dolfin.VertexFunctionSizet.init
     - |todo|
   * - dolfin.VertexFunctionSizet.label
     - |todo|
   * - dolfin.VertexFunctionSizet.leaf_node
     - |todo|
   * - dolfin.VertexFunctionSizet.mesh
     - |todo|
   * - dolfin.VertexFunctionSizet.name
     - |todo|
   * - dolfin.VertexFunctionSizet.parameters
     - |todo|
   * - dolfin.VertexFunctionSizet.parent
     - |todo|
   * - dolfin.VertexFunctionSizet.rename
     - |todo|
   * - dolfin.VertexFunctionSizet.root_node
     - |todo|
   * - dolfin.VertexFunctionSizet.set_all
     - |todo|
   * - dolfin.VertexFunctionSizet.set_child
     - |todo|
   * - dolfin.VertexFunctionSizet.set_parent
     - |todo|
   * - dolfin.VertexFunctionSizet.set_value
     - |todo|
   * - dolfin.VertexFunctionSizet.set_values
     - |todo|
   * - dolfin.VertexFunctionSizet.size
     - |todo|
   * - dolfin.VertexFunctionSizet.str
     - |todo|
   * - dolfin.VertexFunctionSizet.ufl_id
     - |todo|
   * - dolfin.VertexFunctionSizet.where_equal
     - |todo|
   * - dolfin.WARNING
     - |todo|
   * - dolfin.X3DOM.build_x3dom_tree
     - |todo|
   * - dolfin.X3DOMParameters.Representation_surface
     - |todo|
   * - dolfin.X3DOMParameters.Representation_surface_with_edges
     - |todo|
   * - dolfin.X3DOMParameters.Representation_wireframe
     - |todo|
   * - dolfin.X3DOMParameters.get_ambient_intensity
     - |todo|
   * - dolfin.X3DOMParameters.get_background_color
     - |todo|
   * - dolfin.X3DOMParameters.get_color_map
     - |todo|
   * - dolfin.X3DOMParameters.get_emissive_color
     - |todo|
   * - dolfin.X3DOMParameters.get_menu_display
     - |todo|
   * - dolfin.X3DOMParameters.get_representation
     - |todo|
   * - dolfin.X3DOMParameters.get_shininess
     - |todo|
   * - dolfin.X3DOMParameters.get_specular_color
     - |todo|
   * - dolfin.X3DOMParameters.get_transparency
     - |todo|
   * - dolfin.X3DOMParameters.get_viewport_size
     - |todo|
   * - dolfin.X3DOMParameters.get_x3d_stats
     - |todo|
   * - dolfin.X3DOMParameters.set_ambient_intensity
     - |todo|
   * - dolfin.X3DOMParameters.set_background_color
     - |todo|
   * - dolfin.X3DOMParameters.set_color_map
     - |todo|
   * - dolfin.X3DOMParameters.set_emissive_color
     - |todo|
   * - dolfin.X3DOMParameters.set_menu_display
     - |todo|
   * - dolfin.X3DOMParameters.set_representation
     - |todo|
   * - dolfin.X3DOMParameters.set_shininess
     - |todo|
   * - dolfin.X3DOMParameters.set_specular_color
     - |todo|
   * - dolfin.X3DOMParameters.set_transparency
     - |todo|
   * - dolfin.X3DOMParameters.set_x3d_stats
     - |todo|
   * - dolfin.XDMFFile.Encoding
     - |new|
   * - dolfin.XDMFFile.close
     - |todo|
   * - dolfin.XDMFFile.default_encoding
     - |todo|
   * - dolfin.XDMFFile.str
     - |todo|
   * - dolfin.adapt
     - |todo|
   * - dolfin.adapt_markers
     - |todo|
   * - dolfin.add_logfile
     - |todo|
   * - dolfin.as_cell
     - |todo|
   * - dolfin.as_domain
     - |todo|
   * - dolfin.as_ufl
     - |todo|
   * - dolfin.assemble_multimesh
     - |todo|
   * - dolfin.atan_2
     - |todo|
   * - dolfin.begin
     - |todo|
   * - dolfin.block_split
     - |todo|
   * - dolfin.cell_avg
     - |todo|
   * - dolfin.cells.end
     - |todo|
   * - dolfin.cells.end_iterator
     - |todo|
   * - dolfin.cells.next
     - |todo|
   * - dolfin.cells.pos
     - |todo|
   * - dolfin.cofac
     - |todo|
   * - dolfin.compile_cpp_code
     - |new|
   * - dolfin.compile_expressions
     - |todo|
   * - dolfin.compile_extension_module
     - |todo|
   * - dolfin.cosh
     - |todo|
   * - dolfin.custom_integral_types
     - |todo|
   * - dolfin.dI
     - |todo|
   * - dolfin.dO
     - |todo|
   * - dolfin.dS_h
     - |todo|
   * - dolfin.dS_v
     - |todo|
   * - dolfin.dc
     - |todo|
   * - dolfin.debug
     - |todo|
   * - dolfin.diag
     - |todo|
   * - dolfin.diag_vector
     - |todo|
   * - dolfin.dolfin_error
     - |todo|
   * - dolfin.dolfin_terminate
     - |todo|
   * - dolfin.dolfin_version
     - |todo|
   * - dolfin.dorfler_mark
     - |todo|
   * - dolfin.ds_b
     - |todo|
   * - dolfin.ds_t
     - |todo|
   * - dolfin.ds_tb
     - |todo|
   * - dolfin.ds_v
     - |todo|
   * - dolfin.e
     - |todo|
   * - dolfin.edges.end
     - |todo|
   * - dolfin.edges.end_iterator
     - |todo|
   * - dolfin.edges.next
     - |todo|
   * - dolfin.edges.pos
     - |todo|
   * - dolfin.end
     - |todo|
   * - dolfin.energy_norm
     - |todo|
   * - dolfin.entities.end
     - |todo|
   * - dolfin.entities.end_iterator
     - |todo|
   * - dolfin.entities.next
     - |todo|
   * - dolfin.entities.pos
     - |todo|
   * - dolfin.eq
     - |todo|
   * - dolfin.error
     - |todo|
   * - dolfin.exterior_derivative
     - |todo|
   * - dolfin.faces.end
     - |todo|
   * - dolfin.faces.end_iterator
     - |todo|
   * - dolfin.faces.next
     - |todo|
   * - dolfin.faces.pos
     - |todo|
   * - dolfin.facet
     - |todo|
   * - dolfin.facet_avg
     - |todo|
   * - dolfin.facets.end
     - |todo|
   * - dolfin.facets.end_iterator
     - |todo|
   * - dolfin.facets.next
     - |todo|
   * - dolfin.facets.pos
     - |todo|
   * - dolfin.fem_solve
     - |todo|
   * - dolfin.ffc_default_parameters
     - |todo|
   * - dolfin.functional
     - |todo|
   * - dolfin.generate_error_control
     - |todo|
   * - dolfin.generate_error_control_forms
     - |todo|
   * - dolfin.get_global_parameters
     - |todo|
   * - dolfin.get_handler
     - |todo|
   * - dolfin.get_logger
     - |todo|
   * - dolfin.get_tensor_type
     - |todo|
   * - dolfin.has_cholmod
     - |todo|
   * - dolfin.has_lu_solver_method
     - |todo|
   * - dolfin.has_mpi4py
     - |new|
   * - dolfin.has_openmp
     - |todo|
   * - dolfin.has_type
     - |todo|
   * - dolfin.has_umfpack
     - |todo|
   * - dolfin.has_zlib
     - |todo|
   * - dolfin.i
     - |todo|
   * - dolfin.indices
     - |todo|
   * - dolfin.info_blue
     - |todo|
   * - dolfin.info_green
     - |todo|
   * - dolfin.info_red
     - |todo|
   * - dolfin.info_stream
     - |todo|
   * - dolfin.info_underline
     - |todo|
   * - dolfin.init
     - |todo|
   * - dolfin.integral_types
     - |todo|
   * - dolfin.inv
     - |todo|
   * - dolfin.j
     - |todo|
   * - dolfin.jit
     - |todo|
   * - dolfin.k
     - |todo|
   * - dolfin.krylov_solver_methods
     - |todo|
   * - dolfin.krylov_solver_preconditioners
     - |todo|
   * - dolfin.l
     - |todo|
   * - dolfin.la_solve
     - |todo|
   * - dolfin.linear_solver_methods
     - |todo|
   * - dolfin.list_krylov_solver_methods
     - |todo|
   * - dolfin.list_krylov_solver_preconditioners
     - |todo|
   * - dolfin.list_linear_algebra_backends
     - |todo|
   * - dolfin.list_linear_solver_methods
     - |todo|
   * - dolfin.list_lu_solver_methods
     - |todo|
   * - dolfin.log
     - |todo|
   * - dolfin.lu_solver_methods
     - |todo|
   * - dolfin.make_ufc_coordinate_mapping
     - |todo|
   * - dolfin.make_ufc_dofmap
     - |todo|
   * - dolfin.make_ufc_finite_element
     - |todo|
   * - dolfin.make_ufc_form
     - |todo|
   * - dolfin.mark
     - |todo|
   * - dolfin.max_value
     - |todo|
   * - dolfin.memory_usage
     - |todo|
   * - dolfin.min_value
     - |todo|
   * - dolfin.monitor_memory_usage
     - |todo|
   * - dolfin.nabla_div
     - |todo|
   * - dolfin.nabla_grad
     - |todo|
   * - dolfin.ne
     - |todo|
   * - dolfin.not_working_in_parallel
     - |todo|
   * - dolfin.old_init
     - |todo|
   * - dolfin.p
     - |todo|
   * - dolfin.p_refine
     - |todo|
   * - dolfin.perp
     - |todo|
   * - dolfin.product
     - |todo|
   * - dolfin.q
     - |todo|
   * - dolfin.r
     - |todo|
   * - dolfin.rand
     - |todo|
   * - dolfin.rank
     - |todo|
   * - dolfin.register_element
     - |todo|
   * - dolfin.register_integral_type
     - |todo|
   * - dolfin.relabel
     - |todo|
   * - dolfin.replace
     - |todo|
   * - dolfin.replace_integral_domains
     - |todo|
   * - dolfin.residual
     - |todo|
   * - dolfin.rot
     - |todo|
   * - dolfin.s
     - |todo|
   * - dolfin.seed
     - |todo|
   * - dolfin.sensitivity_rhs
     - |todo|
   * - dolfin.set_indentation_level
     - |todo|
   * - dolfin.set_log_active
     - |todo|
   * - dolfin.set_output_stream
     - |todo|
   * - dolfin.shape
     - |todo|
   * - dolfin.show_elements
     - |todo|
   * - dolfin.sign
     - |todo|
   * - dolfin.sinh
     - |todo|
   * - dolfin.sizeof_la_index
     - |todo|
   * - dolfin.sqr
     - |todo|
   * - dolfin.stored_dlopen_flags
     - |new|
   * - dolfin.string_types
     - |todo|
   * - dolfin.supported_elements
     - |todo|
   * - dolfin.supported_elements_for_plotting
     - |todo|
   * - dolfin.tanh
     - |todo|
   * - dolfin.tic
     - |todo|
   * - dolfin.time
     - |todo|
   * - dolfin.toc
     - |todo|
   * - dolfin.transpose
     - |todo|
   * - dolfin.ufc_cell
     - |todo|
   * - dolfin.ufc_cell.cell_shape
     - |todo|
   * - dolfin.ufc_cell.geometric_dimension
     - |todo|
   * - dolfin.ufc_cell.index
     - |todo|
   * - dolfin.ufc_cell.local_facet
     - |todo|
   * - dolfin.ufc_cell.mesh_identifier
     - |todo|
   * - dolfin.ufc_cell.topological_dimension
     - |todo|
   * - dolfin.ufc_coordinate_mapping
     - |todo|
   * - dolfin.ufc_dofmap
     - |todo|
   * - dolfin.ufc_finite_element
     - |todo|
   * - dolfin.ufc_form
     - |todo|
   * - dolfin.ufc_function
     - |todo|
   * - dolfin.ufc_shape_hexahedron
     - |todo|
   * - dolfin.ufc_shape_interval
     - |todo|
   * - dolfin.ufc_shape_quadrilateral
     - |todo|
   * - dolfin.ufc_shape_tetrahedron
     - |todo|
   * - dolfin.ufc_shape_triangle
     - |todo|
   * - dolfin.ufc_signature
     - |todo|
   * - dolfin.unit_matrices
     - |todo|
   * - dolfin.unit_matrix
     - |todo|
   * - dolfin.unit_vector
     - |todo|
   * - dolfin.unit_vectors
     - |todo|
   * - dolfin.vertex
     - |todo|
   * - dolfin.vertices.end
     - |todo|
   * - dolfin.vertices.end_iterator
     - |todo|
   * - dolfin.vertices.next
     - |todo|
   * - dolfin.vertices.pos
     - |todo|
   * - dolfin.warning
     - |todo|
   * - dolfin.zero
     - |todo|

In total there are 3521 differences

Methodology
-----------

Differences between the SWIG wrappers and pybind11 wrappers of the DOLFIN API were made using the following script, which dumps a list of dolfin namespace names for a given dolfin version. Two such lists were postprocessed to make an initial version of the above table.

.. code:: python

  import types, inspect
  import dolfin
  
  for x in sorted(dir(dolfin)):
      if x.startswith('_'):
          continue
      xo = getattr(dolfin, x)
      if isinstance(xo, types.ModuleType):
          continue
  
      print("dolfin.%s" % x)
      xo = getattr(dolfin, x)
  
      if not inspect.isclass(xo):
          continue
  
      for y in sorted(dir(xo)):
          if y.startswith('_'): continue
          print("dolfin.%s.%s" % (x, y))

The postprocessor is:

.. code:: python

  SWIG = set()
  PB11 = set()
  status = {}

  for line in open('api_swig.txt'):
      SWIG.add(line.strip())

  for line in open('api_pybind11.txt'):
      PB11.add(line.strip())

  for name in SWIG.union(PB11):
      if name in SWIG and name in PB11:
          continue
      if name.endswith('.thisown'):
          continue
      if name in SWIG:
          status[name] = '|todo|'
      if name in PB11:
          status[name] = '|new|'

  for name, stat in sorted(status.items()):
      print('   * -', name)
      print('     -', stat)

  print('\n\nIn total there are %d differences' % len(status))

