
import ffc
import ufl
import types
import dolfin_test.cpp as cpp

class FunctionSpace(ufl.FunctionSpace, cpp.function.FunctionSpace):

    def __init__(self, mesh, family, degree):
        """Create finite element function space."""

        print("hello boo")
        print(mesh.ufl_id())

        # Add ufl_cell function to mesh
        def ufl_cell(self):
            return ufl.Cell(self.cell_name(),
                            geometric_dimension=self.geometry().dim())
        mesh.ufl_cell = types.MethodType(ufl_cell, mesh)

        def ufl_coordinate_element(self):
            """Return the finite element of the coordinate vector field of this
            domain."""
            cell = self.ufl_cell()
            degree = self.geometry().degree()
            return ufl.VectorElement("Lagrange", cell, degree,
                                     dim=cell.geometric_dimension())
        mesh.ufl_coordinate_element = types.MethodType(ufl_coordinate_element, mesh)

        def ufl_domain(self):
            """Returns the ufl domain corresponding to the mesh."""
            # Cache object to avoid recreating it a lot
            if not hasattr(self, "_ufl_domain"):
                self._ufl_domain = ufl.Mesh(self.ufl_coordinate_element(),
                                            ufl_id=self.ufl_id(),
                                            cargo=self)
            return self._ufl_domain
        mesh.ufl_domain = types.MethodType(ufl_domain, mesh)

        # Create UFL element
        element = ufl.FiniteElement(family, mesh.ufl_cell(), degree,
                                    form_degree=None)

        # Initialize the ufl.FunctionSpace first to check for good
        # meaning
        ufl.FunctionSpace.__init__(self, mesh.ufl_domain(), element)

        ufc_element, ufc_dofmap = ffc.jit(element, parameters=None)
        print("JIT type:", type(ufc_element), ufc_element)
        ufc_element = cpp.fem.make_ufc_finite_element(ufc_element)

        dolfin_element = cpp.fem.FiniteElement(ufc_element)
        print("Element signature:", dolfin_element.signature())

        ufc_dofmap = cpp.fem.make_ufc_dofmap(ufc_dofmap)
        dolfin_dofmap  = cpp.fem.DofMap(ufc_dofmap, mesh)

        # Initialize the cpp.FunctionSpace
        cpp.function.FunctionSpace.__init__(self, mesh, dolfin_element,
                                            dolfin_dofmap)
