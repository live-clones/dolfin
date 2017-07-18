import ffc
import ufl
import types
import dolfin.cpp as cpp


def _add_ufl_mesh_data(mesh):
    """Add UFL domain data to DOLFIN mesh"""

    def ufl_cell(self):
        return ufl.Cell(self.cell_name(),
                        geometric_dimension=self.geometry().dim())
    mesh.ufl_cell = types.MethodType(ufl_cell, mesh)

    def ufl_coordinate_element(self):
        """Return the finite element of the coordinate vector field of this
            domain.

        """
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


class FunctionSpace(ufl.FunctionSpace, cpp.function.FunctionSpace):

    def __init__(self,  *args, **kwargs):
        """Create finite element function space."""

        if len(args) == 1:
            # Do we relly want to do it this way? Can we get the
            # sub-element from UFL?
            self._init_from_cpp(*args, **kwargs)
        else:
            if len(args) == 0 or not isinstance(args[0], cpp.mesh.Mesh):
                #cpp.dolfin_error("functionspace.py",
                #                 "create function space",
                #                 "Illegal argument, not a mesh: "
                #                 + str(args[0]))
                pass
            elif len(args) == 2:
                self._init_from_ufl(*args, **kwargs)
            else:
                self._init_convenience(*args, **kwargs)

    def _init_from_ufl(self, mesh, element, constrained_domain=None):

        _add_ufl_mesh_data(mesh)

        # Initialize the ufl.FunctionSpace first to check for good
        # meaning
        ufl.FunctionSpace.__init__(self, mesh.ufl_domain(), element)

        # Compile dofmap and element
        ufc_element, ufc_dofmap = ffc.jit(element, parameters=None)
        ufc_element = cpp.fem.make_ufc_finite_element(ufc_element)

        # Create DOLFIN element and dofpa
        dolfin_element = cpp.fem.FiniteElement(ufc_element)
        ufc_dofmap = cpp.fem.make_ufc_dofmap(ufc_dofmap)
        dolfin_dofmap  = cpp.fem.DofMap(ufc_dofmap, mesh)

        # Initialize the cpp.FunctionSpace
        cpp.function.FunctionSpace.__init__(self, mesh, dolfin_element,
                                            dolfin_dofmap)

    def _init_from_cpp(self, cppV, **kwargs):
        """
        if not isinstance(cppV, cpp.FunctionSpace):
            cpp.dolfin_error("functionspace.py",
                             "create function space",
                             "Illegal argument for C++ function space, "
                             "not a cpp.FunctionSpace: " + str(cppV))
        # We don't want to support copy construction. This would
        # indicate internal defficiency in the library
        if isinstance(cppV, FunctionSpace):
            cpp.dolfin_error("functionspace.py",
                             "create function space",
                             "Illegal argument for C++ function space, "
                             "should not be functions.functionspace.FunctionSpace: " + str(cppV))
        if len(kwargs) > 0:
            cpp.dolfin_error("functionspace.py",
                             "create function space",
                             "Illegal arguments, did not expect C++ "
                             "function space and **kwargs: " + str(kwargs))
        """
        # Assign all the members (including 'this' pointer to SWIG wraper)
        # NOTE: This in fact performs assignment of C++ context
        #self.__dict__ = cppV.__dict__


        # Reconstruct UFL element from signature
        #ufl_element = eval(self.element().signature(), ufl.__dict__)
        ufl_element = eval(cppV.element().signature(), ufl.__dict__)

        # Get mesh
        ufl_domain = cppV.mesh().ufl_domain()

        # Initialize the ufl.FunctionSpace (not calling cpp.Function.__init__)
        cpp.function.FunctionSpace.__init__(self, cppV)

        # Initialize the ufl.FunctionSpace (not calling cpp.Function.__init__)
        ufl.FunctionSpace.__init__(self, ufl_domain, ufl_element)

    def _init_convenience(self, mesh, family, degree, form_degree=None,
                          constrained_domain=None, restriction=None):

        _add_ufl_mesh_data(mesh)

        # Create UFL element
        element = ufl.FiniteElement(family, mesh.ufl_cell(), degree,
                                    form_degree=form_degree)

        self._init_from_ufl(mesh, element,
                            constrained_domain=constrained_domain)

    def dolfin_element(self):
        "Return the DOLFIN element."
        return self.element()

    def num_sub_spaces(self):
        "Return the number of sub spaces"
        return self.dolfin_element().num_sub_elements()

    def sub(self, i):
        "Return the i-th sub space"
        # FIXME: Should we have a more extensive check other than
        # whats includeding the cpp code?
        if not isinstance(i, int):
            raise TypeError("expected an int for 'i'")
        if self.num_sub_spaces() == 1:
            raise ValueError("no SubSpaces to extract")
        if i >= self.num_sub_spaces():
            raise ValueError("Can only extract SubSpaces with i = 0 ... %d" % \
                  (self.num_sub_spaces() - 1))
        assert hasattr(self.ufl_element(), "sub_elements")

        # Extend with the python layer
        return FunctionSpace(cpp.function.FunctionSpace.sub(self, i))

    def ufl_function_space(self):
        return self


class VectorFunctionSpace(ufl.FunctionSpace, cpp.function.FunctionSpace):

    def __init__(self, mesh, family, degree, **kwargs):
        """Create finite element function space."""

        print("hello boo")
        print(mesh.ufl_id())

        dim = kwargs.pop("dim", None)

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
        element = ufl.VectorElement(family, mesh.ufl_cell(), degree,
                                    form_degree=None, dim=dim)

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

    def ufl_function_space(self):
        return self
