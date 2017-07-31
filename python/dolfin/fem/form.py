import dolfin.cpp as cpp
import ufl
import ffc

class Form(cpp.fem.Form):
    def __init__(self, form, **kwargs):

        # Check form argument
        if not isinstance(form, ufl.Form):
            raise RuntimeError("Expected a ufl.Form.")

        form_compiler_parameters = kwargs.pop("form_compiler_parameters", None)

        ufc_form = ffc.jit(form, form_compiler_parameters)
        ufc_form = cpp.fem.make_ufc_form(ufc_form[0])

        function_spaces = [func.function_space() for func in form.arguments()]
        cpp.fem.Form.__init__(self, ufc_form, function_spaces)

        original_coefficients = form.coefficients()
        self.coefficients = []
        for i in range(self.num_coefficients()):
            j = self.original_coefficient_position(i)
            self.coefficients.append(original_coefficients[j].cpp_object())

        # Type checking coefficients
        if not all(isinstance(c, (cpp.function.GenericFunction, cpp.function.MultiMeshFunction))
                   for c in self.coefficients):
            # Developer note:
            # The form accepts a MultiMeshFunction but does not set the
            # correct coefficients. This is done in assemble_multimesh
            # at the moment
            coefficient_error = "Error while extracting coefficients. "
            raise TypeError(coefficient_error +
                            "Either provide a dict of cpp.function.GenericFunctions, " +
                            "or use Function to define your form.")

        for i in range(self.num_coefficients()):
            if isinstance(self.coefficients[i], cpp.function.GenericFunction):
                self.set_coefficient(i, self.coefficients[i])
