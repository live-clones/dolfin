
import dolfin_test.cpp as cpp
import ufl
import ffc

class Form(cpp.fem.Form):
    def __init__(self, form, function_spaces):

        form_compiler_parameters = None
        ufc_form = ffc.jit(form, form_compiler_parameters)
        print("jit returns: ", ufc_form)
        ufc_form = cpp.fem.make_ufc_form(ufc_form[0])


#        function_spaces = [func.function_space() for func
#                           in form.arguments()]

        print(function_spaces)

        # Initialize base class
        print("--init--")
        cpp.fem.Form.__init__(self, ufc_form, function_spaces)
        print("--end init--", type(form))

        if isinstance(form, ufl.Form):
            print("* Have a UFL form")
        else:
            print("* Do not have a UFL form")

        print("Handle coefficients")
        original_coefficients = form.coefficients()
        self.coefficients = []
        for i in range(self.num_coefficients()):
            j = self.original_coefficient_position(i)
            print("**** Apend")
            self.coefficients.append(original_coefficients[j])

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
                print("XXXboo")
                self.set_coefficient(i, self.coefficients[i])
