
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
        cpp.fem.Form.__init__(self, ufc_form,
                              function_spaces)
