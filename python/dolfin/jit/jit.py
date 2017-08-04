# -*- coding: utf-8 -*-
import hashlib
import dijitso
import dolfin.cpp as cpp

def compile_class(cpp_data):
    """Compile a user C(++) string or set of statements to a Python object

    cpp_data is a dict

"""
    import pkgconfig
    if not pkgconfig.exists('dolfin'):
        raise RuntimeError("Could not find DOLFIN pkg-config file. Please make sure appropriate paths are set.")

    # Get pkg-config data
    d = pkgconfig.parse('dolfin')

    # Set compiler/build options
    params = dijitso.params.default_params()
    params['build']['include_dirs'] = d["include_dirs"]
    params['build']['libs'] = d["libraries"]
    params['build']['lib_dirs'] = d["library_dirs"]

#    if not isinstance(statements, (string_types, tuple, list)):
#        raise RuntimeError("Expression must be a string, or a list or tuple of strings")

    name = cpp_data['name']
    if name not in ('subdomain', 'expression'):
        raise ValueError("DOLFIN JIT only for SubDomain and Expression")
    statements = cpp_data['statements']
    properties = cpp_data['properties']

    hash_str = str(statements) + str(properties.keys())
    module_hash = hashlib.md5(hash_str.encode('utf-8')).hexdigest()
    module_name = "dolfin_" + name + "_" + module_hash

    try:
        module, signature = dijitso.jit(cpp_data, module_name, params,
                                        generate=cpp_data['jit_generate'])
        submodule = dijitso.extract_factory_function(module, "create_" + module_name)()
    except:
        raise RuntimeError("Unable to compile C++ code with dijitso")

    if name == 'expression':
        python_object = cpp.function.make_dolfin_expression(submodule)
    else:
        python_object = cpp.mesh.make_dolfin_subdomain(submodule)

    # Set properties to initial values (if any)
    for k in properties:
        python_object.set_property(k, properties[k])

    return python_object
