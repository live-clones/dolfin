import cppimport

cpp_code="""
/*
<%
setup_pybind11(cfg)
%>
*/
#include <pybind11/pybind11.h>

namespace py = pybind11;

int square(int x) {
    return x * x;
}

PYBIND11_PLUGIN(testcode) {
    pybind11::module m("testcode", "auto-compiled c++ extension");
    m.def("square", &square);
    return m.ptr();
}
"""

f = open('testcode.cpp', 'w')
f.write(cpp_code)
f.close()

mycode = cppimport.imp("testcode")
print(mycode.square(2))


cpp_code2="""
/*
<%
setup_pybind11(cfg)
%>
*/
#include <pybind11/pybind11.h>

namespace py = pybind11;

class MyExpression
{
public:

  MyExpression(std::size_t i) : _i(i) {}

  std::size_t get() {return _i; }

private:

  std::size_t _i;

};

int square(int x)
{
    return x * x;
}

PYBIND11_PLUGIN(testcode2)
{
  pybind11::module m("testcode", "auto-compiled c++ extension");

  py::class_<MyExpression>(m, "MyExpression")
        .def(py::init<std::size_t>())
        .def("get", &MyExpression::get);

  m.def("square", &square);

  return m.ptr();
}
"""

f = open('testcode2.cpp', 'w')
f.write(cpp_code2)
f.close()

mycode2 = cppimport.imp("testcode2")

my_exp = mycode2.MyExpression(6)
g = my_exp.get()
print(g)



cpp_code3="""
/*
<%
import os
import ffc
setup_pybind11(cfg)
cfg['include_dirs'] = ['../', ffc.get_include_path()]
cfg['libraries'] = ['dolfin']
lib_dirs = []
variables = ('LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH')
for e in variables:
    if e in os.environ:
        lib_dirs += os.environ[e].split(":")
lib_dirs = [x for x in lib_dirs if x != ""]
cfg['library_dirs'] = lib_dirs
%>
*/
#include <pybind11/pybind11.h>
#include <string>
#include <dolfin/function/Expression.h>

namespace py = pybind11;

class TestExpression : public dolfin::Expression
{
public:

  TestExpression(std::string name) : _name(name) {}

  std::string get() {return _name; }

private:

  std::string _name;

};

PYBIND11_PLUGIN(testcode3)
{
  //py::module::import("dolfin_test");
  //py::module::import("cpp");
  py::module::import("cpp.function");
  //expy::module::import("cpp.function.Expression");

    py::module m("example", "pybind11 example plugin");


  py::class_<dolfin::Expression> dolfin_expression(m, "cpp.function.Expression");

  py::class_<TestExpression>(m, "TestExpression", dolfin_expression)
    .def(py::init<std::string>())
    .def("get", &TestExpression::get);

  return m.ptr();
}
"""
cppimport.set_quiet(False)

f = open('testcode3.cpp', 'w')
f.write(cpp_code3)
f.close()


import cpp.function

mycode3 = cppimport.imp("testcode3")


print(type(mycode3))
print(dir(mycode3))


my_exp = mycode3.TestExpression("testing")
print(type(my_exp))
print(dir(my_exp))

#c =  my_exp.cpp.function.Expression()

help(my_exp)
#g = my_exp.get()
#print(g)
