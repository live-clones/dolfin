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
