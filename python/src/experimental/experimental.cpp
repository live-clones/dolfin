#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


namespace py = pybind11;

namespace dolfin_wrappers
{


  class MyMatrix
  {
  public:
    MyMatrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols)
    {
      m_data = new float[rows*cols];
    }
    float *data() { return m_data; }
    size_t rows() const { return m_rows; }
    size_t cols() const { return m_cols; }
    float values(std::size_t i) const { return *(m_data + i); }
  private:
    size_t m_rows, m_cols;
    float *m_data;
  };

  std::shared_ptr<MyMatrix> mfactory()
  { return std::make_shared<MyMatrix>(12, 8); }

  void experimental(py::module& m)
  {


    py::class_<MyMatrix, std::shared_ptr<MyMatrix>>(m, "Matrix", py::buffer_protocol())
      .def_buffer([](MyMatrix &m) -> py::buffer_info
                  {
                    return py::buffer_info(
                      m.data(),                               /* Pointer to buffer */
                      sizeof(float),                          /* Size of one scalar */
                      py::format_descriptor<float>::format(), /* Python struct-style format descriptor */
                      2,                                      /* Number of dimensions */
                      { m.rows(), m.cols() },                 /* Buffer dimensions */
                      { sizeof(float) * m.rows(),             /* Strides (in bytes) for each index */
                          sizeof(float) }
                      );
                  })
      .def("values", &MyMatrix::values, "entry values")
      .def("num_rows", &MyMatrix::rows, "number of rows");


    m.def("mfactory", &mfactory, "Make matrix");
  }




}
