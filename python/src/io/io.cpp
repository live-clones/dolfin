// Copyright (C) 2017 Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

#include <memory>
#include <pybind11/pybind11.h>

#include <dolfin/io/VTKFile.h>
#include <dolfin/io/XDMFFile.h>
#include <dolfin/mesh/Mesh.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void io(py::module& m)
  {
    // dolfin::VTKFile
    py::class_<dolfin::VTKFile, std::shared_ptr<dolfin::VTKFile>>(m, "VTKFile")
      .def(py::init<std::string, std::string>())
      .def("__lshift__",  [](py::object& obj, const dolfin::Mesh& mesh) { dolfin::VTKFile *cls = obj.cast<dolfin::VTKFile*>(); *cls << mesh; })
      .def("write", [](py::object& obj, const dolfin::Mesh& mesh) { dolfin::VTKFile *cls = obj.cast<dolfin::VTKFile*>(); *cls << mesh; });

    // dolfin::XDMFFile
    py::class_<dolfin::XDMFFile, std::shared_ptr<dolfin::XDMFFile>> xdmffile(m, "XDMFFile");
    py::enum_<dolfin::XDMFFile::Encoding>(xdmffile, "Encoding")
    .value("HDF5", dolfin::XDMFFile::Encoding::HDF5)
    .value("ASCII", dolfin::XDMFFile::Encoding::ASCII);

    xdmffile
      .def(py::init<std::string>())
      .def("write", (void (dolfin::XDMFFile::*)(const dolfin::Function&, dolfin::XDMFFile::Encoding))
           &dolfin::XDMFFile::write, py::arg("u"), py::arg("encoding") = dolfin::XDMFFile::Encoding::HDF5)
      .def("write", (void (dolfin::XDMFFile::*)(const dolfin::Mesh&, dolfin::XDMFFile::Encoding))
           &dolfin::XDMFFile::write, py::arg("mesh"), py::arg("encoding") = dolfin::XDMFFile::Encoding::HDF5);

  }

}
