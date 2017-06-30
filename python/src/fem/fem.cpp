// Copyright (C) 2017 Chris Richardson and Garth N. Wells
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

#include <iostream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <ufc/ufc.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void fem(py::module& m)
  {
    //-----------------------------------------------------------------------------
    // dolfin::FiniteElement class
    py::class_<dolfin::FiniteElement, std::shared_ptr<dolfin::FiniteElement>>(m, "FiniteElement", "DOLFIN FiniteElement object")
      .def(py::init<std::shared_ptr<const ufc::finite_element>>());

    //-----------------------------------------------------------------------------
    // dolfin::DofMap class

  }

}
