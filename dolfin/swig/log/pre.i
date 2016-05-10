/* -*- C -*- */
// Copyright (C) 2009 Johan Hake
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
//
// First added:  2009-05-10
// Last changed: 2009-09-23

//=============================================================================
// SWIG directives for the DOLFIN log kernel module (pre)
//
// The directives in this file are applied _before_ the header files of the
// modules has been loaded.
//=============================================================================

//-----------------------------------------------------------------------------
// Due to a SWIG bug when overloading a function that also use elipsis (...)
// argument in C++, we need to ignore other overloaded functions. They are
// reimplemented in log_post.i
//-----------------------------------------------------------------------------
%ignore dolfin::info(const Parameters& parameters, bool verbose=false);
%ignore dolfin::info(const Variable& variable, bool verbose=false);
%rename(_info) dolfin::info;

//-----------------------------------------------------------------------------
// Need to ignore these dues to SWIG confusion of overloaded functions
//-----------------------------------------------------------------------------
%ignore dolfin::Table::set(std::string,std::string,std::size_t);

//-----------------------------------------------------------------------------
// Ignore operators so SWIG stop complaining
//-----------------------------------------------------------------------------
%ignore dolfin::TableEntry::operator std::string;
%ignore dolfin::Progress::operator++;
%ignore dolfin::Progress::operator=;
%ignore dolfin::Table::operator=;
%ignore dolfin::TableEntry::operator=;

//-----------------------------------------------------------------------------
// Ignore DOLFIN C++ stream handling
//-----------------------------------------------------------------------------
%ignore dolfin::LogStream;
%ignore dolfin::cout;
%ignore dolfin::endl;

//-----------------------------------------------------------------------------
// Typemap for passing Python file as output stream
//-----------------------------------------------------------------------------
%{
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#include <unistd.h>
#include <stdio.h>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <dolfin/log/LogManager.h>
#endif
%}

%typemap(in) std::ostream& {
%#if _POSIX_C_SOURCE >= 1 || _XOPEN_SOURCE || _POSIX_SOURCE

  // Get FILE
  FILE* f = PyFile_AsFile($input);
  if (!f)
    SWIG_exception(SWIG_TypeError, "File object expected.");

  // Check the file is writable
  PyObject* mode = PyObject_GetAttrString($input, (char*)"mode");
#if PY_MAJOR_VERSION >= 3
  if(!mode || !PyUnicode_Check(mode) || !PyUnicode_CompareWithASCIIString(mode, "r"))
#else
  if(!mode || !PyString_Check(mode) || !strcmp(PyString_AsString(mode), "r"))
#endif
  {
    Py_XDECREF(mode);
    SWIG_exception(SWIG_TypeError, "File does not seem to be writable. Mode not set or 'r'.");
  }
  Py_DECREF(mode);

  // Flush (possibly buffered) file and build ostream from fd
  PyObject_CallMethod($input, (char*)"flush", NULL);
  auto stream = new boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_sink>
    (fileno(f), boost::iostreams::never_close_handle);
  $1 = new std::ostream(stream);

%#else
  SWIG_exception(SWIG_RuntimeError, "Not on non-POSIX, sorry...");
%#endif
}
