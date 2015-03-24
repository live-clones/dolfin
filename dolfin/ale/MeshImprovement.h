// Copyright (C) 2015 Chris Richardson
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

#ifndef __MESHIMPROVEMENT_H
#define __MESHIMPROVEMENT_H

namespace dolfin
{
  class Mesh;
  template<typename T> class MeshFunction;

  class MeshImprovement
  {
  public:
    /// Return a new mesh, removing some short edges less than local length lideal
    ///
    static std::shared_ptr<Mesh> collapse(std::shared_ptr<const Mesh> mesh,
                                          const MeshFunction<double>& lideal);

    /// Return a new mesh, flipping some long edges greater than 2.8x local length lideal
    static std::shared_ptr<Mesh> flip(std::shared_ptr<const Mesh> mesh,
                                      const MeshFunction<double>& lideal);

    /// Return a new mesh, splitting some long edges greater than 2.8x local length lideal
    /// Basically this is just an interface to 'refine'
    static std::shared_ptr<const Mesh> split(std::shared_ptr<const Mesh> mesh,
                                             const MeshFunction<double>& lideal);

  };
}

#endif
