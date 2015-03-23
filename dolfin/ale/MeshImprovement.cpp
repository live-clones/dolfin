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
//

#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include <dolfin/common/MPI.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/refinement/refine.h>

#include "MeshImprovement.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Mesh MeshImprovement::collapse(const Mesh& mesh, const MeshFunction<double>& lideal)
{
  // L(ideal) should be defined on vertices
  dolfin_assert(lideal.dim() == 0);

  std::size_t tdim = mesh.topology().dim();
  std::size_t gdim = mesh.geometry().dim();

  if (tdim != gdim)
  {
    dolfin_error("MeshImprovement.cpp",
                 "collapse edges on manifold",
                 "Not yet implemented");
  }

  if (MPI::size(mesh.mpi_comm()) > 1)
  {
    dolfin_error("MeshImprovement.cpp",
                 "collapse in parallel",
                 "Not yet implemented");
  }

  // Fix on 2D for now
  dolfin_assert(tdim == 2);

  // Copy topological data
  boost::multi_array<std::size_t, 2> topo(boost::extents[mesh.num_cells()][tdim + 1]);
  std::copy(mesh.cells().begin(), mesh.cells().end(), topo.data());
  // Copy geometry data
  boost::multi_array<double, 2> geom(boost::extents[mesh.num_vertices()][gdim]);
  std::copy(mesh.coordinates().begin(), mesh.coordinates().end(), geom.data());

  mesh.init(1, tdim);

  std::vector<std::pair<double, std::size_t> > candidates;
  std::set<std::size_t> v_ext;

  // Go through all edges, looking for any which meet criteria
  // Also note external edges in v_ext set
  for (EdgeIterator e(mesh); !e.end(); ++e)
  {
    const std::size_t v0 = e->entities(0)[0];
    const std::size_t v1 = e->entities(0)[1];
    const double l_local = 0.5*(lideal[v0] + lideal[v1]);
    const std::pair<double, std::size_t> elen(e->length()/l_local, e->index());
    const std::size_t num_cells = e->num_entities(tdim);
    if (elen.first < 1.0 && num_cells == 2)
      candidates.push_back(elen);
    if (num_cells == 1)
    {
      v_ext.insert(v0);
      v_ext.insert(v1);
    }
  }
  // Sort into ascending order
  std::sort(candidates.begin(), candidates.end());

  mesh.init(0, tdim);
  std::vector<bool> can_collapse(mesh.num_cells(), true);
  std::vector<std::size_t> kill_cells;
  std::vector<std::size_t> kill_verts;

  std::size_t collapse_count = 0;
  for (auto it : candidates)
  {
    std::set<std::size_t> cell_block;
    Edge e(mesh, it.second);

    std::size_t kvert = e.entities(0)[0];
    std::size_t svert = e.entities(0)[1];

    // If vertex to be deleted is external, swap
    std::size_t ext_count = v_ext.count(kvert);
    if (ext_count == 1)
    {
      ext_count += v_ext.count(svert);
      std::swap(kvert, svert);
    }
    else
      ext_count += v_ext.count(svert);

    // If both vertices are external, cannot collapse edge
    if (ext_count == 2)
      continue;

    // Get all affected cells
    Vertex v0(mesh, kvert);
    cell_block.insert(v0.entities(tdim),
                      v0.entities(tdim) + v0.num_entities(tdim));
    Vertex v1(mesh, svert);
    cell_block.insert(v1.entities(tdim),
                      v1.entities(tdim) + v1.num_entities(tdim));

    // Check all cells are not already involved in a collapse
    bool go_collapse = true;
    for (auto idx : cell_block)
      go_collapse &= can_collapse[idx];

    if (go_collapse)
    {
      collapse_count++;
      for (auto idx : cell_block)
        can_collapse[idx] = false;

      for (CellIterator c(e); !c.end(); ++c)
        kill_cells.push_back(c->index());

      kill_verts.push_back(kvert);

      // If neither vertex is external, choose midpoint
      if (ext_count == 0)
      {
        double *midpt = e.midpoint().coordinates();
        std::copy(midpt, midpt + gdim, geom[svert].begin());
      }

      // Update all surrounding cells
      for (CellIterator c(Vertex(mesh, kvert)); !c.end(); ++c)
        std::replace(topo[c->index()].begin(), topo[c->index()].end(), kvert, svert);
    }
  }

  std::cout << "Collapse count = " << collapse_count << std::endl;

  if (collapse_count == 0)
    return mesh;

  // Generate updated mesh, removing killed cells and vertices

  Mesh mesh2;
  MeshEditor ed;
  ed.open(mesh2, tdim, gdim);

  // Vertices
  const std::size_t nverts = mesh.num_vertices() - kill_verts.size();
  std::sort(kill_verts.begin(), kill_verts.end());

  ed.init_vertices(nverts);

  std::size_t kidx = 0;
  const std::size_t kimax = kill_verts.size();
  std::vector<int> new_map(mesh.num_vertices(), -1);
  for (std::size_t i = 0; i != mesh.num_vertices(); ++i)
  {
    if (kidx < kimax and kill_verts[kidx] == i)
      ++kidx;
    else
    {
      std::size_t idx = i - kidx;
      new_map[i] = idx;
      ed.add_vertex(idx, geom[i][0], geom[i][1]);
    }
  }

  // Cells
  std::size_t ncells = mesh.num_cells() - kill_cells.size();
  std::sort(kill_cells.begin(), kill_cells.end());

  ed.init_cells(ncells);

  std::size_t cidx = 0;
  const std::size_t cimax = kill_cells.size();
  for (std::size_t i = 0; i != mesh.num_cells(); ++i)
  {
    if (cidx < cimax and kill_cells[cidx] == i)
      ++cidx;
    else
    {
      std::size_t idx = i - cidx;
      ed.add_cell(idx, new_map[topo[i][0]], new_map[topo[i][1]], new_map[topo[i][2]]);
    }
  }

  ed.close();

  return mesh2;
}

//-----------------------------------------------------------------------------
Eigen::Vector2d calc_crossover(const Mesh& mesh, std::vector<std::size_t> vidx)
{
  Eigen::Vector2d a, b, c, d;
  Point p = Vertex(mesh, vidx[0]).point();
  a << p.x(), p.y();
  p = Vertex(mesh, vidx[1]).point();
  b << p.x(), p.y();
  p = Vertex(mesh, vidx[2]).point();
  c << p.x(), p.y();
  p = Vertex(mesh, vidx[3]).point();
  d << p.x(), p.y();

  Eigen::Matrix2d A;
  A.col(0) = a - b;
  A.col(1) = c - d;

  Eigen::Vector2d rhs = c - b;

  // FIXME - need to catch singular A
  Eigen::Vector2d x = A.colPivHouseholderQr().solve(rhs);

  return x;
}
//-----------------------------------------------------------------------------
Mesh MeshImprovement::flip(const Mesh& mesh, const MeshFunction<double>& lideal)
{
  const std::size_t tdim = mesh.topology().dim();
  const std::size_t gdim = mesh.geometry().dim();

  if (tdim != gdim)
  {
    dolfin_error("MeshImprovement.cpp",
                 "flip edges on manifold",
                 "Not implemented");
  }

  // Fix on 2D for now
  dolfin_assert(tdim == 2);

  // Copy topological data
  boost::multi_array<std::size_t, 2> topo(boost::extents[mesh.num_cells()][tdim + 1]);
  std::copy(mesh.cells().begin(), mesh.cells().end(), topo.data());
  // Make ref to geometry data
  boost::const_multi_array_ref<double, 2> geom(mesh.coordinates().data(),
                                               boost::extents[mesh.num_vertices()][gdim]);

  mesh.init(1, tdim);

  std::vector<std::pair<double, std::size_t> > candidates;
  for (EdgeIterator e(mesh); !e.end(); ++e)
  {
    const std::size_t v0 = e->entities(0)[0];
    const std::size_t v1 = e->entities(0)[1];
    const std::size_t num_cells = e->num_entities(tdim);
    const double l_local = 2.8*0.5*(lideal[v0] + lideal[v1]);
    const std::pair<double, std::size_t> elen(l_local/e->length(), e->index());
    if (elen.first < 1.0 && num_cells == 2)
      candidates.push_back(elen);
  }
  // Sort into ascending order
  std::sort(candidates.begin(), candidates.end());

  // Go through and try to flip
  CellFunction<bool> can_flip(mesh, true);
  std::size_t flip_count = 0;
  std::vector<std::size_t> vidx;
  for (auto it : candidates)
  {
    const Edge e(mesh, it.second);

    // Check if cells have already been flipped
    bool go_flip = true;
    for (CellIterator c(e); !c.end(); ++c)
      go_flip &= can_flip[*c];

    if (go_flip)
    {
      vidx.clear();
      vidx.push_back(e.entities(0)[0]);
      vidx.push_back(e.entities(0)[1]);

      // Get indices of 'off-edge' points
      for (CellIterator c(e); !c.end(); ++c)
      {
        for (VertexIterator v(*c); !v.end(); ++v)
          if (v->index() != vidx[0] and
              v->index() != vidx[1])
          {
            vidx.push_back(v->index());
          }
      }

      Point d = Vertex(mesh, vidx[2]).point()
              - Vertex(mesh, vidx[3]).point();

      // Check flip is advantageous
      if (d.norm() < e.length())
      {
        // Final check that crossover is inside cell, and
        // close to centre of quadrilateral
        Eigen::Vector2d x = calc_crossover(mesh, vidx);

        // Set tolerance (in range 0.0 - 0.5)
        const double al = 0.2;
        if (x[0] > al and x[0] < (1.0 - al) and x[1] > 0.0 and x[1] < 1.0)
        {
          ++flip_count;
          const std::size_t c0 = e.entities(tdim)[0];
          const std::size_t c1 = e.entities(tdim)[1];
          topo[c0][0] = vidx[2];
          topo[c0][1] = vidx[3];
          topo[c0][2] = vidx[0];
          topo[c1][0] = vidx[3];
          topo[c1][1] = vidx[2];
          topo[c1][2] = vidx[1];

          // Prevent cells from being touched again
          // FIXME: do we need this?
          for (CellIterator c(e); !c.end(); ++c)
            can_flip[*c] = false;
        }
      }
    }
  }

  std::cout << "Flip count = " << flip_count << std::endl;

  if (flip_count == 0)
    return mesh;

  Mesh mesh2;
  MeshEditor ed;
  ed.open(mesh2, tdim, gdim);
  ed.init_vertices(mesh.num_vertices());
  for (unsigned int i = 0; i != mesh.num_vertices(); ++i)
    ed.add_vertex(i, geom[i][0], geom[i][1]);

  ed.init_cells(mesh.num_cells());
  for (unsigned int i = 0; i != mesh.num_cells(); ++i)
    ed.add_cell(i, topo[i][0], topo[i][1], topo[i][2]);
  ed.close();

  return mesh2;
}
//-----------------------------------------------------------------------------
Mesh MeshImprovement::split(const Mesh& mesh, const MeshFunction<double>& lideal)
{
  const std::size_t tdim = mesh.topology().dim();

  dolfin_assert(tdim == 2);

  mesh.init(1, tdim);
  EdgeFunction<bool> split(mesh, false);

  std::size_t sct = 0;
  for (EdgeIterator e(mesh); !e.end(); ++e)
  {
    const std::size_t v0 = e->entities(0)[0];
    const std::size_t v1 = e->entities(0)[1];
    double l_local = 0.5*(lideal[v0] + lideal[v1]);
    if (e->length() > 2.8*l_local)
    {
      split[*e] = true;
      ++sct;
    }
  }

  if (sct == 0)
    return mesh;

  // Marker refine
  Mesh mesh2;
  refine(mesh2, mesh, split, false);

  return mesh2;
}
//-----------------------------------------------------------------------------
