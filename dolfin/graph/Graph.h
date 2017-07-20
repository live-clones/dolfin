// Copyright (C) 2010 Garth N. Wells
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
// First added:  2010-11-16
// Last changed: 2010-11-25

#ifndef __GRAPH_TYPES_H
#define __GRAPH_TYPES_H

#include <vector>

#define BOOST_NO_HASH

#include <boost/graph/adjacency_list.hpp>
#include <dolfin/common/Set.h>

namespace dolfin
{

  /// Simple graph data structure

  // Edge storage type
  typedef dolfin::Set<int> graph_set_type;

  class Graph
  {

  public:

    /// Create an empty graph
    Graph() {}

    /// Create a graph with a specified number of nodes
    explicit Graph(std::size_t size) : _graph(size) {}

    std::size_t size() const
    { return _graph.size(); }

    graph_set_type& operator[] (std::size_t i)
    { return _graph[i]; }

    std::vector<graph_set_type>::iterator begin()
    { return _graph.begin(); }
    std::vector<graph_set_type>::iterator end()
    { return _graph.end(); }

    std::vector<graph_set_type>::const_iterator begin() const
    { return _graph.begin(); }
    std::vector<graph_set_type>::const_iterator end() const
    { return _graph.end(); }

    void push_back(const graph_set_type& edges)
    { _graph.push_back(edges); }

  private:

    // Vector of unordered Sets
    std::vector<graph_set_type> _graph;

  };

}

#endif
