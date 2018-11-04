**quotient** is a C++14 header-only implementation of the (Approximate)
Minimum Degree reordering algorithm. Loosely speaking it is applicable for
sparse Symmetric Quasi-SemiDefinite (SQSD) Cholesky factorizations.

### Dependencies
There are no dependencies for installation of the headers.

But running the example drivers requires previous installation of the
header-only command-line processing C++14 library
[specify](https://gitlab.com/jack_poulson/specify).

The build system for the examples uses [meson](http://mesonbuild.com) and
the unit tests use [Catch2](https://github.com/catchorg/Catch2).

### Example usage

Usage through the `quotient::CoordinateGraph` class is fairly straight-forward:
```c++
#include <iostream>
#include "quotient.hpp"

[...]

// Build a symmetric input graph for the (approximate) minimum degree
// reordering.
//
// Alternatively, one could use quotient::CoordinateGraph::FromMatrixMarket
// to read the graph from a Matrix Market file (e.g., from the Davis sparse
// matrix collection). But keep in mind that one often needs to enforce
// explicit symmetry.
quotient::CoordinateGraph graph;
graph.Resize(num_vertices);
graph.ReserveEdgeAdditions(num_edges_upper_bound);
for (quotient::Int index = 0; index < num_edges_to_add; ++index) {
  graph.QueueEdgeAddition(edge_source[index], edge_target[index]);
}
graph.FlushEdgeQueues();

// Run the Approximate Minimum Degree analysis using the usual external degree
// approximation.
quotient::MinimumDegreeControl control;
control.degree_type = quotient::kAmestoyDegree;
const quotient::MinimumDegreeResult analysis = quotient::MinimumDegree(
    graph, control);
```

### Running the unit tests
Assuming that [meson](http://mesonbuild.com) and
[Catch2](https://github.com/catchorg/Catch2) are already installed:
```
mkdir build-debug/
meson build-debug
cd build-debug
ninja
ninja test
```

### License
`quotient` is distributed under the
[Mozilla Public License, v. 2.0](https://www.mozilla.org/media/MPL/2.0/index.815ca599c9df.txt).
