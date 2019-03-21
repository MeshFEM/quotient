.. catamari documentation master file, created by
   sphinx-quickstart on Mon Mar  4 10:29:06 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Quickstart
----------
Lorem ipsum.

Building the examples and tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
While `quotient <https://hodgestar.com/catamari/>`_ is a header-only library,
there are a number of configuration options which are handled through
preprocessor directives and compiler/linker options determined during a
configuration stage. Quotient uses `meson <https://mesonbuild.com>`_, a modern
alternative to `CMake <https://cmake.org/>`_, for this configuration.

One might start with a debug build (the default :samp:`buildtype` in
`meson <https://mesonbuild.com>`_). Assuming that the
`Ninja build system <https://ninja-build.org>`_ was installed alongside
meson (it typically is), one can configure and build catamari with its default
options via:

.. code-block:: bash

  mkdir build-debug/
  meson build-debug
  cd build-debug
  ninja

A release version can be built by specifying the :samp:`buildtype` option as
:samp:`release`:

.. code-block:: bash

  mkdir build-release/
  meson build-release -Dbuild-type=release
  cd build-release
  ninja

By default, quotient is configured with :samp:`quotient::Int` equal to a
64-bit signed integer. But the library can be configured with 32-bit integer
support via the :samp:`-Duse_64bit=false` option.

In any build configuration, the library's unit tests can be run via:

.. code-block:: bash

  ninja test

Using :samp:`quotient::Buffer<T>` instead of :samp:`std::vector<T>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A well-known issue with :samp:`std::vector<T>` is that it cannot be readily used
to allocate data without initializing each entry. In the case of
`catamari's <https://hodgestar.com/catamari/>`_ multithreaded sparse-direct
solver, sequential default initialization overhead was seen to be a significant
bottleneck when running on 16 cores. Due to the interdependence of quotient and
catamari, quotient performs its Minimum Degree reordering on top of a custom
:samp:`quotient::Buffer<T>` template class as an alternative buffer allocation
mechanism.

Both :samp:`std::vector<T>` and :samp:`quotient::Buffer<T>` have the same
:samp:`operator[]` entry access semantics.

The function :samp:`quotient::Buffer<T>::Resize(std::size_t)` is
an alternative to :samp:`std::vector<T>::resize(std::size_t)` which does not
default-initialize members. Likewise,
:samp:`quotient::Buffer<T>::Resize(std::size_t, const T& value)` is an
analogue for :samp:`std::vector<T>::resize(std::size_t, const T& value)`, but
it differs in that it will ensure that **all** members of the result are equal
to the specified value (not just newly allocated ones).

Lastly, the underlying data pointer can be accessed via
:samp:`quotient::Buffer<T>::Data()` instead of
:samp:`std::vector<T>::data()` (the :samp:`begin()` and :samp:`end()` member
functions exist so that range-based for loops function over
:samp:`quotient::Buffer<T>`).

A simple example combining all of these features is:

.. code-block:: cpp

  #include <iostream>
  #include "quotient.hpp"
  const std::size_t num_entries = 5;
  quotient::Buffer<float> entries;
  entries.Resize(num_entries);
  // The five entries are not yet initialized.

  // Initialize the i'th entry as i^2.
  for (std::size_t i = 0; i < num_entries; ++i) {
    entries[i] = i * i;
  }

  // Print the entries.
  std::cout << "entries: ";
  for (const float& entry : entries) { 
    std::cout << entry << " ";
  }
  std::cout << std::endl;

  // Double the length of the buffer and zero-initialize.
  entries.Resize(2 * num_entries, 0.f);

  // Extract a mutable pointer to the entries.
  float* entries_ptr = entries.Data();

Manipulating graphs with :samp:`CoordinateGraph`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The current user-level interface for manipulating graphs is via the
coordinate-format class :samp:`quotient::CoordinateGraph`. Its primary
underlying data is a lexicographically sorted
:samp:`quotient::Buffer<quotient::GraphEdge>` and an associated
:samp:`quotient::Buffer<Int>` of row offsets (which serve the
same role as in a Compressed Sparse Row (CSR) format). Thus, this storage
scheme is a superset of the CSR format that explicitly stores both row and
column indices for each entry.

The :samp:`quotient::GraphEdge` struct is a typedef of a
:samp:`std::pair<quotient::Int, quotient::Int>` containing the row and column
index.

The class is designed so that the sorting and offset computation overhead
can be amortized over batches of edge insertions and removals.

For example, the code block:

.. code-block:: cpp

  #include "quotient.hpp"
  quotient::CoordinateGraph graph;
  graph.Resize(5);
  graph.ReserveEdgeAdditions(6);
  graph.QueueEdgeAddition(3, 4);
  graph.QueueEdgeAddition(2, 3);
  graph.QueueEdgeAddition(2, 0);
  graph.QueueEdgeAddition(4, 2);
  graph.QueueEdgeAddition(4, 4);
  graph.QueueEdgeAddition(3, 2);
  graph.FlushEdgeQueues();
  const quotient::Buffer<quotient::GraphEdge>& edges = graph.Edges();

would return a reference to the underlying
:samp:`quotient::Buffer<quotient::GraphEdge>` of :samp:`graph`,
which should contain the edge sequence:

:samp:`(2, 0), (2, 3), (3, 2), (3, 4), (4, 2), (4, 4)`.

Similarly, subsequently running the code block:

.. code-block:: cpp

  graph.ReserveEdgeRemovals(2);
  graph.QueueEdgeRemoval(2, 3);
  graph.QueueEdgeRemoval(0, 4);
  graph.FlushEdgeQueues();

would modify the Buffer underlying the :samp:`edges` reference to now
contain the edge sequence:

:samp:`(2, 0), (3, 2), (3, 4), (4, 2), (4, 4)`.

(Approximate) Minimum Degree reorderings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Lorem ipsum [AmestoyEtAl-1996]_.

.. [AmestoyEtAl-1996] Patrick R. Amestoy, Timothy A. Davis, and Iain S. Duff, An Approximate Minimum Degree Reordering Algorithm, SIAM J. Matrix Analysis & Applic., 17 (4), pp. 886--905, 1996. DOI: https://doi.org/10.1137/S0895479894278952
