#
# Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
project(
  'quotient',
  'cpp',
  version : '0.1',
  license: 'MPLv2.0',
  default_options : ['cpp_std=c++14'])
include_dir = include_directories('include')

cpp_compiler = meson.get_compiler('cpp')
libgomp_cpp = cpp_compiler.find_library('gomp')

add_global_arguments('-DQUOTIENT_DEBUG', language : 'cpp')

# Unit tests for the CoordinateGraph class.
coordinate_graph_test_exe = executable(
    'coordinate_graph_test',
    ['test/coordinate_graph_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir)
test('CoordinateGraph tests', coordinate_graph_test_exe)

# Unit tests for the QuotientGraph class.
quotient_graph_test_exe = executable(
    'quotient_graph_test',
    ['test/quotient_graph_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir)
test('QuotientGraph tests', quotient_graph_test_exe)

# Unit tests for the RandomAccessHeap class.
random_access_heap_test_exe = executable(
    'random_access_heap_test',
    ['test/random_access_heap_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir)
test('RandomAccessHeap tests', random_access_heap_test_exe)

# Unit tests for the MinimumDegree function.
minimum_degree_test_exe = executable(
    'minimum_degree_test',
    ['test/minimum_degree_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir)
test('MinimumDegree tests', minimum_degree_test_exe)

# A utility for running the MinimumDegree analysis on a Matrix Market graph.
matrix_market_minimum_degree_exe = executable(
    'matrix_market_minimum_degree',
    ['example/matrix_market_minimum_degree.cc', 'include/quotient.hpp'],
    include_directories : include_dir)
