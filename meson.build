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
  license : 'MPLv2.0',
  meson_version : '>=0.46.0',
  default_options : ['cpp_std=c++14'])
include_dir = include_directories('include')

deps = []

# For controlling the tradeoff between quality and speed for the supervariable
# hash function. It seems to typically be preferable to prefer the cheapest-
# possible hash.
if get_option('strong-hashes')
  add_global_arguments('-DQUOTIENT_STRONG_HASHES', language : 'cpp')
endif

# For controlling whether QuotientGraph includes timers for the various stages
# of minimum-degree analysis.
if get_option('enable-timers')
  add_global_arguments('-DQUOTIENT_ENABLE_TIMERS', language : 'cpp')
endif

# For controlling whether (expensive) debugging checks are performed.
if get_option('enable-debug')
  add_global_arguments('-DQUOTIENT_DEBUG', language : 'cpp')
endif

# For controlling whether (still not implemented) OpenMP support is enabled.
if get_option('enable-openmp')
  openmp_dep = dependency('openmp', version : '>=3.0')
  deps.append(openmp_dep)
elif get_option('disable-openmp')
  openmp_dep = false
else
  openmp_dep = dependency('openmp', version : '>=3.0', required : false)
  deps.append(openmp_dep)
endif

# Unit tests for the CoordinateGraph class.
coordinate_graph_test_exe = executable(
    'coordinate_graph_test',
    ['test/coordinate_graph_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir,
    dependencies : deps)
test('CoordinateGraph tests', coordinate_graph_test_exe)

# Unit tests for the QuotientGraph class.
quotient_graph_test_exe = executable(
    'quotient_graph_test',
    ['test/quotient_graph_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir,
    dependencies : deps)
test('QuotientGraph tests', quotient_graph_test_exe)

# Unit tests for the MinimumDegree function.
minimum_degree_test_exe = executable(
    'minimum_degree_test',
    ['test/minimum_degree_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir,
    dependencies : deps)
test('MinimumDegree tests', minimum_degree_test_exe)

# A utility for running the MinimumDegree analysis on a Matrix Market graph.
matrix_market_minimum_degree_exe = executable(
    'matrix_market_minimum_degree',
    ['example/matrix_market_minimum_degree.cc', 'include/quotient.hpp'],
    include_directories : include_dir,
    dependencies : deps)
