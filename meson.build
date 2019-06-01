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
cxx = meson.get_compiler('cpp')

deps = []
test_deps = []
cxx_args = []

# Add Catch
test_deps += subproject('catch2').get_variable('catch2_dep')

# Test if mantis is installed.
found_mantis = cxx.has_header('mantis.hpp')
if not found_mantis
  mantis_proj = subproject('mantis')
  mantis_dep = mantis_proj.get_variable('mantis_dep')
  deps += mantis_dep
endif

# Test if specify is installed.
found_specify = cxx.has_header('specify.hpp')
if not found_specify
  specify_proj = subproject('specify')
  specify_dep = specify_proj.get_variable('specify_dep')
  deps += specify_dep
endif

cpu_family = build_machine.cpu_family()
if cpu_family == 'x86' or cpu_family == 'x86_64'
  cxx_args += '-DX86'
endif

system = build_machine.system()
if system == 'linux'
  cxx_args += '-DLINUX'
endif

# Handle fast (Fortran-like) std::complex arithmetic.
if cxx.get_id() == 'gcc'
  message('Using gcc')
  # Complex arithmetic is greatly accelerated by adopting the Fortran
  # division and NaN-checking conventions.
  cxx_args += '-fcx-fortran-rules'
endif

# Handle the Fused Multiply Add/Subtraction approaches.
have_fma_intrinsic = false
if not have_fma_intrinsic and cxx.get_id() == 'gcc' and cxx.version() >= '4.3'
  if not get_option('disable_gcc_fma')
    have_fma_intrinsic = true
    cxx_args += '-DMANTIS_GCC_FUSED_MULTIPLY_ADD'
  else
    message('Avoiding GCC __builtin_fma due to disable_gcc_fma option.')
  endif
endif
if not have_fma_intrinsic and cxx.get_id() == 'xlc'
  if not get_option('disable_xlc_fma')
    have_fma_intrinsic = true
    cxx_args += '-DMANTIS_XLC_FUSED_MULTIPLY_ADD'
  else
    message('Avoiding XLC __fmadd due to disable_xlc_fma option.')
  endif
endif

if get_option('ieee_sum')
  # Use the more accurate, but slower, summation mechanism.
  cxx_args += '-DMANTIS_IEEE_SUM'
endif

# For controlling whether 32-bit or 64-bit integers are used. By default,
# 64-bit integers are used.
if get_option('use_64bit')
  cxx_args += '-DQUOTIENT_USE_64BIT'
endif

# For controlling the tradeoff between quality and speed for the supervariable
# hash function. It seems to typically be preferable to prefer the cheapest-
# possible hash.
if get_option('strong_hashes')
  cxx_args += '-DQUOTIENT_STRONG_HASHES'
endif

# For controlling whether QuotientGraph includes timers for the various stages
# of minimum-degree analysis.
if get_option('enable_timers')
  cxx_args += '-DQUOTIENT_ENABLE_TIMERS'
endif

# For controlling whether (expensive) debugging checks are performed.
if get_option('enable_debug')
  cxx_args += '-DQUOTIENT_DEBUG'
endif

# For controlling whether (still not implemented) OpenMP support is enabled.
trust_omp_version = true
# Older versions of clang, e.g., 3.8.0-2ubuntu4, falsely report support for
# OpenMP 4.0 despite numerous issues (e.g., with task data sharing attributes).
# As of clang 8.0, only OpenMP 3.1 is supported. We therefore block OpenMP
# if the version of clang is less than 9.0.
if cxx.get_id() == 'clang' and cxx.version() < '9.0'
  message('Ignoring Clang\'s reported OpenMP version')
  trust_omp_version = false
endif
# Older versions of gcc, e.g., 5.4.0-6ubuntu1~16.04.11, falsely report support
# for OpenMP 4.0 despite bugs in task dependencies. However, version
# 5.5.0-12ubuntu1 appears to have fixed this issue.
if cxx.get_id() == 'gcc' and cxx.version() < '5.5'
  message('Ignoring GCC\'s reported OpenMP version')
  trust_omp_version = false
endif
if trust_omp_version
  have_openmp = false
  if get_option('enable_openmp')
    openmp_dep = dependency('openmp', version : '>=4.0')
    have_openmp = true
    deps += openmp_dep
  elif not get_option('disable_openmp')
    openmp_dep = dependency('openmp', version : '>=4.0', required : false)
    if openmp_dep.found()
      have_openmp = true
      deps += openmp_dep
    endif
  endif
  if have_openmp
    cxx_args += '-DQUOTIENT_OPENMP'
  endif
endif

# Unit tests for the CoordinateGraph class.
coordinate_graph_test_exe = executable(
    'coordinate_graph_test',
    ['test/coordinate_graph_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir,
    dependencies : deps + test_deps,
    cpp_args : cxx_args)
test('CoordinateGraph tests', coordinate_graph_test_exe)

# Unit tests for the QuotientGraph class.
quotient_graph_test_exe = executable(
    'quotient_graph_test',
    ['test/quotient_graph_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir,
    dependencies : deps + test_deps,
    cpp_args : cxx_args)
test('QuotientGraph tests', quotient_graph_test_exe)

# Unit tests for the MinimumDegree function.
minimum_degree_test_exe = executable(
    'minimum_degree_test',
    ['test/minimum_degree_test.cc', 'include/quotient.hpp'],
    include_directories : include_dir,
    dependencies : deps + test_deps,
    cpp_args : cxx_args)
test('MinimumDegree tests', minimum_degree_test_exe)

# A utility for running the MinimumDegree analysis on a Matrix Market graph.
matrix_market_minimum_degree_exe = executable(
    'matrix_market_minimum_degree',
    ['example/matrix_market_minimum_degree.cc', 'include/quotient.hpp'],
    include_directories : include_dir,
    dependencies : deps,
    cpp_args : cxx_args)

# For using quotient as a subproject.
quotient_dep = declare_dependency(include_directories : include_dir)

# Install include/quotient.hpp and include/quotient/
install_subdir('include', install_dir : '.')
