.. Quotient documentation master file, created by
   sphinx-quickstart on Thu Mar 21 01:50:41 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Quotient's development documentation
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart

`quotient <https://hodgestar.com/quotient/>`_ is a `C++14 <https://en.wikipedia.org/wiki/C%2B%2B14>`_, header-only implementations of the (Approximate)
`Minimum Degree reordering algorithm <https://en.wikipedia.org/wiki/Minimum_degree_algorithm>`_ (and an alternative to :samp:`std::vector` which avoids default
initialization). Loosely speaking, it is applicable for factoring sparse
`Symmetric Quasi-SemiDefinite <https://epubs.siam.org/doi/10.1137/S0895479897329400>`_ (SQSD) matrices.
