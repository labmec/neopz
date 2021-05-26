The Solver hierarchy
==========================

.. toctree::
   :maxdepth: 2
              
   matrixsolvers.rst
   eigensolvers.rst
   
.. contents:: Table of Contents
   :local:


TPZSolver
---------

The :cpp:expr:`TPZSolver` defines the hierarchy of solvers to be used.

.. note::
   Since the availability of a given solver might depend on the chosen matrix storage format, the choice of solver is directly connected with the choice of :doc:`../structmatrix/structoptions`. Currently this only applies for the class :cpp:expr:`TPZPardisoSolver`, which is only compatible with sparse matrix storage.

The NeoPZ solvers can be divided in two main groups: the :cpp:expr:`TPZMatrixSolver` hierarchy, for solving algebraic equation systems, and the :cpp:expr:`TPZEigenSolver` hierarchy that implements solvers for eigenvalue problems.

.. doxygenclass:: TPZSolver
   :members: