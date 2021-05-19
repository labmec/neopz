The Solver hierarchy
==========================

.. toctree::
   :maxdepth: 2
   
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

TPZMatrixSolver
---------------

The :cpp:expr:`TPZMatrixSolver` represents a solver for algebraic equation systems in which the matrix has entries with the type :cpp:type:`TVar`, where :cpp:expr:`TVar=STATE`, for real problems, and :cpp:expr:`TVar=CSTATE`, for complex problems.

The available solvers are

.. doxygenenum:: TPZMatrixSolver::MSolver
   :no-link:


Further documentation on TPZMatrixSolver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: TPZMatrixSolver
   :members:


TPZStepSolver
^^^^^^^^^^^^^
.. doxygenclass:: TPZStepSolver
   :members:
                  
TPZPardisoSolver
^^^^^^^^^^^^^^^^

The :cpp:expr:`TPZPardisoSolver` class acts as an wrapper for controlling the Intel MKL PARDISO Solver. Therefore, for one to use it, NeoPZ should have been configured using :code:`USING_MKL=ON`. While the PARDISO solver can still be used through :cpp:expr:`TPZStepSolver`, using an instance of :cpp:expr:`TPZPardisoSolver` directly provides a finer control of the PARDISO parameters.

.. note::
   The PARDISO solver is a solver for both symmetric and non-symmetric sparse matrices,
   so it should be used with either :cpp:expr:`TPZSSpStructMatrix` or :cpp:expr:`TPZSpStructMatrix`.

.. note::
   There is still plenty of tuning options to be implemented in this class. Please contact us or submit a Pull Request if you think this class could be improved.

.. doxygenclass:: TPZPardisoSolver
   :members:

TPZEigenSolver
--------------

This class defines the interface for the NeoPZ solvers for eigenvalue problems. It will soon be complemented with an appropriated :cpp:expr:`TPZAnalysis` class.

.. doxygenclass:: TPZEigenSolver
   :members:
   :membergroups: Eigen

TPZLapackEigenSolver
^^^^^^^^^^^^^^^^^^^^
This class acts as a wrapper over LAPACK calls that can be used for solving eigenvalue problems. It supports :cpp:expr:`TPZFMatrix` and :cpp:expr:`TPZSBMatrix`, therefore it can only be used with the structural matrices :cpp:expr:`TPZFStructMatrix` and :cpp:expr:`TPZSBandStructMatrix`.

.. note::
   This class is also used internally by the :cpp:expr:`TPZFMatrix` and :cpp:expr:`TPZSBMatrix` classes, thus the specific (and protected) interfaces.

.. doxygenclass:: TPZLapackEigenSolver
   :members:
   :protected-members:
   :membergroups: Eigen EigenFMatrix EigenSBMatrix
