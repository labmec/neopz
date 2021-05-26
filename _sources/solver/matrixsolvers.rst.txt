TPZMatrixSolver
===============

.. toctree::
   :maxdepth: 2
   
.. contents:: Table of Contents
   :local:


The :cpp:expr:`TPZMatrixSolver` represents a solver for algebraic equation systems in which the matrix has entries with the type :cpp:type:`TVar`, where :cpp:expr:`TVar=STATE`, for real problems, and :cpp:expr:`TVar=CSTATE`, for complex problems.

The available solvers are

.. doxygenenum:: TPZMatrixSolver::MSolver
   :no-link:


Further documentation on TPZMatrixSolver
----------------------------------------

.. doxygenclass:: TPZMatrixSolver
   :members:


TPZStepSolver
-------------
.. doxygenclass:: TPZStepSolver
   :members:
                  
TPZPardisoSolver
----------------

The :cpp:expr:`TPZPardisoSolver` class acts as an wrapper for controlling the Intel MKL PARDISO Solver. Therefore, for one to use it, NeoPZ should have been configured using :code:`USING_MKL=ON`. While the PARDISO solver can still be used through :cpp:expr:`TPZStepSolver`, using an instance of :cpp:expr:`TPZPardisoSolver` directly provides a finer control of the PARDISO parameters.

.. note::
   The PARDISO solver is a solver for both symmetric and non-symmetric sparse matrices,
   so it should be used with either :cpp:expr:`TPZSSpStructMatrix` or :cpp:expr:`TPZSpStructMatrix`.

.. note::
   There is still plenty of tuning options to be implemented in this class. Please contact us or submit a Pull Request if you think this class could be improved.

.. doxygenclass:: TPZPardisoSolver
   :members: