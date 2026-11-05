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

.. _section-pardiso-advanced:
   
Advanced settings
+++++++++++++++++

PARDISO calls are set up primarily through the :code:`iparm` array (more info `here <https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface/pardiso-iparm-parameter.html>`_).
We now provide means to access this array and configure PARDISO to match your needs.

Ater the assembly of the FEM matrices, one can call
:cpp:expr:`TPZMatrixSolver::GetPardisoControl` (or :cpp:expr:`TPZEigenSolver::GetPardisoControlA`, when solving EVPs with the :ref:`section-kryloveigensolver`)
to have access to the corresponding instance of the :cpp:expr:`TPZPardisoSolver` class.

Then, the methods :cpp:expr:`TPZPardisoSolver::GetParam`  and :cpp:expr:`TPZPardisoSolver::SetParam` can be used to customise the :code:`iparm` array.

Upon usage, the parameters that were found to be of most importance are, as follows:

- :code:`iparm[7]`: Maximum number of iterative refinement steps that the solver performs when perturbed pivots are obtained during the numerical factorization.
- :code:`iparm[8]`: Tolerance level for the relative residual in the iterative refinement process. (:math:`10^{-param[8]}`)
- :code:`iparm[9]`: Perturb the pivot elements with :math:`10^{-param[9]}`.
- :code:`iparm[10]`: Use nonsymmetric permutation and scaling MPS
- :code:`iparam[12]`: Maximum weighted matching algorithm is switched-off (default for symmetric).

Further documentation on TPZPardisoSolver
+++++++++++++++++++++++++++++++++++++++++

.. doxygenclass:: TPZPardisoSolver
   :members:
