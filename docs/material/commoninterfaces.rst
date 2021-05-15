Common Interfaces
=================

These interfaces apply to both single and combined spaces materials.

.. _section-common-interfaces-error:

Error
-----

Interface for computing the error of the FEM solution based on an exact solution.

.. doxygenclass:: TPZMatError
   :members:
   :membergroups: Error

.. _section-common-interfaces-loadcases:

.. note:: Materials should use :cpp:expr:`TPZMatErrorSingleSpace` or :cpp:expr:`TPZMatErrorCombinedSpaces`.
   

Load cases
----------

Interface for solving an equation system for multiple right-hand side vectors at once

.. doxygenclass:: TPZMatLoadCasesBase
   :members:
   :membergroups: LoadCases
                  
.. doxygenclass:: TPZMatLoadCases
   :outline:

.. _section-common-interfaces-memory:

Memory
------

Interface for implementing a memory on the integration points

.. doxygenclass:: TPZMatWithMem
   :members:
   :membergroups: Memory

.. _section-common-interfaces-eigen:

Generalised Eigenvalue Problem
------------------------------

.. doxygenclass:: TPZMatGeneralisedEigenVal
   :members:
   :membergroups: Generalised
   :protected-members:
