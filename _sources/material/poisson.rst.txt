Poisson
-------

.. toctree::
   :maxdepth: 1

.. contents:: Table of Contents
   :local:

Material for solving the Poisson Equation

.. math::

   -\text{div}\,\nabla u = f

   
for one, two and three-dimensional domains. The right hand side :math:`f` can be a function of the spatial coordinates, and the class can solve for multiple right hand sides at once, and it can perform error analysis based on a given solution.

TPZMatPoisson
^^^^^^^^^^^^^

.. doxygenclass:: TPZMatPoisson
   :members: