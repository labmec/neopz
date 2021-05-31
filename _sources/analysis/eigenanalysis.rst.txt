FEM Analysis of Eigenvalue Problems
###################################


.. toctree::
   :maxdepth: 1
   
.. contents:: Table of Contents
   :local:

The TPZEigenAnalysis class
***************************

The :cpp:expr:`TPZEigenAnalysis` class should be use for the analysis of (both standard or generalised) eigenvalue problems.

Its usage is fairly simple and the following documentation should suffice.

Users interested in these kind of problems should read the page on :doc:`available eigensolvers in NeoPZ<../solver/eigensolvers>`, specially the sections on :cpp:expr:`TPZKrylovEigenSolver` and :cpp:expr:`TPZSpectralTransform`.


.. doxygenclass:: TPZEigenAnalysis
   :members:
   :membergroups: Constructors FEM