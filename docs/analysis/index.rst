The Analysis hierarchy
######################

.. toctree::

   linearanalysis.rst
   eigenanalysis.rst
   
.. contents:: Table of Contents
   :local:
  
      

The Analysis Group
******************

The :cpp:expr:`TPZAnalysis` hierarchy is where whole workflow for performing Finite Element Analysis in NeoPZ is defined.


This hierarchy is responsible for coordinating the interaction between all the other classes involved in a simulation, such as :cpp:expr:`TPZCompMesh`, :cpp:expr:`TPZStructMatrix` and :cpp:expr:`TPZSolver`.

For the FEM analysis of equation systems, one can use the :cpp:expr:`TPZLinearAnalysis` class, described :doc:`here<linearanalysis>`. For the analysis of eigenvalue problems, one can use the :cpp:expr:`TPZEigenAnalysis`, described :doc:`here<eigenanalysis>`.

The TPZAnalysis class
=====================

The :cpp:expr:`TPZAnalysis` class is the base of its hierarchy. Its attributions include (but are not limited to)

- Renumbering the equations in order to optimize the bandwidth
- Choosing the storage format of the global matrix (:cpp:expr:`TPZStructMatrix`)
- Choosing the solver to be used (:cpp:expr:`TPZSolver`)
- Coordinating the assemble process
- Coordinating the solving process
- Coordinating numerical post-processing on the calculated solution
- Coordinating graphical post-processing of the results

.. note::
   The :cpp:expr:`TPZAnalysis` class has a :cpp:expr:`fTime` attribute which is mainly used for post processing. However, anyone interested in creating a derived class for a transient problem may find it usefl.

Further TPZAnalysis documentation
---------------------------------

.. doxygenclass:: TPZAnalysis
   :members:
   :membergroups: Constructors MainFEM PostFEM GettersSetters Utils Graphical

Basic steps for running a FEM analysis
======================================

This section is under construction.


Customising your FEM analysis
-----------------------------

This section is under construction.

Numerical post-processing
=========================

This section is under construction.

Graphical post-processing
=========================

This section is under construction.