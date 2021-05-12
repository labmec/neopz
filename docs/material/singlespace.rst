.. _section-single-mat:

The Single Space Material Group
===============================

.. toctree::
   :maxdepth: 1

.. contents:: Table of Contents
   :local:

Single Space Materials
----------------------
      
The :cpp:expr:`TPZMatSingleSpaceT` interface is dedicated for implementing weak formulations using only one approximation space (regardless of which). Any given material using only one approximation space should then inherit from :cpp:expr:`TPZMatBase<TPZMatSingleSpaceT<TVar>,AnyInterfaces>`. As an example of single space materials, check :cpp:expr:`TPZElasticity3D` or :cpp:expr:`TPZMatPoisson`.


.. doxygenclass:: TPZMatSingleSpace
   :members:
      
.. doxygenclass:: TPZMatSingleSpaceT
   :members:

Single Space Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As discussed in :ref:`section-interfaces`, each interface shall also define a boundary condition counterpart. In this case, it is the :cpp:expr:`TPZMatSingleSpaceBC`. This class mostly forwards any :cpp:expr:`TPZMatSingleSpaceBC::Contribute` calls to the appropriate :cpp:expr:`TPZMatSingleSpaceT::ContributeBC`.

.. doxygenclass:: TPZMatSingleSpaceBC
   :members:
   :protected-members:

.. _section-single-interfaces:

Available Interfaces
--------------------


Interface
^^^^^^^^^

Interface for discontinuous materials

.. doxygenclass:: TPZMatInterfaceSingleSpace
   :members:
   :membergroups: Interface

Error
^^^^^

Interface for computing the error of the FEM solution based on an exact solution.

.. doxygenclass:: TPZMatError
   :members:
   :membergroups: Error
                  
.. doxygenclass:: TPZMatErrorSingleSpace
   :members:
   :membergroups: Error

.. include:: commoninterfaces.rst