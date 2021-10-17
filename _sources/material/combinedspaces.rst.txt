.. _section-combined-mat:

The Combined Spaces Material Group
==================================

.. toctree::
   :maxdepth: 1

.. contents:: Table of Contents
   :local:

Combined Spaces Materials
-------------------------
      
The :cpp:expr:`TPZMatCombinedSpacesT` interface is dedicated for implementing weak formulations using a combination approximation spaces. 
Any given material using only one approximation space should then inherit from :cpp:expr:`TPZMatBase<TPZMatCombinedSpacesT<TVar>,AnyInterfaces>`. 
As an example of a combined spaces material, check :cpp:expr:`TPZWaveguideModalAnalysis`.


.. doxygenclass:: TPZMatCombinedSpaces
   :members:
      
.. doxygenclass:: TPZMatCombinedSpacesT
   :members:

Combined Spaces Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As discussed in :ref:`section-interfaces`, each interface shall also define a boundary condition counterpart. In this case, it is the :cpp:expr:`TPZMatCombinedSpacesBC`. This class mostly forwards any :cpp:expr:`TPZMatCombinedSpacesBC::Contribute` calls to the appropriate :cpp:expr:`TPZMatCombinedSpacesT::ContributeBC`.

.. doxygenclass:: TPZMatCombinedSpacesBC
   :members:
   :protected-members:

.. _section-combined-interfaces:

Available Interfaces
--------------------

.. _section-combined-interfaces-interface:

Interface
^^^^^^^^^

Interface for discontinuous materials

.. doxygenclass:: TPZMatInterfaceCombinedSpaces
   :members:
   :membergroups: Interface

.. _section-combined-interfaces-error:

Error
^^^^^

Interface for computing the error of the FEM solution based on an exact solution.
                  
.. doxygenclass:: TPZMatErrorCombinedSpaces
   :members:
   :membergroups: Error

.. _section-combined-interfaces-common:

Common
^^^^^^

See also the :doc:`commoninterfaces`.