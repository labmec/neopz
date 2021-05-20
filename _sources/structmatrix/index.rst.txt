The StructMatrix hierarchy
==========================

.. toctree::

   parallelinterfaces
   structoptions


.. contents:: Table of Contents

The :cpp:expr:`TPZStructMatrix` hierarchy is responsible for creating the interface between the Linear Algebra and the Finite Element. It is the set of classes responsible for creating the Finite Element matrix and coordinating its assembly. The available Structural Matrices are listed in :doc:`structoptions`, and the Parallel Schemes in :doc:`parallelinterfaces`.

The Structural Matrix Template
------------------------------

This is the class from which the structural matrices should derive. It defines how to implement a structural matrix for an arbitrary matrix storage format. See its subclasses in :doc:`structoptions`

.. doxygenclass:: TPZStructMatrixT
   :members:

The Parallel Layer Interface
----------------------------

This class defines the interface for a custom parallel strategy for the assembly of a structural matrices. Users that wish to implement a custom parallel strategy should do so in a derived class. See its available derived classes in :doc:`parallelinterfaces`.

.. doxygenclass:: TPZStrMatParInterface
   :members:
                  

The Structural Matrix Base Class
--------------------------------

In order to provide an insight of its members and functionalities that will be transmitted to the derived classes, here we present :cpp:expr:`TPZStructMatrix`.

.. doxygenclass:: TPZStructMatrix
   :members:
