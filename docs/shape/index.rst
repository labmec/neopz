The Shape hierarchy
=========================

.. toctree::
   :maxdepth: 2

.. contents:: Table of Contents
   :local:

Approximation spaces
--------------------

An essential ingredient of Galerkin finite element approximations if the definition of a finite dimensional
subspace whose basis functions will be used to build the algebraic system of equations. The NeoPZ environment
allows to user to define a variety of approximation spaces : :math:`H^1, H(div), H(curl)` or discontinuous spaces

In NeoPZ the different flavors of approximation spaces are based on polynomial spaces associated with
a master element. These basic "building blocks" of the polynomial spaces are each implemented in a class dedicated
to each element topology. The shape template classes are :

- :cpp:expr:`pzshape::TPZShapePoint`

- :cpp:expr:`pzshape::TPZShapeLinear`

- :cpp:expr:`pzshape::TPZShapeTriang`

- :cpp:expr:`pzshape::TPZShapeQuad`

- :cpp:expr:`pzshape::TPZShapeCube`

- :cpp:expr:`pzshape::TPZShapeTetra`

- :cpp:expr:`pzshape::TPZShapePrism`

- :cpp:expr:`pzshape::TPZShapePiram`

..
	Or referencing by :ref:`TPZShapePointTTE<ShapePoint>` works?

The full (doxygen generated) documentation of the classes can be found at :ref:`listofshapeclasses`.

Based on these classes a generic classes define :math:`H^1,\;H(div)\;H(curl)` functions associated with each of the shape building blocks.
 
TPZShapeH1
^^^^^^^^^^

the :cpp:expr:`TPZShapeH1` class computes H1 compatible shape functions. The method *Initialize* takes as argument the polynomial order associated with
each side of the element, an integer identifier associated with each node and initializes datastructures in the :cpp:expr:`TPZShapeData` data structure.

The method ``Shape`` fills in the data structure fPhi and fDPhi of the :cpp:expr:`TPZShapeData` object.

.. doxygenstruct:: TPZShapeH1
   :members:

TPZShapeHDiv
^^^^^^^^^^^^

.. doxygenstruct:: TPZShapeHDiv
   :members:

TPZShapeHDivKernel
^^^^^^^^^^^^^^^^^^

.. doxygenstruct:: TPZShapeHDivKernel
   :members:

TPZShapeHDivConstant
^^^^^^^^^^^^^^^^^^^^

.. doxygenstruct:: TPZShapeHDivConstant
   :members:

TPZShapeHCurl
^^^^^^^^^^^^^

.. doxygenstruct:: TPZShapeHCurl
   :members:
   
.. _listofshapeclasses:

TPZShapeData
^^^^^^^^^^^^^

.. doxygenclass:: TPZShapeData
   :members:


The list of shape template classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ShapePoint:

.. doxygenclass:: pzshape::TPZShapePoint
.. doxygenclass:: pzshape::TPZShapeLinear
.. doxygenclass:: pzshape::TPZShapeTriang
.. doxygenclass:: pzshape::TPZShapeQuad
.. doxygenclass:: pzshape::TPZShapeCube
.. doxygenclass:: pzshape::TPZShapeTetra
.. doxygenclass:: pzshape::TPZShapePrism
.. doxygenclass:: pzshape::TPZShapePiram

   
