TPZNullMaterial
---------------

.. toctree::
   :maxdepth: 1

.. contents:: Table of Contents
   :local:

When dealing with combined approximation spaces, there are atomic FEM meshes associated with a multiphysic mesh. In these cases, one needs a dummy material in the atomic meshes, since no computation needs to be done. This can be easily done by using the following class

.. doxygenclass:: TPZNullMaterial
   :outline: