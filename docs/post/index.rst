Post-processing results
==========================

.. toctree::
   :maxdepth: 2
   
.. contents:: Table of Contents
   :local:


TPZVTKGenerator (usage)
---------------------------

The :cpp:expr:`TPZVTKGenerator` is a utility class for post-processing quantities
obtained by a FEM simulation in a straightforward way.
It efficiently generates files in the `VTK Legacy format <http://www.princeton.edu/~efeibush/viscourse/vtk.pdf>`__ .

In the following examples, it is assumed that a solution was computed for a given
computatinal mesh ``cmesh`` and that the header
``#include <TPZVTKGenerator.h>`` was included.

For further information on how to setup the post-processing quantities, see documentation
for the :cpp:expr:`TPZMaterial` class, specially :cpp:expr:`TPZMaterial::VariableIndex` and :cpp:expr:`TPZMaterial::NSolutionVariables`.


Post-processing all materials of a given dimension
++++++++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: cpp

   const std::string plotfile = "myfile";//no .vtk at the end
   constexpr int vtkRes{2};//vtk post-processing resolution

   //list of variables to be post-processed. It can be scalar, vectors, tensors, etc.
   TPZVec<std::string> fields = {
   "scalar_var",
   "vec_var",
   "another_scalar_var",
   "yet_another_scalar_var"};

   const int postprocdim{2};
   //if postprocdim is not set, it will be set as cmesh->Dimension()
   auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes,postprocdim); 
   

   //generates the .vtk file
   vtk.Do();
   //compute another solution, etc, etc
   vtk.Do();
   //if the elements are refined after a call to TPZVTKGenerator::Do,
   //one must then call

   vtk.ResetArrays();
   vtk.Do();

Post-processing a given set of materials
++++++++++++++++++++++++++++++++++++++++

.. note::
   All materials must have the same dimension
        
.. code-block:: cpp

   const std::string plotfile = "myfile";//no .vtk at the end
   constexpr int vtkRes{2};//vtk post-processing resolution

   //list of variables to be post-processed. It can be scalar, vectors, tensors, etc.
   TPZVec<std::string> fields = {
   "scalar_var",
   "vec_var",
   "another_scalar_var",
   "yet_another_scalar_var"};

   std::set<int> postprocmats = {1,2,3};
   
   auto vtk = TPZVTKGenerator(cmesh,postprocmats,fields, plotfile, vtkRes); 
   //proceed as in the above example


Full documentation is available below

TPZVTKGenerator (full docs)
+++++++++++++++++++++++++++

.. doxygenclass:: TPZVTKGenerator
   :members: