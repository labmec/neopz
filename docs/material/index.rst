The Material hierarchy
==========================

.. toctree::
   :maxdepth: 2
              
   singlespace.rst
   combinedspaces.rst
   availablemats.rst
   commoninterfaces.rst
   
.. contents:: Table of Contents
   :local:
  
      

The Material Group
------------------

The :cpp:expr:`TPZMaterial` hierarchy is where the weak statement of a differential
equation is defined in NeoPZ.

For a single integration point, the main role of the material is to compute
the contribution to the FEM matrix and to the right hand side vector and,
after the linear system is solved, the post-processing variables.

Different regions of the domain can have different materials, or maybe
multiple instances of a same material with different properties (e.g. in
scenarios where the domain is homogeneous by parts).

For a list of materials (weak statements) available in NeoPZ, please check :doc:`availablemats`.

The Material Hierarchy Explained
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

During a simulation, the computational element knows which approximation space (or which *combination* of approximation spaces) is being used. The element will then call its associated material for computations on the integration points. These calls vary wildly depending on the simulation being performed: using discontinuous approximation spaces one often has to perform calculations on the interface between elements, which is not needed for continuous spaces, for instance.

In order not to pollute the base class of the :cpp:expr:`TPZMaterial` hierarchy with every possible interface, the following inheritance scheme was devised:


.. cpp:alias:: TPZMaterial
      
This class has all the common type-agnostic interface that is available in a given material.

.. cpp:alias:: TPZMaterialT

This class complements the interfaces of :cpp:expr:`TPZMaterial` with some type-related interface, which depends if the material is solving a problem resulting in a real matrix :cpp:expr:`TVar=STATE` or a complex matrix :cpp:expr:`TVar=CSTATE`. It also introduces the Forcing Function (more about that later).
      
.. cpp:alias:: TPZMatBase

The :cpp:expr:`TPZMatBase` is a variadic class template, meaning that it can be instantiated with an arbitrary number of the :cpp:expr:`Interfaces` parameters. This will be used so one can create interfaces according to their needs. 
The class TPZMatBase is a subclass of :cpp:expr:`TPZMaterialT<TVar>` and of each :cpp:expr:`Interfaces` parameter.
You can read more about these interfaces in the section :ref:`section-interfaces`.

.. note::
   :cpp:expr:`TPZMatBase` is a variadic class template. Therefore
        
   .. code-block:: cpp
                   
      auto *mat = dynamic_cast<TPZMatBase<STATE>*>(someMat);

   is expected to fail, since :cpp:expr:`TPZMatBase<STATE>` and :cpp:expr:`TPZMatBase<STATE,Interface1,Interface2>` are different types. Cast to :cpp:expr:`TPZMaterial` is always possible, to :cpp:expr:`TPZMaterialT<TVar>` with the correct :cpp:expr:`TVar` type or to any of the interface classes.

There are two special interfaces that have their own detailed documentation, :cpp:expr:`TPZMatSingleSpaceT` and :cpp:expr:`TPZMatCombinedSpacesT`. These are the interfaces that define, respectively, the behaviour of materials implementing a formulation with only one approximation space and with combined approximation spaces. See :doc:`singlespace` and :doc:`combinedspaces` for further details on these interfaces.

Forcing Function
^^^^^^^^^^^^^^^^

The :cpp:expr:`TPZMaterialT::ForcingFunction` is a mechanism for position-dependent functions to be evaluated at an integration point. Its type is defined as:

.. doxygentypedef:: ForcingFunctionType

Lambdas can be used for easily creating forcing functions as

.. code-block:: cpp

   /*let us assume that mat inherits from a TPZMaterialT<STATE>*/

   
   auto forcingFunction = [](const TPZVec<REAL>&x, TPZVec<STATE>&u){
     u[0] = x[0]*x[0];
     u[1] = x[1]*x[1];
     u[2] = x[2]*x[2];
   };

   //pOrder is the suggested integration rule for the forcing function
   constexpr int pOrder{2};
   mat->SetForcingFunction(forcingFunction,pOrder);

The class :cpp:expr:`TPZMaterialT` contains the forcing function-related methods.



Further documentation on TPZMaterial Hierarchy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. doxygenclass:: TPZMaterial
   :members:
   :membergroups: Constructors Utilities PostProcess

.. doxygenclass:: TPZMaterialT
   :members:
      
.. doxygenclass:: TPZMatBase
   :members:
   :membergroups: MandatoryOverride OptionalOverride


Material Data
-------------

At an integration point, the information available to the material is given by the
:cpp:expr:`TPZMaterialData` hierarchy.


.. note:: See :doc:`TPZShapeData<../shape/index>` for information regarding
          the reference element


These classes relate to the information that is needed to compute the element matrices.
They contain the basis functions and its derivatives in the deformed element,
as well as informations regarding the geometric transform that relates the master
and the deformed elements.

.. doxygenclass:: TPZMaterialData
   :members:
   :membergroups: Flags Data

.. doxygenclass:: TPZMaterialDataT
   :members:
   :membergroups: Data

.. _section-boundaryconditions:

Boundary Conditions
-------------------

The Boundary Condition Hierarchy explained
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Usually there is a different contribution to be made by the materials associated with a boundary. For this purpose, there is a special hierarchy for boundary conditions, the :cpp:expr:`TPZBndCond` hierarchy. The boundary conditions are usually created from a given :cpp:expr:`TPZMaterial` and it will be of the type :cpp:expr:`TPZBndCondBase`, a variadic class template with as many :cpp:expr:`InterfaceBC` parameters as are present in original material. Using the CreateBC interface therefore ensures that the boundary condition will contain all the interfaces present in the original material.

The programmer will note that there is a unique boundary interface associated with each interface. The interface type is
declared by a "using InterfaceBC" statement declared in each material interface

.. cpp:alias:: TPZBndCond
      
.. cpp:alias:: TPZBndCondT
      
.. cpp:alias:: TPZBndCondBase

Creating Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following snippet illustrates how to create a boundary condition

.. code-block:: cpp
                
   /*
   * mat: pointer to a given material
   *      representing a weak form using only real variables
   * cmesh: pointer to the computational mesh
   */

   /*
   val1 is the matrix of values that goes on the FEM matrix.
   In this example, it is a 1x1 matrix
   */
   TPZFMatrix<STATE> v1(1,1,0.);
   //val2 is the vector of values that goes on the rhs
   TPZVec<STATE> v2 = {0};
   /*
   boundary identifier: each identifier should be unique in a given mesh
   */
   const int bcIdentifier{-1};
   /*
   boundary condition type: In the available NeoPZ materials,
                            Dirichlet=0,
                            Neumann=1,
                            Robin=2
                            but in your own materials any value can be used
   */
   const int bcType = 0;
   
   //Creates a boundary condition associated with mat
   TPZBndCond * bnd = mat->CreateBC(mat,bcIdentifier,boundType,val1,val2);

   //inserts into the computational mesh
   cmesh->InsertMaterialObject(bnd)


Read the next section for a better comprehension of the :cpp:expr:`val1` and :cpp:expr:`val2` variables.

.. note:: The :cpp:expr:`bcIdentifier` should correspond to a (geometrical) mesh region (each element in the mesh are associated with a region through an identifier).

.. note:: The created boundary condition will automatically inherit any interface that :cpp:expr:`mat` has.

Implementing Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Examples of boundary conditions can be seen in the :doc:`availablemats`, but for clarity, a quick example follows.

Taking the :ref:`TPZMatSingleSpaceT<section-single-mat>` class as an example, the boundary condition contribution is performed on


.. cpp:alias:: TPZMatSingleSpaceT::ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight, TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc) 
      
The function parameter :cpp:expr:`bc` has all the information regarding the boundary condition. We will give special attention to four methods  that are specially useful for computing contributions to boundaries:

.. cpp:alias:: TPZBndCondT::HasForcingFunctionBC
      
.. cpp:alias:: TPZBndCondT::Val1
      
.. cpp:alias:: TPZBndCondT::Val2
      
.. cpp:alias:: TPZMaterial::BigNumber
      
With :cpp:expr:`HasForcingFunctionBC`, one can check if the boundary condition has an associated Forcing Function, which is a position-dependent function that might be used to specify boundary conditions. Its signature is:

.. doxygentypedef:: ForcingFunctionBCType
      
When forcing functions are not being used, than :cpp:expr:`Val1` and :cpp:expr:`Val2` should be used for extracting the user-prescribed values for the boundary conditions. A simple example to illustrate their usage follows:

.. code-block:: cpp
                
   /*
   * This is inside the ContributeBC method
   * we suppose that the material is expecting
   * a vec of size 3 for val 2
   * and a 3x3 matrix for val1
   */

   //TPZManVector has static and dynamic storage
   TPZManVector<TVar,3> v2(3,0.);
   //TPZFNMatrix is a matrix equivalent of TPZManVector
   TFNMatrix<9,TVar> v1 (3,3,0.);
   
   if(bc.HasForcingFunctionBC()){
     bc.ForcingFunctionBC()(data.x,v2,v1);
   } else {
     v2 = bc.Val2();
     v1 = bc.Val1();
   }

   /*example of dirichlet BC using big number*/
   const auto &big = TPZMaterial::BigNumber();
   switch (bc.Type()){		
   // Dirichlet condition
   case 0 : {      
     for(auto in = 0 ; in < phr; in++) {
       ef(in,0) += big * v2[0] * phi.GetVal(in,0) * weight;
         for (auto jn = 0 ; jn < phr; jn++) {
           ek(in,jb) += big * phi(in,0) * phi(jn,0) * weight;
         }//jn
       }//in
     break;
  }
  //more code...

.. note::
   Homogeneous Dirichlet boundary conditions can also be implemented by making use of the method

   .. cpp:alias:: TPZStructMatrix::FilterEquations
         
   This is a slightly more advanced usage. As soon as there is a tutorial or quick example on that, this section will be updated.


Further documentation on TPZBndCond hierarchy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: TPZBndCond
   :members:
   

.. doxygenclass:: TPZBndCondT
   :members:
      
.. doxygenclass:: TPZBndCondBase
   :members:
      
.. _section-interfaces:

Material Interfaces
-------------------

Additional functionality may be incorporated to a material through an interface.

The most important interfaces are, as already mentioned, :cpp:expr:`TPZMatSingleSpaceT` and :cpp:expr:`TPZMatCombinedSpacesT`. These classes provide the interface for materials to perform computations using only one approximation space or a combination of approximation spaces, respectively. These classes share many similarities, as can be seen on :doc:`singlespace` and :doc:`combinedspaces`.

However, a material can have an arbitrary number of interfaces. Examples of available interfaces are:

- :cpp:expr:`TPZMatErrorSingleSpace`/:cpp:expr:`TPZMatErrorCombinedSpaces`: to perform error analysis
- :cpp:expr:`TPZMatWithMem`: for materials with memory in the integration point
- :cpp:expr:`TPZMatLoadCases`: for materials with multiple load cases
- etc

Interfaces are specific to the number of spaces of the material.
For details on the interfaces available for each material group,
see :doc:`singlespace` and :doc:`combinedspaces`. The common interfaces are in :doc:`commoninterfaces`.

..
	I believe the :doc: directive indicates a link to the singlespace.rst file
	
Creating custom material interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Material interfaces can be easily created to adapt a material to one's needs.

However, special care must be taken due to the materials associated with boundary conditions: In some cases the boundary condition material might need to forward function calls to its associated material, or sometimes it needs to be cast to the interface type.

For this to happen, an auxiliary class for the interface is needed, as it can be seen in the following example:

.. code-block:: cpp
                
   //forward declare bc auxiliary class
   class TPZMyInterfaceBC;

   class TPZMyInterface : public TPZSavable{
   public:
     using TInterfaceBC = TPZMyInterfaceBC;
     //declare all the virtual methods and etc as usual
     //@{
     //! Read and Write methods
     int ClassId() const override;
     void Read(TPZStream&,int) override;
     void Write(TPZStream&,void*) override;
     //@}
   };

   //maybe this forward declaration is needed
   class TPZMaterial;
   
   class TPZMyInterfaceBC : public TPZMyInterface{
   protected:
     TPZMyInterface *fMatMyInterface{nullptr};
     void SetMaterialImpl(TPZMaterial* mat){
       fMatMyInterface = dynamic_cast<TPZMyInterface*>(mat);
     }
   public:
   //overload the methods as needed.
   };

In this example, an interface is given in the class :cpp:expr:`TPZMyInterface`. The rules for :cpp:expr:`TPZMyInterface` to be a valid interface are, as follows:

- It must derive from :cpp:expr:`TPZSavable`.
- It must have an associated boundary condition class, even if there is nothing to be done in this class.
- The associated boundary condition class must be publicly avaliable under the alias :cpp:expr:`TPZMyInterface::TInterfaceBC`.
- The associated boundary condition class must have a (protected) method :cpp:func:`void SetMaterialImpl(TPZMaterial *)`. This method can be useful if a casted pointer to the associated material is needed when forwarding calls.
- The read and write methods must be overriden.