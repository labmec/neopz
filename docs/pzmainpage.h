/**
 \mainpage The NeoPZ environment
 
 \author Philippe Remy Bernard Devloo <a href="http://lattes.cnpq.br/6051486998967925">Lattes</a>
 \author Jorge Lizardo Diaz Calle <a href="http://lattes.cnpq.br/2049910703027682">Lattes</a>
 \author Edimar Cesar Rylo <a href="http://lattes.cnpq.br/7462096912445959">Lattes</a>
 \author Gustavo Camargo Longhin <a href="http://lattes.cnpq.br/9121612523149859">Lattes</a>
 \author Erick Raggio Slis dos Santos <a href="http://lattes.cnpq.br/6586851137916033">Lattes</a>
 \author Tiago Luis Duarte Forti <a href="http://lattes.cnpq.br/9586074227742751">Lattes</a>
 \author Paulo Cesar de Alvarenga Lucci <a href="http://lattes.cnpq.br/5381087404504911">Lattes</a>
 \author Denise de Siqueira <a href="http://lattes.cnpq.br/8437756334087793">Lattes</a>
 \author Agnaldo Monteiro Farias <a href="http://lattes.cnpq.br/2401725550781559">Lattes</a>
 \author Joao Luis Gonçalves <a href="http://lattes.cnpq.br/2719190119956611">Lattes</a>
 \author Diogo Lira Cecilio <a href="http://lattes.cnpq.br/2594284000782489">Lattes</a>
 \author Nathan Shauer <a href="http://lattes.cnpq.br/5762871737832497">Lattes</a>
 \author Cedric Marcelo Augusto Ayala Bravo <a href="http://lattes.cnpq.br/3642648349492905">Lattes</a>
 \author Renato Gomes Damas <a href="http://lattes.cnpq.br/9705909592533525">Lattes</a>
 \author Misael Luis Santana Mandujano
 \author Others
 
 The NeoPZ environment is a object oriented environment for the development finite element simulations.
 
 The NeoPZ environment (in the future quoted as simply NeoPZ) incorporates several advanced finite element
 technologies in a single coherent structure, allowing its user to develop sophisticated simulations
 in a short period of time.
 
 \section sec_motivation Motivation: Why develop a finite element library?
 
 During my PhD work (late 1980's) I developed hp-adaptive finite element algorithms applied to the
 simulation of compressible fluid flow. The first version of the adaptive mesh datastructure dates
 back to 1984.
 
 I soon noticed that adaptivity is a universal concept which can be applied to virtually any finite
 element simulation. During the time I studied in Texas, adaptivity was applied to the Stokes equations,
 to plasticity, to thermal problems, convection problems etc.
 
 On the other hand, It was obvious that writing an hp-adaptive code requires a
 major investment. It takes at least two years to write and validate a three dimensional adaptive
 finite element code.
 
 At that time I imagined it would be possible to write a finite element framework that would be
 allow its user to apply hp-adaptive strategies to different systems of differential equations in
 a single framework.
 
 More recently, the concept of generality has been extended in that the NeoPZ library allows its user
 to choose the approximation space as well. One can approximate a differential equation with continuous
 or discontinuous approximation spaces. Denise Siqueira implemented two dimensional HDiv spaces in the library.
 Douglas Castro is working on its three dimensional extension.
 We are working on incorporating HCurl spaces as well.
 
 
 \section sec_obective Objectives
 
 The objective of the NeoPZ environment is to provide its user access to advanced finite element
 technologies within a coherent framework. Wherever possible those technologies should be able
 to interact with each other.
 
 What is meant by "advanced technologies" is documented in the section \ref adv_technologies
 
 \section sec_doc_structure Structure of the Documentation
 There are many ways to define a library of classes. A global view of the NeoPZ environment is
 found in \ref page_structure. This same structure is "more or less" recognized in the
 <a href="modules.html">Modules</a> section.
 The section \ref page_finite_element_different is dedicated to describing which algorithms within the NeoPZ
 environment are different from regular finite element codes
 
 
 \page init Initial Information
 
 The information over utilitaries to work, update and get documentation is in the following pages:
 \li Utilities needed to configure the NeoPZ environment \ref utilitaries
 \li External Libraries used in NeoPZ \ref externlibs
 
 \page utilitaries - Utilities needed to configure the NeoPZ environment
 
 \section svn Getting NeoPZ code
 
 It is recommended to use <a href="http://www.syntevo.com/smartsvn/download.html?all=true">SmartSVN</a> to get the NeoPZ code. Actually we are using SmartSVN 6.6.9 .
 
 \section cmake Creating project depending on the system user
 
 It is recommended to use <a href="http://www.cmake.org/cmake/resources/software.html">CMake</a> to generate the neopz project depending on your system. \n
 Actually we are using CMake 2.8.5 .
 
 \section doxygen Generating documentation
 
 It is recommended to use <a href="http://www.stack.nl/~dimitri/doxygen/download.html#latestsrc">Doxygen</a> to generate the neopz documentation. \n
 Actually we are using Doxygen 1.7.5.1 . To right compiling using doxygen you must to have the following executables:
 
 \li <a href="http://www.cs.utah.edu/dept/old/texinfo/dvips/dvips.html">dvips</a> or
 <a href="http://rpmfind.net/linux/rpm2html/search.php?query=ghostscript-dvipdf">dvipdf</a> - Convert tex to ps (post script) or pdf format.
 
 \section manuals Manuals
 
 To get or access the manuals clik on following links:
 \li <a href="http://www.syntevo.com/download/smartsvn/smartsvn-reference.pdf">SmartSVN</a>
 \li <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html">CMake</a>
 \li <a href="ftp://ftp.stack.nl/pub/users/dimitri/doxygen_manual-1.7.5.1.pdf.zip">Doxygen</a>
 
 
 \page externlibs - External Libraries used in NeoPZ
 
 It is recommended to create a directory libs or externallibs at the same level as neopz project. \n
 External libraries which are not necessarily installed at root lib directory or root include directory \n
 can to be installed or copied into the indicated directory.
 
 NeoPZ uses until five external libraries:
 \li boost - It provides free portable peer-reviewed C++ libraries. We mainly was used as unit test framework.
 \li log4cxx - For efficient log messages.
 \li fad - For automatic differentiation.
 \li metis - To partitioning finite element meshes, and producing fill reducing orderings for sparse matrices
 
 \section metis Metis library
 If you set USING_METIS as TRUE, you must to install Metis library. \n
 Metis library is avaliable
 from <a href="http://www.labmec.org.br/pz/libexternal/metis">Metis 5.0</a>
 
 The least metis version actually suported by PZ is 5.0.x \n
 METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes, \n
 and producing fill reducing orderings for sparse matrices.
 
 \section boost Boost library
 If you set USING_BOOST as TRUE it is necessary to install the Boost library. \n
 Get the latest version of BOOST library at download from <a href="http://sourceforge.net/projects/boost/files/boost">Sourceforge</a>. \n
 
 It is recommended to use version <a href="http://sourceforge.net/projects/boost/files/boost/1.49.0">1.49.0</a> or newer.
 
 Use boost_1_49_0.tar.gz or boost_1_49_0.tar.bz2 for unix or mac systems. \n
 Use boost_1_49_0.7z or boost_1_49_0.zip for windows system.
 
 To install following next steps:\n
 Uncompress the version downloaded. \n
 Using command line change into the uncompress directory. \n
 For mac or unix systems, type the following commands and execute: \n
 \li sudo ./bootstrap.sh
 \li sudo ./bjam install
 
 For Windows systems, execute: \n
 \li .\\bootstrap.bat
 \li .\\b2 install or .\\bjam install
 
 See <a href="http://www.boost.org/doc/libs/1_49_0/doc/html/bbv2/installation.html">Installation</a>
 
 \section fad Fad library
 
 The neopz project uses a library that implements automatic differentiation. This library is distributed toghether with the source code of NeoPZ
 under the subdirectory External
 
 \section log4cxx Log4cxx library
 
 Apache <a href="http://www.labmec.org.br/pz/libexternal/log4cxx/">log4cxx</a> is a logging framework for C++ pattern. It has three main components: loggers, appenders and layouts.
 These three types of components work together to enable developers to log messages according to message type and level,
 and to control at runtime how these messages are formatted and where they are reported.
 
 To install for unix or mac systems, <a href="http://www.labmec.org.br/pz/libexternal/log4cxx/">download</a> apache-log4cxx-0.10.0.zip for windows system, or apache-log4cxx-0.10.0.tar.gz for unix or mac systems. \n
 Then uncompress the archive. Using command line change into uncompress directory. \n
 Type de following commands: \n
 \li ./configure
 \li make check
 \li sudo make install
 
 \page neopz About the PZ Library
 
 The description of the neopzlib is in the following pages:
 \ref page_finite_element_different
 \ref page_structure
 \ref adv_technologies
 \ref teoria
 
 \page page_finite_element_different - Differences from Regular Finite Element Computations
 
 \section sec_finite_element_different Differences from Regular Finite Element Computations
 NeoPZ integrates zero, one, two and three dimensional simulations into a single finite element library.
 It also incorporates non linear geometric maps, hp adaptive meshes and runs a large variety of finite
 element simulations. It should therefore come as no surprise that its structure is somewhat different
 from textbook finite element structures.
 
 In this section we describe which finite element concepts were modified or extended in the NeoPZ
 environment and how these concepts translated in an object oriented framework
 
 \subsection sec_neighbour Neigbouring Information
 Within the geometric mesh, all geometric elements keep track of their neighbours along all the sides
 (see \ref sec_topological) of the element. This information is usually not encountered in other finite element programs.
 
 Within NeoPZ a neighbour is represented by objects of type TPZGeoElSide. An object of type TPZGeoElSide is \em almost a geometrical element:
 \li it has a set of nodes
 \li it has a parametric space associated with it
 \li it has a geometric map, implemented by a method which computes the coordinates at a parametric point and a jacobian
 \li it has a set of sons of type TPZGeoElSide (if the associated geometric element is divided)
 \li it has a father element
 
 \subsection sec_jacobian Jacobian Matrix
 
 The concept of Jacobian matrix in regular finite elements is the \(n \otimes n\) matrix which maps a parameter in the master element to a point
 in \(R^3\). In NeoPZ all elements are in \(R^3\). As such a line element is mapped to a segment in three dimensional space, etc. Therefore the gradient of the
 map can not be represented in a traditional way. The Jacobian matrix is computed as the QR decomposition of the gradient of \(\bf x(\xi)\). The upper triangular
 square matrix /(R/) is \em called the Jacobian matrix and the orthogonal matrix \(Q\) is called the \em axes througout the software. This has evoluated this way
 because at an earlier stage of the NeoPZ development it seemed to P. Devloo he had developed a new approach to computing the Jacobian matrix.
 
 The \(QR\) decomposition is computed using a Gram Schmidt orthogonalization. It is worth mentioning that this approach makes the NeoPZ library ready to describe
 geometric maps in any space dimension.
 
 \subsection sec_topological Topological Concepts associated with an Element
 Within NeoPZ an element is considered as the union of open sets of points.
 These sets of points are named sides. As such:
 \li linear element (pztopology::TPZLinear) has 3 sides (2 points and one line)
 \li quadralaterial element (pztopology::TPZQuad) has 9 sides (4 points 4 lines and one area)
 \li triangular element (pztopology::TPZTriangle) has 7 sides (3 points 3 lines and one area)
 \li hexahedral element (pztopology::TPZCube) has 27 sides (8 points 12 lines 6 quadrilaterials and one volume)
 \li prism element (pztopology::TPZPrism) has 21 sides (6 points 9 lines 2 triangles 3 quadrilaterials and one volume)
 \li pyramid element (pztopology::TPZPyramid) has 19 sides (5 points 8 lines 5 triangles and one volume)
 \li tetrahedral element (pztopology::TPZTetrahedra) has 15 sides (4 points 6 lines 4 triangles and one volume)
 \li point element (pztopology::TPZPoint) has one side: the point itself
 
 
 All geometries are grouped in the namespace \ref pzgeom. The topology themselves are defined in the namespace \ref pztopology.
 
 Each topology is associated with an area within the dimension associated with the topology.
 For example the one dimensional line element is associated with the line segment \f$]-1,1[\subset R\f$.
 A quadrilateral element is associated with the area \f$]-1,1[\times]-1,1[\subset R^2\f$.
 The area associated with a topology is named parameter space.
 In finite element textbooks the parameter space is associated with the space of the master element.
 Theoretically each finite element code can define its own parameter space.
 In the NeoPZ environment the parameter space is defined and/or can be modified by specifying other topologies.
 
 Each sides of an element associated with a topology (point, line, quadrilateral, etc).
 The closure of a side (remember that a side is an open set of points) includes its neighbouring topologies.
 For instance the closure of the line includes two point topologies, the closure of a quadrilateral topology
 includes the four lines and four points.
 
 The topology associated with a side of a topology is returned in the method Type(int side).
 This method exists in all classes of the \ref pztopology namespace
 
 The sides included in the closure of a given side are returned in the method LowerDimensionSides.
 
 As each side has its own parameter space, an affine parameter transformation can be defined between
 the lower dimension sides and the side itself.
 This affine transformation is returned in the SideToSideTransform method
 
 \subsection sec_template_elements Elements based on templates
 
 A more intuitive approach to object oriented finite element programming is to associate a distinct class with each element topology. This was the first approach
 to programming the NeoPZ library. After incorporating the three dimensional elements it became apparent that very similar code was replicated in the seven
 traditional topologies. This code replication generated a huge overhead in code maintenance : improvement and/or corrections for one topology had be
 duplicated to all topologies. These replicated code remained frequently untested.
 
 Migrating the code to a class template concept was a maior improvement in terms of code maintenance. The downside is that template oriented code is much more difficult
 to read and more distant from traditional (i.e. fortran) coding concepts.
 
 The template oriented implementation of geometric elements allows us to incorporate nonstandard maps as a geometric element specialization. As such we have implemented
 maps from a line to circel segment (in 3D space), a quadratic map, a \em wavyline map and a map to a NACA profile. All nonstandard maps are put in the Special Maps directory.
 
 \subsection sec_matrix Matrix concept as a Linear Transformation
 
 Finite element programming is essentially matrix analysis. Matrix objects are an essential ingredient of finite element programming. Within the NeoPZ environment a matrix
 object is a linear transformation capable of mapping a vector (matrix) to a vector (matrix). Some \em types of matrices implement decomposition procedures but the availability
 of such procedure is not essential. The storage pattern of each matrix class should not be of concern to the user.
 
 Matrices can be of any type defined by a template parameter. All matrix classes are instantiated with float, double, long double and the corresponding complex types.
 
 There is one matrix class which is \em storage \em declared : TPZFMatrix<T> implements a full matrix and stores its elements in a column-maior order (http://en.wikipedia.org/wiki/Row-major_order ). This allows the other matrix classes to implement matrix vector multiplications efficiently by using pointer arithmetic.
 
 \subsection sec_solver A Matrix inversion procedure as an object
 
 One of the characteristics of Finite element research is the quest of the more efficient procedure for inverting the global system of equations. Inumerous possibilities
 present themselves: direct inversion (LU, Cholesky or LDLt), iterative inversion using a Krylov method (conjugate gradient, GMRES, CGStab, Bi-CGStab), preconditioning
 based on approximate inverses, multigrid etc. It is impractical to provide the user of NeoPZ with a comprehensive list of inversion procedures.
 
 Instead, a TPZMatrixSolver<T> class (where T is the data type e.g. float, double, etc) was idealized which represents a transformation process applied to a matrix object
 and a right hand side. Its derived class
 TPZMatrixSolver<T> stands for the same concept, but has a matrix object associated with it. TPZStepSolver<T> presents an interface allowing to choose between a
 direct solver a regular iterative solver (e.g. Jacobi, SOR, SSOR) or de preconditioned Krylov solver. The preconditioner is represented by a Solver<T> object.
 
 Most known solution procedures can and have been implemented using this class structure.
 
 \subsection sec_restraints Shape function restraints
 
 Ever since his first developments in adaptivity, one of the leading authors of NeoPZ has used the concept of /em hangin /em nodes. The name hanging node refers to
 the fact that the element mesh is nonconforming and that it can be seen that a node /em hangs on the side between two other nodes (in 2D). The concept was generalized
 to hp-meshes by the same author in 1987 and extended to three dimensions in the same year. Withing NeoPZ, hanging nodes are referred to as /em coefficient /em restraints
 because this is how they are effectively implemented. Moreover, shape functions are generally not associated with nodes.
 
 Coefficient restraints convey the idea that, in order to obtain conforming approximation space (i.e. continuous), the multiplying coeficients of certain shape functions
 can not be choosen independently: they are linearly dependent on the value of the multiplying coefficients of other shape functions. This then introduces the concept of
 dependent and independent connects.
 
 The values which define the depdendency of the coefficient is stored in the TPZConnect objects. The TPZConnect class implements a method HasDependency() that verifies
 whether the TPZConnect object has a datastructure which defines the relationship between its coefficients and the coefficients of other TPZConnect objects. It is
 conceivable that a TPZConnect object depends on a TPZConnect object that also has dependency.
 
 The TPZConnect objects which are dependent on other objects have no correspondence in the global system of equations. There value is updated in the TPZCompMesh::LoadSolution method.
 
 \subsection sec_connect Grouping Multiplier Coefficients in a TPZConnect object
 The TPZConnect class represents multiplier coefficients associated with a set of shape functions.
 When two elements share a vertex, the continuity of the solution is obtained by establishing that the multiplying
 coeficients associated with all elements which share the node are identical.
 
 Withing PZ, H1 shapefunctions are associated with the sides of the elements (see \ref sec_neighbour). As a consequence a TPZConnect object is created for each
 side of H1 elements. In a discontinuous Galerkin approximation a unique TPZConnect object is associated with each element. In this case no other element
 will be associated with this connect leading a discontinuous approximation.
 
 The TPZConnect object contains information related to the numbering of the global system of equations, the maximum order of the shape functions associated
 with the connect, the number of state variables associated with each connect, the number of shape functions associated with the connect, whether the connect is restrained, whether the connect is associated with a Lagrange
 multiplier and whether the connect has been condensed at the element level.
 
 \subsubsection subsec_sequencenumber Sequence number
 The structure of the global system of equations generated by the finite element approximation is determined by the global equation number associated with each shape function.
 (The set of shape functions determine the approximation space). By modifying the global equation number associated with the shape function, the finite element approximation
 remains the same. The computational effort needed to invert the global system of equations, on the other hand, is strongly dependent on the ordering of the shape functions.
 Within PZ, equations are renumbered by the TPZRenumbering class. The TPZRenumbering class relies on the Boost library (http://www.boost.org ) to optimize the bandwidth and/or fillin
 of the global system of equations.
 
 Rather than associating a global equation number with each shape function, the global equation number are associated with a group of shape functions (which is a TPZConnect object).
 This is why each connect contains a sequence number as the fSequenceNumber variable.
 
 If the TPZConnect object isn't associated with any element, its sequence number will be set to -1. A TPZConnect object with sequence number equal to -1 is called a "free" connect.
 
 \subsubsection subsec_orderofapproximation The order of the shape functions
 Within PZ the maximum polynomial order of the shapefunctions associated with each side of the elements can be choosen independently. The data structure which represents the
 polynomial order of approximation is the fOrder variable of the TPZConnect object. This variable is of type unsigned char. This means that the order can vary from 0 to 255.
 
 \subsubsection subsec_numberstate The number of state variables associated with each connect
 In finite element approximations of systems of partial differential equations it is customary to associate the same shape function with each state variable. For instance, in two
 dimensional elasticity each shape function is associated with the horizontal and vertical displacement. In most finite element approximations the number of state variables associated
 with each shapefunction is constant. In these cases the number of state variables can be a value associated with a mesh. In multiphysics problems, the number of state variables can
 vary acording to the physical quantity being represented. For instance, in numerical approximations of flow through porous media, two (or three) state variables are associated with
 the shapefunctions which approximate the displacements of the porous matrix and a single shape function is associated with the pressure variable. This is the reason why each connect
 keeps track of the number of state variables associated with its shape functions
 
 
 
 \page page_structure - Structure of the Environment
 
 \page adv_technologies - Advanced Finite element Technologies
 
 As advanced finite element technologies we denominate finite element techniques which are
 generally not available in textbook finite element codes.
 NeoPZ is able to generate adaptive meshes, interpolation between meshes, nonlinear geometric maps,
 multigrid iterations, continuous and discontinuous approximation spaces, among others.
 
 \section sec_nonlinear Nonlinear Geometric Maps
 
 Commercial finite element libraries allow the user to choose between linear or quadratic elements. In most finite element
 textbooks, the map between the master element and its actual location in space is computed using the same shape functions
 that are used to approximate the solution. These elements are consequently called \em isoparametric elements.
 
 The appproximation spaces implemented in NeoPZ use hierarchical shape functions that are not easily used for mapping the
 geometry of the elements. Therefore, in NeoPZ the definition of the geometric map is strongly separated from the definition of the approximation space.
 This approach allows us to implement any geometric map between a master element and a deformed element.
 
 In his master's thesis Cesar Lucci implemented an extension of the concept of transfinite blending functions originally conceived by Gordon and Hall in 1973
 to the complete family of known element topologies. This effort, developed under the sponsorship of ANP/Petrobras, allows to develop high order approximations
 of differential equations on curved domains.
 
 \section sec_uniformh Geometric Element Uniform Refinement
 
 Considering one the first motivations for developing NeoPZ was to develop a library which applies adaptive finite elements to different families of partial
 differential equations, it is natural that h-adaptivity is \em built \em in its structure.
 
 In early versions each type of element (linear, quadrilateral, triangular, etc) was implemented using its proper class, derived from the abstract TPZGeoEl class.
 This approach was very intuitive, but very difficult to maintain. Each change in the interface of the geometric element had to be repeated in each element topology.
 Frequently a feature implemented in a topology that was frequently used did not work for other topologies. This was one of the main motivations for reimplementing
 the geometric maps within a template structure. In its current version, a single c++ class, specialized by a template implements the complete family of finite elements.
 
 A geometric element is responsible for implementing/computing the mapping between the master element and deformed element. Within an adaptive context, the geometric
 element keeps track of its father and children.
 
 One possibility of refining geometric elements is to divide the elements in a given pattern: a one dimensional element is divided in two elements, a quadrilateral element
 is divided into four elements, a hexahedral element is divided into eight elements, etc. Within the NeoPZ context this refinement pattern is called \em uniform, because
 each element is divided using elements.
 
 The refinement of the individual elements is also specialized by a template class, allowing linear or nonlinear geometric maps to be divided alike.
 
 In order to understand the concepts of the geometric map, the advanced user should consult the paper published by Calle, Devloo, Bravo and Gomes. There, the essential
 concepts which guide the development are explained.
 
 One key feature of the implementation of the geometric refinement is that the elements division is performed without dynamic memory allocation.
 
 A data structure item which distinguishes NeoPZ from traditional finite element codes is that within NeoPZ each element keeps track of its neighbours. The size of
 this data structure is also constant for each element type.
 
 \section sec_patternref Geometric Element and Refinement Patterns
 
 Adaptive geometric refinement using uniform refinement patterns is capable of improving the efficiency of the computations, but in most cases the singularity of
 the functions that are being approximated are one dimensional in nature. A typical example is a problem with boundary layers.
 
 This was the motivation for developing a more advanced procedure for dividing geometric elements. In the concept of refinement patterns, each element can be divided
 in an arbitrary number of sub elements, as long as the sub elements form a partition of the original element. Although this concept is relatively simple, it
 creates issues with respect to the compatibility of the refinements of neighbouring elements.
 
 This research effort was sponsored by FAPESP in a collaborative project with Embraer as a PICTA project.
 
 Using refinement patterns very specialized meshes can be created starting from an essential geometric description of the computational domain.
 
 On the downside is that creating refined meshes now becomes a programming issue instead of a decision of which element to refine. One technique which has worked well
 is to dicide on the most appropriate refinement pattern based on identifying the edges which need to be halfed.
 
 The technology of refinement patterns was further extended by Cesar Lucci when applying refinement patterns to the numerical simulation of hydraulic fracturing. In
 this effort, refinement patterns are dynamically created to represent the edge of the fracture.
 
 \section sec_prefinement Shape Functions of Arbitrary Order
 
 Within NeoPZ the approximation space is implemented in a separate class structure. One of the main reasons is that different approximation spaces can be generated
 based on a single geometry. At the geometry level the full refinement tree is kept. The approximation space at the other hand is a partition of the computational domain.
 
 The definition of the approximation space, similarly to the geometry, is implemented as a specialization of computational element class. The geometry and approximation space
 class are derived from a single topology class. As such the pzgeom::TPZGeoQuad class and the pzshape::TPZShapeQuad classes are derived from pztopology::TPZQuadrilateral class
 
 \section sec_restraints Shape Function Restraints
 \section sec_discontinous Discontinous Approximation Spaces
 
 \page teoria - Theorical concepts implemented in NeoPZ
 
 \section topology Topology
 
 \subsection symm_quads Symmetric quadrature rules
 
 \section shapes Shape functions
 
 \section analysis Analysis: Solving process
 
 \page pg_app_projects Documentation of the application projects
 
 The application projects are projects created by students to support their graduate work. The documentation depends
 on the effort of the student.
 
 \page pg_tut_projects Documentation of the tutorial projects
 
 The tutorial projects demonstrate the functioning of the different modules of the NeoPZ library
 
 \section tut_testintegration How to create and use integration rules
 
 NeoPZ includes integration rules of arbitrary order for all element topologies. This program illustrates
 how to create an integration rule for a given topology and validates the rule for a simple integral
 
 This tutorial also verifies if the objects are able to integrate a polynom of order 500 and prints the number of necessary points
 
 It is nice introduction in the use c++ classes
 
 \section tut_testgeom To create geometric mesh
 
 The tutorial discShape implements a bi-dimensional problem using as mesh a disc with a central hole.
 
 It initializes given the coordinates of the 10 points on the \f$ x^2 + y^2 = 4 \f$ circunference and another 10 points on the \f$ x^2 + y^2 = r^2 \f$, where \f$ r \f$ is the radio of the hole.
 Then it constructs a quadrilateral elements jointing two next nodes on external boundary with another two nodes on internal boundary (hole). Also creates
 20 boundary (one-dimensional) elements on external and internal boundary.
 
 After creating the geometrical mesh with nodes and elements, it is necessary to construct the connectivity between them calling BuildConnectivity().
 
 \section tut_compmesh Computational mesh construction
 
 \section tut_material Creating material from differential equation
 
 \section tut_testonedim To solve one-dimensional differential equation
 \section tut_testtwodim To solve two-dimensional differential equation
 \section tut_analysis Solving differential equation
 
 \page pg_unit_projects Unit tests used for project validation
 
 In this projects we are using Boost framework for unit test boost_unit_test_framework. See information in \ref boost Boost .
 
 \section unit_integral To test numerical integration module
 See <a href="group__integral.html">Numerical integration</a> module
 \section unit_matrix To test matrix module
 See <a href="group__matrix.html">Matrix</a> module
 \section unit_solvers To test linear solvers for systems module
 \section unit_topology To test topology module
 See <a href="group_geom.html">Geometry</a> module.
 \section unit_linmaterial To test linear material
 \section unit_util Using util classes: vector, chunk vector ...
 See <a href="group__util.html">Utility</a> module
 
 \page special_pages SPECIAL PAGES
 
 This pages contains suplementary information over the description and implementation of the PZ modules that can be called separately.
 
 \ref numinteg
 
 
 \page numinteg - Over numeric integration
 
 
 \page deprecated OBSOLETE
 
 \li \ref TutorialGeometry
 
 */
