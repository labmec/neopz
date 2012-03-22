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
or discontinuous approximation spaces. We are working on incorporating HDiv and HCurl spaces as well.

\section sec_obective Objectives

The objective of the NeoPZ environment is to provide its user access to advanced finite element
technologies within a coherent framework. Wherever possible those technologies should be able
to interact with each other.

What is meant by "advanced technologies" is documented in the section \ref sec_advanced

\section sec_doc_structure Structure of the Documentation
There are many ways to define a library of classes. A global view of the NeoPZ environment is
found in \ref page_structure. This same structure is "more or less" recognized in the 
<a href="modules.html">Modules</a> section.
The section \ref page_finite_element_different is dedicated to describing which algorithms within the NeoPZ 
environment are different from regular finite element codes


 \page init INITIAL INFORMATION
 
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
 \li pthread - For the POSIX pthread threading.
 \li boost - It provides free portable peer-reviewed C++ libraries. We mainly was used as unit test framework.
 \li log4cxx - For efficient log messages.
 \li fad - For automatic differentiation.
 \li metis - To partitioning finite element meshes, and producing fill reducing orderings for sparse matrices
 
 \section metis Metis library
 If you set USING_METIS as TRUE, you must to install Metis library. \n
 Metis library is avaliable 
 from <a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/download">Karypis Lab</a>
 
 The least metis version actually suported by PZ is 5.0.x \n
 METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes, \n
 and producing fill reducing orderings for sparse matrices.
 
 \section boost Boost library
 If you set USING_BOOST as TRUE it is necessary to install the Boost library. \n
 Get the latest version of BOOST library at download from <a href="http://sourceforge.net/projects/boost/files/boost">Sourceforge</a>. \n
 
 It is recommended to use version <a href="http://sourceforge.net/projects/boost/files/boost/1.47.0">1.47.0</a>
 
 Use boost_1_47_0.tar.gz or boost_1_47_0.tar.bz2 for unix or mac systems. \n
 Use boost_1_47_0.7z or boost_1_47_0.zip for windows system.
 
 To install following next steps:\n
 Uncompress the version downloaded. \n
 Using command line change into the uncompress directory. \n
 For mac or unix systems, type the following commands and execute: \n
 \li sudo ./bootstrap.sh
 \li sudo ./bjam install
 
 For Windows systems, execute: \n
 \li ./bootstrap.bat
 \li ./bjam install
 
 See <a href="http://www.boost.org/doc/libs/1_47_0/doc/html/bbv2/installation.html">Installation</a>
 
 \section fad Fad library
 
 The neopz project uses a old library to fast automatic differentiation <a href=" http://www.ann.jussieu.fr/~pironnea/">FAD</a>, but it seems inactive. 
 You can to claim a copy of the source code sending e-mail to phil@fec.unicamp.br.
 
 We are testing now using another library <a href="http://www.fadbad.com/fadbad.html">FADBAD++</a>. (2007)\n
 This library implements the forward, backward and Taylor methods utilizing C++ templates and operator overloading. 
 
 Also we are testing using <a href="http://admb-project.org/downloads">AUTODIF</a> which is a library for automatic differentiation used as the building block for AD Model Builder.
 This library is current and it has versions to Windows, unix and Max systems.
 
 \section log4cxx Log4cxx library
 
 Apache <a href="http://logging.apache.org/log4cxx/">log4cxx</a> is a logging framework for C++ pattern. It has three main components: loggers, appenders and layouts.
 These three types of components work together to enable developers to log messages according to message type and level, 
 and to control at runtime how these messages are formatted and where they are reported.
 
 To install, <a href="http://logging.apache.org/log4cxx/download.html">download</a> apache-log4cxx-0.10.0.zip for windows system, or apache-log4cxx-0.10.0.tar.gz for unix or mac systems. \n
 Then uncompress the archive. Using command line change into uncompress directory. \n
 Type de following commands: \n
 \li ./configure
 \li make check
 \li sudo make install
 
 \section pthread Pthread library
 
 The neopz project uses <a href="http://staff.science.uva.nl/~bterwijn/Projects/PThread/">PThread</a> library for the POSIX pthread threading. 
 To install for unix or mac systems, make <a href="http://staff.science.uva.nl/~bterwijn/Projects/PThread/PThread.tar.gz">download</a> of the source code. Uncompress the archive and using command line change into the uncompress directory with PThread.
 Use the following commands:
 \li ./configure
 \li make
 \li sudo make install
 
 For windows system, get the header files and release libraries (dll, lib) from <a href="http://sourceware.org/pthreads-win32/>sourceware</a>. You can to connect with server as guest and then copy the include and lib directories into your system.

 
 \page neopz OVER PZ LIBRARY

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
(see \ref sec_topological) of the element
\subsection sec_jacobian Jacobian Matrix
\subsection sec_topological Topological Concepts associated with an Element
Within NeoPZ a geometric element is considered as the union of open sets of points. 
These sets of points are named sides. As such:
- linear element (pzgeom::TPZGeoLinear) has 3 sides (2 points and one line)
- quadralaterial element (pzgeom::TPZGeoQuad) has 9 sides (4 points 4 lines and one area)
- triangular element (pzgeom::TPZGeoTriangle) has 7 sides (3 points 3 lines and one area)
- hexahedral element (pzgeom::TPZGeoCube) has 27 sides (8 points 12 lines 6 quadrilaterials and one volume)
- prism element (pzgeom::TPZGeoPrism) has 21 sides (6 points 9 lines 2 triangles 3 quadrilaterials and one volume)
- pyramid element (pzgeom::TPZGeoPyramid) has 19 sides (5 points 8 lines 5 triangles and one volume)
- tetrahedral element (pzgeom::TPZGeoTetrahedra) has 15 sides (4 points 6 lines 4 triangles and one volume)
- point element (pzgeom::TPZGeoPoint) has one side: the point itself

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
\subsection sec_matrix Matrix concept as a Linear Transformation
\subsection sec_solver A Matrix inversion procedure as an object
\subsection sec_restraints Shape function restraints
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

\page adv_technologies - Advance Finite element Technologies

\section sec_advanced Advanced Finite Element Technologies
As advanced finite element technologies we denominate finite element techniques which are 
generally not available in textbook finite element codes. 
NeoPZ is able to generate adaptive meshes, interpolation between meshes, nonlinear geometric maps,
multigrid iterations, continuous and discontinuous approximation spaces, among others.

\subsection sec_nonlinear Nonlinear Geometric Maps
\subsection sec_uniformh Geometric Element Uniform Refinement
\subsection sec_patternref Goemetric Element and Refinement Patterns
\subsection sec_prefinement Shape Functions of Arbitrary Order
\subsection sec_restraints Shape Function Restraints
\subsection sec_discontinous Discontinous Approximation Spaces

 \page teoria - Theorical concepts implemented in NeoPZ
 
 \section topology Topology
 
 \subsection symm_quads Symmetric quadrature rules
 
 \section shapes Shape functions
 
 \section analysis Analysis: Solving process
 

 \page app_pz TUTORIAL AND APPLICATION PROJECTS
 
 The description of the examples, application projects and unit test projects are in the following pages:
 \ref tutorial
 \ref projects
 \ref unit_test
 
 \page tutorial - Tutorial examples

 The examples as tutorials are related by module of the PZ.
  
 \section tut_matrix Creating, Filling and Operating matrizes and Solving Linear Systems
 
 \section use_matrix Using matrices classes
See <a href="group__matrix.html">Matrix</a> module
 
\section use_integral Using numerical integration classes
See <a href="group__integral.html">Numerical integration</a> module
 
\section use_util Using util classes: vector, chunk vector ...
See <a href="group__util.html">Utility</a> module

 \section use_geometric Creating geometric objects
 
 See <a href="group_geom.html">Geometry</a> module.
 
 The tutorial discShape implements a bi-dimensional problem using as mesh a disc with a central hole.
 
 It initializes given the coordinates of the 10 points on the \f$ x^2 + y^2 = 4 \f$ circunference and another 10 points on the \f$ x^2 + y^2 = r^2 \f$, where \f$ r \f$ is the radio of the hole. 
 Then it constructs a quadrilateral elements jointing two next nodes on external boundary with another two nodes on internal boundary (hole). Also creates 
 20 boundary (one-dimensional) elements on external and internal boundary.
 
 After creating the geometrical mesh with nodes and elements, it is necessary to construct the connectivity between them calling BuildConnectivity().
 
\section use_material Creating material from differential equation
 
\section tut_compmesh Computational mesh construction

\section use_analysis Solving differential equation
 

 
\page projects - Projects with NeoPZ
 
 \section steam_injection Steam Injection in Reservoir
 
 \section cons_law Conservation Laws
 
 \section highperform High performance
 
 \section adaptive hp-Adaptivity
  
\page unit_test - Unit Test Projects for PZ Modules 
 
 In this projects we are using Boost framework for unit test boost_unit_test_framework. See information in \ref boost Boost .
 
 \page special_pages SPECIAL PAGES
 
 This pages contains suplementary information over the description and implementation of the PZ modules that can be called separately.
 
 \ref numinteg
 
 
\page numinteg - Over numeric integration

 
 \page deprecated OBSOLETE
 
 */
/*! \li \ref TutorialGeometry

*/
