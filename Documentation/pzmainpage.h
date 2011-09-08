/**
\mainpage The PZ environment

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

The PZ environment is a object oriented environment for the development finite element simulations.

The PZ environment (in the future quoted as simply PZ) incorporates several advanced finite element
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

More recently, the concept of generality has been extended in that the PZ library allows its user 
to choose the approximation space as well. One can approximate a differential equation with continuous
or discontinuous approximation spaces. We are working on incorporating HDiv and HCurl spaces as well.

\section sec_obective Objectives

The objective of the PZ environment is to provide its user access to advanced finite element
technologies within a coherent framework. Wherever possible those technologies should be able
to interact with each other.

What is meant by "advanced technologies" is documented in the section \ref sec_advanced

\section sec_doc_structure Structure of the Documentation
There are many ways to define a library of classes. A global view of the PZ environment is
found in \ref page_structure. This same structure is "more or less" recognized in the 
<a href="modules.html">Modules</a> section.
The section \ref page_finite_element_different is dedicated to describing which algorithms within the PZ 
environment are different from regular finite element codes


\page page_finite_element_different Differences from Regular Finite Element Computations
\section sec_finite_element_different Differences from Regular Finite Element Computations
PZ integrates zero, one, two and three dimensional simulations into a single finite element library.
It also incorporates non linear geometric maps, hp adaptive meshes and runs a large variety of finite
element simulations. It should therefore come as no surprise that its structure is somewhat different
from textbook finite element structures.

In this section we describe which finite element concepts were modified or extended in the PZ 
environment and how these concepts translated in an object oriented framework

\subsection sec_neighbour Neigbouring Information
Within the geometric mesh, all geometric elements keep track of their neighbours along all the sides 
(see \ref sec_topological) of the element
\subsection sec_jacobian Jacobian Matrix
\subsection sec_topological Topological Concepts associated with an Element
Within PZ a geometric element is considered as the union of open sets of points. 
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
In the PZ environment the parameter space is defined and/or can be modified by specifying other topologies.

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
\subsection sec_solver A Matrix invertion procedure as an object
\subsection sec_restraints Shape function restraints
\subsection sec_connect Grouping Multiplier Coefficients in an object

\page page_structure Structure of the Environment

\page adv_technologies Advance Finite element Technologies

\section sec_advanced Advanced Finite Element Technologies
As advanced finite element technologies we denominate finite element techniques which are 
generally not available in textbook finite element codes. 
PZ is able to generate adaptive meshes, interpolation between meshes, nonlinear geometric maps,
multigrid iterations, continuous and discontinuous approximation spaces, among others.

\subsection sec_nonlinear Nonlinear Geometric Maps
\subsection sec_uniformh Geometric Element Uniform Refinement
\subsection sec_patternref Goemetric Element and Refinement Patterns
\subsection sec_prefinement Shape Functions of Arbitrary Order
\subsection sec_restraints Shape Function Restraints
\subsection sec_discontinous Discontinous Approximation Spaces

\page externlibs External Libraries used in PZ
 
\section metis Metis library
If you set USING_METIS as TRUE, you must to install Metis library. Metis library is avaliable 
from <a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/download">Karypis Lab</a>
 
\section boost Boost library
 
\section fad Fad library
 
\section log4cxx Log4cxx library
 
\section pthread Pthread library
 
\page projects Projects with PZ
 
\section steam_injection Steam Injection in Reservoir
 
\section cons_law Conservation Laws
 
\section highperform High performance
 
\section adaptive hp-Adaptivity
 
\page teoria Theorical fundaments implemented in PZ
 
\section topology Topology
\section integral Numerical integration
 
\section shapes Shape functions
 
\section analysis Analysis: Solving process
 
\page tutorial Tutorial examples
\section use_matrix Using matrices classes
See <a href="group__matrix.html">Matrix</a> module
 
\section use_integral Using numerical integration classes
See <a href="group__integral.html">Numerical integration</a> module
 
\section use_util Using util classes: vector, chunk vector ...
See <a href="group__util.html">Utility</a> module

\section use_material Creating material from differential equation
\section use_analysis Solving differential equation
*/
