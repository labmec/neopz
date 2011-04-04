/**
\mainpage The PZ environment

\author Philippe R B Devloo
\author Misael M Santana
\author Jorge Lizardo Diaz Calle
\author Edimar Cesar Rylo
\author Gustavo Camargo Longhin
\author Cedric Marcelo Ayalo Bravo
\author Erick Raggio Slis dos Santos
\author Tiago Luis Duarte Forti
\author Cesar Lucci
\author Denise A Siqueira
\author Agnaldo Monteiro Farias
\author Diogo Cecilio
\author Nathan Shauer
\author Others

The PZ environment is a object oriented environment for the deveopment finite element simulations.

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
found in \ref sec_structure. This same structure is "more or less" recognized in the 
<a href="modules.html">Modules</a> section.
The section \ref sec_finite_element_different is dedicated to describing which algorithms within the PZ 
environment are different from regular finite element codes


\page page_finite_element_different
\section sec_finite_element_different Differences from Regular Finite Element Computations
PZ integrates zero, one, two and three dimensional simulations into a single finite element library.
It also incorporates non linear geometric maps, hp adaptive meshes and runs a large variety of finite
element simulations. It should therefore come as no surprise that its structure is somewhat different
from textbook finite element structures.

In this section we describe which finite element concepts were modified or extended in the PZ 
environment and how these concepts translated in an object oriented framework

\subsection sec_neighbour Neigbouring Information
\subsection sec_jacobian Jacobian Matrix
\subsection sec_topological Topological Concepts associated with an Element
\subsection sec_template_elements Elements based on templates
\subsection sec_matrix Matrix concept as a Linear Transformation
\subsection sec_solver A Matrix invertion procedure as an object
\subsection sec_restraints Shape function restraints
\subsection sec_connect Grouping Multiplier Coefficients in an object


\page page_structure
\section sec_structure Structure of the Environment



\page adv_technologies Advance Finite element Technologies

\section sec_advanced Advanced Finite Element Technologies
As advanced finite element technologies we denominate finite element techniques which are 
generally not available in textbook finite element codes. PZ is able to generate adaptive meshes,
interpolation between meshes, nonlinear geometric maps, multigrid iterations, continuous and
discontinuous approximation spaces, among others.

\subsection sec_nonlinear Nonlinear Geometric Maps
\subsection sec_uniformh Geometric Element Uniform Refinement
\subsection sec_patternref Goemetric Element and Refinement Patterns
\subsection sec_prefinement Shape Functions of Arbitrary Order
\subsection sec_restraints Shape Function Restraints
\subsection sec_discontinous Discontinous Approximation Spaces
*/
