.. NeoPZ documentation master file, created by
   sphinx-quickstart on Wed Mar 24 17:39:38 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   within the documentation one will find
   general restructured text info : https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html
   sphinx syntax  https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#basic-markup
   breathe syntax https://breathe.readthedocs.io/en/latest/directives.html
   doxygen directives https://www.doxygen.nl/index.html
   blog on installing doxygen breathe sphinx https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/
   

Welcome to NeoPZ's documentation!
=================================

.. toctree::
   :maxdepth: 1
   :caption: Contents:
             
   analysis/index.rst
   material/index.rst
   structmatrix/index.rst
   solver/index.rst
   shape/index.rst
   util/index.rst
   
   
Introduction
------------

This documentation is an ongoing work. The following sections have been written:

- :doc:`analysis/index` How a FEM analysis is performed in NeoPZ
- :doc:`material/index` How weak formulations are implemented
- :doc:`structmatrix/index` How different matrix storage formats can be used
- :doc:`solver/index` Available solvers in NeoPZ
- :doc:`shape/index` Shape functions available in NeoPZ
- :doc:`util/index` Utility classes available in NeoPZ


.. note::
   In this documentation the most relevant methods and classes are presented and discussed. Users are encouraged to use Doxygen to generate a complete documentation.

About
-----

The NeoPZ environment is a object oriented environment for the development
finite element simulations.
 
The NeoPZ environment (in the future quoted as simply NeoPZ) incorporates
several advanced finite element technologies in a single coherent structure,
allowing its user to develop sophisticated simulations in a short period of
time.


Motivation: Why develop a finite element library?
-------------------------------------------------

During my PhD work (late 1980's) I developed hp-adaptive finite element 
algorithms applied to the simulation of compressible fluid flow. The first 
version of the adaptive mesh datastructure dates back to 1984.

I soon noticed that adaptivity is a universal concept which can be applied to
virtually any finite element simulation. During the time I studied in Texas,
adaptivity was applied to the Stokes equations,to plasticity, to thermal
problems, convection problems etc.

On the other hand, It was obvious that writing an hp-adaptive code requires a
major investment. It takes at least two years to write and validate a three
dimensional adaptive finite element code.

At that time I imagined it would be possible to write a finite element
framework that would be allow its user to apply hp-adaptive strategies to
different systems of differential equations in a single framework.

More recently, the concept of generality has been extended in that the NeoPZ
library allows its user to choose the approximation space as well. One can
approximate a differential equation with continuous or discontinuous
approximation spaces. Denise Siqueira implemented two dimensional HDiv spaces
in the library.
 
Objectives
----------

The objective of the NeoPZ environment is to provide its user access to
advanced finite element technologies within a coherent framework. Wherever
possible those technologies should be able to interact with each other.

.. What is meant by "advanced technologies" is documented in the section ADD_SEC

Authors
-------

.. list-table::
   :widths: 50 50
   :header-rows: 0

   * - Philippe Remy Bernard Devloo (`Lattes <http://lattes.cnpq.br/6051486998967925>`__)
     - Jorge Lizardo Diaz Calle (`Lattes <http://lattes.cnpq.br/2049910703027682>`__)
   * - Edimar Cesar Rylo (`Lattes <http://lattes.cnpq.br/7462096912445959>`__)
     - Gustavo Camargo Longhin (`Lattes <http://lattes.cnpq.br/9121612523149859>`__)
   * - Erick Raggio Slis dos Santos (`Lattes <http://lattes.cnpq.br/6586851137916033>`__)
     - Tiago Luis Duarte Forti (`Lattes <http://lattes.cnpq.br/9586074227742751>`__)
   * - Paulo Cesar de Alvarenga Lucci (`Lattes <http://lattes.cnpq.br/5381087404504911>`__)
     - Denise de Siqueira (`Lattes <http://lattes.cnpq.br/8437756334087793>`__)
   * - Agnaldo Monteiro Farias (`Lattes <http://lattes.cnpq.br/2401725550781559>`__)
     - Joao Luis Gonçalves (`Lattes <http://lattes.cnpq.br/2719190119956611>`__)
   * - Diogo Lira Cecilio (`Lattes <http://lattes.cnpq.br/2594284000782489>`__)
     - Nathan Shauer (`Lattes <http://lattes.cnpq.br/5762871737832497>`__)
   * - Cedric Marcelo Augusto Ayala Bravo (`Lattes <http://lattes.cnpq.br/3642648349492905>`__)
     - Renato Gomes Damas (`Lattes <http://lattes.cnpq.br/9705909592533525>`__)
   * - Misael Luis Santana Mandujano
     - Others

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
