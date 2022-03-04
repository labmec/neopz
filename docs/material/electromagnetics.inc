Electromagnetics Group
----------------------

.. toctree::
   :maxdepth: 2

.. contents:: Table of Contents
   :local:
      
This group has materials related to electromagnetic problems. For usage of these
materials, the reader is referred to the `wgma library <http://github.com/labmec/wgma>`_


Modal Analysis of Waveguides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following classes implement a :math:`H(\text{curl},\Omega)\times H^1(\Omega)` formulation for the analysis of waveguides. More precisely, they implement:

Find non-trivial :math:`\left(\beta^2,{e_t}, {e_z}\right) \in (\mathbb{C} \times [\mathbb{C}]^N \times [\mathbb{C}]^M)` such that:

.. math::
   \begin{bmatrix}
   A_{tt} & 0 \\
   0      & 0
   \end{bmatrix}
   \begin{Bmatrix}
   e_t \\
   e_z
   \end{Bmatrix}
   = - \beta ^2
   \begin{bmatrix}
   B_{tt} & B_{tz} \\
   B_{zt} & B_{zz}
   \end{bmatrix}
   \begin{Bmatrix}
   e_t \\
   e_z
   \end{Bmatrix}
   \text{,}

where

.. math::

   [A_{tt}]_{ij} &= \int_{\Omega}
   \left[
   \mu_{zz}^{-1}\left( \nabla_t\times \boldsymbol{\varphi}_i \right) \cdot
   \left( \nabla_t\times \boldsymbol{\varphi}_j \right)
   - k_0^2 \boldsymbol{\epsilon}_{xy} \boldsymbol{\varphi}_i\cdot \boldsymbol{\varphi}_j\right] d\Omega\text{,}\label{eq:fem-mat-1}\\
     [B_{tt}]_{ij} &= \int_{\Omega}  \boldsymbol{\mu}_{xy^S}^{-1}\boldsymbol{\varphi}_i  \cdot
     \boldsymbol{\varphi}_j\, d\Omega\text{,}\label{eq:fem-mat-2}\\
     [B_{tz}]_{ij} &= \int_{\Omega}  \boldsymbol{\mu}_{xy^S}^{-1}\boldsymbol{\varphi}_i  \cdot
     \nabla_t\, \varphi_j\, d\Omega\text{,}\label{eq:fem-mat-3}\\
     [B_{zt}]_{ij} &= \int_{\Omega}  \boldsymbol{\mu}_{xy^S}^{-1}\nabla_t\, \varphi_i  \cdot
     \boldsymbol{\varphi}_j
     \, d\Omega\text{,}\label{eq:fem-mat-4}\\
     [B_{zz}]_{ij} &= \int_{\Omega}
     \left[
     \boldsymbol{\mu}_r^{-1}\nabla_t\, \varphi_i  \cdot
     \nabla_t\, \varphi_j
     - k_0^2 \epsilon_{zz} \varphi_i \varphi_j\right] d\Omega\text{.}\label{eq:fem-mat-5}

This formulation was presented originally in this work :footcite:`lee1991`. In this work :footcite:`orlandini2018` a very detailed derivation along with results obtained using NeoPZ are presented.

TPZWaveguideModalAnalysis
"""""""""""""""""""""""""
This class supports materials with the following form of anisotropy:

.. math::
   \epsilon_r =
   \begin{bmatrix}
   \epsilon_{xx} & 0 & 0\\
   0 & \epsilon_{yy}  & 0\\
   0 & 0 & \epsilon_{zz}\\
   \end{bmatrix}\,\,\,\,\,
   \mu_r =
   \begin{bmatrix}
   \mu_{xx} & 0 & 0\\
   0 & \mu_{yy}  & 0\\
   0 & 0 & \mu_{zz}\\
   \end{bmatrix}

.. doxygenclass:: TPZWaveguideModalAnalysis
   :members:
   :protected-members:
   :membergroups: ContributeMethods InternalContributeMethods SolutionMethods GeneralisedMethods InternalGeneralisedTypes InternalGeneralisedAttributes

.. note::
   
   The methods :cpp:expr:`TPZWaveguideModalAnalysis::SetMatrixA` and :cpp:expr:`TPZWaveguideModalAnalysis::SetMatrixB` are responsible for ensuring that the right call to the contribute method is done depending on which matrix is being assembled.

TPZWaveguideModalAnalysisPML
""""""""""""""""""""""""""""

This class implements a rectangular PML to be used with TPZWaveguideModalAnalysis.


.. doxygenclass:: TPZWaveguideModalAnalysisPML
   :members:


Scattering Analysis of Planar Waveguides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following classes implement a :math:`H^1(\Omega)` formulation for the scattering analysis of planar waveguides. 

Assuming that

- the computational domain is in the :math:`yz` plane
- the direction of no variation of the waveguide is :math:`x`
- at the incidence plane :math:`z=z'`, denoted by :math:`\Gamma`
  there is a source :math:`f_s: f_s(y)` that propagates with :math:`e^{-j\beta(z-z')}`
  in the :math:`+z` direction and with :math:`e^{j\beta(z-z')}`
  in the :math:`-z` direction (such as a mode field obtained by a
  previous one-dimensional modal analysis, ffor instance),

this material implements:

Find non-trivial :math:`\{\phi\} \in [\mathbb{C}]^N` such that:

.. math::
   [A]\{\phi\} = [f]\text{,}

where

.. math::
   [A]_{ij} &= \int_{\Omega}
   \left[
   c_z \frac{\partial \phi_i}{\partial z}\frac{\partial \phi_j}{\partial z}
   +
   c_y \frac{\partial \phi_i}{\partial y}\frac{\partial \phi_j}{\partial y}  
   - k_0 c_s \phi_i \phi_j\right] \,d\Omega\\
   [f]_i &= -2 \int_{\Gamma}c_z \phi_i -j \beta f_s(y)\,d\Gamma\,

where we have, for both TE and TM modes, the following meanings:

+---------------+-------------------------+------------------------------+
|               | TE                      | TM                           |
+===============+=========================+==============================+
| :math:`\phi`  | :math:`E_x`             | :math:`H_x`                  |
+---------------+-------------------------+------------------------------+
| :math:`c_z`   | :math:`\mu_{r,yy}^{-1}` | :math:`\epsilon_{r,yy}^{-1}` |
+---------------+-------------------------+------------------------------+
| :math:`c_y`   | :math:`\mu_{r,zz}^{-1}` | :math:`\epsilon_{r,zz}^{-1}` |
+---------------+-------------------------+------------------------------+
| :math:`c_s`   | :math:`\mu_{r,xx}^{-1}` | :math:`\epsilon_{r,xx}^{-1}` |
+---------------+-------------------------+------------------------------+

Note that the media are isotropic everywhere except in the PML regions,
which justify the subscripts in the table above.

This formulation was presented originally in this work :footcite:`koshiba02`, and it
is specially interesting when analysing waveguide discontinuities.


.. note::
   For simplicity, in the implementation, the computational domain is assumed
   to be in the :math:`xy` plane (source propagates in the :math:`x` direction).

TPZPlanarWGScattering
"""""""""""""""""""""
Implements the scattering analysis of a planar waveguide.

.. doxygenclass:: TPZPlanarWGScattering
   :members:
      
TPZPlanarWGScatteringPML
""""""""""""""""""""""""
Implements a PML for the planar waveguide scattering analysis.

.. doxygenclass:: TPZPlanarWGScatteringPML
   :members:
.. footbibliography::