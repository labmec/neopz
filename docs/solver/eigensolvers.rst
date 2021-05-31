TPZEigenSolver
##############

.. toctree::
   :maxdepth: 2
   
.. contents:: Table of Contents
   :local:


The :cpp:expr:`TPZEigenSolver` represents a solver for Eigenvalue problems (EVP) in which the matrix has entries with the type :cpp:type:`TVar`, where :cpp:expr:`TVar=STATE`, for real problems, and :cpp:expr:`TVar=CSTATE`, for complex problems. It should be used with the :cpp:expr:`TPZEigenAnalysis` class. See :ref:`section-lapackeigensolver` and :ref:`section-kryloveigensolver` for details on the available solvers in NeoPZ.

Sorting eigenpairs
******************

The returned eigenvalues can be sorted in the following ways:

.. doxygenenum:: TPZEigenSort

Further documentation on TPZEigenSolver
***************************************

.. doxygenclass:: TPZEigenSolver
   :members:
   :membergroups: Eigen

.. _section-lapackeigensolver:
                  
TPZLapackEigenSolver
********************

This class acts as a wrapper over LAPACK calls that can be used for solving eigenvalue problems. It supports :cpp:expr:`TPZFMatrix` and :cpp:expr:`TPZSBMatrix`, therefore it can only be used with the structural matrices :cpp:expr:`TPZFStructMatrix` and :cpp:expr:`TPZSBandStructMatrix`. It provides no additional public interfaces.

.. note::
   This class is also used internally by the :cpp:expr:`TPZFMatrix` and :cpp:expr:`TPZSBMatrix` classes, thus the specific (and protected) interfaces.

.. doxygenclass:: TPZLapackEigenSolver

.. _section-kryloveigensolver:
                  
TPZKrylovEigenSolver
********************

This class implements a set of Krylov-based eigensolvers for solving standard and generalised EVP problems with no restraints regarding the storage format of the matrix. This class of eigensolvers do not solve for all eigenvalues, but for a small set of the original eigenvalues.

.. note:: For an excellent and detailed discussion on this class of eigensolvers, one can refer to the Chapters 2 and 3 of `SLEPc User Manual <https://slepc.upv.es/documentation/slepc.pdf>`_ and their numerous technical reports.
          
The Krylov-based eigensolvers search for a set of eigenvalues from a given matrix :math:`A` on the Krylov subspace:

.. math::
   \text{span}\left\{ v, Av, A^2v, \ldots, A^{n-1}v\right\}.

Currently :cpp:expr:`TPZKrylovEigenSolver` implements an Arnoldi solver, which is a good solver for computing external eigenvalues.

In this class of solvers, an orthonormal basis for the Krylov subspace is obtained by means of an Arnoldi iteration:


.. math::
   :nowrap:

   \begin{align*}
   &\text{Take an input vector $v$ and set $v_0=v$}\\
   &\text{for }k=0,1,\ldots,m\\
   &\qquad w_k = A v_k\\
   &\qquad \text{for }j =0,\ldots,k\\
   &\qquad \qquad H(j,k) = w_k.v_k\\
   &\qquad \qquad w_k -= H(j,k) v_k\\
   &\qquad \text{if } k < m-1\\
   &\qquad \qquad H(k+1,k) = \lVert w_k \rVert\\
   &\qquad \qquad v_{k+1} = w_k / H(k+1,k)
   \end{align*}


It is interesting to note that in this iteration, only the product :math:`A v_k` has to be performed, thus the sparsity of the matrix is left unchanged.

.. note::
   The matrix :math:`A` does not necessarily correspond to the matrix of a standard EVP. See :ref:`section-spectraltransform` for details on the matrix that is passed to the Arnoldi iteration.
   
After the Arnoldi iteration is performed, an EVP for :math:`H` is then solved using LAPACK. Given that the dimensions of the Krylov subspace is :math:`m\ll n_{eq}`, where :math:`n_{eq}` is the original number of equations, using a dense matrix for :math:`H` does not represent a big impact.

Further documentation on TPZKrylovEigenSolver
=============================================

.. doxygenclass:: TPZKrylovEigenSolver
   :members:
   :protected-members:
   :membergroups: BasicUsage


                  
.. _section-spectraltransform:
                  
TPZSpectralTransform
********************

Spectral transformations are used to create a transformed EVP in which the exterior eigenvalues correspond to the sought eigenvalues of the original problem, and the eigenvectors are left unchanged.

The following table summarises the transformed problem for both standard and generalised EVPs:

+-------------------------------+----------------------------+----------------------------+
| Transformation                | Standard                   | Generalised                |
+===============================+============================+============================+
| None                          | :math:`A`                  | :math:`B^{-1}A`            |
+-------------------------------+----------------------------+----------------------------+
| :code:`TPZSTShiftOrigin`      | :math:`A-\sigma I`         | :math:`B^{-1}A -\sigma I`  |
+-------------------------------+----------------------------+----------------------------+
| :code:`TPZSTShiftAndInvert`   | :math:`(A-\sigma I)^{-1}`  | :math:`(A-\sigma B)^{-1}B` |
+-------------------------------+----------------------------+----------------------------+



.. doxygenclass:: TPZSpectralTransform
   :members:
   :protected-members:
   :membergroups: BasicUsage Protected

TPZSTShiftOrigin
================

.. doxygenclass:: TPZSTShiftOrigin
   :members:
   :protected-members:
   :membergroups: BasicUsage Protected


TPZSTShiftAndInvert
===================

.. doxygenclass:: TPZSTShiftAndInvert
   :members:
   :protected-members:
   :membergroups: BasicUsage Protected