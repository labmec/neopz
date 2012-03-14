/**
 * @file
 * @brief Creates integral group for Doxygen documentation.
 */

/** 
 * @defgroup integral The Numerical Integration classes.
 * @brief Defines integration rules for lines, quadrilaterals, triangles,
 * cubes, tetrahedra, pyramids and prisms.
 * 
 * Allow the user to define the order of the polynomial, the rule should be able to integrate exactly.
 *
 * All the computations into this module use and return long double values. Because the precision of the points and weights will be right and are independent of the 
 * typedef REAL definition.
 */

/**
 \page numinteg Over numeric integration
 
 \section integral Quadrature rules

 \subsection symm_quads Symmetric quadrature rules for triangles and tetrahedra (see [1])
 
 In this module, we present a quadrature rules used for computation of numeric integrations. This is essential because the main computational work
 is on weak formulation, generally integration form. There are many implementations for quadrature rules in the literature, some of them finding the zeros
 of the special polinomials (p.e. [2]) or resolving minimization of high-order multivariate polynomials (p.e. 2).
 
 We use just symmetric quadrature rules with quadrature points lying within the master element and the 
 weights always positives. The symmetry of the quadrature points will can to help reducing the computational cost for construction of the stiffness matrices
 by the symmetry of finite element basis functions. A quadrature rule is said to be symmetric if it is invariant under permutations of the barycentric coordinates.
 Given a simplex \f$E\f$ with \f$n\f$ vertices \f${ v_1, ..., v_n}\f$, the \f$n\f$ barycentric coordinates \f$ (\xi_1, ..., \xi_n) \f$ of a point \f$P=(x,y,z)\f$
 with respect to \f$E\f$ is determined by \f$P=\Sigma_{i=1}^n \xi_i v_i\f$ and \f$\Sigma_{i=1}^n \xi_i = 1 \f$.
 Then whether a quadrature point (in barycentric coordinates) associated with weight \f$w\f$ has its barycentric coordinates permuted, the resultant point is also a quadrature point
 of the same quadrature rule and associated with same weight \f$w\f$.
 
 
 \subsection Some References
 
 [1] Linbo Zhang, Tao Cui and Hui Liu, "Set of symmetric quadrature rules of triangles and tetrahedra." J. Comput. Mathematics, v 27, 89-96, 2009.\n
 [2] H. Li, J. Sun and Y. Xu, Discrete fourier analysis, cubature and interpolation on a hexagon and
 a triangle, SIAM J. Numer. Anal., to appear. \n
 [3] M.A. Taylor, B.A. Wingate and L.P. Bos, A cardinal function algorithm for computing multivari-
 ate quadrature points, SIAM J. Numer. Anal., 45:1 (2007), 193–205. \n
 [4] R. Cools, An encyclopaedia of cubature formulas, J. Complex., 19:3 (2003), 445–453. \n
 [5] S. Heo and Y. Xu, Constructing fully symmetric cubature formulae for the sphere, Math. Comput., 70:233 (2001), 269–279. \n
 [6] P. Hammer and A. Stroud, Numerical integration over simplexes, Mathematical Tables and Aids
 to Computation, 10 (1956), 137–139. \n
 [7] A. Stroud, Approximate Calculation of Multiple Integrals, Prentice-Hall, Englewood Cliffs, N.J., 1971. \n
 [8] A. Stroud and D. Secrest, Gaussian Quadrature Formulas, Prentice-Hall, Englewood Cliffs, N.J., 1966. \n
 [9] P. Solin, K. Segeth and I. Dolezel, Higher-Order Finite ElementMethods, Chapman and Hall/CRC Press, 2003. \n
 [10] S. Wandzura and H. Xiao, Symmetric quadrature rules on a triangle, Comput. Math. Appl., 45 (2003), 1829–1840. \n
 [11] C.A. Felippa, A compendium of FEM integration rules for finite element work, Eng. Computation, 21 (2004), 867–890. \n
 
 */