/** @defgroup interpolation Generation of an aproximation space
 * classes which define the aproximation space combine the geometric aproximation, definition of the element shapefunctions
 * integration rules, and variational statement to build the element stiffness matrices
 * the finite elements within this tree are also responsible to generate C0 compatible interpolation spaces
 */

/** @defgroup geometry The Geometric approximation classes
 * The geometry classes are responsible for the aproximation of the geometry of the problem and
 * for keeping track of the topology of the mesh
 */

/** @defgroup CompMesh Comutational mesh classes
 * @ingroup interpolation
 * The computational mesh classes are responsible for the organise mesh data structure
 */

/** @defgroup CompElement  Comutational element classes
 * @ingroup interpolation
 * The computational element classes are responsible for the organise element data structure
 */

/** @defgroup CompElements-1D One dimensional computational elements
 * @ingroup CompElement
 * One dimensional computational elements
 */

/** @defgroup CompElements-2D Bi dimensional computational elements
 * @ingroup CompElement
 * Bi dimensional computational elements
 */

/** @defgroup CompElements-3D Tri dimensional computational elements
 * @ingroup CompElement
 * Tri dimensional computational elements
 */
