/**
 * \file
 * @brief Creates material group for doxygen documentation.
 */


/** @defgroup material The Material Classes
 *
 * @brief Classes which implement the weak statement of the differential equation over a domain region within the NeoPZ environment.
 *
 * TPZMaterial is the Base Class, but any given Material should inherit from TPZMatBaseallows ditional interfa
 * @note It is noteworthy to observe that this definition does not depend on the definition of the interpolation space.
 *
 **/

/** @defgroup matsinglespace The Single Approximation Space Material Classes
 * @brief Classes which implement interfaces for representing weak statements of differential equations involving only one approximation space.
 * @ingroup material
 **/

/** @defgroup singlespaceinterface Interfaces for the Single Approximation Space Material Classes
 * @brief Classes which implement interfaces for enhancing single approximation spaces materials.
 * @ingroup matsinglespace
 **/


/** @defgroup matcombinedspaces The Combined Approximation Spaces Material Classes
 * @brief Classes which implement interfaces for representing weak statements of differential equations combining different approximation spaces.
 * @ingroup material
 **/

/** @defgroup combinedspacesinterface Interfaces for the Combined Approximation Spaces Material Classes
 * @brief Classes which implement interfaces for enhancing combined approximation spaces materials.
 * @ingroup matsinglespace
 **/