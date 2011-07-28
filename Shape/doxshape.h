
/** @defgroup shape Computation of shape functions
 * 
 * The shape functions are implemented in static classes, i.e. classes whose methods are not associated with
 * a given object. This should facilitate their reuse in existing finite element codes
 * In the future, all element features which are object independent will be transferred to these classes
 * There is no base class of the Shape classes in order to improve efficiency (no virtual call mechanisms)
 */
