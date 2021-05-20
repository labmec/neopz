/**
 * @file
 * @brief Creates analysis group for Doxygen documentation.
 */

/**
 * @defgroup analysis The Analysis classes.
 * @brief Objects of this classes implement the analysis procedures
 */
/**
 
 \page pg_analysis The TPZLinearAnalysis class and its attributes
 

In order to perform a finite element analysis using a given computational mesh, the following steps need to be performed
\li The equations need to renumbered to optimize the bandwidth
\li The storage format of the global matrix needs to be choosen
\li The system resolution procedure needs to be choosen (e.g. direct solver, iterative solver, preconditioner)
\li The postprocessing variables and output file format needs to be specified

These steps are coordinated by the TPZLinearAnalysis class

\section an_classes Renumbering the system of equations
\section an_storage Choices of global matrix storage patterns
\section an_solver Choices of system inversion procedures
\section an_postprocessing Post processing the results
*/