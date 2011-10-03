/**
 * @file
 * @brief Doxygen file for geometric tutorial example.
 */

/*! \page TutorialGeometry Example of the Creation of a Geometric Mesh
 * \dontinclude geom.cpp
 * This example can be found in the Project directory TestGeom and file geom.cpp
 * <b>MAIN Program</b>
 * \skipline main
 * Initialize the Log4cxx system, depending on the compiler directive
 * \until endif
 * Create an object of type TPZGeoMesh
 * \skipline TPZGeoMesh
 * Call the procedure which will read the mesh
 * \skipline LerMalhaGeom
 * Print the mesh to the file "output.dat"
 * \until Print
 * Create a .vtk file representing the geometric mesh
 * \until }
 * <b>Procedure to read a file and build a geometric mesh</b>
 * \skipline LerMalhaGeom
 * open the file indicated by name
 * \until ifstream
 * read control parameters from the file
 * \until >>
 * \until >>
 * \until getline
 * \until getline
 * Resize the vector of geometric nodes of the geometric mesh to hold the geometric nodes of the mesh
 * \until Resize
 * Read the elements of the mesh
 * \until getline
 * \until getline
 * Read the nodes of the mesh
 * \until getline
 * \until getline
 * Read the boundary conditions of the mesh
 * \until }
 * Initialize the connectivity of the elements
 * \until return
 * \skipline }
 */