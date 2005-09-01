#include <pzvisualmatrix.h>
#include <pzfmatrix.h>
using namespace std;

/** This function creats a Data Explorer file that allow to visualization of the value of a matrix passed as parameter */
void VisualMatrix(TPZFMatrix & matrix, char *outfilename)
{
	const int nelx = matrix.Cols();
	const int nely = matrix.Rows();
	const int neltotal = nelx * nely;
	int i,j;
	ofstream out(outfilename);
	out << "# Graphical Visualization of Matrix." << endl;
	out << "# Positions as the indexes of the matrix, beginning by column." << endl;
	out << "# The number of elements in x direction correspond to the number of the columns of the matrix." << endl;
	out << "# The number of elements in y direction correspond to the number of the rows of the matrix." << endl;

	out  << "object 1 class gridpositions counts " << nelx+1 << " " << nely +1 << endl;
	out << "origin 0. 0." << endl;
	out << "delta 1. 0." << endl;
	out << "delta 0. 1." << endl;
	out << "attribute \"dep\" string \"positions\"" << endl;
	out << endl;


	out << "object 2 class gridconnections counts " << nelx+1 << " " << nely +1 << endl;

 	out << "attribute \"element type\" string \"quads\"" << endl;
	out << "attribute \"ref\" string \"positions\"" << endl;

	out.precision(5);
	out  << "object 3 class array type float rank 0 items " << neltotal << " data follows" << endl;
	for (i = 0; i < nelx; i++) {
		for(j=0; j< nely ; j++) out << matrix(i,j) << endl;
	}
	out << "attribute \"dep\" string \"connections\" " << endl;
	out << endl;

	out << "object 4 class field" << endl;
	out << "component \"data\" value 3" << endl;
	out << "component \"positions\" value 1" << endl;
	out << "component \"connections\" value 2" << endl;
	out << "attribute \"name\" string \"Matrix\"" << endl;

	out << endl;
	out << "end" << endl;

	out.close();

	cout << "Data Explorer file " << outfilename << " was created with success!\n";
}



