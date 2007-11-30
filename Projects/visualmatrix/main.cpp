#include <pzvisualmatrix.h>
#include <fstream>

using namespace std;

int main()
{
  cout << "visual matrix." << endl;
  TPZFMatrix matrix(3,3);
  int i, j;
  for(i=0; i< matrix.Cols(); i++){
    for(j=0; j<matrix.Rows(); j++){
      matrix(i,j)= i*3 +j;
    }
  }
  VisualMatrix(matrix, "visualmatrix.dx");

}
