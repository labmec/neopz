// Sorting TPZVec and friends.

#include <iostream>

#include "pzvec.h"
#include "pzstack.h"
#include "pzvec_extras.h"

#define S 10

using namespace std;

int main()
{
   cout << "Vector sorting test using STL (std::sort)..." << endl << endl;

   typedef TPZVec< int > Vector_t;
   typedef TPZStack< int > Stack_t;

   Vector_t v( S );
   Stack_t s;

   for( int ii = 0; ii < S; ii++ )
   {
      v[ ii ] = S - ii;

      s.Push( S - ii );
   }

   cout.precision( 2 );

   cout << "Original vector:" << endl << setw( 6 ) << v << endl;
   cout << "Sorted vector:"   << endl << setw( 6 ) << Sort( v ) << endl;

   cout << "Original stack:" << endl << setw( 6 ) << s << endl;
   cout << "Sorted stack:"   << endl << setw( 6 ) << Sort( s ) << endl;

   return 0;
};

//--| PZ |----------------------------------------------------------------------
