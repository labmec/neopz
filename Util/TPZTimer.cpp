/** @file TPZTimer.cpp */

// $Id: TPZTimer.cpp,v 1.1 2005-01-10 16:11:56 phil Exp $

#include <sstream>
#include <algorithm>

#include "TPZTimer.h"

using namespace std;

#define DIGITS2( STREAM, TIME )			\
if( TIME == 0 )					\
{						\
   STREAM << "00";				\
}						\
else						\
{						\
   if( ( TIME > 0 ) && ( TIME  < 10 ) )		\
   {						\
      STREAM << "0";				\
   }						\
   						\
   STREAM << TIME;				\
}

// Starts the timer.

void TPZTimer::start()
{
#if HAVE_WAIT3
   gettimeofday( &resources.start, ( struct timezone * ) 0 );
#else
   long value;
   struct tms tms;

   value = times( &tms );
   resources.start.tv_sec = value / HZ;
   resources.start.tv_usec = value % HZ * ( 1000000 / HZ );
#endif
}

// Stops the timer and accumulate.

void TPZTimer::stop()
{
#if HAVE_WAIT3
   gettimeofday( &resources.elapsed, ( struct timezone * ) 0 );
#else  // !HAVE_WAIT3
   long value;
   struct tms tms;

   value = times( &tms );

   resources.ru.ru_utime.tv_sec = tms.tms_cutime / HZ;
   resources.ru.ru_stime.tv_sec = tms.tms_cstime / HZ;

#if HAVE_SYS_RUSAGE_H
   resources.ru.ru_utime.tv_nsec = tms.tms_cutime % HZ * (1000000000 / HZ);
   resources.ru.ru_stime.tv_nsec = tms.tms_cstime % HZ * (1000000000 / HZ);
#else
   resources.ru.ru_utime.tv_usec = tms.tms_cutime % HZ * (1000000 / HZ);
   resources.ru.ru_stime.tv_usec = tms.tms_cstime % HZ * (1000000 / HZ);
#endif

   resources.elapsed.tv_sec = value / HZ;
   resources.elapsed.tv_usec = value % HZ * (1000000 / HZ);
#endif  // !HAVE_WAIT3

   resources.elapsed.tv_sec -= resources.start.tv_sec;

   if( resources.elapsed.tv_usec < resources.start.tv_usec )
   {
      // Manually carry a one from the seconds field.
      resources.elapsed.tv_usec += 1000000;
      --resources.elapsed.tv_sec;
   }

   resources.elapsed.tv_usec -= resources.start.tv_usec;

   AccumSec += ( resources.elapsed.tv_sec +
		 resources.elapsed.tv_usec / 1000000.0 );
}

// Prints the time nicely formated.
ostream& operator<<( ostream& Out, const TPZTimer& t )
{
   double AcS = t.AccumSec;

   int elapsedHours = static_cast< int >( AcS / 3600.0 );

   AcS -= static_cast< double >( elapsedHours * 3600 );

   int elapsedMins = static_cast< int >( AcS / 60.0 );

   AcS -= static_cast< double >( elapsedMins * 60 );

   int elapsedSecs = static_cast< int >( AcS );

   AcS -= static_cast< double >( elapsedSecs );

   double elapsedMSecs = AcS * 1000.0;

   ostringstream o;

   DIGITS2( o, elapsedHours );

   o << ":";

   DIGITS2( o, elapsedMins );

   o << ":";

   DIGITS2( o, elapsedSecs );

   o << " :: " << elapsedMSecs;

   Out << o.str();

   return Out;
}

// Prints the time nicely formated.

ostream& operator<<( ostream& Out, const TPZMultiTimer& t )
{
   // First, we find the bigger process name, to align things.

   int MaxNameSize = 0;

   for( int ii = 0; ii < t.nTimers(); ii++ )
   {
      MaxNameSize = max( MaxNameSize,
			 static_cast< int >( t.processName( ii ).size() ) );
   }

   int Temp;

   for( int ii = 0; ii < t.nTimers(); ii++ )
   {
      Temp = MaxNameSize - static_cast< int >( t.processName( ii ).size() ) + 1;

      Out << t.processName( ii ) << string( Temp, ' ' )
	  << ": " << t.getTimer( ii ) << endl;
   }

   return Out;
}

//--| MB |----------------------------------------------------------------------
