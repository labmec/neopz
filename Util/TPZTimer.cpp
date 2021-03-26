/** 
 * @file
 * @brief Contains the implementation of the methods to TPZTimer class.
 */

#include <sstream>
#include <algorithm>

#include "TPZTimer.h"

/// Define to print two digits (increment zeros)
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
	clock_t value = clock( );
	resources.start = value;
	
}

// Stops the timer and accumulate.
void TPZTimer::stop()
{
	clock_t value;
	
	value = clock( );

	resources.elapsed = value - resources.start;
	
	AccumSec = ((double) resources.elapsed)/CLOCKS_PER_SEC;
}

// Prints the time nicely formated.
std::ostream& operator<<( std::ostream& Out, const TPZTimer& t )
{
	double AcS = t.AccumSec;
	
	int elapsedHours = static_cast< int >( AcS / 3600.0 );
	AcS -= static_cast< double >( elapsedHours * 3600 );
	
	int elapsedMins = static_cast< int >( AcS / 60.0 );
	AcS -= static_cast< double >( elapsedMins * 60 );
	
	int elapsedSecs = static_cast< int >( AcS );
	AcS -= static_cast< double >( elapsedSecs );
	
	double elapsedMSecs = AcS * 1000.0;
	std::ostringstream o;
	
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
std::ostream& operator<<( std::ostream& Out, const TPZMultiTimer& t )
{
	// First, we find the bigger process name, to align things.
	int MaxNameSize = 0;
	
	for( int ii = 0; ii < t.nTimers(); ii++ )
	{
		MaxNameSize = std::max( MaxNameSize,
							   static_cast< int >( t.processName( ii ).size() ) );
	}
	
	int Temp;
	
	for( int ii = 0; ii < t.nTimers(); ii++ )
	{
		Temp = MaxNameSize - static_cast< int >( t.processName( ii ).size() ) + 1;
		
		Out << t.processName( ii ) << std::string( Temp, ' ' )
		<< ": " << t.getTimer( ii ) << std::endl;
	}
	
	return Out;
}
