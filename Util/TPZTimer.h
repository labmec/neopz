/**
 * @file TPZTimer.h
 * @brief Timing class. Absolutely copied from GNU time. Take a look at <br>
 * @see http://www.gnu.org/software/time/time.html 
 */

#ifndef TPZTIMER_H
#define TPZTIMER_H

#include <string>
#include <vector>
#include <iostream>

#include <time.h>

/**
 * \addtogroup util
 * @{
 */

/** @brief Information on the resources used by a child process. \ref util "Utility" */
struct PZResourceUsage
{
	/** @brief Waiting status */
	int waitstatus;

	/** @brief Wallclock time of process (start) */
	clock_t start;
	/** @brief Wallclock time of process (elapsed) */
	clock_t elapsed;
};   // Change from the original "RESUSE".


//--| TPZTimer |----------------------------------------------------------------

/** 
 * @brief The timer class. \ref util "Utility"
 * @author Cantao!
 * @warning Check its behavior on non-GNU systems.
 */
/**
 *  The purpose of this class is to provide a simple, mono-threaded
 *  timer. It is based entirely on GNU timer, so I don't know its
 *  behavior on non-GNU systems.
 */
class TPZTimer
{
public:
	/* @brief Empty constructor */
	TPZTimer();
	
	/**
	 * @brief Default constructor.
	 * @param pn The user can give a "process name" to the timer. \n
	 * This name will be printed when the extraction operator is called.
	 */
	TPZTimer( std::string pn );
	
	/** @brief Default Destructor */
	~TPZTimer();
	
	/// Gets the process name (for reporting purposes).
	std::string& processName();
	
	/// Gets the process name (for reporting purposes).
	const std::string& processName() const;
	
	/// Turns the timer on.
	void start();
	
	/// Turns the timer off, and computes the elapsed time.
	void stop();
	
	/// Zeroes the timer.
	void reset();
	
	/// Returns the elapsed time in seconds.
	double seconds() const;
	
	/// Prints the time nicely formated.
	friend std::ostream& operator<<( std::ostream& Out, const TPZTimer& t );
	
private:
	/** @brief Information on the resources used */
	PZResourceUsage resources;
	
	/// Total accumulated time in seconds.
	double AccumSec;
	
	/// Name of the process being timed.
	std::string ProcessName;
};

//--| TPZMultiTimer |-----------------------------------------------------------

/** 
 * @brief Controls several timers at once. \ref util "Utility"
 * @warning Check its behavior on non-GNU systems.
 */
class TPZMultiTimer
{
public:
	/**
	 * @brief Default constructor.
	 * @param nT Number of different timers we want to use.
	 */
	TPZMultiTimer( int nT );
	
	/// Destructor.
	~TPZMultiTimer();
	
	/// Number of active timers.
	int nTimers() const;
	
	/// Returns a specific timer.
	TPZTimer& getTimer( int i );
	
	/// Returns a specific timer.
	const TPZTimer& getTimer( int i ) const;
	
	/// Gets the process name (for reporting purposes).
	std::string& processName( int i );
	
	/// Gets the process name (for reporting purposes).
	const std::string& processName( int i ) const;
	
	/// Turns the timer on.
	void start( int i );
	
	/// Turns the timer on.
	void start();
	
	/// Turns the timer off, and computes the elapsed time.
	void stop( int i );
	
	/// Turns the timer off, and computes the elapsed time.
	void stop();
	
	/// Zeroes the timer.
	void reset( int i );
	
	/// Zeroes the timer.
	void reset();
	
	/// Returns the elapsed time in seconds.
	double seconds( int i ) const;
	
	/// Prints the time nicely formated.
	friend std::ostream& operator<<( std::ostream& Out, const TPZMultiTimer& t );
	
private:
	/// Vector of timers.
	std::vector< TPZTimer > timers;
};

//--| IMPLEMENTATION :: TPZTimer |----------------------------------------------

// Empty constructor.
inline TPZTimer::TPZTimer() : AccumSec( 0.0 ), ProcessName( "" ) {}

// Default constructor.
inline TPZTimer::TPZTimer( std::string pn ) : AccumSec( 0.0 ), ProcessName( pn ) {}

// Destructor.
inline TPZTimer::~TPZTimer() {}

// Gets the process name (for reporting purposes).
inline std::string& TPZTimer::processName()
{
	return ProcessName;
}

// Gets the process name (for reporting purposes).
inline const std::string& TPZTimer::processName() const
{
	return ProcessName;
}

// Resets the timer.
inline void TPZTimer::reset()
{
	AccumSec = 0.0;
}

// Returns the total accumulated time in seconds.
inline double TPZTimer::seconds() const
{
	return AccumSec;
}

//--| IMPLEMENTATION :: TPZMultiTimer |-----------------------------------------

// Default constructor.
inline TPZMultiTimer::TPZMultiTimer( int nT ) : timers()
{
	timers = std::vector< TPZTimer >( nT );
}

// Destructor.
inline TPZMultiTimer::~TPZMultiTimer()
{
	// Nothing to do here!
}

// Number of active timers.
inline int TPZMultiTimer::nTimers() const
{
	return timers.size();
}

// Returns a specific timer.
inline TPZTimer& TPZMultiTimer::getTimer( int i )
{
	return timers[ i ];
}

// Returns a specific timer.
inline const TPZTimer& TPZMultiTimer::getTimer( int i ) const
{
	return timers[ i ];
}

// Gets the process name (for reporting purposes).
inline std::string& TPZMultiTimer::processName( int i )
{
	return timers[ i ].processName();
}

// Gets the process name (for reporting purposes).
inline const std::string& TPZMultiTimer::processName( int i ) const
{
	return timers[ i ].processName();
}

// Turns the timer on.
inline void TPZMultiTimer::start( int i )
{
	timers[ i ].start();
}

// Turns the timer on.
inline void TPZMultiTimer::start()
{
	for( unsigned int ii = 0; ii < timers.size(); ii++ )
	{
		timers[ ii ].start();
	}
}

// Turns the timer off, and computes the elapsed time.
inline void TPZMultiTimer::stop( int i )
{
	timers[ i ].stop();
}

// Turns the timer off, and computes the elapsed time.
inline void TPZMultiTimer::stop()
{
	for( unsigned int ii = 0; ii < timers.size(); ii++ )
	{
		timers[ ii ].stop();
	}
}

// Zeroes the timer.
inline void TPZMultiTimer::reset( int i )
{
	timers[ i ].reset();
}

// Zeroes the timer.
inline void TPZMultiTimer::reset()
{
	for( unsigned int ii = 0; ii < timers.size(); ii++ )
	{
		timers[ ii ].reset();
	}
}

// Returns the elapsed time in seconds.
inline double TPZMultiTimer::seconds( int i ) const
{
	return timers[ i ].seconds();
}

/** @} */

#endif // TPZTIMER_H
