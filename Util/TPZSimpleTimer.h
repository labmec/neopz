#ifndef TPZSimpleTimer_hpp
#define TPZSimpleTimer_hpp

#include <chrono>
#include <string>
/** 
 * @brief Simple timer \ref util "Utility"
 * This timer allows to measure execution time with nested timers support.
 * The time range corresponds to the lifetime of the TPZSimpleTimer instance.
 * For nested timers, one can configure if the information will be printed
 * only when the outermost timer is complete or not.
 * @ingroup util
 */
class TPZSimpleTimer
{
public:
    /** @brief Default constructor. Starts the timer.*/
	TPZSimpleTimer();
    /** @brief Start the timer and gives it a name.*/
	TPZSimpleTimer(const std::string name, bool alwaysPrint=false);
	/** @brief Default destructor (stops the timer) */
	~TPZSimpleTimer();
	/** @brief When called, returns the time since the creation of the object in a string */
	std::string ReturnTimeString();
	/** @brief When called, returns the time since the creation of the object in a double */
	double ReturnTimeDouble();

private:
    //! Init function, common to all constructors.
    void Init();
    //! Timer Name
    const std::string fName{"No name"};
    /**@{*/
	//! Initial and final time calculates.
    std::chrono::time_point<std::chrono::high_resolution_clock> fBegin;
    std::chrono::time_point<std::chrono::high_resolution_clock> fEnd;
    /**@}*/
    //! Pointer to last timer to be created
    inline static TPZSimpleTimer *gLastTimerCreated{nullptr};
    //! Pointer to the parent timer
    TPZSimpleTimer *fParentTimer{nullptr};
    //! String for children timers info
    std::string fNestedTimersInfo{};
    //! Whether to print info as soon as finished or wait for parent timer
    const bool fAlwaysPrint{false};
};

#endif