#ifndef TPZSimpleTimer_hpp
#define TPZSimpleTimer_hpp

#include <chrono>
#include <string>
/** 
 * @brief Simple timer \ref util "Utility"
 * @ingroup util
 */
class TPZSimpleTimer
{
public:
	/** @brief Start the timer when the object is created */
	TPZSimpleTimer();
    /** @brief Start the timer and gives it a name */
	TPZSimpleTimer(const std::string name);
	/** @brief Default destructor (stops the timer) */
	~TPZSimpleTimer();
	/** @brief When called, returns the time since the creation of the object in a string */
	std::string ReturnTimeString();
	/** @brief When called, returns the time since the creation of the object in a double */
	double ReturnTimeDouble();

private:
    //! Timer Name
    std::string fName{"No name"};
//@{
	//! Initial and final time calculates.
    std::chrono::time_point<std::chrono::high_resolution_clock> fBegin;
    std::chrono::time_point<std::chrono::high_resolution_clock> fEnd;
//@}
};

#endif