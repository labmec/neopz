/**
 * @file
 * @brief Contains declaration of the TPZSemaphore class which implements semaphore to threads. Check Dijkstra's semaphore for more information.
 */

#ifndef TPZSEMAPHOREHPP
#define TPZSEMAPHOREHPP

#include <mutex>
#include <condition_variable>

/**
 * @ingroup util
 * @brief Implements semaphore to threads. \ref util "Utility"
 */
class TPZSemaphore
{
private:
	/** @brief Counter of the times the semaphore is locked */
	int fCounter{0};
	/** @brief Mutex for the thread */
	mutable std::mutex fMutex;
	/** @brief Condition for the thread must to be waiting */
	mutable std::condition_variable fCond;
	
public:
	/** @brief Default constructor */
	TPZSemaphore() = default;
	/** @brief Default destructor */
	~TPZSemaphore() = default;
  /** @brief Copy constructor */
  TPZSemaphore(const TPZSemaphore &);
  /** @brief Move constructor */
  TPZSemaphore(const TPZSemaphore &&) = delete;//std::mutex has no move operator
  /** @brief Assignment (copy) operator */
  TPZSemaphore &operator=(const TPZSemaphore &);
  /** @brief Assignment (move) operator */
  TPZSemaphore &operator=(TPZSemaphore &&) = delete;//std::mutex has no move operator
	/** @brief Constructor with initial number for the counter */
	TPZSemaphore(int initCount);
	/** @brief Consumes one resource (operation P in Dijkstra's definition) */
	void Wait();
  /** @brief Frees one resource (operation V in Dijkstra's definition)*/
	void Post();
};

#endif
