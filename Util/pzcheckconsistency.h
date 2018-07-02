/*
 * @file
 * @brief Interface to check the consistency of two implementations. To parallelism using.
 */
#ifndef PZCHECKCONSISTENCY
#define PZCHECKCONSISTENCY
#include "TPZSavable.h"

#include <string>

/**
 * @brief Implements an interface to check the consistency of two implementations. \ref util "Utility"
 * @ingroup util
 * @author Philippe Devloo
 * @since 27/12/2009
 */
class TPZCheckConsistency
{
	/** @brief A counter which will be used to compose the file name */
	int fCounter;
	/** 
	 * @brief Path where the file will be stored. If CHECKPATH is defined then this will be used as the file path, 
	 * else the current directory will be the path
	 */
	std::string fPath;
	/** @brief base file name. This name has to be unique to avoid overwrites */
	std::string fFileName;
	/** @brief Boolean indicating whether the objects should be overwritten or not */
	/** This flag is set to false by default */
	bool fOverWrite;
	
	/** @brief Boolean indicating whether the object will write the objects to disk or read them */
	/**
	 * The value of this variable is set default to true or false depending on the compiler directive READCHECK or WRITECHECK
	 */
	bool fWriteFlag;
	
public:
	/** @brief Constructor indicating the root filename */
	TPZCheckConsistency(const std::string &filename);
	
	/** @brief Set the overwrite flag */
	void SetOverWrite(bool flag = true);
	
	/** @brief Set the read mode. No write. */
	void SetReadMode()
	{
		fWriteFlag = false;
	}
	/** @brief Set the write mode. */
	void SetWriteMode()
	{
		fWriteFlag = true;
	}
	/** @brief Reads or writes depending on the compiler directive */
	/**
	 * when writing the method will write a binary copy of the object to a binary file \n
	 * when reading the method will read an object from the binary file and compare both copies
	 */
	bool CheckObject(TPZSavable &obj);
};

#endif
