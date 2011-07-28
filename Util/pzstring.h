/**
 * @file pzstring.h
 * @brief String implementation.
 */
// $Id: pzstring.h,v 1.6 2009-07-07 18:21:47 longhin Exp $

#ifndef PZSTRING_H
#define PZSTRING_H

#include "pzstack.h"
#include "pzmanvector.h"

#include <string.h>

/**
 * @ingroup util
 * @brief Implements strings as stack 
 */
class TPZString : public TPZStack < char >
{
public:
	
	/** @brief Default Constructor. */
	TPZString();
	
	/** Constructs a TPZString object based on a null char ended string. */
	//      inline TPZString(const char * source);
	
	/** @brief Constructs a TPZString object based on a null char ended string. */
	TPZString(char const * source);
	
	/** @brief Initializes the TPZString object with a predefined size. */
	TPZString(const int size);
	
	/** @brief Single char constructor. */
	TPZString(const char chr);
	
	/** @brief Destructor. */
	~TPZString() { }
	
	/**
	 * @brief Concatenates the 'this' string with another character string.
	 *
	 * Original strings are left unchanged.
	 */
	TPZString operator + (const char * increment) const;
	
	/**
	 * @brief Concatenates the 'this' string with a character.
	 *
	 * Original string and character are left unchanged.
	 */
	TPZString operator + (const TPZString & increment) const;
	
	/** @brief Appends a string at the tail. Resizes the TPZString if necessary. */
	void operator += (const char * increment);
	
	/** @brief Appends a character at the end. Resizes if necessary. */
	void operator += (const char increment);
	
	/** @brief operator =. Resizes if necessary. */
	bool operator == (const TPZString cmp)
	{
		if(!strcmp(fStore, cmp.fStore)) return true;
		else return false;
		
	}
	
	/** @brief operator =. Resizes if necessary. */
	void operator = (const char * source)
	{
		int len = strlen(source);
		
		if (NElements() < (len + 1))
		{
			Resize(len + 1);
		}
		
		strcpy(fStore, source);
	}
	
	/** @brief Appends a character at the end. Resizes if necessary. */
	void Append(const char TailIncrement);
	
	/** @brief Appends a string at the tail. Resizes the TPZString if necessary. */
	void Append(const char * TailIncrement);
	
	/**
	 * @brief Explicitly convertes a TPZString into a const null ended
	 * char string.
	 */
	const char * Str() const;
	
	/** @brief Implicit conversion */
	operator const char * () const;
	
	/**
	 * @brief Similar to strlen(string). Also returns the number of
	 * non-null characters.
	 */
	int Length() const;
	
	/**
	 * @brief Returns a subset string.
	 *
	 * If start and end are the same, a substring containing one
	 * character is returned.
	 * <li>
	 * <ul> if start <=0, start is assumed to be 0;
	 * <ul> if end < start, null pointer is returned;
	 * <ul> if end > length, end is assumed to be equal to length.
	 * </li>
	 *
	 * @param start start of the string (zero-based)
	 * @param end end of the string (zero-based)
	 */
	TPZString SubStr(const int start, const int end) const;
	
	/** @brief Empties the string. */
	void Empty();
	
	/**
	 * @brief Internally allocates the exact string size (length + 1)
	 * to store it.
	 */
	void Optimize();
	
	/** @brief Remove the repeat white spaces */
	void SimplifyWhiteSpace();
	
	/** @brief Replace the subset of string. Return the times of replacement. */
	int Replace(const char * old_str, const char * new_str);
	
	/** @brief Find the positions of the first occurence of the find string */
	int Find(const char * find_str);
};


typedef TPZVec < TPZString > TPZText;

#endif

//--| PZ |----------------------------------------------------------------------
