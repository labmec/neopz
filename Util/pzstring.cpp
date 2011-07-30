/** 
 * @file
 * @brief Contains the implementation of the methods to TPZString class.
 */
// $Id: pzstring.cpp,v 1.6 2005-04-25 02:55:52 phil Exp $

#include "pzstring.h"

#include <stdio.h>
using namespace std;

TPZString::TPZString()
{
	// NOTHING TO DO HERE!
}

TPZString::TPZString(char const * source)
{
	int len = strlen(source);
	Resize(len + 1);
	strcpy(fStore, source);
	
	// although the TPZStack class already stores the length, the null
	// character is necessary to export string as a null character
	// ended string
	fStore[len] = '\0';
}


TPZString::TPZString(const int size)
{
	Resize(size);
}

TPZString::TPZString(const char chr)
{
	Resize(2);
	fStore[0] = chr;
	fStore[1] = '\0';
}

TPZString TPZString::operator + (const char * increment) const
{
	TPZString newstring(* this);
	newstring.Append(increment);
	return newstring;
}

TPZString TPZString::operator + (const TPZString & increment) const
{
	TPZString newstring(* this);
	newstring.Append(increment);
	return newstring;
}

void TPZString::operator += (const char increment)
{
	this->Append(increment);
}

void TPZString::operator += (const char * increment)
{
	this->Append(increment);
}

/*
 void TPZString::operator=(const char * source)
 {
 int len=strlen(source);
 if(NElements()<(len+1))Resize(len+1);
 strcpy(fStore, source);
 }
 */

TPZString::operator const char * () const
{
	return Str();
}

const char * TPZString::Str() const
{
	if (fNElements > 0)
	{
		return fStore;
	}
	
	return NULL;
}

int TPZString::Length() const
{
	if (fNElements == 0) return 0;
	return strlen(fStore);
}

void TPZString::Append(const char TailIncrement)
{
	int len = Length();
	// if the allocated memory is less than the length +1 (null char)
	// +1 (new char), reallocates it
	if (fNElements < len + 2)
	{
		Push('\0');
	}
	else
	{
		fStore[len + 1] = '\0';
	}
	
	fStore[len] = TailIncrement;
}

void TPZString::Append(const char * TailIncrement)
{
	int OldLength = Length();
	int len = strlen(TailIncrement);
	
	if (fNElements < OldLength + len + 1) Resize(OldLength + len + 1); // the 1 stands for the null char
	
	strcpy(fStore + OldLength, TailIncrement);
	fStore[Length()] = '\0'; // just to ensure the string contains a null character
}

TPZString TPZString::SubStr(const int start, const int end) const
{
	int len = Length();
	if (start > len || start > end)
	{
		TPZString newstring;
		return newstring;
	}
	
	int startpos = start, endpos = end, i;
	
	if (startpos < 0) startpos = 0;
	
	if (endpos > len) endpos = len - 1; // null ending character index
	
	TPZString newstring(endpos - startpos + 2);
	for (i = startpos; i <= endpos; i++) newstring[i - startpos] = fStore[i];
	
	newstring[newstring.NElements() - 1] = '\0';
	return newstring;
}

void TPZString::Empty()
{
	Resize(0);
}

void TPZString::Optimize()
{
	int len = Length();
	
	if (len + 1 < fNElements) Resize(len + 1);
}

void TPZString::SimplifyWhiteSpace()
{
	unsigned int i, newpos = 0;
	const unsigned int length = Length();
	char current = '0', next = '0';
	for (i = 0; i < length - 1; i++)
	{
		current = fStore[i];
		if (current == ' ')
		{
			if (i)
			{
				fStore[newpos] = ' ';
				newpos++;
			}
			next = fStore[i + 1];
			while (next == ' ')
			{
				i++;
				next = fStore[i + 1];
			}
			fStore[newpos] = next;
			newpos++;
		}
		else
		{
			if (newpos != i) fStore[newpos] = current;
			newpos++;
		}
	}
	fStore[newpos] = '\0';
}

int TPZString::Replace( const char * replace_str, const char *new_substr){
	string newstring(fStore);
	const int replace_len = strlen(replace_str);
	int  count = 0;
	string::size_type pos = newstring.find(replace_str, 0);
	while( pos!= string::npos) {
		newstring.replace(pos, replace_len, new_substr);
		pos = newstring.find(replace_str, 0);
		count++;
	}
	if(count) strcpy(fStore, newstring.c_str());
	return(count);
}

int TPZString::Find(const char * find_str){
	string fchr(fStore);
	return fchr.find(find_str, 0);
}
//--| PZ |----------------------------------------------------------------------
