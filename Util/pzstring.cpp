/** @file pzstring.cc */
// $Id: pzstring.cpp,v 1.2 2003-09-26 18:30:26 erick Exp $

#include "pzstring.h"

#include <stdio.h>

TPZString::TPZString()
{
   // NOTHING TO DO HERE!
}
/*
TPZString::TPZString(const char * source)
{
   int len=strlen(source);
   Resize(len+1);
   strcpy(fStore, source);

   // although the TPZStack class already stores the length, the null
   // character is necessary to export string as a null character
   // ended string
   fStore[len]='\0';
}
*/
TPZString::TPZString(char const * source)
{
   int len=strlen(source);
   Resize(len+1);
   strcpy(fStore, source);

   // although the TPZStack class already stores the length, the null
   // character is necessary to export string as a null character
   // ended string
   fStore[len]='\0';
}


TPZString::TPZString(const int size)
{
   Resize(size);
}

TPZString::TPZString(const char chr)
{
   Resize(2);
   fStore[0]=chr;
   fStore[1]='\0';
}

TPZString TPZString::operator+(const char * increment)const
{
   TPZString newstring(*this);
   newstring.Append(increment);
   return newstring;
}

TPZString TPZString::operator+(const TPZString & increment)const
{
   TPZString newstring(*this);
   newstring.Append(increment);
   return newstring;
}

void TPZString::operator+=(const char increment)
{
   this->Append(increment);
}

void TPZString::operator+=(const char * increment)
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

TPZString::operator const char *()const
{
   return Str();
}

const char* TPZString::Str()const
{
   if(fNElements>0)
   {
      return fStore;
   }

   return NULL;
}

int TPZString::Length()const
{
   if(fNElements==0)return 0;
   return strlen(fStore);
}

void TPZString::Append(const char TailIncrement)
{	
   int len=Length();
   // if the allocated memory is less than the length +1 (null char)
   // +1 (new char), reallocates it
   if(fNElements<len+2)
   {
      Push('\0');
   }
   else
   {
      fStore[len+1]='\0';
   }
   
   fStore[len]=TailIncrement;
}
	
void TPZString::Append(const char * TailIncrement)
{
   int OldLength=Length();
   int len=strlen(TailIncrement);

   if(fNElements<OldLength+len+1)Resize(OldLength+len+1); // the 1 stands for the null char
	
   strcpy(fStore+OldLength, TailIncrement);
   fStore[Length()]='\0'; // just to ensure the string contains a null character
}

TPZString TPZString::SubStr(const int start,const int end)const
{
   int len=Length();
   if(start>len || start>end)
   {
      TPZString newstring;
      return newstring;
   }

   int startpos = start, endpos = end, i;

   if(startpos<0)startpos=0;

   if(endpos>len)endpos=len-1; // null ending character index

   TPZString newstring(endpos-startpos+2);
   for(i=startpos; i<=endpos; i++)newstring[i-startpos]=fStore[i];

   newstring[newstring.NElements()-1]='\0';
   return newstring;
}

void TPZString::Empty()
{
   Resize(0);
}

void TPZString::Optimize()
{
   int len = Length();

   if(len+1<fNElements)
      Resize(len+1);
}

//--| PZ |----------------------------------------------------------------------
