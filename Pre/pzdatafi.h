/**
 * @file
 * @brief Contains the TDatafile class which defines the methods used to read a grid.
 */

#ifndef PZDATAFIHPP
#define PZDATAFIHPP

#include "pzerror.h"

#include <stdlib>
#include <fstream>


class TGeoGrid;
class TCompGrid;

/* **************************************************************************
 
 CLASS NAME :
 
 datafile
 
 PURPOSE OF THE CLASS :
 
 datafile is a virtual class from which classes with different data formats can
 be derived. The purpose of the datafile class is to facilitate the opening of an
 input stream and define the methods used to read a grid
 
 DEPENDENCIES :
 
 datafile.h :
 
 pzcondef.h
 pzerror.h
 pzcondef.h
 
 ERRORS :
 
 the length of the parameter fn in the constructor is not checked before
 being copied into filename. This may lead to system crashes...
 
 ANOMALIES :
 
 
 WARNING :
 
 If the input stream could not be opened, the file_error parameter is set to TRUE.
 It is up to the programmer of the derived class to check this data when writing
 the read method
 
 IMPROVEMENTS :
 
 The filename variable could be allocated dynamically to fit any filename length
 
 PROGRAMMERS :
 
 Jose Sergio Rodrigues Alves Filho
 Philippe Remy Bernard Devloo
 
 ************************************************************************** */

/** 
 * @ingroup pre
 * @brief datafile is a virtual class from which classes with different data formats can be derived. \ref pre "Getting Data"
 */
/**
 * The purpose of the datafile class is to facilitate the opening of an
 * input stream and defines the methods used to read a grid.
 */
class TDatafile 
{
protected:
	filebuf  fBuffer;
	istream  fData;
	char  fFilename[256];
	int	fFileError;
	//------------------------------------------------------------
public:
	
	TDatafile(char *fn) : fData(&fBuffer) {
		
		fFileError = 0;
		strcpy (fFilename, fn);
		if ( fBuffer.open( fn, ios::in) == 0) {
			pzerror << "\nERROR(TDatafile::const)-> Cannot open " << fn << "\n";
			//		pzerror.show();
			fFileError = 1;
		}
	}
	int Error() {return fFileError;}
	
	virtual ~TDatafile(){ fBuffer.close();}
	
	//virtual short Read (grid & malha) = 0;
	
	virtual short Read (TCompGrid & malha) = 0;
	
};

//------------------------------------------------------------

#endif
