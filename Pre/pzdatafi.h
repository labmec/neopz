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

/** 
 * @ingroup pre
 * @brief TDatafile is a virtual class from which classes with different data formats can be derived. \ref pre "Getting Data"
 * @author Jose Sergio Rodrigues Alves Filho
 * @author Philippe Remy Bernard Devloo
 */
/**
 * The purpose of the datafile class is to facilitate the opening of an
 * input stream and defines the methods used to read a grid.
 *
 * \par ERRORS:
 * The length of the parameter fn in the constructor is not checked before
 * being copied into filename. This may lead to system crashes...
 *
 * \warning
 * If the input stream could not be opened, the file_error parameter is set to TRUE. \n
 * It is up to the programmer of the derived class to check this data when writing the read method.
 *
 * \par IMPROVEMENTS:
 * The filename variable could be allocated dynamically to fit any filename length (Who? When?)
 */
class TDatafile 
{
protected:
	/// File from import data
	filebuf  fBuffer;
	/// Data
	istream  fData;
	/// File name 
	char  fFilename[256];
	/// Flag to report if the file is opening
	int	fFileError;

public:
	
	/// Constructor with filename of the data to import
	TDatafile(char *fn) : fData(&fBuffer) {
		fFileError = 0;
		strcpy (fFilename, fn);
		if ( fBuffer.open( fn, ios::in) == 0) {
			pzerror << "\nERROR(TDatafile::const)-> Cannot open " << fn << "\n";
			//		pzerror.show();
			fFileError = 1;
		}
	}
	/// Returns if could opening or not the data file
	int Error() {return fFileError;}
	/// Destructor
	virtual ~TDatafile(){ fBuffer.close();}

	/// To import data into computational mesh
	virtual short Read (TCompGrid & malha) = 0;

};

#endif
