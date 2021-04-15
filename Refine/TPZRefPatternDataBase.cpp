/**
 * @file
 * @brief Contains the implementation of the TPZRefPatternDataBase methods. 
 */
/*
 *  Created by caju on 12/17/09.
 *  Copyright 2009 LabMeC. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include "pz_config.h"
#endif

#include "TPZRefPatternDataBase.h"
#include "pzgeoelside.h"
#include "pzlog.h"
#include <fstream>

#ifdef PZ_LOG
static TPZLogger logger("pz.refpattern.TPZRefPatternDataBase");
#endif

#ifdef BORLAND
#include <io.h>
#include <fcntl.h>
#endif

using namespace std;

TPZRefPatternDataBase gRefDBase;

//.........................................................................................................................................
TPZRefPatternDataBase::TPZRefPatternDataBase()
{

}

//.........................................................................................................................................
TPZRefPatternDataBase::~TPZRefPatternDataBase()
{
    this->clear();
}

void TPZRefPatternDataBase::clear()
{
    fElTypeRefPatterns.clear();
    fIdRefPatterns.clear();
}


//.........................................................................................................................................
int TPZRefPatternDataBase::ReturnUniqueId()
{
	int uniqueId = TPZRefPattern::fNonInitializedId;
	std::set<int> Ids;
	Ids.clear();
	
	std::map< MElementType , std::list< TPZAutoPointer<TPZRefPattern> > >::iterator mapIt;
	std::list< TPZAutoPointer<TPZRefPattern> >::iterator listIt;
	for(mapIt = fElTypeRefPatterns.begin(); mapIt != fElTypeRefPatterns.end(); mapIt++)
	{
		for(listIt = mapIt->second.begin(); listIt != mapIt->second.end(); listIt++)
		{
			Ids.insert((*listIt)->Id());
		}
	}
	
	std::set<int>::iterator it;
	
	int nRefPatterns = gRefDBase.NRefPatterns();
	for(int id = 0; id <= nRefPatterns; id++)
	{
		it = Ids.find(id);
		if(it == Ids.end())
		{
			uniqueId = id;
			break;
		}
	}
	
	return uniqueId;
}

//.........................................................................................................................................
void TPZRefPatternDataBase::ReadRefPatternDBase(const std::string &filename)
{
	std::ifstream in(filename.c_str());
	ReadRefPatternDBase(in);
}	
//.........................................................................................................................................
void TPZRefPatternDataBase::ReadRefPatternDBase(std::ifstream &filename)
{
    fElTypeRefPatterns.clear();
    fIdRefPatterns.clear();
    
	int nRefpatterns;
	filename >> nRefpatterns;
	for(int i = 0; i < nRefpatterns; i++)
	{
		TPZAutoPointer<TPZRefPattern> refP = new TPZRefPattern;
		refP->ReadAndCreateRefinementPattern(filename);
		
        MElementType eltype = refP->Element(0)->Type();
        fElTypeRefPatterns[eltype].push_back(refP);
        fIdRefPatterns[refP->Id()] = refP;
	}
}

//.........................................................................................................................................
void TPZRefPatternDataBase::WriteRefPatternDBase(std::ofstream &filename)
{
	filename << NRefPatterns() << std::endl ;
	
	std::map< MElementType, std::list< TPZAutoPointer<TPZRefPattern> > >::iterator iterMap;
	for (iterMap = fElTypeRefPatterns.begin(); iterMap != fElTypeRefPatterns.end(); iterMap++)
	{
		std::list<TPZAutoPointer<TPZRefPattern> > &refpList = (*iterMap).second;
		std::list< TPZAutoPointer<TPZRefPattern> >::iterator refp;
		for(refp = refpList.begin(); refp != refpList.end(); refp++)
		{
			TPZAutoPointer<TPZRefPattern> &refpattern = *refp;
			refpattern->WritePattern(filename);
		}
	}
}

//.........................................................................................................................................
int TPZRefPatternDataBase::ImportRefPatterns(int maxdim)
{
	std::string DefaulPath = PZSOURCEDIR;
	
	DefaulPath += "/Refine/RefPatterns";
//#define StartPathDefined 1;
	
	return ImportRefPatterns(DefaulPath, maxdim);
}

//.........................................................................................................................................
int TPZRefPatternDataBase::ImportRefPatterns(std::string &Path, int maxdim)
{
	std::string bar = "/";
	
	std::string FileTypes ("*.rpt");
	std::string Command = std::string("ls -1 ") + Path + bar + FileTypes;
	std::cout << "Generated command: " << Command.c_str() << std::endl;
	FILE   *fp;
	
#ifdef VC
	fp = _popen(Command.c_str(), "r");
#else
	#ifndef BORLAND
		fp = popen(Command.c_str(), "r");
	#else
		fp = (FILE *)open(Command.c_str(), O_RDONLY );
	#endif
#endif
	
	if (!fp)
	{
		return -1;
	}
	
	int count = 0;
	char psBuffer[1024];
	
	while( !feof( fp ) )
	{
		if( fgets(psBuffer, sizeof(psBuffer), fp ) != NULL )
		{
			if (psBuffer[strlen(psBuffer)-1] == '\n') 
			{
				psBuffer[strlen(psBuffer)-1] = 0;
			}
			std::cout << "Reading refinement patern file : " << psBuffer << std::endl;
			std::string filref(psBuffer);
			
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filref);
            
            if (refpat->fRefPatternMesh.Dimension() > maxdim) {
                std::cout << "skipped\n";
                continue;
            }
			
			if(!this->FindRefPattern(refpat))
			{
				this->InsertRefPattern(refpat);
			}
			refpat->InsertPermuted();
			
			count++;
		}
	}

#ifndef BORLAND
#ifdef VC
	_pclose(fp);
#else
	pclose(fp);
#endif
#else
	close((int)fp);
#endif
	
	return count;
}

//.........................................................................................................................................
TPZAutoPointer<TPZRefPattern> TPZRefPatternDataBase::GetUniformRefPattern(MElementType type)
{
	std::list< TPZAutoPointer<TPZRefPattern> > ElTypeList = RefPatternList(type);
	std::list< TPZAutoPointer<TPZRefPattern> >::iterator it;
	
	std::string unifRefName;
	switch (type)
	{
		case (EOned) :
		{
			unifRefName = "UnifLin";
			break;
		}
		case (ETriangle) :
		{
			unifRefName = "UnifTri";			
			break;
		}
		case (EQuadrilateral) :
		{
			unifRefName = "UnifQua";
			break;
		}
		case (ETetraedro) :
		{
			unifRefName = "UnifTet";
			break;
		}
		case (EPiramide) :
		{
			unifRefName = "UnifPyr";
			break;
		}
		case (EPrisma) :
		{
			unifRefName = "UnifPri";
			break;
		}
		case (ECube) :
		{
			unifRefName = "UnifHex";
			break;
		}
		default:
		{
			//none
			break;
		}
	}
	TPZAutoPointer < TPZRefPattern > UnifRefPat = NULL;
	for(it = ElTypeList.begin(); it != ElTypeList.end(); it++)
	{
		std::string compareRefrName = (*it)->Name();
		if( compareRefrName == unifRefName )
		{
			UnifRefPat = (*it);
			break;
		}
	}
	
#ifdef PZDEBUG
	if(!UnifRefPat && type != EPoint)
	{
		std::cout << "Uniform refpattern " << unifRefName << " was not initialized!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
	}
#endif
	
	return UnifRefPat;
}

//.........................................................................................................................................
void TPZRefPatternDataBase::InitializeUniformRefPattern(MElementType elType)
{
	switch (elType)
	{
		case 0://EPoint
		{
			break;
		}
		case 1://EOned
		{
			std::cout << "inserting uniform refpattern: line\n";
			char buf[] =
			"3     3  "
			"-50       UnifLin	"
			"-1.     0.     0. "
			" 1.     0.     0. "
			" 0.     0.     0. "
			" 1     2     0     1 "
			" 1     2     0     2 "
			" 1     2     2     1 ";
			std::istringstream str(buf);
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
			TPZAutoPointer<TPZRefPattern> refpatFound = FindRefPattern(refpat);
			if(!refpatFound)
			{
				InsertRefPattern(refpat);
			}
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
			
			break;
		}
		case 2://ETriangle
		{
			std::cout << "inserting uniform refpattern: triangle\n";
			char buf[] =
			"6     5  "
			"-50       UnifTri	"
			"0.     0.     0. "
			"1.     0.     0. "
			"0.     1.     0. "
			"0.5    0.     0. "
			"0.5    0.5    0. "
			"0.     0.5     0. "
			"2     3     0     1     2 "
			"2     3     0     3     5 "
			"2     3     3     1     4 "
			"2     3     5     4     2 "
			"2     3     4     5     3 ";
			std::istringstream str(buf);
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
			TPZAutoPointer<TPZRefPattern> refpatFound = FindRefPattern(refpat);
			if(!refpatFound)
			{
				InsertRefPattern(refpat);
			}
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
			
			break;
		}
		case 3://EQuadrilateral
		{
			std::cout << "inserting uniform refpattern: quadrilateral\n";
			char buf[] =
			"9     5  "
			"-50       UnifQua	"
			"-1.    -1.     0. "
			" 1.    -1.     0. "
			" 1.     1.     0. "
			"-1.     1.     0.	"
			" 0.    -1.     0. "
			" 1.     0.     0. "
			" 0.     1.     0. "
			"-1.     0.     0. "
			" 0.     0.     0. "
			" 3     4     0     1     2     3 "
			" 3     4     0     4     8     7 "
			" 3     4     4     1     5     8 "
			" 3     4     8     5     2     6 "
			" 3     4     7     8     6     3 ";
			std::istringstream str(buf);
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
			TPZAutoPointer<TPZRefPattern> refpatFound = FindRefPattern(refpat);
			if(!refpatFound)
			{
				InsertRefPattern(refpat);
			}
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
			
			break;
		}
		case 4://ETetraedro
		{
			std::cout << "inserting uniform refpattern: tetrahedra\n";
			char buf[] =
			"10     7"
			"-50       UnifTet	"
			"0.     0.     0. "
			"1.     0.     0. "
			"0.     1.     0. "
			"0.     0.     1. "
			"0.5    0.     0. "
			"0.5    0.5    0. "
			"0.     0.5    0. "
			"0.     0.     0.5 "
			"0.5    0.     0.5 "
			"0.     0.5    0.5 "
			"4     4     0     1     2     3 "
			"4     4     0     4     6     7 "
			"4     4     4     1     5     8 "
			"4     4     6     5     2     9 "
			"4     4     7     8     9     3 "
			"5     5     4     8     9     6     7 "
			"5     5     8     4     6     9     5 ";
			std::istringstream str(buf);
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
			TPZAutoPointer<TPZRefPattern> refpatFound = FindRefPattern(refpat);
			if(!refpatFound)
			{
				InsertRefPattern(refpat);
			}
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
			
			break;
		}
		case 5://EPiramide
		{
			std::cout << "inserting uniform refpattern: pyramid\n";
			char buf[] =
			"14     11	"
			"-50       UnifPyr	"
			"-1.    -1.     0. "
			" 1.    -1.     0. "
			" 1.     1.     0. "
			"-1.     1.     0. "
			" 0.     0.     1. "
			" 0.    -1.     0. "
			" 1.     0.     0. "
			" 0.     1.     0. "
			"-1.     0.     0. "
			"-0.5   -0.5    0.5 "
			" 0.5   -0.5    0.5 "
			" 0.5    0.5    0.5 "
			"-0.5    0.5    0.5 "
			" 0.     0.     0. "
			" 5     5     0     1     2     3     4 "
			" 5     5     0     5    13     8     9	"
			" 5     5     5     1     6    13    10	"
			" 5     5    13     6     2     7    11	"
			" 5     5     8    13     7     3    12	"
			" 5     5     9    10    11    12     4	"
			" 5     5     9    12    11    10    13	"
			" 4     4     9     5    13    10	"
			" 4     4    10     6    13    11	"
			" 4     4    12    13     7    11	"
			" 4     4     8    13    12     9	";
			std::istringstream str(buf);
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
			TPZAutoPointer<TPZRefPattern> refpatFound = FindRefPattern(refpat);
			if(!refpatFound)
			{
				InsertRefPattern(refpat);
			}
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
			
			break;
		}
		case 6://EPrisma
		{
			std::cout << "inserting uniform refpattern: prism\n";
			char buf[] =
			"18     9  "
			"-50       UnifPri 	"
			" 0.     0.    -1.	"
			" 1.     0.    -1.	"
			" 0.     1.    -1.	"
			" 0.     0.     1.	"
			" 1.     0.     1.	"
			" 0.     1.     1.	"
			" 0.5    0.    -1.	"
			" 0.5    0.5   -1.	"
			" 0.     0.5   -1.	"
			" 0.     0.     0.	"
			" 1.     0.     0.	"
			" 0.     1.     0.	"
			" 0.5    0.     1.	"
			" 0.5    0.5    1.	"
			" 0.     0.5    1.	"
			" 0.5    0.     0.	"
			" 0.5    0.5    0.	"
			" 0.     0.5    0.	"
			" 6     6     0     1     2     3     4     5	"
			" 6     6     0     6     8     9    15    17	"
			" 6     6     6     1     7    15    10    16	"
			" 6     6     8     7     2    17    16    11	"
			" 6     6    17    16    15     8     7     6	"
			" 6     6     9    15    17     3    12    14	"
			" 6     6    15    10    16    12     4    13	"
			" 6     6    17    16    11    14    13     5	"
			" 6     6    14    13    12    17    16    15	";
			std::istringstream str(buf);
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
			TPZAutoPointer<TPZRefPattern> refpatFound = FindRefPattern(refpat);
			if(!refpatFound)
			{
				InsertRefPattern(refpat);
			}
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
			
			break;
		}
		case 7://ECube
		{
			std::cout << "inserting uniform refpattern: hexahedre\n";
			char buf[] =
			"27     9  "
			"-50     UnifHex  "
			"-1.    -1.    -1.  "
			" 1.    -1.    -1.  "
			" 1.     1.    -1.  "
			"-1.     1.    -1.  "
			"-1.    -1.     1.  "
			" 1.    -1.     1.  "
			" 1.     1.     1.  "
			"-1.     1.     1.  "
			" 0.    -1.    -1.  "
			" 1.     0.    -1.  "
			" 0.     1.    -1.  "
			"-1.     0.    -1.  "
			"-1.    -1.     0.  "
			" 1.    -1.     0.  "
			" 1.     1.     0.  "
			"-1.     1.     0.  "
			" 0.    -1.     1.  "
			" 1.     0.     1.  "
			" 0.     1.     1.  "
			"-1.     0.     1.  "
			" 0.     0.    -1.  "
			" 0.    -1.     0.  "
			" 1.     0.     0.  "
			" 0.     1.     0.  "
			"-1.     0.     0.  "
			" 0.     0.     1.  "
			" 0.     0.     0.  "
			"7     8     0     1     2     3     4     5     6     7  "
			"7     8     0     8    20    11    12    21    26    24  "
			"7     8     8     1     9    20    21    13    22    26  "
			"7     8    20     9     2    10    26    22    14    23  "
			"7     8    11    20    10     3    24    26    23    15  "
			"7     8    12    21    26    24     4    16    25    19  "
			"7     8    21    13    22    26    16     5    17    25  "
			"7     8    26    22    14    23    25    17     6    18  "
			"7     8    24    26    23    15    19    25    18     7  ";
			std::istringstream str(buf);
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
			TPZAutoPointer<TPZRefPattern> refpatFound = FindRefPattern(refpat);
			if(!refpatFound)
			{
				InsertRefPattern(refpat);
			}
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
			
			break;
		}
		default:
		{
			cout << "Cant generate uniform refpattern because MElementType was not found on " << __PRETTY_FUNCTION__ << endl;
			DebugStop();
		}
	}
}

//.........................................................................................................................................
void TPZRefPatternDataBase::InitializeAllUniformRefPatterns()
{
	InitializeUniformRefPattern(EOned);
	InitializeUniformRefPattern(ETriangle);
	InitializeUniformRefPattern(EQuadrilateral);
	InitializeUniformRefPattern(ETetraedro);
	InitializeUniformRefPattern(EPiramide);
	InitializeUniformRefPattern(EPrisma);
	InitializeUniformRefPattern(ECube);
}

//.........................................................................................................................................
void TPZRefPatternDataBase::InsertRefPattern(TPZAutoPointer<TPZRefPattern> &refpat)
{
	if (!refpat)
	{
		PZError << "TPZGeoMesh::InsertRefPattern ERROR : NULL refinement pattern! " << std::endl;
		return;
	}
	
	MElementType eltype = refpat->Element(0)->Type();
	
	if(!refpat->IdInitialized())
	{
		int id = ReturnUniqueId();
		refpat->SetId(id);
	}
	
	else//the refpattern already have an Id initialized
	{
		std::map< int , TPZAutoPointer<TPZRefPattern> >::iterator itid = fIdRefPatterns.find(refpat->Id());
		if(itid != fIdRefPatterns.end())//the Id is already in use by another refpattern
		{
			int id = ReturnUniqueId();
			refpat->SetId(id);
		}
	}
	
	fElTypeRefPatterns[eltype].push_back(refpat);
	fIdRefPatterns[refpat->Id()] = refpat;
}

//.........................................................................................................................................
TPZAutoPointer<TPZRefPattern> TPZRefPatternDataBase::FindRefPattern(TPZAutoPointer<TPZRefPattern> &refpat)
{
	TPZAutoPointer<TPZRefPattern> NullRefPat;
	if(!refpat )
	{
		return NullRefPat;
	}
	
	std::map< int, TPZAutoPointer<TPZRefPattern> >::iterator it;
	for(it = fIdRefPatterns.begin(); it != fIdRefPatterns.end(); it++)
	{
		if( (TPZRefPattern &)refpat == it->second)
		{
			return it->second;
		}
	}
	
	//else
	return NullRefPat;
}

//.........................................................................................................................................
TPZAutoPointer<TPZRefPattern>  TPZRefPatternDataBase::FindRefPattern(int id)
{
	TPZAutoPointer<TPZRefPattern> refPatReturned;
	std::map< int , TPZAutoPointer<TPZRefPattern> >::iterator it;
	
	it = fIdRefPatterns.find(id);
	
	if(it != fIdRefPatterns.end())
	{
		refPatReturned = it->second;
	}
	
	return refPatReturned;
}

//.........................................................................................................................................
TPZAutoPointer<TPZRefPattern> TPZRefPatternDataBase::FindRefPattern(std::string name)
{
	TPZAutoPointer<TPZRefPattern> refp = NULL;
	
	std::map< int , TPZAutoPointer<TPZRefPattern> >::iterator it;
	for(it = 	fIdRefPatterns.begin(); it != fIdRefPatterns.end(); it++)
	{
		if( (it->second)->fName == name)
		{
			refp = it->second;
			break;
		}
	}
	
	return refp;
}

//.........................................................................................................................................
const std::list< TPZAutoPointer<TPZRefPattern> > &TPZRefPatternDataBase::RefPatternList(MElementType eltype)
{
	return fElTypeRefPatterns[eltype];
}

//.........................................................................................................................................
int TPZRefPatternDataBase::NRefPatterns()
{
	int nRefPatterns = fIdRefPatterns.size();
	
	return nRefPatterns;
}

//.........................................................................................................................................
void TPZRefPatternDataBase::Print(std::ostream &out)
{
	std::map< MElementType , std::list< TPZAutoPointer<TPZRefPattern> > >::iterator it;
	
	for(it = fElTypeRefPatterns.begin(); it != fElTypeRefPatterns.end(); it++)
	{
		out << "======================================================\n";
		out << "**************************************************\n";
		out << "Element type: " << MElementType_Name(it->first) << std::endl;
		out << "**************************************************";
		
		std::list< TPZAutoPointer<TPZRefPattern> >::iterator itt;
		
		for(itt = it->second.begin(); itt != it->second.end(); itt++)
		{
			(*itt)->Print(out);
		}
	}
}
