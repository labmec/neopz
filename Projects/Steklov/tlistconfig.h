#ifndef TLISTCONFIGH
#define TLISTCONFIGH

/*
 *  tlistconfig.h
 *  SubstructEigen
 *
 *  Created by Philippe Devloo on 11/26/08.
 *  Copyright 2008 UNICAMP. All rights reserved.
 *
 */
#include "tconfig.h"
#include "pzgmesh.h"

TPZGeoMesh *LShape(TPZGeoMesh &moodel, TConfig &cf);

struct TListConfig
{
	enum configtype {LShapeC = 0, SquareC = 1, SquareSingularC = 2, SquareAdaptive = 3, SquareQuarterpoint = 4, CurvedTriangleSingular = 5, CurvedQuadrilateral = 6};
	
	std::list<TConfig> fList;
	configtype fConfigType;
	std::string fName;
	
	TListConfig(const std::string &groupname, configtype iconfig): fName(groupname)
	{
		switch(iconfig)
		{
			case 0:
				fConfigType = LShapeC;
				break;
			case 1:
				fConfigType = SquareC;
				break;
			case 2:
				fConfigType = SquareSingularC;
				break;
			case 3:
				fConfigType = SquareAdaptive;
				break;
			case 4:
				fConfigType = SquareQuarterpoint;
				break;
			case 5:
				fConfigType = CurvedTriangleSingular;
				break;
			case 6:
				fConfigType = CurvedQuadrilateral;
				break;
			default:
				fConfigType = LShapeC;
				std::cout << "Wrong configuration parameter iconfig " << iconfig << std::endl;
				break;
		}
	}
	
	std::string Name()
	{
		switch (fConfigType)
		{
			case LShapeC:
				return "LShape";
			case SquareC:
				return "Square";
			case SquareSingularC:
				return "SquareSingular";
			case SquareAdaptive:
				return "SquareAdaptive";
			case SquareQuarterpoint:
				return "SquareQuarterpoint";
			case CurvedTriangleSingular:
				return "CurvedTriangleSingular";
			case CurvedQuadrilateral:
				return "CurvedQuadrilateral";
			default:
				return "Undefined";
		}
	}
	
	TConfig &GenerateConfig(int data[4])
	{
		TConfig gen(Name(),data);
		fList.push_back(gen);
		return fList.back();
	}
	
	TPZGeoMesh *GenerateMesh(TPZGeoMesh &cp, TConfig &cf)
	{
		TPZGeoMesh *mesh = 0;
		switch (fConfigType) {
			case LShapeC:
				mesh = LShape(cp, cf);
				mesh->BuildConnectivity();
				break;
			case SquareC:
				mesh = Square(cp, cf);
				mesh->BuildConnectivity();
				break;
			case SquareSingularC:
				mesh = SquareSingular(cp,cf);
				mesh->BuildConnectivity();
			case SquareAdaptive:
			{
				TConfig cf2(cf);
				cf2.nref = cf2.nrefborder;
				cf2.porder = cf2.pborder;
				mesh = SquareSingular(cp,cf2);
				mesh->BuildConnectivity();
				break;
			}
			case SquareQuarterpoint:
			{
				mesh = SquareQuad(cp,cf);
				mesh->BuildConnectivity();
				break;
			}
			case CurvedTriangleSingular:
			{
				mesh = CurvedTriangleMesh(cp,cf);
//				mesh->BuildConnectivity();
				break;
			}
			case CurvedQuadrilateral:
			{
				mesh = QuadrilateralMesh(cp,cf);
				mesh->BuildConnectivity();
				break;
			}
			default:
				break;
		}
		return mesh;
	}
	
	std::string PrintFinal()
	{
		// Print the values
		std::stringstream sout;
		sout << fName << "VAL = {";
		std::list<TConfig>::iterator it = fList.begin();
		if(it != fList.end())
		{
			sout << it->AsList();
			it++;
		}
		for(;it != fList.end(); it++)
		{
			sout << "," << it->AsList();
		}
		sout << "};\n";
		/*
		// Create a list of all eigenvalues
		sout << fName << "EIG = {";
		it = fList.begin();
		if(it != fList.end())
		{
			sout << it->Eig();
			it++;
		}
		for(;it != fList.end(); it++)
		{
			sout << "," << it->Eig();
		}
		sout << "};\n";
		 */
		// Create a list of all stiffnesses
		sout << fName << "EK = {";
		it = fList.begin();
		if(it != fList.end())
		{
			sout << it->ek();
			it++;
		}
		for(;it != fList.end(); it++)
		{
			sout << "," << it->ek();
		}
		sout << "};\n";
		// Create a list of all masses
		sout << fName << "EM = {";
		it = fList.begin();
		if(it != fList.end())
		{
			sout << it->em();
			it++;
		}
		for(;it != fList.end(); it++)
		{
			sout << "," << it->em();
		}
		sout << "};\n";
		// Create a list of all rhs
		sout << fName << "RHS = {"; 
		it = fList.begin();
		if(it != fList.end())
		{
			sout << it->rhs();
			it++;
		}
		for(;it != fList.end(); it++)
		{
			sout << "," << it->rhs();
		}
		sout << "};\n";
#ifdef EXTCOMP
		sout << fName << "EIGEXT = {";
		it = fList.begin();
		if(it != fList.end())
		{
			sout << it->EigEXT();
			it++;
		}
		for(;it != fList.end(); it++)
		{
			sout << "," << it->EigEXT();
		}
		sout << "};";
#endif
		// Create a list of condensation
		sout << fName << "INDEX = {";
		it = fList.begin();
		if(it != fList.end())
		{
			sout << it->condense();
			it++;
		}
		for(;it != fList.end(); it++)
		{
			sout << "," << it->condense();
		}
		sout << "};\n";
		
		sout << fName << "EIG" << " = ListEigenvalues[" << fName << "EK" << "," << fName << "EM" << "];\n";
		sout << fName << "EIGCOND" << " = ListEigenvaluesCondense[" << fName << "EK" << "," << fName << "EM" << "," << fName << "INDEX" << "];\n";
		sout << fName << "EIGVectorsCOND" << " = ListEigenvectorsCondense[" << fName << "EK" << "," << fName << "EM" << "," << fName << "INDEX" << "];\n";
		sout << fName << "RHSDecompose" << " = ListRhsContribute[" << fName << "RHS" << "," << fName << "EIGVectorsCOND" << "," << fName << "EIGCOND" << "," << fName << "INDEX" << "];\n";
		// Export the eigenvalues file
		sout << "CreateDirectory[NotebookDirectory[] <> \"../Eigenvectors\"];\n";
		sout << "SetDirectory[NotebookDirectory[] <> \"../Eigenvectors\"];\n";

		it = fList.begin();
		int num = 0;
		for(;it != fList.end(); it++)
		{
			sout << "Export[\"" << it->ExportFileName(fName) << "\"," << fName << "EIGVectorsCOND[[" << num+1 << "]]];\n";
			num++;
		}
		return sout.str();
	}
};


#endif