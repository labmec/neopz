#ifndef TCONFIGH
#define TCONFIGH

#include <string>
#include <sstream>
/*
 *  tconfig.h
 *  SubstructEigen
 *
 *  Created by Philippe Devloo on 11/26/08.
 *  Copyright 2008 UNICAMP. All rights reserved.
 *
 */
struct TConfig
{
	int materialid;
	int materialborder1;
	int materialborder2;
	int materialborder3;
	int materialcorner;
	int materialcondense;
	int nref;
	int nrefborder;
	int porder;
	int pborder;
	std::string fMesh;
	TConfig(const std::string &meshname, int nrefinementsborder, int borderorder, int nrefinements, int pinternalorder)
	{
		materialid = 1;
		materialborder1 = -1;
		materialborder2 = -2;
		materialborder2 = -3;
		materialcorner = -4;
		materialcondense = -5;
		nref = nrefinements;
		nrefborder = nrefinementsborder;
		porder = pinternalorder;
		pborder = borderorder;
		fMesh = meshname;
	}
	TConfig(const std::string &meshname, int val[4])
	{
		materialid = 1;
		materialborder1 = -1;
		materialborder2 = -2;
		materialborder3 = -3;
		materialcorner = -4;
		materialcondense = -5;
		nrefborder = val[0];
		pborder = val[1];
		nref = val[2];
		porder = val[3];
		fMesh = meshname;
	}
	
	TConfig(const TConfig &cf) : materialid(cf.materialid), materialborder1(cf.materialborder1),materialborder2(cf.materialborder2),materialborder3(cf.materialborder3),materialcorner(cf.materialcorner), materialcondense(cf.materialcondense), nrefborder(cf.nrefborder),
	pborder(cf.pborder),nref(cf.nref), porder(cf.porder), fMesh(cf.fMesh)
	{
		
	}
	TConfig &operator=(const TConfig &cf)
	{
		materialid = cf.materialid;
		materialborder1 = cf.materialborder1;
		materialborder2 = cf.materialborder2;
		materialborder3 = cf.materialborder3;
		materialcorner = cf.materialcorner;
		materialcondense = cf.materialcondense;
		nrefborder = cf.nrefborder;
		pborder = cf.pborder;
		nref = cf.nref;
		porder = cf.porder;
		fMesh = cf.fMesh;
		return *this;
	}
	std::string cfg()
	{
		std::stringstream sout;
		sout << nrefborder << pborder << nref << porder ;
		return sout.str();
	}
	std::string Eig()
	{
		std::stringstream sout;
		sout << fMesh << "Eig" << cfg();
		return sout.str();
	}
	std::string ek()
	{
		std::stringstream sout;
		sout << fMesh << "ek" << cfg();
		return sout.str();
	}
	std::string em()
	{
		std::stringstream sout;
		sout << fMesh << "em" << cfg();
		return sout.str();
	}
	std::string rhs()
	{
		std::stringstream sout;
		sout << fMesh << "rhs" << cfg();
		return sout.str();
	}
	std::string condense()
	{
		std::stringstream sout;
		sout << fMesh << "Index" << cfg();
		return sout.str();
	}
	std::string emEXT()
	{
		std::stringstream sout;
		sout << fMesh << "emEXT" << cfg();
		return sout.str();
	}
	std::string EigEXT()
	{
		std::stringstream sout;
		sout << fMesh << "EigEXT" << cfg();
		return sout.str();
	}
	std::string ExportFileName(const std::string &root)
	{
		std::stringstream sout2;
		sout2 << "Eig_" << root << cfg() << ".dat";
		return sout2.str();
		
	}
	std::string EigValDecomposition()
	{
		std::stringstream sout;
		sout << fMesh << "EigCoef" << cfg();
		return sout.str();
	}
	std::string EigenvectorFileName(const std::string &root)
	{
		std::stringstream sout2;
		sout2 << "VTK/Eig_" << root << cfg() << ".vtk";
		return sout2.str();
	}
	std::string SolutionFileName(const std::string &root)
	{
		std::stringstream sout2;
		sout2 << "VTK/Sol_" << root << cfg() << ".vtk";
//		sout2 << "VTK/Sol_" << root << cfg() << ".dx";
		return sout2.str();
	}
	
	std::string Eigenvalues()
	{
		std::stringstream sout;
		sout << Eig() << " = Reverse[Eigenvalues[{" << ek() << "," << em() << "}]];\n";
#ifdef EXTCOMP		
		sout << EigEXT() << " = Reverse[Eigenvalues[{" << ek() << "," << emEXT() << "}]];\n";
#endif
		return sout.str();
		
	}
	std::string AsList()
	{
		std::stringstream sout;
		sout << "{" << nrefborder << "," << pborder << "," << nref << "," << porder << "}";
		return sout.str();
	}
};

#endif