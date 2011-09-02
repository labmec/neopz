/**
 * @file
 * @brief Contains the TPZRefPatternDataBase class which defines data base of patterns.
 */
/*
 *  TPZRefPatternDataBase.h
 *  NeoPZ
 *
 *  Created by caju on 12/17/09.
 *  Copyright 2009 LabMeC. All rights reserved.
 */

#ifndef TPZREFPATTERNDATABASEH
#define TPZREFPATTERNDATABASEH

#include <iostream>
#include <string>
#include <map>
#include <list>

#include "pzeltype.h"
#include "tpzautopointer.h"
#include "TPZRefPattern.h"

/**
 * @ingroup refine
 * @brief Defines data base of patterns. \ref refine "Refine"
 */
class TPZRefPatternDataBase
{
public:
	
	TPZRefPatternDataBase();
	~TPZRefPatternDataBase();
	
	int ReturnUniqueId();
	
	/** @brief Read all refpatterns available in the given file */
	void ReadRefPatternDBase(const std::string &filename);
	
	void ReadRefPatternDBase(std::ifstream &filein);
	
	void WriteRefPatternDBase(std::ofstream &fileout);
	
	/**
	 * @brief Import a library of refinement patterns from the install directory
	 * @return Return the number of refpatterns imported
	 */
	int ImportRefPatterns();
	
	/**
	 * @brief Import a library of refinement patterns from the given directory
	 * @return Return the number of refpatterns imported
	 */
	int ImportRefPatterns(std::string &Path);
	
	/** @brief Retrieves the uniform refinement pattern for given element type */ 
	TPZAutoPointer<TPZRefPattern> GetUniformRefPattern(MElementType type);
	
	/** @brief Initialize the uniform refinement pattern from hard coaded data for an specific geometric element */
	void InitializeUniformRefPattern(MElementType elType);
	
	void InitializeRefPatterns();
	
	/** @brief Initialize the uniform refinement pattern from hard coaded data for all linear geometric elements */
	void InitializeAllUniformRefPatterns();
	
	/** @brief Insert the refinement pattern in the list of availabe refinement patterns assigns an Id to refPattern */
	void InsertRefPattern(TPZAutoPointer<TPZRefPattern> & refpat);
	
	/** @brief Check whether the refinement pattern already exists */
	TPZAutoPointer<TPZRefPattern> FindRefPattern(TPZAutoPointer<TPZRefPattern> & refpat);
	
	TPZAutoPointer<TPZRefPattern> FindRefPattern(int id);
	
	TPZAutoPointer<TPZRefPattern> FindRefPattern(std::string name);
	
	/** @brief Return the complete set of refinement patterns availabe */
	const std::list< TPZAutoPointer<TPZRefPattern> > &RefPatternList(MElementType eltype);
	
	int NRefPatterns();
	
	void Print(std::ostream &out = std::cout);	
	
	void clear()
	{
		fElTypeRefPatterns.clear();
		fIdRefPatterns.clear();
	}
	
protected:
	
	/** @brief Maps all refinement pattern objects in the mesh, indexed by refpattern element type */
	std::map< MElementType , std::list< TPZAutoPointer<TPZRefPattern> > > fElTypeRefPatterns;
	
	/** @brief Maps all refinement pattern objects in the mesh, indexed by refpattern Id */
	std::map< int , TPZAutoPointer<TPZRefPattern> > fIdRefPatterns;

};

/// External variable to data base of patterns
extern TPZRefPatternDataBase gRefDBase;

#endif
