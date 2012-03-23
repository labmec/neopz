/*
 *  TPZRefPatternDataBase.h
 *  NeoPZ
 *
 *  Created by caju on 12/17/09.
 *  Copyright 2009 LabMeC. All rights reserved.
 *
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

using namespace std;

class TPZRefPatternDataBase
{
	public:
	
		TPZRefPatternDataBase();
	    ~TPZRefPatternDataBase();
	
		int ReturnUniqueId();
	
		/**
		 * Read all refpatterns available in the given file
		 */
	void ReadRefPatternDBase(const std::string &filename);
	
		void ReadRefPatternDBase(std::ifstream &filein);
	
		void WriteRefPatternDBase(std::ofstream &fileout);
	
		/**
		 * Import a library of refinement patterns from the install directory
		 * Return the number of refpatterns imported
		 */
		int ImportRefPatterns();
	
		/**
		 * Import a library of refinement patterns from the given directory
		 * Return the number of refpatterns imported
		 */
		int ImportRefPatterns(std::string &Path);
	
		/**
		 * Retrieves the uniform refinement pattern for given element type
		 */ 
		TPZAutoPointer<TPZRefPattern> GetUniformRefPattern(MElementType type);
	
	  /**
		* Initialize the uniform refinement pattern from hard coaded data for an specific geometric element
		*/
		void InitializeUniformRefPattern(MElementType elType);
	
		void InitializeRefPatterns();

		/**
		 * Initialize the uniform refinement pattern from hard coaded data for all linear geometric elements
		 */
		void InitializeAllUniformRefPatterns();

		/**
		 * Insert the refinement pattern in the list of availabe refinement patterns
		 * assigns an Id to refPattern
		 */
		void InsertRefPattern(TPZAutoPointer<TPZRefPattern> & refpat);
	
		/**
		 * Check whether the refinement pattern already exists
		 */
		TPZAutoPointer<TPZRefPattern> FindRefPattern(TPZAutoPointer<TPZRefPattern> & refpat);
	
		TPZAutoPointer<TPZRefPattern> FindRefPattern(int id);
	
		TPZAutoPointer<TPZRefPattern> FindRefPattern(std::string name);
	
		/**
		 * Return the complete set of refinement patterns availabe
		 */
		const std::list< TPZAutoPointer<TPZRefPattern> > &RefPatternList(MElementType eltype);

		int NRefPatterns();
	
		void Print(std::ostream &out = std::cout);	
	
		void clear()
		{
			fElTypeRefPatterns.clear();
			fIdRefPatterns.clear();
		}
	
		protected:
	
		/** 
		 * Maps all refinement pattern objects in the mesh, indexed by refpattern element type
		 */
		std::map< MElementType , std::list< TPZAutoPointer<TPZRefPattern> > > fElTypeRefPatterns;
	
		/** 
		 * Maps all refinement pattern objects in the mesh, indexed by refpattern Id
		 */
		std::map< int , TPZAutoPointer<TPZRefPattern> > fIdRefPatterns;
	
		public:
};

extern TPZRefPatternDataBase gRefDBase;

#endif