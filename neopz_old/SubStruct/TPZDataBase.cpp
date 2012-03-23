/*
 *  TPZDataBase.cpp
 *  SubStruct
 *
 *  Created by Bandit on 8/9/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZDataBase.h"

TPZDataBase::TPZDataBase()
{
	
}

void TPZDataBase::Read(std::string &FileName)
{
	
	SetFileName(FileName);
	std::ifstream ReadFile;
	ReadFile.open(fFileName.c_str());
	
	int nlinhas = NumberOfLines();
	
	fTimes.Resize(nlinhas-1);

	char buf[1024];	
	ReadFile.getline(buf, 1024);							// passando pelo cabecalho
	
	for (int i = 0; i < nlinhas - 1; i++)				// armazena sem contar o cabecalho
	{				
		fTimes[i].ReadLine(ReadFile);
	}
	
	//std::cout << fTimes[0].ft0sub << "\n" << fTimes[1].ft0sub << "\n" << fTimes[1].ft1comput << "\n"<< fTimes[2].ft0sub << std::endl;
	
}

void TPZDataBase::SetFileName(std::string &FileName)
{
	fFileName = FileName;
}

int TPZDataBase::NumberOfLines()
{
	std::ifstream Read;
	Read.open(fFileName.c_str());
	int nlinhas = -1;											// ele sempre conta uma linha a mais (por causa do endl no final)
	
	{
		while (Read)											// conta o numero de linhas
		{
			char buf[1024];
			Read.getline(buf, 1024);		
			nlinhas++;
		}
	}
	return nlinhas;
	

}

