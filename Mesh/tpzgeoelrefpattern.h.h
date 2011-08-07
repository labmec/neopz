/**
 * @file
 * @brief Contains the implementation of the TPZGeoElRefPattern methods.
 */
//
// C++ Interface: tpzgeoelrefpattern
//
// Description:
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef TPZGEOELREFPATTERN_H_H
#define TPZGEOELREFPATTERN_H_H

#include "pzlog.h"

#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"

template <class TGeo>
void TPZGeoElRefPattern<TGeo>::Read(TPZStream &str, void *context)
{
	TPZGeoElRefLess<TGeo>::Read(str, context);
	//TPZGeoMesh *gmesh = (TPZGeoMesh *) context;
	int refpatternindex;
	str.Read(&refpatternindex, 1);
	if(refpatternindex != -1)
	{
		const std::list< TPZAutoPointer<TPZRefPattern> > &RefPatternList = gRefDBase.RefPatternList(this->Type());
		std::list< TPZAutoPointer<TPZRefPattern> >::const_iterator it;
		
		for(it = RefPatternList.begin(); it != RefPatternList.end(); it++)
		{
			if((*it)->Id() == refpatternindex)
			{
				break;
			}
		}
		
		if(it != RefPatternList.end()) fRefPattern =*it;
	}
	TPZSaveable::ReadObjects(str, this->fSubEl);
}

template <class TGeo>
void TPZGeoElRefPattern<TGeo>::Write(TPZStream &str, int withclassid){
	TPZGeoElRefLess<TGeo>::Write(str, withclassid);
	int refpatternindex = -1;
	if(fRefPattern) refpatternindex = fRefPattern->Id();
	str.Write(&refpatternindex, 1);
	TPZSaveable::WriteObjects(str, this->fSubEl);
}

// template<class TGeo>
// TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index) :
// TPZGeoElRefLess<TGeo>(nodeindices,matind,mesh,index)
// {
// }

template <class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZGeoMesh &DestMesh, const TPZGeoElRefPattern<TGeo> &cp):TPZGeoElRefLess<TGeo>(DestMesh,cp),
fRefPattern(cp.fRefPattern) {
	this->fSubEl = cp.fSubEl;
}

template <class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TGeo>::Clone(TPZGeoMesh &DestMesh) const{
	return new TPZGeoElRefPattern<TGeo>(DestMesh, *this);
}


template <class TGeo>
TPZGeoElRefPattern<TGeo>::TPZGeoElRefPattern(TPZGeoMesh &DestMesh,
											 const TPZGeoElRefPattern<TGeo> &cp,
											 std::map<int,int> &gl2lcNdMap,
											 std::map<int,int> &gl2lcElMap):
TPZGeoElRefLess<TGeo>(DestMesh,cp,gl2lcNdMap,gl2lcElMap),
fRefPattern ( cp.fRefPattern )
{
	int i;
	for (i=0;i<cp.fSubEl.NElements();i++)
	{
		if (cp.fSubEl[i] == -1)
		{
			this->fSubEl[i] = -1;
			continue;
		}
		if (gl2lcElMap.find(cp.fSubEl[i]) == gl2lcElMap.end())
		{
			std::stringstream sout;
			sout << "ERROR in - " << __PRETTY_FUNCTION__
			<< " subelement " << i << " index = " << cp.fSubEl[i] << " is not in the map.";
			LOGPZ_ERROR (loggerrefpattern,sout.str().c_str());
			exit(-1);
		}
		this->fSubEl[i] = gl2lcElMap[cp.fSubEl[i]];
	}
}


template <class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TGeo>::ClonePatchEl(TPZGeoMesh &DestMesh,
												  std::map<int,int> &gl2lcNdMap,
												  std::map<int,int> &gl2lcElMap) const{
	return new TPZGeoElRefPattern<TGeo>(DestMesh, *this, gl2lcNdMap, gl2lcElMap);
}

#endif
