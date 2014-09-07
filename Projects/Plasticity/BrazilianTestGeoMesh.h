#ifndef QUARTERWELLBOREMESH_H
#define QUARTERWELLBOREMESH_H


#include <iostream>
#include <cstdlib>
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZTensor.h"
#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"
#include "TPZYCVonMises.h"
#include "TPZVonMises.h"
#include "pzfstrmatrix.h"
#include "pzbndmat.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "tpzgeoelrefpattern.h"
#include "pzbndcond.h"
#include "pzstepsolver.h"
#include "TPZTensor.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZMohrCoulomb.h"
#include "TPZDruckerPrager.h"
#include "GeoMeshClass.h"
#include "pzelastoplastic2D.h"
#include "tpzycvonmisescombtresca.h"
#include "TPZMohrCoulombNeto.h"
#include "TPZSandlerDimaggio.h"

#include "pzstring.h"

inline REAL propPG(REAL ratio, int nrad, int ncirc, int i)
{
	REAL q = 1.5;
	
	//REAL a0 = (ratio - 1.) * (q - 1.)/(pow(q, nrad+1) -1.);
	
	REAL SPi = (pow(q,i+1) - 1.) / (q - 1.) -1.;
	REAL SPn = (pow(q,nrad+1) - 1.) / (q - 1.) -1.;
	
	return SPi / SPn;
}

const REAL pi = 3.1415927;

template <class T>
inline void AddVec(T & target, const T & v1, const T & v2, const REAL multipl1 = 1., const REAL multipl2 = 1.)
{
    long size  = v1.NElements();
	long size2 = v2.NElements();
	if(size2 < size) size = size2;
	
	target.Resize(size);
	
	for( long i = 0; i < size; i++)
	{
		target[i] = v1[i] * multipl1 +
        v2[i] * multipl2;
	}
}

template <class T>
inline void SubtrVec(T & target, const T & v1, const T & v2)
{
	AddVec(target, v1, v2, 1., -1.);
}

inline TPZGeoMesh * CreateGMesh(TPZVec< TPZVec< REAL > > & pt,
					     TPZVec< TPZVec< long > > & el,
					     MElementType ElType,
					     int matId)
{

	long npt = pt.NElements();
	long nel = el.NElements();
	long i;
	
	TPZGeoMesh * pGMesh = new TPZGeoMesh;
	
	pGMesh->NodeVec().Resize(npt);
	
	for(i = 0; i < npt; i++)
        pGMesh->NodeVec()[i].Initialize(pt[i], *pGMesh);
	
	for(i = 0; i < nel; i++)
        pGMesh->CreateGeoElement(ElType, el[i], matId, i);
	
	pGMesh->BuildConnectivity();
	
	return pGMesh;
}



inline REAL propPA(REAL ratio, int nrad, int ncirc, int i)
{
	REAL rncirc = ncirc;
	REAL rnrad  = nrad;
	REAL ri = i;
	
	REAL d = ( 2. * (ratio - 1.) / (rnrad) - pi / rncirc ) / (rnrad - 1.);
	
	REAL SAi = (pi/rncirc + (ri - 1.) * d) * ri;
	REAL SAn = (pi/rncirc + (rnrad - 1.) * d) * rnrad;
	
	return SAi / SAn;
}

#define PLASTICITY_CLEAN_OUT

inline void PrepareInitialMat(TPZPlasticBase & mat, TPZTensor<REAL> &initialStress, TPZTensor<REAL> &endStress, int steps)
{

	REAL multipl;
	int i;
	TPZTensor<REAL> strain, localLoad, diffStress;
	
	diffStress = endStress;
	diffStress.Add(initialStress, -1.);
	
	for(i = 1; i <= steps; i++)
	{
#ifndef PLASTICITY_CLEAN_OUT
		cout << "Starting step " << i << " of " << steps << endl;
#endif
		if(i == 0)
		{
			multipl = 0;
		}else
		{
			multipl = (REAL)i / (REAL)steps;
		}
		
		localLoad = initialStress;
		localLoad.Add(diffStress, multipl);
        
		
		mat.ApplyLoad(localLoad, strain);
//        cout << "\n load  = "<< localLoad <<endl;
//        cout << "\n strain  = "<< strain <<endl;
        if(i==1||i==10)
        {
            //cout << "\n State = "<< mat.GetState() <<endl;
        }
		
	}
	
}



static void QuarterWellboreGeom(int ncirc,
						 REAL ioratio,
						 TPZVec< TPZVec<REAL> > &pt,
						 TPZVec< TPZVec<long> > &el,
						 TPZVec< MElementType > &elType,
						 TPZVec< TPZString > &elName,
						 int & nrad)
{

	
	REAL pi = M_PI;
	int i, j;
	ncirc = 2 * (ncirc / 2); //ensuring it's even
	int ncircpt = ncirc + 1;
	REAL q = 1.4; // 1.5
	nrad = static_cast<int>(( log(2.*((REAL) ncirc)*(ioratio - 1.)/pi)/log(q) -1. )*1.4 ); //1.6
	
	int nradpt = nrad + 1;
	long nel3d = ncirc * nrad;
	
	long nel = 3*nel3d + 2*ncirc + 2*nrad;
	
	///nel = nel3d;////
	
	int nlayerpt = ncircpt * nradpt;
	
	pt.Resize(2 * nlayerpt);
	el.Resize(nel);
	elType.Resize(nel);
	elName.Resize(nel);
	
	//creating the vertices
	// wellbore vertices
	for( i = 0 ; i < ncircpt; i++)
	{
		REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
		TPZVec< REAL > pti(3,0.);
		pti[0] = cos(theta);
		pti[1] = sin(theta);
		pt[i] = pti;
	}
	// external x = ratio points
	for( i = 0; i <= ncirc / 2; i ++)
	{
		REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
		//REAL scale = ioratio / cos(theta);
		TPZVec< REAL > pti(3,0.);
		//pti[0] = ioratio;
		//pti[1] = sin(theta) * scale;
		pti[0] = ioratio * cos(theta);
		pti[1] = ioratio * sin(theta);
		pt[i + ncircpt * nrad] = pti;
	}
	// external y = ratio points
	for( i = ncirc / 2; i < ncircpt ; i ++)
	{
		REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
		//REAL scale = ioratio / sin(theta);
		TPZVec< REAL > pti(3,0.);
		//pti[0] = cos(theta) * scale;
		//pti[1] = ioratio;
		pti[0] = ioratio * cos(theta);
		pti[1] = ioratio * sin(theta);
		pt[i + ncircpt * nrad] = pti;
	}
	// filling the inbetween points
	for( i = 1; i < nrad; i++)
	{
		REAL prop = propPG(ioratio, nrad, ncirc, i);
		for( j = 0; j < ncircpt; j++)
		{
			TPZVec< REAL > pti(3,0.);
			TPZVec< REAL > ptExt(pt[ncircpt * nrad + j]);
			REAL sizeExt = ptExt[0] * ptExt[0] + ptExt[1] * ptExt[1] + ptExt[2] * ptExt[2];
			sizeExt = sqrt(sizeExt);
			pti[0] = pt[j][0] * (1. - prop) + ptExt[0] * prop * sqrt(ioratio / sizeExt);
			pti[1] = pt[j][1] * (1. - prop) + ptExt[1] * prop * sqrt(ioratio / sizeExt);
			pti[2] = pt[j][2] * (1. - prop) + ptExt[2] * prop * sqrt(ioratio / sizeExt);
			pt[i*ncircpt + j] = pti;
		}
	}
	// creating the z-layer points
	for( i = 0 ; i < nlayerpt; i++)
	{
		TPZVec< REAL > pti(3,0.);
		pti[0] = pt[i][0];
		pti[1] = pt[i][1];
		pti[2] = 1.;
		pt[i + nlayerpt] = pti;
	}
	//creating the 3D elements
	for(i = 0; i < nrad; i++)
	{
		for(j = 0; j < ncirc; j++)
		{
			long index = i * ncirc + j;
			TPZVec< long > eli(8,0);
			eli[0] = ncircpt * i + j;
			eli[1] = ncircpt * (i + 1) + j;
			eli[2] = eli[1] + 1;
			eli[3] = eli[0] + 1;
			eli[4] = eli[0] + nlayerpt;
			eli[5] = eli[1] + nlayerpt;
			eli[6] = eli[2] + nlayerpt;
			eli[7] = eli[3] + nlayerpt;
			el[index] = eli;
			elType[index] = ECube;
			elName[index] = "Interior";
		}
	}
	long lastIndex = nel3d;
	
	// Wellbore Faces
	for(i = 0; i < ncirc; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		
		eli[0] = el[i][0];
		eli[1] = el[i][3];
		eli[2] = el[i][7];
		eli[3] = el[i][4];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Wellbore";
	}
	lastIndex += ncirc;
	
	//creating the 2D border elements
	for(i = 0; i < nel3d; i++)
	{
		TPZVec< long > elTop(4), elBot(4);
		long indexBot = i + lastIndex;
		long indexTop = i + lastIndex + nel3d;
		
		for(j = 0 ; j < 4; j++)
		{
			elBot[j] = el[i][j];
			elTop[j] = el[i][j+4];
		}
		el[indexBot] = elBot;
		elType[indexBot] = EQuadrilateral;
		elName[indexBot] = "Z- Plane Strain";
        
		el[indexTop] = elTop;
		elType[indexTop] = EQuadrilateral;
		elName[indexTop] = "Z+ Plane Strain";
	}
	lastIndex += 2* nel3d;
	
	// Lower Symmetry
	for(i = 0; i < nrad; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		long refIndex = ncirc*i;
		
		eli[0] = el[refIndex][0];
		eli[1] = el[refIndex][1];
		eli[2] = el[refIndex][5];
		eli[3] = el[refIndex][4];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Lower Symmetry";
	}
	lastIndex += nrad;
    
	// Left Symmetry
	for(i = 0; i < nrad; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		long refIndex = ncirc*(i+1) - 1;
		
		eli[0] = el[refIndex][3];
		eli[1] = el[refIndex][2];
		eli[2] = el[refIndex][6];
		eli[3] = el[refIndex][7];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Left Symmetry";
	}
	lastIndex += nrad;
	
	// Right Farfield
	for(i = 0; i < ncirc / 2; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		long refIndex = (nrad-1)*ncirc + i;
		
		eli[0] = el[refIndex][1];
		eli[1] = el[refIndex][2];
		eli[2] = el[refIndex][6];
		eli[3] = el[refIndex][5];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Right Farfield";
	}
	//lastIndex += ncirc / 2;
	
	// Top Farfield
	for(i = ncirc / 2; i < ncirc; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		long refIndex = (nrad-1)*ncirc + i;
		
		eli[0] = el[refIndex][1];
		eli[1] = el[refIndex][2];
		eli[2] = el[refIndex][6];
		eli[3] = el[refIndex][5];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Top Farfield";
	}
	

    
}

#include "pzpostprocanalysis.h"
#include "pzelastoplasticanalysis.h"

static TPZCompMesh * CreateQuarterWellboreMesh( int gOrder,
										int ncirc,
										REAL ioratio,
						  				TPZMaterial * pMat,
										TPZFMatrix<REAL> & BCStressState,
										TPZFMatrix<REAL> & WellboreStressState,
										int allNeumannBC = 0)
{

	
	int matId = pMat->Id();
	int nstate = pMat->NStateVariables();
	long i;
	int nrad;//int ncirc = 4, nrad;
	//REAL ioratio = 10.;
	
	TPZVec< TPZVec< REAL > > pt;
	TPZVec< TPZVec< long > > el;
	TPZVec< MElementType > elType;
	TPZVec< TPZString > elName;
	TPZFMatrix<REAL> val1(nstate,nstate), val2(nstate,1);
	TPZFMatrix<REAL> bcNormal(nstate, 1, 0.);
	TPZFMatrix<REAL> bcStressState(nstate, nstate, 0.),
    wellboreStressState(nstate, nstate, 0.);
    
	QuarterWellboreGeom(ncirc, ioratio, pt, el, elType, elName, nrad);
	long nel = el.NElements();
	TPZGeoMesh * pGMesh = new TPZGeoMesh;
    
    //TPZPoroElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pGMesh); // self explanatory
    
	
	// preparing nodes and elements
	long npt = pt.NElements();
	pGMesh->NodeVec().Resize(npt);
	for(i = 0; i < npt; i++)
        pGMesh->NodeVec()[i].Initialize(pt[i], *pGMesh);
	
	for(i = 0; i < nel; i++)
	{
		matId = 0;
		if(!strcmp(elName[i].Str(),"Interior"))       matId =  1;
		if(!strcmp(elName[i].Str(),"Z- Plane Strain"))matId = -1;
		if(!strcmp(elName[i].Str(),"Z+ Plane Strain"))matId = -2;
		if(!strcmp(elName[i].Str(),"Lower Symmetry")) matId = -3;
		if(!strcmp(elName[i].Str(),"Left Symmetry"))  matId = -4;
		if(!strcmp(elName[i].Str(),"Right Farfield")) matId = -5;
		if(!strcmp(elName[i].Str(),"Top Farfield"))   matId = -5;
		if(!strcmp(elName[i].Str(),"Wellbore"))		  matId = -7;
        
		if(matId != 0)
		{
			pGMesh->CreateGeoElement(elType[i], el[i], matId, i);
		}else{
			PZError << "\nQuarterWellboreMesh error - element " << i << " without material assignement.\n";
		}
	}
    
	pGMesh->BuildConnectivity();
    
	TPZCompEl::SetgOrder(gOrder);
	
	// Creating Computation Mesh
	TPZCompMesh * pCMesh = new TPZCompMesh(pGMesh);
    
     TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
   
    
	//Preparing Material
	pMat -> SetForcingFunction(NULL);
	TPZMaterial  *mat(pMat);
	pCMesh->InsertMaterialObject(mat);
	
	//preparing BCs
	
	if(allNeumannBC)
	{
		for(i = -1; i > -7; i--)
		{
            //TPZMaterial *bc1 = mat->CreateBC(mat,-1,0,k1,f1);
			TPZMaterial * bc;
			bc = mat->CreateBC(mat, i, 4 /*StressField Neumann BCType*/, BCStressState, val2);
			pCMesh->InsertMaterialObject(bc);
		}
	}else
	{
		{//"Z- Plane Strain" bc -1
			bcNormal.Zero();
			bcNormal(2,0) = -1.;
			val2.Zero();
			TPZMaterial * bc;
			val2(2,0) = 1.;
			bc = mat->CreateBC(mat, -1, 3 /*Directional Dirichlet BCType*/, val1, val2);
			pCMesh->InsertMaterialObject(bc);
		}
		{//"Z+ Plane Strain" bc -2
			bcNormal.Zero();
			bcNormal(2,0) = 1.;
			val2.Zero();
			TPZMaterial * bc;
			val2(2,0) = 1.;
			bc = mat->CreateBC(mat, -2, 3 /*Directional Dirichlet BCType*/, val1, val2);
			pCMesh->InsertMaterialObject(bc);
		}
		{//"Lower Symmetry"  bc -3
			bcNormal.Zero();
			bcNormal(1,0) = -1.;
			val2.Zero();
			TPZMaterial * bc;
			val2(1,0) = 1.;
			bc = mat->CreateBC(mat, -3, 3 /*Directional Dirichlet BCType*/, val1, val2);
			pCMesh->InsertMaterialObject(bc);
		}
		{//"Left Symmetry"   bc -4
			bcNormal.Zero();
			bcNormal(0,0) = -1.;
			val2.Zero();
			TPZMaterial * bc;
			val2(0,0) = 1.;
			bc = mat->CreateBC(mat, -4, 3 /*Directional Dirichlet BCType*/, val1, val2);
			pCMesh->InsertMaterialObject(bc);
		}
		{//"Right & Top Farfield"  bc -5
			val2.Zero();
			val2(0,0) = 1.;
			val2(1,0) = 1.;
			TPZMaterial * bc;
			bc = mat->CreateBC(mat, -5, 3, val1, val2);	//Directional Dirichlet BCType
			pCMesh->InsertMaterialObject(bc);
		}
		
	}
    
	{//"Wellbore" bc -7
		val2.Zero();
		TPZMaterial * bc;
		bc = mat->CreateBC(mat, -7, 4 , WellboreStressState, val2);//StressField Neumann
		pCMesh->InsertMaterialObject(bc);
	}
	
	// building mesh connections

   	pCMesh->AutoBuild();

	return pCMesh;
}

static void QuarterWellboreGeom2d(int ncirc,
						 REAL ioratio,
						 TPZVec< TPZVec<REAL> > &pt,
						 TPZVec< TPZVec<long> > &el,
						 TPZVec< MElementType > &elType,
						 TPZVec< TPZString > &elName,
						 int & nrad)
{
    
	
	REAL pi = 3.1415927;
	int i, j;
	ncirc = 2 * (ncirc / 2); //ensuring it's even
	int ncircpt = ncirc + 1;
	REAL q = 1.4; // 1.5
	nrad = static_cast<int>(( log(2.*((REAL) ncirc)*(ioratio - 1.)/pi)/log(q) -1. )*1.4 ); //1.6
	
	int nradpt = nrad + 1;
	long nel3d = ncirc * nrad;
	
	long nel = 3*nel3d + 2*ncirc + 2*nrad;
	
	///nel = nel3d;////
	
	int nlayerpt = ncircpt * nradpt;
	
	pt.Resize(2 * nlayerpt);
	el.Resize(nel);
	elType.Resize(nel);
	elName.Resize(nel);
	
	//creating the vertices
	// wellbore vertices
	for( i = 0 ; i < ncircpt; i++)
	{
		REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
		TPZVec< REAL > pti(3,0.);
		pti[0] = cos(theta);
		pti[1] = sin(theta);
		pt[i] = pti;
	}
	// external x = ratio points
	for( i = 0; i <= ncirc / 2; i ++)
	{
		REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
		//REAL scale = ioratio / cos(theta);
		TPZVec< REAL > pti(3,0.);
		//pti[0] = ioratio;
		//pti[1] = sin(theta) * scale;
		pti[0] = ioratio * cos(theta);
		pti[1] = ioratio * sin(theta);
		pt[i + ncircpt * nrad] = pti;
	}
	// external y = ratio points
	for( i = ncirc / 2; i < ncircpt ; i ++)
	{
		REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
		//REAL scale = ioratio / sin(theta);
		TPZVec< REAL > pti(3,0.);
		//pti[0] = cos(theta) * scale;
		//pti[1] = ioratio;
		pti[0] = ioratio * cos(theta);
		pti[1] = ioratio * sin(theta);
		pt[i + ncircpt * nrad] = pti;
	}
	// filling the inbetween points
	for( i = 1; i < nrad; i++)
	{
		REAL prop = propPG(ioratio, nrad, ncirc, i);
		for( j = 0; j < ncircpt; j++)
		{
			TPZVec< REAL > pti(3,0.);
			TPZVec< REAL > ptExt(pt[ncircpt * nrad + j]);
			REAL sizeExt = ptExt[0] * ptExt[0] + ptExt[1] * ptExt[1] + ptExt[2] * ptExt[2];
			sizeExt = sqrt(sizeExt);
			pti[0] = pt[j][0] * (1. - prop) + ptExt[0] * prop * sqrt(ioratio / sizeExt);
			pti[1] = pt[j][1] * (1. - prop) + ptExt[1] * prop * sqrt(ioratio / sizeExt);
			pti[2] = pt[j][2] * (1. - prop) + ptExt[2] * prop * sqrt(ioratio / sizeExt);
			pt[i*ncircpt + j] = pti;
		}
	}
	// creating the z-layer points
	for( i = 0 ; i < nlayerpt; i++)
	{
		TPZVec< REAL > pti(3,0.);
		pti[0] = pt[i][0];
		pti[1] = pt[i][1];
		pti[2] = 1.;
		pt[i + nlayerpt] = pti;
	}
	//creating the 3D elements
	for(i = 0; i < nrad; i++)
	{
		for(j = 0; j < ncirc; j++)
		{
			long index = i * ncirc + j;
			TPZVec< long > eli(8,0);
			eli[0] = ncircpt * i + j;
			eli[1] = ncircpt * (i + 1) + j;
			eli[2] = eli[1] + 1;
			eli[3] = eli[0] + 1;
			eli[4] = eli[0] + nlayerpt;
			eli[5] = eli[1] + nlayerpt;
			eli[6] = eli[2] + nlayerpt;
			eli[7] = eli[3] + nlayerpt;
			el[index] = eli;
			elType[index] = ECube;
			elName[index] = "Interior";
		}
	}
	int lastIndex = nel3d;
	
	// Wellbore Faces
	for(i = 0; i < ncirc; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		
		eli[0] = el[i][0];
		eli[1] = el[i][3];
		eli[2] = el[i][7];
		eli[3] = el[i][4];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Wellbore";
	}
	lastIndex += ncirc;
	
	//creating the 2D border elements
	for(i = 0; i < nel3d; i++)
	{
		TPZVec< long > elTop(4), elBot(4);
		long indexBot = i + lastIndex;
		long indexTop = i + lastIndex + nel3d;
		
		for(j = 0 ; j < 4; j++)
		{
			elBot[j] = el[i][j];
			elTop[j] = el[i][j+4];
		}
		el[indexBot] = elBot;
		elType[indexBot] = EQuadrilateral;
		elName[indexBot] = "Z- Plane Strain";
        
		el[indexTop] = elTop;
		elType[indexTop] = EQuadrilateral;
		elName[indexTop] = "Z+ Plane Strain";
	}
	lastIndex += 2* nel3d;
	
	// Lower Symmetry
	for(i = 0; i < nrad; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		long refIndex = ncirc*i;
		
		eli[0] = el[refIndex][0];
		eli[1] = el[refIndex][1];
		eli[2] = el[refIndex][5];
		eli[3] = el[refIndex][4];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Lower Symmetry";
	}
	lastIndex += nrad;
    
	// Left Symmetry
	for(i = 0; i < nrad; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		long refIndex = ncirc*(i+1) - 1;
		
		eli[0] = el[refIndex][3];
		eli[1] = el[refIndex][2];
		eli[2] = el[refIndex][6];
		eli[3] = el[refIndex][7];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Left Symmetry";
	}
	lastIndex += nrad;
	
	// Right Farfield
	for(i = 0; i < ncirc / 2; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		int refIndex = (nrad-1)*ncirc + i;
		
		eli[0] = el[refIndex][1];
		eli[1] = el[refIndex][2];
		eli[2] = el[refIndex][6];
		eli[3] = el[refIndex][5];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Right Farfield";
	}
	//lastIndex += ncirc / 2;
	
	// Top Farfield
	for(i = ncirc / 2; i < ncirc; i++)
	{
		TPZVec< long > eli(4);
		long index = i + lastIndex;
		long refIndex = (nrad-1)*ncirc + i;
		
		eli[0] = el[refIndex][1];
		eli[1] = el[refIndex][2];
		eli[2] = el[refIndex][6];
		eli[3] = el[refIndex][5];
		
		el[index] = eli;
		elType[index] = EQuadrilateral;
		elName[index] = "Top Farfield";
	}
}


static TPZCompMesh * CreateQuarterWellboreMesh2d( int gOrder,
										int ncirc,
										REAL ioratio,
						  				TPZMaterial * pMat,
										TPZFMatrix<REAL> & BCStressState,
										TPZFMatrix<REAL> & WellboreStressState,
										int allNeumannBC = 0)
{

	int matId = pMat->Id();
	int nstate = pMat->NStateVariables();
	long i;
	int nrad;//int ncirc = 4, nrad;
	//REAL ioratio = 10.;
	
	TPZVec< TPZVec< REAL > > pt;
	TPZVec< TPZVec< long > > el;
	TPZVec< MElementType > elType;
	TPZVec< TPZString > elName;
	TPZFMatrix<REAL> val1(nstate,nstate), val2(nstate,1);
	TPZFMatrix<REAL> bcNormal(nstate, 1, 0.);
	TPZFMatrix<REAL> bcStressState(nstate, nstate, 0.),
    wellboreStressState(nstate, nstate, 0.);
    
	QuarterWellboreGeom(ncirc, ioratio, pt, el, elType, elName, nrad);
	long nel = el.NElements();
	TPZGeoMesh * pGMesh = new TPZGeoMesh;
    
    //TPZPoroElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pGMesh); // self explanatory
    
	
	// preparing nodes and elements
	long npt = pt.NElements();
	pGMesh->NodeVec().Resize(npt);
	for(i = 0; i < npt; i++)
        pGMesh->NodeVec()[i].Initialize(pt[i], *pGMesh);
	
	for(i = 0; i < nel; i++)
	{
		matId = 0;
		if(!strcmp(elName[i].Str(),"Interior"))       matId =  1;
		if(!strcmp(elName[i].Str(),"Z- Plane Strain"))matId = -1;
		if(!strcmp(elName[i].Str(),"Z+ Plane Strain"))matId = -2;
		if(!strcmp(elName[i].Str(),"Lower Symmetry")) matId = -3;
		if(!strcmp(elName[i].Str(),"Left Symmetry"))  matId = -4;
		if(!strcmp(elName[i].Str(),"Right Farfield")) matId = -5;
		if(!strcmp(elName[i].Str(),"Top Farfield"))   matId = -5;
		if(!strcmp(elName[i].Str(),"Wellbore"))		  matId = -7;
        
		if(matId != 0)
		{
			pGMesh->CreateGeoElement(elType[i], el[i], matId, i);
		}else{
			PZError << "\nQuarterWellboreMesh error - element " << i << " without material assignement.\n";
		}
	}
    
	pGMesh->BuildConnectivity();
    
	TPZCompEl::SetgOrder(gOrder);
	
	// Creating Computation Mesh
	TPZCompMesh * pCMesh = new TPZCompMesh(pGMesh);
    
    TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
    
    
	//Preparing Material
	pMat -> SetForcingFunction(NULL);
	TPZMaterial  *mat(pMat);
	pCMesh->InsertMaterialObject(mat);
	
	//preparing BCs
	
	if(allNeumannBC)
	{
		for(i = -1; i > -7; i--)
		{
            //TPZMaterial *bc1 = mat->CreateBC(mat,-1,0,k1,f1);
			TPZMaterial * bc;
			bc = mat->CreateBC(mat, i, 4 /*StressField Neumann BCType*/, BCStressState, val2);
			pCMesh->InsertMaterialObject(bc);
		}
	}else
	{
		{//"Z- Plane Strain" bc -1
			bcNormal.Zero();
			bcNormal(2,0) = -1.;
			val2.Zero();
			TPZMaterial * bc;
			val2(2,0) = 1.;
			bc = mat->CreateBC(mat, -1, 3 /*Directional Dirichlet BCType*/, val1, val2);
			pCMesh->InsertMaterialObject(bc);
		}
		{//"Z+ Plane Strain" bc -2
			bcNormal.Zero();
			bcNormal(2,0) = 1.;
			val2.Zero();
			TPZMaterial * bc;
			val2(2,0) = 1.;
			bc = mat->CreateBC(mat, -2, 3 /*Directional Dirichlet BCType*/, val1, val2);
			pCMesh->InsertMaterialObject(bc);
		}
		{//"Lower Symmetry"  bc -3
			bcNormal.Zero();
			bcNormal(1,0) = -1.;
			val2.Zero();
			TPZMaterial * bc;
			val2(1,0) = 1.;
			bc = mat->CreateBC(mat, -3, 3 /*Directional Dirichlet BCType*/, val1, val2);
			pCMesh->InsertMaterialObject(bc);
		}
		{//"Left Symmetry"   bc -4
			bcNormal.Zero();
			bcNormal(0,0) = -1.;
			val2.Zero();
			TPZMaterial * bc;
			val2(0,0) = 1.;
			bc = mat->CreateBC(mat, -4, 3 /*Directional Dirichlet BCType*/, val1, val2);
			pCMesh->InsertMaterialObject(bc);
		}
		{//"Right & Top Farfield"  bc -5
			val2.Zero();
			val2(0,0) = 1.;
			val2(1,0) = 1.;
			TPZMaterial * bc;
			bc = mat->CreateBC(mat, -5, 3, val1, val2);	//Directional Dirichlet BCType
			pCMesh->InsertMaterialObject(bc);
		}
		
	}
    
	{//"Wellbore" bc -7
		val2.Zero();
		TPZMaterial * bc;
		bc = mat->CreateBC(mat, -7, 4 , WellboreStressState, val2);//StressField Neumann
		pCMesh->InsertMaterialObject(bc);
	}
	
	// building mesh connections
    
   	pCMesh->AutoBuild();
    
	return pCMesh;
}



#endif
