//
//  mymeshes.h
//  PZ
//
//  Created by Agnaldo Farias on 29/04/13.
//
//

#ifndef __PZ__mymeshes__
#define __PZ__mymeshes__

#include <iostream>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzinterpolationspace.h"

#include "pzporoelasticmf2d.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"

#include "pzfunction.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>

using namespace std;


class DadosMalhas{

public:

    DadosMalhas();
    
    ~DadosMalhas();
    
    DadosMalhas(const DadosMalhas &copy);
    
    DadosMalhas &operator=(const DadosMalhas &copy);
    
    TPZGeoMesh *GMesh(bool triang_elements, REAL L, REAL w);
    
 
    TPZGeoMesh *GMesh2(REAL L, REAL w, REAL La);
    TPZGeoMesh *GMesh3(REAL L, REAL w);
    TPZGeoMesh *GMesh4(REAL L, REAL w,int h, int nrefdir);
    
    void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
    
    TPZCompMesh *MalhaCompElast(TPZGeoMesh * gmesh,int pOrder, bool twomaterial, bool stripload);
    
    TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, bool twomaterial, bool stripload);
    
    TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder,bool triang, bool twomaterial);
    
    TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZPoroElasticMF2d * &mymaterial,TPZAutoPointer<TPZFunction<STATE> > solExata);
    
    TPZCompMesh *MalhaCompTerzaghi(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZPoroElasticMF2d * &mymaterial,TPZAutoPointer<TPZFunction<STATE> > solExata);
    
    TPZCompMesh *MalhaCompBarryMercer(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZAutoPointer<TPZFunction<STATE> > sourceterm, TPZAutoPointer<TPZFunction<STATE> > solExata);
    
    TPZCompMesh *MalhaCompBarryMercerPressureSolution(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZPoroElasticMF2d * &mymaterial, TPZAutoPointer<TPZFunction<STATE> > BCterm, TPZAutoPointer<TPZFunction<STATE> > solExata);

    
    void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);


    TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZPoroElasticMF2d *mymateria, TPZCompMesh *mphysics, int nthreads);
    
    TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZCompMesh *mphysics, int nthreads);

    void StiffMatrixLoadVec(TPZPoroElasticMF2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<STATE> &matK1, TPZFMatrix<STATE> &fvec, int nthreads);
    
    void StiffMatrixLoadVec(TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<STATE> &matK1, TPZFMatrix<STATE> &fvec, int nthreads);

    void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);

    void SolveSistTransient(REAL deltaT,REAL maxTime, TPZPoroElasticMF2d * &mymaterial,
                        TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, int ntimestep, REAL &timeatual);
    
    void SolveSistBarryMercert(REAL deltaT,REAL maxTime, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, int ntimestep, REAL &timeatual);
    
    void SolveSistWithError(REAL deltaT,REAL maxTime,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, int ntimestep, REAL &timeatual, TPZAutoPointer<TPZFunction<STATE> > solExata1, TPZAutoPointer<TPZFunction<STATE> > solExata2, int h,  ofstream saidaerro);


    TPZCompMesh *CMeshPressureL2(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini, bool triang);
    
    void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);
    
    void RefiningNearLine(int dim,TPZGeoMesh *gmesh,int nref);
    
    void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance);
    
    void AjustarContorno(TPZGeoMesh *gmesh);
    
    
    void SetParameters(REAL mod_young, REAL mod_poisson, REAL coef_alpha, REAL coef_Se, REAL permeabil_fluido, REAL visc_fluido, REAL fx, REAL fy,REAL sign){
        
        fEyoung = mod_young;
        fpoisson= mod_poisson;
        falpha = coef_alpha;
        fSe = coef_Se;
        fperm = permeabil_fluido;
        fvisc = visc_fluido;
        ffx = fx;
        ffy = fy;
        fsign = sign;
    }
    
    void SetPrefLref(REAL pref, REAL Lref,REAL kovervisc){
        fpref = pref;
        fLref = Lref;
        fkovervisc = kovervisc;
    }
    
    int GetIdSourceTerm(){
        return fbcSourceTerm;
    }
    
    void SetValSourceTerm(REAL sf){
        
        fvalsourceterm = sf;
    }
    
    int GetMatId(){
        return fmatId;
    }
    
protected:
    
    //dados do material
    REAL fEyoung;
    REAL fpoisson;
    REAL falpha;
    REAL fSe;
    REAL fperm;
    REAL fvisc;
    REAL ffx;
    REAL ffy;
    REAL fsign;
    
    REAL fpref;
    REAL fLref;
    REAL fkovervisc;
    
    REAL fvalsourceterm;
    
    //dados da malha geometrica
    int fmatId;
    int fbcBottom;
    int fbcRight;
    int fbcTop;
    int fbcTopStripLoad;
    int fbcLeft;
    int fbcSourceTerm;
    
    int  fdirichlet;
    int  fneumann;
    int  fneumdir;
    int  fdirfreey_neum;
    int  fdirneum;
    int  fmixedneum;
    int  fmixeddirich;
    
    int fmixedFreeXYdirich;
    int fmixedFreeYXdirich;
};

#endif /* defined(__PZ__mymeshes__) */
