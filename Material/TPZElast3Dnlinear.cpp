//
//  TPZElast3Dnlinear.cpp
//  PZ
//
//  Created by Cesar Lucci on 22/10/13.
//
//

#include "TPZElast3Dnlinear.h"

TPZElast3Dnlinear::TPZElast3Dnlinear() : TPZElasticity3D()
{
    
}

TPZElast3Dnlinear::TPZElast3Dnlinear(int nummat, STATE E, STATE poisson, TPZVec<STATE> &force,
                                     STATE preStressXX, STATE preStressYY, STATE preStressZZ) :
TPZElasticity3D(nummat, E, poisson, force, preStressXX, preStressYY, preStressZZ)
{
    
}

TPZElast3Dnlinear::~TPZElast3Dnlinear()
{
    
}

void TPZElast3Dnlinear::Contribute(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<STATE> ek(ef.Rows(),ef.Rows());
    Contribute(data, weight, ek, ef);
}

void TPZElast3Dnlinear::Contribute(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef)
{
#ifdef DEBUG
    if(ek.Rows() != ef.Rows())
    {
        std::cout << "\n\n" << "ek and ef should have same number of rows!!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
    
    if(ef.Cols() != this->fNumLoadCases)
    {
        std::cout << "\n\n" << "ef should have fNumLoadCases equals to its NCols()!!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
#endif
    
    TPZFMatrix<STATE> ekTemp(ek.Rows(),ek.Cols());
    ekTemp.Zero();
    
    TPZFMatrix<STATE> efTemp(ef.Rows(),ef.Cols());
    efTemp.Zero();
    
    TPZElasticity3D::Contribute(data, weight, ekTemp, efTemp);
    
    for(int i = 0; i < ek.Rows(); i++)
    {
        for(int loadCase = 0; loadCase < fNumLoadCases; loadCase++)
        {
            ef(i,loadCase) += efTemp(i,loadCase);//Residuo (etapa +f)
        }
        for(int j = 0; j < ek.Cols(); j++)
        {
            ek(i,j) += ekTemp(i,j);//Matriz tangente
            
            for(int loadCase = 0; loadCase < fNumLoadCases; loadCase++)
            {
                ef(i,loadCase) -= ekTemp(i,j)*efTemp(j,loadCase);//Residuo (etapa -K.u)
            }
        }
    }
}

void TPZElast3Dnlinear::ContributeBC(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef,
                                     TPZBndCond &bc)
{
#ifdef DEBUG
    if(ek.Rows() != ef.Rows())
    {
        std::cout << "\n\n" << "ek and ef should have same number of rows!!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
    
    if(ef.Cols() != this->fNumLoadCases)
    {
        std::cout << "\n\n" << "ef should have fNumLoadCases equals to its NCols()!!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
#endif
    
    TPZFMatrix<STATE> ekTemp(ek.Rows(),ek.Cols());
    ekTemp.Zero();
    
    TPZFMatrix<STATE> efTemp(ef.Rows(),ef.Cols());
    efTemp.Zero();
    
    TPZElasticity3D::ContributeBC(data, weight, ekTemp, efTemp, bc);
    
    for(int i = 0; i < ek.Rows(); i++)
    {
        for(int loadCase = 0; loadCase < fNumLoadCases; loadCase++)
        {
            ef(i,loadCase) += efTemp(i,loadCase);//Residuo (etapa +f)
        }
        for(int j = 0; j < ek.Cols(); j++)
        {
            ek(i,j) += ekTemp(i,j);//Matriz tangente
            
            for(int loadCase = 0; loadCase < fNumLoadCases; loadCase++)
            {
                ef(i,loadCase) -= ekTemp(i,j)*efTemp(j,loadCase);//Residuo (etapa -K.u)
            }
        }
    }
}

int TPZElast3Dnlinear::ClassId() const
{
	return TPZELASTICITY3DNLINEARMATERIALID;
}


