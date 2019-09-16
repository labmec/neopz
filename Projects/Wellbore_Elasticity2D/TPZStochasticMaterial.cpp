//
//  TPZMatElasticity2D.cpp
//  PZ
//
//  Created by Omar on 10/27/14.
//
//


#include "TPZStochasticMaterial.h"


TPZStochasticMaterial::TPZStochasticMaterial(){

    fnSquareElements = 0.;
    fdirection = 0.;
    finclination = 0.;
    frw = 0.;
    frext = 0.;
    fscale = 0.;
}

TPZStochasticMaterial::TPZStochasticMaterial(TPZGeoMesh* geometricMesh, int numSquareElems, int stochasticInclined, REAL direction,
                       REAL inclination, REAL rw, REAL rext, const TPZFMatrix<REAL> &M, REAL scale, int funcE,
                                                 int funcnu, int distribE, int distribnu){
    //This should be given by the user
    fgmesh = geometricMesh;
    
    //Any way to get this information from the geometric mesh?
    fnSquareElements = numSquareElems;
    fstochasticInclined = stochasticInclined;
    fdirection = direction;
    finclination = inclination;
    frw = rw;
    frext = rext;
    
    // this should not be here, delete after SVD implementation
    fM = M;
    
    // this should be given by the user
    fE_dist = distribE;
    fnu_dist = distribnu;
    
    // this should be given by the user
    fE_funct = funcE;
    fnu_funct = funcnu;
    
    // exponential function scale
    fscale = scale;

}

TPZStochasticMaterial::TPZStochasticMaterial (const TPZStochasticMaterial &cp){

    fgmesh = cp.fgmesh;
    fnSquareElements = cp.fnSquareElements;  // number of Square Elements
    fstochasticInclined = cp.fstochasticInclined;
    fdirection = cp.fdirection;
    finclination = cp.finclination;

    // this should not be here, try to get from fgmesh
    frw = cp.frw;
    frext = cp.frext;
    fM = cp.fM;

    return *this;
}

TPZStochasticMaterial::~TPZStochasticMaterial()
{
}

void TPZStochasticMaterial::SetYoungField(int distribution, int function){
    fE_dist = distribution;
    fE_funct = function;
}

void TPZStochasticMaterial::SetPoissonField(int distribution, int function){
    fnu_dist = distribution;
    fnu_funct = function;
}

void TPZStochasticMaterial::SetReadMatrix(const TPZFMatrix<REAL> &M){
    fM = M;
}


void TPZStochasticMaterial::SetInclinedField(int stochasticInclined,REAL direction, REAL inclination)
{
    fstochasticInclined = stochasticInclined;
    fdirection = direction;
    finclination = inclination;
}

void TPZStochasticMaterial::SetFieldGeometry()
{
    if (fstochasticInclined==true){
        InclinedFieldGeometry();
    }
}

void TPZStochasticMaterial::InclinedFieldGeometry(){
    int nLayers = 8;
    fH = 2 * frext;
    fh = fH / nLayers;
    fmatsize = fnSquareElements * (fH/fh) + fnSquareElements;
}

TPZFMatrix<REAL>  TPZStochasticMaterial::EvaluateCorrelation(int function)
{
    if (fstochasticInclined == 1) {
        
        fK = calcCorrelationMatrixInclined(function);
        
        //PrintCorrelation();
    }
    else{
        
        fK = calcCorrelationMatrix(function);
        
        //PrintCorrelation();
    }
    return fK;
}


TPZFMatrix<REAL> TPZStochasticMaterial::GetCorrelatedVector(int distribution)
{
    if (fstochasticInclined == 1) {
        fU = GetDistribution(fmatsize, distribution);
    }
    else {
        fU = GetDistribution(fnSquareElements, distribution);
    }
    return fU;
}

TPZFMatrix<REAL> TPZStochasticMaterial::GetDistribution(int matrixSize, int distribution)
{
    //normal distribution
    if(distribution==1) {
        
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution(0.,1.0);
        
        // Random Vector U
        TPZFMatrix<REAL> Rand_U (matrixSize, 1, 0.);
        
        for (int i = 0; i < matrixSize; i++) {
            Rand_U(i,0) = distribution(generator);
            distribution.reset();
            fRand_U = Rand_U;
        }
    }
    
    //lognormal distribution
    else if(distribution==2){
        
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::lognormal_distribution<double> distribution(0.,1.0);
        
        // Random Vector U
        TPZFMatrix<REAL> Rand_U (matrixSize, 1, 0.);
        
        for (int i = 0; i < matrixSize; i++) {
            Rand_U(i,0) = distribution(generator);
            distribution.reset();
            fRand_U = Rand_U;
        }
    }
    return fRand_U;
}

void TPZStochasticMaterial::GetStochasticField( TPZFMatrix<REAL> f_E, TPZFMatrix<REAL> f_nu)
{
    // Use global parameter for U and sqrt(S) for the field
    // fE = U*sqrt(S)*fU_E
    // fnu = U*sqrt(S)*fU_nu
}

TPZFMatrix<REAL> TPZStochasticMaterial::calcCorrelationMatrix(int function) {
    
    std::cout << "\nCria matriz da norma entre os centroides (para a matriz de correlacao)" << std::endl;
    
    TPZFMatrix<REAL> CenterNorm(fnSquareElements, fnSquareElements, 0.0);
    
    TPZManVector<REAL, 3> CenterPoint1, CenterPoint2;
    
    // Elemento analizado
    TPZGeoEl *gel1;
    TPZVec<TPZGeoEl *> sub1;
    TPZManVector<REAL> centerpsi1(3), center1(3);
    
    // Outros elementos
    TPZGeoEl *gel2;
    TPZVec<TPZGeoEl *> sub2;
    TPZManVector<REAL> centerpsi2(3), center2(3);
    
    // Matriz de correlacao
    TPZFMatrix<REAL> KCorr(fnSquareElements, fnSquareElements, 0.0);
    
    // Matriz da distancia entre os centroides
    for (int i = 0; i < fnSquareElements; i++) {
        for (int j = 0; j < fnSquareElements; j++) {
            gel1 = fgmesh->ElementVec()[i];
            gel1->CenterPoint(8, centerpsi1);
            gel1->X(centerpsi1, center1);
            
            CenterPoint1 = center1;
            
            gel2 = fgmesh->ElementVec()[j];
            gel2->CenterPoint(8, centerpsi2);
            gel2->X(centerpsi2, center2);
            
            CenterPoint2 = center2;
            
            //    /*3*/    EQuadrilateral
            if (gel1->Type() == 3 && gel2->Type() == 3) {
                
                REAL dx = pow((CenterPoint2[0]-CenterPoint1[0]), 2);
                REAL dy = pow((CenterPoint2[1]-CenterPoint1[1]), 2);
                REAL dz = pow((CenterPoint2[2]-CenterPoint1[2]), 2);
                
                CenterNorm(i,j) = sqrt(dx + dy + dz);
                
                REAL r = CenterNorm(i,j);
                REAL r2 = pow(r, 2);
                
                //exponential function
                if (function==1){
                    KCorr(i,j) = pow(M_E, -((r2*r2)/(fscale*fscale)));
                }
                
                //spherical function
                else if (function==2){
                    //insert function
                }
            }
            
            else {
                // Verifica se el atual eh quadrilatero
                std::cout<< "Element Type Error" << std::endl;
            }
        }
    }
    return KCorr;
}


TPZFMatrix<REAL> TPZStochasticMaterial::calcCorrelationMatrixInclined(int function) {
    
    std::cout << "\nCria matriz dos centroides dos elementos " << std::endl;
    
    TPZFMatrix<REAL> CenterNorm(fmatsize, fmatsize, 0.0);
    
    // Matriz de correlacao
    TPZFMatrix<REAL> KCorr(fmatsize, fmatsize, 0.0);
    
    REAL Pi = M_PI;
    
    //******* angulos COLOCADOS A MAO para fazer teste *********
    REAL alpha = 0.; // azimuth
    REAL beta = 0.; // inclination
    alpha = (fdirection*(Pi/180)); // azimuth
    beta = (finclination*(Pi/180)); // inclination
    
    //  Geeting all coordinates
    TPZGeoEl *gel;
    TPZFMatrix<REAL> Coordinates(fmatsize, 4, 0.0); //nanana
    TPZFMatrix<REAL> rotCoordinates(fnSquareElements, 4, 0.0);
    TPZManVector<REAL> centerpsi(3), center(3);
    TPZManVector<REAL, 3> CenterPoint;
    
    for (int i = 0; i < fnSquareElements; i++) {
        gel = fgmesh->ElementVec()[i];
        gel->CenterPoint(8, centerpsi);
        gel->X(centerpsi, center);
        
        CenterPoint = center;
        
        //    /*3*/    EQuadrilateral
        if (gel->Type() == 3) {
            //Coordinates
            REAL xx = CenterPoint[0];
            REAL yy = CenterPoint[1];
            REAL zz = CenterPoint[2];
            
            Coordinates(i, 0) = i;
            Coordinates(i, 1) = xx;
            Coordinates(i, 2) = yy;
            Coordinates(i, 3) = zz;
        }
    }
    
    int z = 0; // z <= (rext/h);
    int signal = 1;
    REAL altura = (z+(z-1))*(fh/2);
    for (int k = fnSquareElements; k < fmatsize; k += fnSquareElements) {
        if (k >= fmatsize/2 && signal > 0) {
            z = 0;
            signal = -1;
        }
        
        if (k % fnSquareElements == 0) {
            z++;
            altura = signal * (z + (z-1)) * (fh/2);
        }
        
        //std::cout << k << std::endl;
        for (int j = 0; j < fnSquareElements; j++) {
            Coordinates(k+j, 0) = k+j;
            Coordinates(k+j, 1) = Coordinates(j, 1);
            Coordinates(k+j, 2) = Coordinates(j, 2);
            Coordinates(k+j, 3) = altura;
        }
    }
    
    
    // Rotate fnSquareElements Coordinates and alocate in rotCoordinates
    for (int i = 0; i < fnSquareElements; i++) {
        rotCoordinates(i, 0) = i;
        rotCoordinates(i, 1) = Coordinates(i,1)*cos(alpha)*cos(beta) + Coordinates(i,2)*
        cos(beta)*sin(alpha) - Coordinates(i,3)*sin(beta);
        rotCoordinates(i, 2) = Coordinates(i,2)*cos(alpha) - Coordinates(i,1)*sin(alpha);
        rotCoordinates(i, 3) = Coordinates(i,3)*cos(beta) + Coordinates(i,1)*cos(alpha)*
        sin(beta) + Coordinates(i,2)*sin(alpha)*sin(beta);
    }
    
    // Getting rotCoordinates in Coordinates
    for (int i = 0; i < fnSquareElements; i++) {
        Coordinates(i, 0) = i;
        Coordinates(i, 1) = rotCoordinates(i, 1);
        Coordinates(i, 2) = rotCoordinates(i, 2);
        Coordinates(i, 3) = rotCoordinates(i, 3);
    }
    
    
    //std::cout << Coordinates << std::endl;
    std::ofstream out_Coordinates("Coordinates.txt");
    Coordinates.Print("XYZ = ",out_Coordinates,EMathematicaInput);
    
    //std::cout << rotCoordinates << std::endl;
    std::ofstream out_rotCoordinates("rotCoordinates.txt");
    rotCoordinates.Print("XYZ = ",out_rotCoordinates,EMathematicaInput);
    
    std::cout << "\nCria matriz da norma entre os centroides e Matriz de Correlacao" << std::endl;
    
    // Matriz da distancia entre os centroides
    for (int i = 0; i < fmatsize; i++) {
        for (int j = 0; j < fmatsize; j++) {
            
            REAL dx = pow((Coordinates(i,1)-Coordinates(j,1)), 2);
            REAL dy = pow((Coordinates(i,2)-Coordinates(j,2)), 2);
            REAL dz = pow((Coordinates(i,3)-Coordinates(j,3)), 2);
            
            CenterNorm(i,j) = sqrt(dx + dy + dz);
            
            REAL r = CenterNorm(i,j);
            REAL r2 = pow(r, 2);
            
            //exponential function
            if (function==1){
                KCorr(i,j) = pow(M_E, -((r2*r2)/(fscale*fscale)));
            }
            
            //spherical function
            else if (function==2){
                //insert function
            }
        }
    }
    
    //        std::cout << "Numero colunas: " << KCorr.Cols() << std::endl;
    //        std::cout << "Numero linhas: " << KCorr.Rows() << std::endl;
    //        std::cout << "Penultimo valor " << KCorr(fmatsize-1,fmatsize-2) << std::endl;
    //        std::cout << "Ultimo valor " << KCorr(fmatsize-1,fmatsize-1) << std::endl;
    
    return KCorr;
}


void TPZStochasticMaterial::PrintCorrelation() {
    std::ofstream out_kmatrix("KCorr.txt");
    fK.Print("KCorr = ",out_kmatrix,EMathematicaInput);
}


TPZVec<REAL> TPZStochasticMaterial::calcStochasticField(){
    
//    // Stochastic Field
//    TPZFMatrix<REAL> K = calcCorrelationMatrix();
//
    return NULL;
}

