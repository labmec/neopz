//
//  TPZMatElasticity2D.cpp
//  PZ
//
//  Created by Omar on 10/27/14.
//
//

#include <math.h>

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

    //** Calculation of stochastic field - This should be defined in "main or computational mesh *****/////
    //This can be defined in "main" or computational mesh
    if (fstochasticInclined==1) {
        SetInclinedField(frw, frext, fstochasticInclined, fdirection, finclination);
    }
    calcStochasticField();
    //******* End
    
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
}

TPZStochasticMaterial::~TPZStochasticMaterial()
{
}

TPZFMatrix<STATE> TPZStochasticMaterial::calcStochasticField(){
    
    //This should be removed after using SVD decomposition
    SetReadMatrix(fM);
    
    /* Start */
    
    //Get fK
    // if fE_func==nu_func, fK will be the same for both variables
    if(fE_funct==fnu_funct){
        fKE = EvaluateCorrelation(fE_funct);
        fKnu = fKE;
        
        // Create function to decompose fK using SVD decomposition
        // fK = U sqt(S) V'
        // create glogal parameter for U and sqrt(S)
    }
    
    else {
        //Evaluate for E
        fKE = EvaluateCorrelation(fE_funct);
        //Evaluate for nu
        fKnu = EvaluateCorrelation(fnu_funct);
        
        // Create function to decompose fK using SVD decomposition
        // fKE = U sqt(S) V'
        // fKnu = U sqt(S) V'
        // create global parameter for U and sqrt(S)
    }
    
    //Get a random and not correlated distribution for each variable, they might be different or the same
    fU_E = GetRandomDistribution(fE_dist);
    fU_nu = GetRandomDistribution(fnu_dist);
    
    //* Needs E and nu mean values and coef of variation or standard deviation */
    // f_E = U*sqrt(S)*fU_E * std_E + mean_E
    // f_nu = U*sqrt(S)*fU_nu * std_nu + mean_nu
    
    //***** Should use this one to get the random correlated and scaled field, in TPZMaterial use fE and fnu
    //GetStochasticField(fE, fnu);
    
    /* End */
    
    //Overwrtie fU just to test - This should be deleted
    // Multiplying decomposed Matrix M (U*Sqrt(S)) and random normal vector fRand_U
    fU = fM * fU_E; // Get correlated random distribution
    // In this function fM should be replaced by the left singular vetor U and the square root of the diagonal matrix S, then multiply by fRand_U
    
    return NULL;
}

void TPZStochasticMaterial::SetYoungField(int distribution, int function){
    fE_dist = distribution;
    fE_funct = function;
}

void TPZStochasticMaterial::SetPoissonField(int distribution, int function){
    fnu_dist = distribution;
    fnu_funct = function;
}

void TPZStochasticMaterial::SetReadMatrix(const TPZFMatrix<STATE> &M){
    fM = M;
}

void TPZStochasticMaterial::SetInclinedField(REAL rw, REAL rext,int stochasticInclined,REAL direction, REAL inclination)
{
    fstochasticInclined = stochasticInclined;
    fdirection = direction;
    finclination = inclination;
    frw = rw;
    frext = rext;
    InclinedFieldGeometry();
}

void TPZStochasticMaterial::InclinedFieldGeometry(){
    int nLayers = 8;
    fH = 2 * frext;
    fh = fH / nLayers;
    fmatsize = fnSquareElements * (fH/fh) + fnSquareElements;
}

TPZFMatrix<STATE>  TPZStochasticMaterial::EvaluateCorrelation(int function)
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


TPZFMatrix<STATE> TPZStochasticMaterial::GetRandomDistribution(int distribution)
{
    if (fstochasticInclined == 1) {
        fU = GetDistribution(fmatsize, distribution);
    }
    else {
        fU = GetDistribution(fnSquareElements, distribution);
    }
    return fU;
}

TPZFMatrix<STATE> TPZStochasticMaterial::GetDistribution(int matrixSize, int distribution)
{
    //normal distribution
    if(distribution==1) {
        
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution(0.,1.0);
        
        // Random Vector U
        TPZFMatrix<STATE> Rand_U (matrixSize, 1, 0.);
        
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
        TPZFMatrix<STATE> Rand_U (matrixSize, 1, 0.);
        
        for (int i = 0; i < matrixSize; i++) {
            Rand_U(i,0) = distribution(generator);
            distribution.reset();
            fRand_U = Rand_U;
        }
    }
    return fRand_U;
}

void TPZStochasticMaterial::GetStochasticField( TPZFMatrix<STATE> f_E, TPZFMatrix<STATE> f_nu)
{
    // Use global parameter for U and sqrt(S) for the field
    // This method needs the mean values of E and nu, as well as their coefficient of variatin or standard deviation
    // f_E = U*sqrt(S)*fU_E * std_E + mean_E
    // f_nu = U*sqrt(S)*fU_nu * std_nu + mean_nu
}

TPZFMatrix<STATE> TPZStochasticMaterial::calcCorrelationMatrix(int function) {
    
    std::cout << "\nCria matriz da norma entre os centroides (para a matriz de correlacao)" << std::endl;
    
    TPZFMatrix<STATE> CenterNorm(fnSquareElements, fnSquareElements, 0.0);
    
    TPZManVector<REAL, 3> CenterPoint1, CenterPoint2;
    
    // Element of reference
    TPZGeoEl *gel1;
    TPZVec<TPZGeoEl *> sub1;
    TPZManVector<REAL> centerpsi1(3), center1(3);
    
    // Other elements
    TPZGeoEl *gel2;
    TPZVec<TPZGeoEl *> sub2;
    TPZManVector<REAL> centerpsi2(3), center2(3);
    
    // Correlation Matrix
    TPZFMatrix<REAL> KCorr(fnSquareElements, fnSquareElements, 0.0);
    
    // Matrix of distance between centroids
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
            
            //    /*3*/    EQuadrilateral (element type)
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
                // Error in element type (not quadrilateral)
                std::cout<< "Element Type Error" << std::endl;
                DebugStop();
            }
        }
    }
    return KCorr;
}


TPZFMatrix<STATE> TPZStochasticMaterial::calcCorrelationMatrixInclined(int function) {
    
    std::cout << "\nCria matriz dos centroides dos elementos " << std::endl;
    
    TPZFMatrix<REAL> CenterNorm(fmatsize, fmatsize, 0.0);
    
    // Correlation Matrix
    TPZFMatrix<STATE> KCorr(fmatsize, fmatsize, 0.0);
    
    REAL Pi = M_PI;
    
    //Wellbore inclination
    REAL alpha = 0.; // azimuth
    REAL beta = 0.; // inclination
    alpha = (fdirection*(Pi/180));
    beta = (finclination*(Pi/180));
    
    //  Geeting all coordinates
    TPZGeoEl *gel;
    TPZFMatrix<STATE> Coordinates(fmatsize, 4, 0.0);
    TPZFMatrix<STATE> rotCoordinates(fnSquareElements, 4, 0.0);
    TPZManVector<REAL> centerpsi(3), center(3);
    TPZManVector<REAL, 3> CenterPoint;
    
    for (int i = 0; i < fnSquareElements; i++) {
        gel = fgmesh->ElementVec()[i];
        gel->CenterPoint(8, centerpsi);
        gel->X(centerpsi, center);
        
        CenterPoint = center;
        
        //    /*3*/    EQuadrilateral (element type)
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
    
    // Matrix of distance between centroids
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

