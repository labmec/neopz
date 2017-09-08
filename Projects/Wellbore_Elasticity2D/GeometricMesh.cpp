//
//  CircularGeoMesh.cpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#include "GeometricMesh.hpp"

//********************************** Cria malha Geometrica Circular (360 graus) *********************************************************/

// Para gerar a malha eh necessario:
// rwb -> raio do poco
// re -> raio do contorno
// ncirc -> nro elem ao longo 1/4 do poco
// nrad -> nro elem da parede do poco ate contorno externo
// DrDcirc -> proporcao dos elementos da parede do poco

TPZGeoMesh *CircularGeoMesh (REAL rwb, REAL re, int ncirc, int nrad,
                             REAL DrDcirc, REAL alpha, REAL beta,
                             TPZFMatrix<REAL> GetKCorr) {
    
    
    // calcula comprimento radial do primeiro elemento
    REAL szmin;
    REAL Pi = M_PI;
    szmin = (2*Pi)*(rwb/ncirc)*(DrDcirc); //NNNNNNNNN
    
    // calcula comprimento radial da parede do poco ate contorno
    REAL radiallength;
    radiallength = re - rwb;
    
    // definindo variacao do angulo theta ao redor do poco
    // em rads!!!!!
    TPZVec<REAL> theta;
    theta.Resize(ncirc+1);
    REAL firsttheta = (2*Pi) / (ncirc);
    for (int k = 0; k<ncirc+1; k++) {
        REAL sumtheta = 0.;
        sumtheta += firsttheta * k;
        theta[k] = sumtheta;
    }
    
    
    //    // *******Imprime variacao dos angulos (em rads)
    //    std::cout<< "elementos de theta: " << endl;
    //    for (int t=0; t<ncirc+1; t++) {
    //        std::cout<< "Theta[" << t << "] :" << theta[t] << endl;
    //    }
    //    std::cout<< "Print theta " << endl;
    //    theta.Print();
    //    std::cout << endl;
    
    
    // nx = number of nodes in x direction
    // ny = number of nodes in y direction
    int nx,ny;
    nx = nrad+1;
    ny = ncirc+1;
    
    
    // Geometric Progression of the elements
    REAL q;
    if(nrad >1) {
        q = TPZGenGrid::GeometricProgression(szmin, radiallength, nrad);
    }
    else {
        q=radiallength;
    }
    
    // q = 1; // CALCULAR TAXAS DE CONVERGENCIA, malha com mesmo tamanho de elementos
    // std::cout<< "valor de q " << q << endl; // imprime razao da PG
    
    
    /* Creates the geometric mesh... The nodes and elements
     * will be inserted into mesh object during initilize process
     */
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    long i,j;
    long id, index;
    
    
    /****************** Malha Circunferencial (360 graus) *******************/
    
    //vector to store a coordinate
    TPZVec <REAL> coord (3,0.);
    TPZVec <REAL> coordT (3,0.);
    
    
    // aloca valor acumulado dos raios
    REAL rsum = 0.;
    REAL sz = 0.;
    
    if (q == 1) {
        sz = radiallength/nrad;
    }
    else {
        sz = szmin; // para taxas
    }
    
    //Nodes initialization
    for(i = 1; i < nx+1; i++){
        for(j = 1; j < ny+1; j++){
            // aloca coordenadas em cada no
            coord[0] = (rwb + rsum)* cos(theta[j-1]);
            coord[1] = (rwb + rsum)* sin(theta[j-1]);
            coord[2] = 0.;
            
            // id do elemento
            id = (i) * ny + (j);
            
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
            
            
            ////               Print
            // std::cout << "*****Iteracao nro: " << j << endl;
            // std::cout << "rsum: " << rsum << endl;
            // std::cout << "cos" << "[" << theta[j-1] << "]" <<": " ;
            // std::cout << cos(theta[j-1]);
            // std::cout << endl;
            // std::cout << "sin" << "[" << theta[j-1] << "]" << ": ";
            // std::cout << sin(theta[j-1]);
            // std::cout << endl;
            // std::cout << "Coord x: " << coord[0] << ";" << " Coord y: " << coord[1] << ";" << " Coord z: " << coord[2] << endl;
            // std::cout << endl;
            // std::cout << "alpha: " << alpha << ";" << " beta: " << beta  << endl;
            // std::cout << "CoordT x: " << coordT[0] << ";" << " CoordT y: " << coordT[1] << ";" << " CoordT z: " << coordT[2] << ";" << endl;
            // std::cout << endl;
        }
        
        if (q == 1) {
            rsum = sz  * i * q; // para taxas
        }
        else {
            rsum += sz; //valor acumulado dos raios
            sz *= q;
        }
    }
    
    //vector to store a element connectivities
    TPZVec <long> connect(4,0);
    
    //Element connectivities
    for(i = 0; i < (nx -1); i++){
        for(j = 0; j < (ny -1); j++){
            // index = (i)*(ny - 1)+ (j);
            connect[0] = (i) * ny + (j);
            connect[1] = connect[0]+(ny);
            connect[2] = connect[1]+1;
            connect[3] = connect[0]+1;
            //Allocates and define the geometric element
            gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
            gmesh->ElementVec()[id];
        }
    }
    
    // Conecta os nos dos ultimos elementos da circunferencia com os primeiros
    for (int k=0; k < nrad; k++) {
        TPZGeoEl *gel = gmesh->Element(((k+1)*(ncirc-1))+k);
        TPZGeoEl *gelAbove = gmesh->Element((((k+1)*(ncirc-1))+k)-(ncirc-1));
        TPZVec <int> gelAboveIndex(4,0);
        gelAbove->GetNodeIndices(gelAboveIndex);
        gel->SetNodeIndex(3, gelAboveIndex[0]);
        gel->SetNodeIndex(2, gelAboveIndex[1]);
    }
    
    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
    //********** Criando Geo de BC ********//
    
    // bc = -1 -> Normal Pressure condition na parede do poco
    for (int i = 0; i<ncirc; i++ ) {
        // TPZGeoElBC gbc(gmesh->ElementVec()[i],7,-1);
        gmesh->ElementVec()[i]->CreateBCGeoEl(7, -1);
    }
    
    // bc = -2 -> StressField condition circunferencia externa do farfield
    for (int i = 1; i<ncirc+1; i++ ) {
        gmesh->ElementVec()[(ncirc*nrad)-i]->CreateBCGeoEl(5, -2);
    }
    
    // bc -3 -> Mixed, ponto fixo canto externo do poco (farfield) bottom
    TPZGeoElBC gbc1(gmesh->ElementVec()[(ncirc*nrad)-1],2,-3);
    
    // bc -4 -> Mixed, ponto fixo canto externo do poco (farfield) lateral direita
    TPZGeoElBC gbc2(gmesh->ElementVec()[(ncirc*nrad)-(ncirc/4)],1,-4);
    
    // Cria matriz da norma entre os centroides (para a matriz de correlacao)
    
    TPZManVector<REAL> centerpsi(3), center(3);
    // Refinamento de elementos selecionados
    TPZGeoEl *gel;
    TPZVec<TPZGeoEl *> sub;
    
    REAL nelemtsr = nrad*ncirc; //Nro de elementos - quadrilateros
    
    TPZFMatrix<REAL> CenterNorm(nelemtsr,nelemtsr,0.0);
    
    TPZManVector<REAL,3> CenterPoint1;
    TPZManVector<REAL,3> CenterPoint2;
    
    // Matriz da distancia entre os centroides
    for (i = 0; i < nelemtsr; i++) {
        for (j = 0; j < nelemtsr; j++) {
            
            gel = gmesh->ElementVec()[i];
            gel->CenterPoint(8,centerpsi);
            gel->X(centerpsi,center);
            
            CenterPoint1 = center;
            
            gel = gmesh->ElementVec()[j];
            gel->CenterPoint(8,centerpsi);
            gel->X(centerpsi,center);
            
            CenterPoint2 = center;
            
            CenterNorm(i,j) = sqrt(pow((CenterPoint2[0]-CenterPoint1[0]),2)+pow((CenterPoint2[1]-CenterPoint1[1]),2)+pow((CenterPoint2[2]-CenterPoint1[2]),2));
        }
    }
    
    
    // Matriz de correlacao
    TPZFMatrix<REAL> KCorr(nelemtsr,nelemtsr,0.0);
    REAL r = 0.;        // define distancia r entre cada centroide
    REAL e = M_E;       //Numero de Euler
    REAL scale = 0.;    //Valor de alpha, escala normalizada
    
    for (i = 0; i < nelemtsr; i++) {
        for (j = 0; j < nelemtsr; j++) {
            r = CenterNorm(i,j);
            KCorr(i,j) = pow(e, (-scale*(pow(r, 2))));
        }
    }
    
    GetKCorr = KCorr;
    // std::cout << center[0] << ";" << center[1] << ";" << center[2] << ";" << std::endl;
    return gmesh;
}

// ******************************************** Cria malha Geometrica 1/4 do Poco *************************************************************/

// Para gerar a malha eh necessario:
// rwb -> raio do poco
// re -> raio do contorno
// ncirc -> nro elem ao longo 1/4 do poco
// nrad -> nro elem da parede do poco ate contorno externo
// DrDcirc -> proporcao dos elementos da parede do poco

TPZGeoMesh *QuarterGeoMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc) {
    
    
    // calcula comprimento radial do primeiro elemento
    REAL szmin;
    REAL Pi = 3.14159;
    szmin = (Pi/2)*(rwb/ncirc)*(DrDcirc);
    
    // calcula comprimento radial da parede do poco ate contorno
    REAL radiallength;
    radiallength = re - rwb;
    
    // definindo variacao do angulo theta ao redor do poco
    // em rads!!!!!
    TPZVec<REAL> theta;
    theta.Resize(ncirc+1);
    
    REAL firsttheta = ((Pi/2)/5) / (ncirc);
    for (int k = 0; k<ncirc+1; k++) {
        REAL sumtheta = 0.;
        sumtheta += firsttheta * k;
        theta[k] = sumtheta;
    }
    
    //
    //    // *******Imprime variacao dos angulos (em rads)
    //    std::cout<< "elementos de theta: " << endl;
    //        for (int t=0; t<ncirc+1; t++) {
    //            std::cout<< "Theta[" << t << "] :" << theta[t] << endl;
    //            }
    //            std::cout<< "Print theta " << endl;
    //            theta.Print();
    //            std::cout << endl;
    
    
    // nx = number of nodes in x direction
    // ny = number of nodes in y direction
    int nx,ny;
    nx = nrad+1;
    ny = ncirc+1;
    
    
    // Geometric Progression of the elements
    REAL q;
    if(nrad >1)
    {
        q = TPZGenGrid::GeometricProgression(szmin, radiallength, nrad);
    }
    else
    {
        q=radiallength;
    }
    // std::cout<< "valor de q " << q << endl; // imprime razao da PG
    
    //Creates the geometric mesh... The nodes and elements
    //will be inserted into mesh object during initilize process
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    long i,j;
    long id, index;
    
    
    //*********************** Malha Circunferencial (1/4) *********************//
    
    //vector to store a coordinate
    TPZVec <REAL> coord (2,0.);
    
    // aloca valor acumulado dos raios
    REAL rsum = 0.;
    REAL sz = szmin;
    
    //Nodes initialization
    for(i = 1; i < nx+1; i++){
        for(j = 1; j < ny+1; j++){
            
            // aloca coordenadas em cada no
            coord[0] = (rwb + rsum)* cos(theta[j-1]);
            coord[1] = (rwb + rsum)* sin(theta[j-1]);
            
            // id do elemento
            id = (i) * ny + (j);
            
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
            
            
            //            // print
            //
            //            std::cout << "*****Iteracao nro: " << j << endl;
            //            std::cout << "rsum: " << rsum << endl;
            //            std::cout << "cos" << "[" << theta[j-1] << "]" <<": " ;
            //            std::cout << cos(theta[j-1]);
            //            std::cout << endl;
            //            std::cout << "sin" << "[" << theta[j-1] << "]" << ": ";
            //            std::cout << sin(theta[j-1]);
            //            std::cout << endl;
            //            std::cout << "Coord x: " << coord[0] << ";" << " Coord y: " << coord[1] << ";" << " Rad: " << theta[j-1] << endl;
            //            std::cout << endl;
            
            
        }
        
        rsum += sz; //valor acumulado dos raios
        sz *= q;
    }
    
    
    
    //vector to store a element connectivities
    TPZVec <long> connect(4,0);
    
    
    //Element connectivities
    for(i = 0; i < (nx -1); i++){
        for(j = 0; j < (ny -1); j++){
            // index = (i)*(ny - 1)+ (j);
            connect[0] = (i) * ny + (j);
            connect[1] = connect[0]+(ny);
            connect[2] = connect[1]+1;
            connect[3] = connect[0]+1;
            //Allocates and define the geometric element
            gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
            //std::cout << "connect: " << connect << endl;
            
            //std::cout << "id: " << id << endl;
            
            gmesh->ElementVec()[id];
        }
    }
    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
    
    
    //********** Criando Geo de BC,********//
    
    // bc = -1 -> Normal Pressure condition
    for (int i = 0; i<ncirc; i++ ) {
        //        TPZGeoElBC gbc(gmesh->ElementVec()[i],7,-1);
        gmesh->ElementVec()[i]->CreateBCGeoEl(7, -1);
        
    }
    
    // bc = -2 -> Neumann condition contorno externo bottom
    for (int i = 0; i<nrad; i++ ) {
        gmesh->ElementVec()[ncirc*i]->CreateBCGeoEl(4, -2);
        
    }
    
    // bc = -3 -> Neumann condition contorno externo upper
    for (int i = 1; i<nrad+1; i++ ) {
        gmesh->ElementVec()[(ncirc*i)-1]->CreateBCGeoEl(6, -3);
        
    }
    
    // bc = -4 -> Neumann condition arco externo do farfield
    for (int i = 1; i<ncirc+1; i++ ) {
        gmesh->ElementVec()[(ncirc*nrad)-i]->CreateBCGeoEl(5, -4);
        
    }
    
    
    // bc -5 -> Mixed, ponto fixo canto interno parede do poco bottom
    TPZGeoElBC gbc1(gmesh->ElementVec()[0],0,-5);
    
    // bc -6 -> Mixed, ponto fixo externo ao 1/4 de poco bottom
    TPZGeoElBC gbc2(gmesh->ElementVec()[0],1,-6);
    
    // bc -5 -> Mixed, ponto fixo canto interno parede do poco bottom
    TPZGeoElBC gbc3(gmesh->ElementVec()[0],2,-7);
    
    // bc -6 -> Mixed, ponto fixo externo ao 1/4 de poco bottom
    TPZGeoElBC gbc4(gmesh->ElementVec()[0],3,-8);
    
    
    //    {
    //        std::ofstream malha("../malhageo.vtk");
    //        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malha);
    //    }
    
    return gmesh;
    
}
