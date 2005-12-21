//
// C++ Interface: tpzdxfdraw
//
// Description: 
//
//
// Author:  Luis Guilherme Mello Decourt <luis@labmec.fec.unicamp.br> , (C) 2004
//          Edimar Cesar Rylo <cesar@labmec.fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZ_DXFTPZDXFDRAW_H
#define PZ_DXFTPZDXFDRAW_H

#include <string>
#include <fstream>
class TPZGeoMesh;
class TPZGeoEl;


namespace pz_dxf {

/**
Writes a geometric mesh into autocad dxf format.

@author Luís Guilherme Mello Décourt
*/
class TPZDXFDraw{

  public:
    
    TPZDXFDraw(std::string &file, TPZGeoMesh *gmesh);

    ~TPZDXFDraw();
    
    /**
     * Write a geometric PZ Mesh into a autocad dxf file
     */
    void DXFDraw ();
    
    /**
     * Scale son elements for best visualization
     */
    void DXFDrawSep ();
    
  protected:
    std::ofstream fDXFFile;
    
    TPZGeoMesh *fMesh;
    
  protected :
    
    void PolyFaceMesh (int meshindex);
    
    void FaceMaker (int faces ,TPZGeoEl *gel , int el);
};

};

#endif
