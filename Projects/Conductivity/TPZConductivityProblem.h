//
//  TPZConductivityProblem.h
//  PZ
//
//  Created by Philippe Devloo on 5/14/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//
#ifndef CONDUCTIVITYPROBLEM
#define CONDUCTIVITYPROBLEM

#include "pzreal.h"
#include "pzsave.h"
#include "pzcmesh.h"
#include <string>

const int TPZCONDUCTIVITYID = 600;
/// Class to set up a conductivity problem in a simple way
/**
 * the idea of a problem set is to simplify the setup of a standard problem
 * and exemplify how to run a particular simulation
 * ALL parameters can be choosen default
 * the default configuration is used as unit test configuration
 */
class TPZConductivityProblem : public TPZSaveable
{
public:
    
    /// create a problem object with default configuration
    TPZConductivityProblem() : fDelx(2),fNx(2),fFluidFlux(0.)
    {
        SetDomainSize();
        SetMesh();
        SetBridgeVoidRatio();
        SetConductivity();
        SetPressureDifferential();
        SetGraphicsFile();
    }
    
    TPZConductivityProblem(const TPZConductivityProblem &cp) : fDelx(cp.fDelx), fNx(cp.fNx),
        fBridgeVoidRatio(cp.fBridgeVoidRatio), fConductivity(cp.fConductivity), 
        fDelPressure(cp.fDelPressure), fGraphicsFile(cp.fGraphicsFile), fFluidFlux(cp.fFluidFlux)
    {
    
    }
    
    /// set the size of the computational domain
    void SetDomainSize(REAL delx = 1., REAL dely = 1.)
    {
        fDelx[0] = delx;
        fDelx[1] = dely;
    }
    
    /// set the number of elements in the x and y direction
    void SetMesh(int nx = 10, int ny = 5)
    {
        fNx[0] = nx;
        fNx[1] = ny;
    }
    
    /// set the ratio between the area in contact and voids
    void SetBridgeVoidRatio(REAL ratio = 0.01)
    {
        fBridgeVoidRatio = ratio;
    }
    
    /// set the conductivity of the porous media
    void SetConductivity(REAL conductivity = 1)
    {
        fConductivity = conductivity;
    }
    
    /// set the pressure diferential of the stimulated fracture
    void SetPressureDifferential(REAL delp = 1.)
    {
        fDelPressure = delp;
    }
    
    /// set the file for post processing
    void SetGraphicsFile(std::string filename = "pressure.vtk")
    {
        fGraphicsFile = filename;
    }

    /// rturn the size of the computational domain in x
    void GetDomainSize(TPZVec<REAL> &delx)
    {
        delx = fDelx;
    }
        
    /// return the number of elements in the xdirection
    void GetMesh(TPZVec<int> &nx)
    {
        nx = fNx;
    }
    
    /// return the ratio between the area in contact and voids
    REAL GetBridgeVoidRatio()
    {
        return fBridgeVoidRatio;
    }
    
    /// return the conductivity of the porous media
    REAL GetConductivity()
    {
        return fConductivity;
    }
    
    /// return the pressure diferential of the stimulated fracture
    REAL GetPressureDifferential()
    {
        return fDelPressure;
    }
    
    /// set the file for post processing
    std::string GetGraphicsFile()
    {
        return fGraphicsFile;
    }

    /// set up the finite element mesh and compute the flux
    REAL ComputeFlux();
    
    /// Method to compare the current object with a copy
    virtual bool Compare (TPZSaveable *copy, bool override=false);

    /// Method to compare the current object with a copy
    virtual bool Compare (TPZSaveable *copy, bool override=false) const;

    virtual int ClassId () const;

    /// write this object to the TPZStream buffer. Include the classid if withclassid = true
    virtual void Write(TPZStream &buf, int withclassid);
    
    /// read objects from the stream
    virtual void Read(TPZStream &buf, void *context);
    
    /// Generate a computational mesh according to the specification of the problem
    TPZAutoPointer<TPZCompMesh> GenerateCompMesh();
    
    /// compute the perimeter of all two dimensional elements
static REAL Perimeter(TPZGeoMesh &gmesh);
    
    /// compute the perimeter of all two dimensional elements
static REAL DomainArea(TPZGeoMesh &gmesh);
    

    
private:
    
    /// compute the flux over the right side
    REAL MeshFlux(TPZAutoPointer<TPZCompMesh> cmesh);
    
    /// size of the domain
    TPZManVector<REAL,2> fDelx;
    
    /// number of elements in the x and y direction (the fluid flows in the x direction
    TPZManVector<int,2> fNx;
    
    /// bridge void ratio
    REAL fBridgeVoidRatio;
    
    /// conductivity of the porous media
    REAL fConductivity;
    
    /// pressure differential
    REAL fDelPressure;
    
    /// name of the graphics post processing file
    std::string fGraphicsFile;
    
    /// computed fluid flux
    REAL fFluidFlux;
};



#endif