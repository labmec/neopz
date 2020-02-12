//
//  TPZWavyLine.h
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#ifndef __PZ__TPZWavyLine__
#define __PZ__TPZWavyLine__

#include <iostream>
#include "TPZGeoLinear.h"

namespace pzgeom {
    
    class TPZWavyLine : public TPZGeoLinear
    {
        int fNumWaves;
        TPZManVector<REAL,3> fWaveDir;

    public:
			
        int ClassId() const override;

       
        /** @brief Constructor with list of nodes */
		TPZWavyLine(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZWavyLine::ClassId),
        TPZGeoLinear(nodeindexes), fNumWaves(0), fWaveDir()
		{
		}
		
		/** @brief Empty constructor */
		TPZWavyLine() : TPZRegisterClassId(&TPZWavyLine::ClassId),
        TPZGeoLinear(), fNumWaves(0), fWaveDir()
		{
		}
		
		/** @brief Constructor with node map */
		TPZWavyLine(const TPZWavyLine &cp,
				   std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZWavyLine::ClassId),
        TPZGeoLinear(cp,gl2lcNdMap), fNumWaves(cp.fNumWaves), fWaveDir(cp.fWaveDir)
		{
		}
		
		/** @brief Copy constructor */
		TPZWavyLine(const TPZWavyLine &cp) : TPZRegisterClassId(&TPZWavyLine::ClassId),
        TPZGeoLinear(cp), fNumWaves(cp.fNumWaves), fWaveDir(cp.fWaveDir)
		{
		}
		
		/** @brief Copy constructor */
		TPZWavyLine(const TPZWavyLine &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZWavyLine::ClassId),
        TPZGeoLinear(cp), fNumWaves(cp.fNumWaves), fWaveDir(cp.fWaveDir)
		{
		}
        
        void SetData(TPZVec<REAL> &wavedir, int numwaves)
        {
#ifdef PZDEBUG
            if (wavedir.size() != 3) {
                DebugStop();
            }
#endif
            fWaveDir = wavedir;
            fNumWaves = numwaves;
        }
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }

		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Wavy";}
		
		/* @brief Computes the coordinate of a point given in parameter space */
//        template<class T>
//        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            X(coord,loc,result);
//        }
        
//        template<class T>
//        void GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
//        {
//            gradx.Resize(3,1);
//            gradx.Zero();
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            T cosval = cos(fNumWaves*M_PI*par[0]);
//            for (int i=0; i<3; i++) {
//                gradx(i) = (coord(i,1)-coord(i,0))/2. + fNumWaves*M_PI*cosval*fWaveDir[i];
//            }
//        }
        template<class T>
        void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx) const {

            int nrow = nodes.Rows();
            int ncol = nodes.Cols();

            gradx.Resize(nrow,1);
            gradx.Zero();

            T cosval = cos(fNumWaves*M_PI*loc[0]);
            for (int i=0; i<3; i++) {
                gradx(i) = (nodes.GetVal(i,1)-nodes.GetVal(i,0))/2. + fNumWaves*M_PI*cosval*fWaveDir[i];
            }

        }
        /* @brief Computes the jacobian of the map between the master element and deformed element */
//        void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
//        {
//            jacobian.Resize(1,1); axes.Resize(1,3); jacinv.Resize(1,1);
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            TPZManVector<REAL,3> gradx(3);
//            CornerCoordinates(gel, coord);
//            REAL cosval = cos(fNumWaves*M_PI*param[0]);
//            for (int i=0; i<3; i++) {
//                gradx[i] = (coord(i,1)-coord(i,0))/2. + fNumWaves*M_PI*cosval*fWaveDir[i];
//            }
//            REAL normv = 0.;
//            for (int i=0; i<3; i++) {
//                normv += gradx[i]*gradx[i];
//            }
//            normv = sqrt(normv);
//            jacobian(0,0) = normv;
//            detjac = normv;
//            jacinv(0,0) = 1./normv;
//            for (int i=0; i<3; i++) {
//                axes(0,i) = gradx[i]/normv;
//            }
//        }
        
        template<class T>
        void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZGeoLinear::X(nodes,loc,result);
            T sinval = sin(this->fNumWaves*M_PI*loc[0]);

            for (int i=0; i<3; i++) {
                result[i] += this->fWaveDir[i]*sinval;
            }

        }
		
		
		// static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

		// /** @brief Creates a geometric element according to the type of the father element */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);
		
        void Read(TPZStream& buf, void* context) override {
            pzgeom::TPZGeoLinear::Read(buf,0);
            buf.Read(&fNumWaves,1);
            buf.Read<3>( fWaveDir);
        }
        
        void Write(TPZStream &buf, int withclassid) const override{
            pzgeom::TPZGeoLinear::Write(buf, withclassid);
            buf.Write(&fNumWaves,1);
            buf.Write( fWaveDir);
		}

        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        

	};

    
}

#endif /* defined(__PZ__TPZWavyLine__) */
