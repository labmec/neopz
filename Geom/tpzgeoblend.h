/**
 * @file
 * @brief Contains the TPZGeoBlend class which implements a blending map from curved boundaries to the interior of the element.
 */

#ifndef TPZGEOBLEND_H
#define TPZGEOBLEND_H

#include "pzgeotriangle.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeoel.h"
#include "pznoderep.h"
#include "pzgeoelside.h"

#include <iostream>

namespace pzgeom {

    /**
     * @ingroup geometry
     * @brief Implements a blending map from curved boundaries to the interior of the element. \ref geometry "Geometry"
     * @author Paulo Cesar de Alvarenga Lucci
     * @since 2007
     */
    template<class TGeo>
    class TPZGeoBlend : public TGeo {

    public:

        int ClassId() const override;

        bool IsLinearMapping(int side) const;

        bool IsGeoBlendEl() const {
            return true;
        }
        typedef typename TGeo::Top Top;
        /** @brief Constructor with list of nodes */
        TPZGeoBlend(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoBlend::ClassId),
                                                    TGeo(nodeindexes) {
        }

        /** @brief Empty constructor */
        TPZGeoBlend() : TPZRegisterClassId(&TPZGeoBlend::ClassId),
                        TGeo() {
        }

        /** @brief Constructor with node map */
        TPZGeoBlend(const TPZGeoBlend &cp, std::map<int64_t, int64_t> &gl2lcNdMap) : TPZRegisterClassId(
                &TPZGeoBlend::ClassId), TGeo(cp, gl2lcNdMap) {
             std::cout << "Please implement me\n";
             DebugStop();
        }

        /** @brief Copy constructor */
        TPZGeoBlend(const TPZGeoBlend &cp) : TPZRegisterClassId(&TPZGeoBlend::ClassId), TGeo(cp) {
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fNeighbours[is] =  cp.fNeighbours[is];
                fTrans[is] =  cp.fTrans[is];
            }
            fGeoEl = 0;
            DebugStop();
        }

        /** @brief Copy constructor */
        TPZGeoBlend(const TPZGeoBlend &cp, TPZGeoMesh &destmesh);

        void Read(TPZStream &buf, void *context) override{
            TGeo::Read(buf, context);
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fNeighbours[is].Read(buf, context);
            }
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fTrans[is].Read(buf, context);
            }
        }

        void Write(TPZStream &buf, int withclassid) const override{
            TGeo::Write(buf, withclassid);
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fNeighbours[is].Write(buf, withclassid);
            }
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fTrans[is].Write(buf, withclassid);
            }
        }

        void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform<> &trans) override;

        bool ResetBlendConnectivity(const int64_t &side, const int64_t &index);

        TPZGeoElSide Neighbour(int side, TPZGeoMesh *gmesh) const {
            if (side < TGeo::NNodes) {
                DebugStop();
            }
            return TPZGeoElSide(fNeighbours[side - TGeo::NNodes], gmesh);
        }

        template<class T>
        void TransfBetweenNeigh(int side, TPZTransform<T> &tr) const {
            tr.CopyFrom(fTrans[side - TGeo::NNodes]);
        }

        /** @brief Returns the type name of the element */
        static std::string TypeName() { return TGeo::TypeName(); }

        /** @brief Get the coordinates of the point at geometric elements from coordinates of the parametric point at the master element */
        template<class T>
        void X(TPZFMatrix<REAL> &cornerco, TPZVec<T> &par, TPZVec<T> &result) const;

        /** @brief Computes the Jacobian for parametric point at master element */
        void Jacobian(TPZFMatrix<REAL> &cornerco, TPZVec<REAL> &par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                      REAL &detjac, TPZFMatrix<REAL> &jacinv) const;

        /** @brief Computes the gradient of the transformation for parametric point at master element */
        template<class T>
        void GradX(TPZFMatrix<REAL> &cornerco, TPZVec<T> &par, TPZFMatrix<T> &gradx) const;

        /** @brief Print all relevant data of the element to cout*/
        void Print(std::ostream &out = std::cout) const;

        /**
         * @brief Initialize the element checking connectivity on all sides.
         *
         / Este método percorre todos os lados do elemento fornecido (el) e, para cada lado,
         / percorre todos os vizinhos até retornar a si próprio. Se neste processo encontrar algum elemento que
         / seja quadrático e não seja do tipo TPZGeoBlend, recorre ao método SetNeighbourInfo() para estabelecer que este
         / elemento encontrado será seu vizinho pelo respectivo lado.
         */
        void Initialize(TPZGeoEl *refel);

        //void Initialize(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh);

        // /**
        //  * @brief Method which creates a geometric boundary condition
        //  * element based on the current geometric element,
        //  * a side and a boundary condition number
        //  */
        // TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig, int side, int bc);

        // TPZGeoEl *CreateBCGeoBlendEl(TPZGeoEl *orig, int side, int bc);

//		TPZGeoEl *CreateGeoBlend(TPZGeoMesh &mesh, MElementType type, TPZVec<int64_t>& nodeindexes, int matid, int64_t& index);



    public:

        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

        /**
         * @brief Creates a geometric element according to the type of the father element
         */
        // static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
        //                                   TPZVec<int64_t> &nodeindexes,
        //                                   int matid,
        //                                   int64_t &index);

    protected:

        /// Project the InternalPar parameter to the parameter of the neighbour along side. Return true if the map is nonsingular
        template<class T>
        bool MapToNeighSide(int side, int sidedim, TPZVec<T> &InternalPar, TPZVec<T> &NeighPar,
                            TPZFMatrix<T> &JacNeighSide) const;
        /** @brief Vector of indexes of the neighbours */
        //TPZGeoElSideIndex fNeighbours[1+TGeo::NSides - TGeo::NNodes];
        TPZGeoEl *fGeoEl;
        TPZGeoElSideIndex fNeighbours[1 + TGeo::NSides - TGeo::NNodes];
        TPZTransform<> fTrans[1 + TGeo::NSides - TGeo::NNodes];
    };


    template<class TGeo>
    template<class T>
    inline bool
    pzgeom::TPZGeoBlend<TGeo>::MapToNeighSide(int side, int SideDim, TPZVec<T> &InternalPar, TPZVec<T> &NeighPar,
                                              TPZFMatrix<T> &JacNeighSide) const {
        TPZFNMatrix<9, T> JacSide;

        TPZManVector<T, 3> SidePar(SideDim);
        const bool check = TGeo::CheckProjectionForSingularity(side, InternalPar);

        this->MapToSide(side, InternalPar, SidePar, JacSide);
        if (!check) {
            #ifdef PZ_LOG2
            td::stringstream sout;
            sout << "side " << side << std::endl;
            sout << "InternalPar: ";
            for(int i = 0; i < InternalPar.NElements(); i++) sout << InternalPar[i] << "\t";
            sout << "\n";
            sout << "\tmapping is not regular"<<std::endl;
            LOGPZ_DEBUG(logger,sout.str())
            #endif
            return false;
        }
#ifdef PZ_LOG2
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "side " << side << std::endl;
            sout << "InternalPar: ";
            for(int i = 0; i < InternalPar.NElements(); i++) sout << InternalPar[i] << "\t";
            sout << "\n";

            sout << "SidePar: ";
            for(int i = 0; i < SidePar.NElements(); i++) sout << SidePar[i] << "\t";
            sout << "\n";

            JacSide.Print("JacSide = ",sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif

        NeighPar.Resize(SideDim);
        TPZTransform<T> tr;
        TransfBetweenNeigh(side, tr);
        tr.Apply(SidePar, NeighPar);

        JacNeighSide.Resize(0, 0);
        TransfBetweenNeigh(side, tr);
        JacNeighSide.Resize(tr.Mult().Rows(), JacSide.Cols());
        tr.Mult().Multiply(JacSide, JacNeighSide);

#ifdef PZ_LOG2
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            
            sout << "NeighPar: ";
            for(int i = 0; i < NeighPar.NElements(); i++) sout << NeighPar[i] << "\t";
            sout << "\n";
            
            JacNeighSide.Print("JacNeighSide = ",sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif

        return true;
    }


    template<class TGeo>
    int TPZGeoBlend<TGeo>::ClassId() const {
        return Hash("TPZGeoBlend") ^ TGeo::ClassId() << 1;
    }

};
#endif
