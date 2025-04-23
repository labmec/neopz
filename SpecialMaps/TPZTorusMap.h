#ifndef TPZTORUSMAP
#define TPZTORUSMAP

#include <pzfmatrix.h>
#include <pzmanvector.h>

namespace pzgeom{

  /**
     @brief Implements elements in a toroidal shell.
    It is based on storing information in the parametric coordinates
    phi, theta, along with the relevant parameters of the torus
    (small radius, large radius, origin).
    Currently, it only supports toruses whose axis are aligned with the z-axis.
    Supported values of TGeo are pzgeom::TPZGeoQuad and pzgeom::TPZGeoTriangle
  */
  template<class TGeo>
  class TPZTorusMap : public TGeo{
    //! Parametric coordinates of vertices
    TPZFNMatrix<TGeo::NNodes*2,REAL> fPhiTheta;
    //! Center of the torus
    TPZManVector<REAL,3> fOrigin;
    //! Large radius
    REAL fR{0};
    //! Small radius
    REAL fr{0};

  public:
    int ClassId() const override;

    static bool IsLinearMapping(int side)
    {
      return false;
    }
    /** @brief Constructor with list of nodes */
		TPZTorusMap(TPZVec<int64_t> &nodeindexes) :
      TGeo(nodeindexes), fOrigin(3,0.), fPhiTheta(2,TGeo::NNodes,0.)
		{
		}
		
		/** @brief Empty constructor */
		TPZTorusMap() : TGeo(), fOrigin(3,0.), fPhiTheta(2,TGeo::NNodes,0.)
		{
		}
		
		/** @brief Constructor with node map */
		TPZTorusMap(const TPZTorusMap &cp, std::map<int64_t,int64_t> & gl2lcNdMap) :
      TGeo(cp,gl2lcNdMap), fR(cp.fR), fr(cp.fr),
      fOrigin(cp.fOrigin), fPhiTheta(cp.fPhiTheta)
		{
		}
		
		/** @brief Copy constructor */
		TPZTorusMap(const TPZTorusMap &cp) = default;

    /** @brief Copy constructor for another mesh*/
    TPZTorusMap(const TPZTorusMap &cp, TPZGeoMesh &destmesh) :
      TGeo(cp,destmesh), fR(cp.fR), fr(cp.fr),
      fOrigin(cp.fOrigin), fPhiTheta(cp.fPhiTheta)
    {}
		    
    TPZTorusMap &operator=(const TPZTorusMap &cp) = default;
		
    void SetRadii(const REAL &R, const REAL &r)
    {
#ifdef PZDEBUG
      if (R < r) 
			{
        DebugStop();
      }
#endif
      fR = R;
      fr = r;			
    }
        
    void SetOrigin(const TPZVec<REAL> &origin)
    {
      fOrigin = origin;
    }

    /// compute the parametric coordinates of the corner nodes
    void ComputeCornerCoordinates(TPZGeoMesh &gmesh);

		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "TPZTorusMap";}
		
		/* @brief Computes the coordinate of a point given in parameter space */
        
    template<class T>
    void GradX(TPZFMatrix<REAL> &cornerco, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
    {
      TPZFNMatrix<6,T> DxDphi(3,2,0.), gradphi(2,2);
      TPZManVector<T,3> ft(3,0.);
      TGeo::X(fPhiTheta,par,ft);
      TGeo::GradX(fPhiTheta, par, gradphi);

      const auto phi = ft[0];
      const auto  theta = ft[1];
      DxDphi(0,0) = -sin(phi)*(fR + fr*cos(theta));
      DxDphi(0,1) = -fr*cos(phi)*sin(theta);
      DxDphi(1,0) = cos(phi)*(fR + fr*cos(theta));
      DxDphi(1,1) = -fr*sin(phi)*sin(theta);
      DxDphi(2,0) = 0;
      DxDphi(2,1) = fr*cos(theta);
      DxDphi.Multiply(gradphi, gradx);
    }

    template<class T>
		void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &result) const
    {
      TPZManVector<T,2> resloc(2);
      TGeo::X(this->fPhiTheta,loc,resloc);

      const auto phi = resloc[0];
      const auto theta = resloc[1];
      result[0] = fOrigin[0]+(fR + fr*cos(theta))*cos(phi);
      result[1] = fOrigin[1]+(fR + fr*cos(theta))*sin(phi);
      result[2] = fOrigin[2]+fr*sin(theta);
            
    }

		
    void Read(TPZStream& buf, void* context) override
    {
      DebugStop();
    }
        
    void Write(TPZStream &buf, int withclassid) const override
    {
      DebugStop();
		}

    static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
  };
};
#endif