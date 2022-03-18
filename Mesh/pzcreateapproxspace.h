/**
 * @file
 * @brief Contains the functions to create different computational elements (one- two- three-dimensional). The creation of SBFem elements depends on the NeoPZ library configuration USING_MKL=ON.
 */

#ifndef CREATECONTINUOUSHPP
#define CREATECONTINUOUSHPP

class TPZCompEl;
class TPZCompMesh;
class TPZGeoEl;
class TPZCompEl;
class TPZCompMesh;
#include <set>
#include <functional>
#include "pzvec.h"
#include "TPZSavable.h"
#include "TPZEnumApproxFamily.h"

typedef std::function<TPZCompEl* (TPZGeoEl* el, TPZCompMesh &mesh)> TCreateFunction;
/*
 * @brief Administer the creation of approximation spaces
 * @author Philippe Devloo
 * @since 2009
 * @ingroup interpolation
 */
class TPZCreateApproximationSpace : public TPZSavable {
    /** @brief Function pointer which determines what type of computational element will be created */
    TCreateFunction fp[8];
    
    /// @brief boolean indicating if each element should be created disconnected from the others
    /**
     * this flag allows to create hybrid meshes (default is false)
     */
    bool fCreateHybridMesh;
    
    /// flag indicating whether each element should have an aditional lagrange multiplier
    bool fCreateLagrangeMultiplier;
    
    /// flag indicating that the elements need to be created with memory
    bool fCreateWithMemory;
    
    /// flags indicating "flavor" ofthe approximation space
    HDivFamily fhdivfam = DefaultFamily::fHDivDefaultValue;
    H1Family fh1fam = DefaultFamily::fH1DefaultValue;
    HCurlFamily fhcurlfam = DefaultFamily::fHCurlDefaultValue;
    
public:
    
    enum MApproximationStyle {ENone,EContinuous,EDiscontinuous,EHDiv,EHCurl, EMultiphysics, EMultiphysicsSBFem, ESBFem, ECustom};
private:
    /// approximation space style last used
    MApproximationStyle fStyle = ENone;

public:
    
    
    TPZCreateApproximationSpace() : fCreateHybridMesh(false), fCreateLagrangeMultiplier(false), fCreateWithMemory(false)
    {
        SetAllCreateFunctionsContinuous();
    }
    
    TPZCreateApproximationSpace(const TPZCreateApproximationSpace &copy) : fCreateHybridMesh(copy.fCreateHybridMesh), fCreateLagrangeMultiplier(copy.fCreateLagrangeMultiplier)
    ,fCreateWithMemory(copy.fCreateWithMemory)
    {
        for (int i=0; i<8; i++) {
            fp[i] = copy.fp[i];
        }
    }
    
    TPZCreateApproximationSpace &operator=(const TPZCreateApproximationSpace &copy)
    {
        for (int i=0; i<8; i++) {
            fp[i] = copy.fp[i];
        }
        fCreateHybridMesh = copy.fCreateHybridMesh;
        fCreateLagrangeMultiplier = copy.fCreateLagrangeMultiplier;
        fCreateWithMemory = copy.fCreateWithMemory;
        return *this;
    }
        int ClassId() const override;
    
    void Read(TPZStream &buf, void *context) override;
    
    void Write(TPZStream &buf, int withclassid) const override;
    
    void SetCreateLagrange(bool flag)
    {
        fCreateLagrangeMultiplier = flag;
    }
    
    void CreateWithMemory(bool flag)
    {
        fCreateWithMemory = flag;
    }

    // Get set methods for space families
    const HDivFamily &HDivFam() const {return fhdivfam;}
    const HDivFamily &HDivFam() {return fhdivfam;}
    const void SetHDivFamily(HDivFamily fam){fhdivfam = fam;}

    const H1Family &H1Fam() const {return fh1fam;}
    const H1Family &H1Fam() {return fh1fam;}
    const void SetH1Family(H1Family fam){fh1fam = fam;}

    const HCurlFamily &HCurlFam() const {return fhcurlfam;}
    const HCurlFamily &HCurlFam() {return fhcurlfam;}
    const void SetHCurlFamily(HCurlFamily fam){fhcurlfam = fam;}
    
    /** @brief Create discontinuous approximation spaces */
    void SetAllCreateFunctionsDiscontinuous();
    /** @brief Create continuous approximation spaces */
	void SetAllCreateFunctionsContinuous();
    /** @brief Create an approximation space with HDiv elements */
	void SetAllCreateFunctionsHDiv(int meshdim);
    /** @brief Create an approximation space with HCurl elements */
    void SetAllCreateFunctionsHCurl(int meshdim);
	/** @brief Create an approximation space with HDiv elements and full basis for quadrilateral element */
//    void SetAllCreateFunctionsHDivFull(int meshdim);
    
    /** @brief Create SBFem approximation space. Depends on MKL. */
    void SetAllCreateFunctionsSBFem(int meshdim);
    /** @brief Create an approximation space based on SBFem multiphysics elements */
    void SetAllCreateFunctionsSBFemMultiphysics(int meshdim);

#ifndef STATE_COMPLEX
    /** @brief Create an approximation space with HDivxL2 elements */
	void SetAllCreateFunctionsHDivPressure(int meshdim);
#endif
    /** @brief Create approximation spaces corresponding to the space defined by cel */
	void SetAllCreateFunctions(TPZCompEl &cel, TPZCompMesh *mesh);
    /** @brief Create an approximation space based on multiphysics elements */
	void SetAllCreateFunctionsMultiphysicElem();
    /** @brief Create an approximation space based on multiphysics elements with memory*/	
	void SetAllCreateFunctionsMultiphysicElemWithMem();
    /** @brief Create an approximation space with continous elements with memory. Only dimension 3 elements quem have memory in viscoelastic materials
		 @ param dimension dimension of the mesh
		 */	
    void SetAllCreateFunctionsContinuousWithMem();
    
    /** @brief Set custom function pointers */
    void SetCreateFunctions(TPZVec<TCreateFunction> &createfuncs);

    MApproximationStyle& Style()
    {
        return fStyle;
    }

    const MApproximationStyle& Style() const
    {
        return fStyle;
    }
    
    /** @brief Create a computational element using the function pointer for the topology */
    TPZCompEl *CreateCompEl(TPZGeoEl *gel, TPZCompMesh &mesh) const;
    
	/** @brief Creates the computational elements, and the degree of freedom nodes */ 
	/** Only element of material id in the set<int> will be created */
	void BuildMesh(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs) const;
	
	/** @brief Creates the computational elements, and the degree of freedom nodes */
	void BuildMesh(TPZCompMesh &cmesh) const;
    
    /** @brief Creates the computational elements, and the degree of freedom nodes */
    void BuildMesh(TPZCompMesh &cmesh, const TPZVec<int64_t> &gelindexes) const;
    
    /** @brief Creates the interface elements */
	/** Only element of material id in the set<int> will be created */
	static void CreateInterfaces(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs);
	
	/** @brief Creates the interface elements */
	static void CreateInterfaces(TPZCompMesh &cmesh);

	/** @brief Creates the computational elements, and the degree of freedom nodes */
	/**
	 * Elements created may be TPZInterpolatedElement or TPZCompElDisc. \n
	 * indices contains the type of the element. Element type are given by the enumerate MCreationType.
	 */
	static void AutoBuildContDisc(const TPZVec<TPZGeoEl*> &continuous, const TPZVec<TPZGeoEl*> &discontinuous);
    
    /** @brief Encapsulate the elements in condensed computational elements */
    static void CondenseLocalEquations(TPZCompMesh &cmesh);
    
    /** @brief Undo the encapsulate elements */
    static void UndoCondenseLocalEquations(TPZCompMesh &cmesh);
    
    /** @brief transform in low order Raviar Tomas */
    static void MakeRaviartThomas(TPZCompMesh &cmesh);
    
    /** @brief transform in low order Raviar Tomas */
    static void UndoMakeRaviartThomas(TPZCompMesh &cmesh);
    
    /** @brief Create interface elements between the computational elements */
    static void CreateInterfaceElements(TPZCompMesh *mesh, bool onlydiscontinuous = true, bool multiphysics = false);
    
    /** @brief Determine if the mesh will be created with disconnected elements
     * After the mesh is created, interface elements need to be created "by hand"
     */
    void CreateDisconnectedElements(bool create)
    {
        fCreateHybridMesh = create;
    }
    
    bool NeedsMemory()
    {
        return fCreateWithMemory;
    }
    
    /// this method will substitute all interface elements with materialid within the set by three elements : one H1 element and two interface elements
    static void Hybridize(TPZCompMesh &cmesh,const std::set<int> &matids, bool isconnectedElem = false);
};

#endif
