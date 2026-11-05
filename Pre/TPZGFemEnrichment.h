#ifndef TPZGFEMENRICHMENT_H
#define TPZGFEMENRICHMENT_H

#include "pzvec.h"
#include "pzmanvector.h"
#include "TPZSavable.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include <array>
#include <string>
#include <functional>

#include "TPZBuildSBFem.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;


/**
 * @brief Generate local SBFEM crack-tip modes for use as GFEM enrichments.
 * Evaluate tip-centered displacement modes and derivatives with optional basis changes.
 * @ingroup analysis
 */
class TPZGFemEnrichment {

public:
    /// Requested crack-side family for sector selection near the angular cut.
    /// EUpper/ELower override ambiguous side selection near a crack face;
    /// EAny/ENone leave selection to the point's angle.
    enum MSectorFamily {ENone, EUpper, ELower, EAny};

private:
    /// Radius of the local circular SBFEM domain, used in ksi = r/fRadius.
    /// This is not the radius selecting enriched nodes in the global GFEM mesh.
    REAL fRadius = 0.;
    
    /// Geometric mesh of the local circular SBFEM domain.
    TPZAutoPointer<TPZGeoMesh> fGeoMesh = 0;
    
    /// Local computational mesh containing the SBFEM approximation.
    TPZCompMesh *fCompMesh = 0;

    /// SBFEM element group storing the radial spectrum and boundary modes.
    TPZSBFemElementGroup *fElementGroup = 0;

    /// How a selected eigenmode contributes a real enrichment: a real mode,
    /// the real part of a complex mode, or its imaginary part.
    enum EPartOfEigenvector {ENonInitialized,EOnlyReal, ERealPart, EImagPart};
    /// Construct real enrichments from one complex representative or from
    /// explicitly combined conjugate eigenmodes.
    enum EComplexModeReconstruction {ERepresentativeParts, EConjugateSum};
    /// Eigenvector index for each output enrichment. A complex representative
    /// is listed twice to supply its two real-valued functions.
    TPZVec<int> fEigenvectorIndexes;
    /// Conjugate partner for a complex mode, or -1 for a real mode.
    TPZVec<int> fConjugateEigenvectorIndexes;
    /// Real/imaginary extraction tag for each entry of fEigenvectorIndexes.
    TPZVec<EPartOfEigenvector> fPartOfEigenvector;
    /// Selected complex-pair reconstruction policy; representative parts
    /// are used unless the JSON configuration requests conjugate sums.
    EComplexModeReconstruction fComplexModeReconstruction = ERepresentativeParts;

    /// Optional basis change among selected real-valued enrichments.
    /// Rows define output modes; columns multiply the original modes.
    std::vector<std::vector<REAL>> fModeTransform;
    /// Enable the stored basis transformation for evaluated mode values
    /// and Cartesian derivatives; false keeps the original modal basis.
    bool fUseModeTransform = false;

    /// Material ID for each input angular sector, before skeleton refinement.
    TPZVec<int> fMaterialIds;

    /// Geometric/material ID identifying the circular boundary arcs.
    int fArcMatId = -1;

    /// Material ID assigned to skeleton elements by the SBFEM builder.
    int fSkeletonMatId = 4;

    /// Point-element ID identifying the crack tip and singular partition.
    int fPointMatId = 5;


    /// Offset identifying auxiliary triangular-sector material IDs before
    /// the builder maps them to the physical SBFEM materials.
    int fMaterialIdOffset = 1000;

    /// Ordered sector-endpoint angles in radians, spanning exactly 2*pi.
    /// Initially read from JSON; expanded when the skeleton is refined.
    TPZManVector<REAL,10> fAngles;

    /// Sector index defining the upper-family angular bound; this is not
    /// a material ID and does not imply an upper/lower material assignment.
    int fMaxUpper;

    /// Sector index defining the lower-family angular bound used to resolve
    /// crack-side selection near the angular cut.
    int fMinLower;

    /// Boundary interpolation elements, stored in the same order as sectors.
    TPZVec<TPZInterpolatedElement *> fSkeletonElements;

    /// Locate an angle in the ordered sector intervals. The requested
    /// family resolves crack-face ambiguity; endpoint indices are clamped.
    int Sector(REAL angle, MSectorFamily sector) const {
        int nangles = fAngles.size();
        auto it = std::upper_bound(&fAngles[0], &fAngles[0]+nangles, angle);
        int pos = it-&fAngles[0]-1;
        if(sector == ELower && angle <= MinLowerAngle()) {
            pos = nangles-2;
        }
        if(sector == EUpper && angle >= MaxUpperAngle()) {
            pos = 0;
        }
        // std::atan2 may return either endpoint of the closed angular
        // interval. Both endpoints belong to a valid adjacent sector.
        pos = std::max(0, std::min(pos, nangles-2));
        return pos;
    }
    
    /// Map an angle in the chosen sector to its boundary coordinate eta.
    /// The sector endpoints correspond to eta = -1 and eta = +1.
    REAL Eta(int sector, REAL angle) const {
#ifdef PZDEBUG
        int nangles = fAngles.size();
        if(sector < 0 || sector >= nangles-1) {
            DebugStop();
        }
#endif
        REAL angle1 = fAngles[sector];
        REAL angle2 = fAngles[(sector+1)];
        return 2.*(angle - angle1)/(angle2 - angle1)-1.;
    }

    /// Compute the in-plane distance from the scaling center at (0,0).
    /// The supplied coordinates must already be in the generator's frame.
    REAL R(const std::vector<REAL> &x) const {
        REAL r = std::sqrt(x[0]*x[0]+x[1]*x[1]);
        return r;
    }

    /// Convert physical radius to the dimensionless scaling coordinate
    /// ksi = r/fRadius, with ksi = 1 on the circular boundary.
    REAL Ksi(REAL r) const {
        return r/fRadius;
    }

    /// Fill the Cartesian mapping matrix at a sector point (ksi,eta).
    /// Rows are x,y; columns are radial and angular mapping derivatives.
    /// The radial column includes the reference-coordinate factor 1/2.
    void gradX(REAL ksi, REAL eta, int sector, TPZFMatrix<REAL> &gradxieta) const {
        REAL dangle = (fAngles[sector+1]-fAngles[sector]);
        REAL angle = fAngles[sector] + (eta+1.)*dangle/2.;
        REAL cosa = std::cos(angle);
        REAL sina = std::sin(angle);
        REAL r = ksi*fRadius;
        gradxieta(0,0) = fRadius*cosa/2.;
        gradxieta(1,0) = fRadius*sina/2.;
        gradxieta(0,1) = -r*sina*dangle/2.;
        gradxieta(1,1) = r*cosa*dangle/2.;
    }
    
public:

    /// Callable interface for enrichment evaluation at a Cartesian point.
    /// Returns two displacement components and four Cartesian derivatives
    /// per real-valued mode, using the requested crack-side family.
    typedef std::function<void(TPZGFemEnrichment::MSectorFamily family, const std::vector<REAL> &x, 
                                   std::vector<std::array<REAL, 2>> &phi, 
                                   std::vector<std::array<REAL ,4>> &dphixy)> GFemEnrichmentFunction;
    /// Default constructor
    TPZGFemEnrichment() = default;
    
    /// Copy constructor
    TPZGFemEnrichment(const TPZGFemEnrichment &copy) = default;
    
    /// Assignment operator
    TPZGFemEnrichment &operator=(const TPZGFemEnrichment &copy) = default;
    
    /// Destructor
    virtual ~TPZGFemEnrichment();
    
    /// Initialize geometry, materials, SBFEM spectrum and selected modes
    /// from an already parsed JSON object; optionally build a modal transform.
    void LoadFromJSON(json &input);
    
    
    /// Evaluate selected real modes at x in the tip-centered generator frame.
    /// Overwrites phi with {ux,uy} and dphi with {dux/dx,dux/dy,duy/dx,duy/dy};
    /// applies the configured complex reconstruction and modal transform.
    void ComputeEnrichmentWithDerivatives(TPZGFemEnrichment::MSectorFamily family, const std::vector<REAL> &x, 
                                          std::vector<std::array<REAL, 2>> &phi, 
                                          std::vector<std::array<REAL ,4>> &dphi);
    
    /// Return a callback to ComputeEnrichmentWithDerivatives, capturing this.
    /// The enrichment object must outlive all uses of the callback.
    GFemEnrichmentFunction GetEnrichmentFunction() {
        return [this](TPZGFemEnrichment::MSectorFamily family, const std::vector<REAL> &x, 
                                   std::vector<std::array<REAL, 2>> &phi, 
                                   std::vector<std::array<REAL ,4>> &dphixy) {
            this->ComputeEnrichmentWithDerivatives(family, x, phi, dphixy);
        };
    }

    /// Alternative evaluator using loaded mesh eigenvectors and volume
    /// solution recovery. Appends results; GetEnrichmentFunction instead
    /// calls ComputeEnrichmentWithDerivatives.
    void EnrichmentFunction(MSectorFamily family, const std::vector<REAL> &x, 
                                   std::vector<std::array<REAL, 2>> &phi, 
                                   std::vector<std::array<REAL ,4>> &dphixy);
    /// Set the radius used for radial scaling; does not rebuild the meshes.
    void SetRadius(REAL radius) { fRadius = radius; }
    
    /// Return the radius of the local circular SBFEM domain.
    REAL Radius() const { return fRadius; }

    /// Forward crack-axis angle in radians, opposite the first crack-face ray.
    REAL CrackAxisAngle() const { return fAngles[0]+M_PI; }
    
    /// Return the angular threshold used by upper-family sector selection.
    REAL MaxUpperAngle() const {
        return fAngles[fMaxUpper+1];
    }

    /// Return the angular threshold used by lower-family sector selection.
    REAL MinLowerAngle() const {
        return fAngles[fMinLower];
    }
    /// Return the local geometric mesh without transferring ownership.
    TPZGeoMesh *GetGeoMesh() const { return fGeoMesh.operator->(); }
    
    /// Return the local computational mesh without transferring ownership.
    TPZCompMesh *GetCompMesh() const { return fCompMesh; }

    /// Export selected mesh eigenvectors as VTK steps under input["directory"].
    /// Fields and sampling resolution come from the postprocessing settings.
    void PlotEnrichmentFunction(const std::string &filename, json &input);
    
    /// Fill the 2x2 counterclockwise rotation matrix for an angle in radians.
    static void ComputeRotationMatrix(REAL angle, TPZFNMatrix<4,REAL> &rotmat);
    /// Rotate the Cartesian components of a vector in place by the given
    /// counterclockwise angle in radians.
    static void RotateVector(REAL angle, std::array<REAL,2> &vec);
    /// Rotate a 2x2 tensor in place as R*T*R^T, with the angle in radians.
    /// Components are stored in row-major order: Txx, Txy, Tyx, Tyy.
    static void RotateTensor(REAL angle, std::array<REAL,4> &tensor);

    /// Number of selected real-valued enrichment modes, not the number of
    /// distinct eigenvalues; a complex pair contributes two real modes.
    const int getNEigenvalues() const { return fEigenvectorIndexes.size(); }

    /// Return eigenvalues indexed by the original selected enrichments,
    /// before any modal transform. A complex representative appears twice.
    TPZVec<CSTATE> SelectedEigenvalues() const;
    
private:

    /// Construct the circular, tip-centered geometric mesh from JSON sectors,
    /// including circular arcs and separate crack-face endpoints.
    TPZGeoMesh *CreateGeometry(json &input);


    /// Build the computational mesh using the configured SBFEM builder.
    /// Refine the skeleton and record the element group, sector angles
    /// and boundary interpolation elements needed for mode evaluation.
    TPZCompMesh *CreateSBFemSpace(TPZBuildSBFem &builder, json &input);

    /// Set the builder's skeleton order, refinement level, partition and
    /// material-ID translations from the JSON configuration.
    void ConfigureSBFemBuilder (TPZBuildSBFem &builder, json &input);
    /// Insert physical materials and auxiliary point, arc and skeleton
    /// materials into the local computational mesh according to JSON.
    void InsertMaterialObjectsH1(TPZCompMesh &cmesh, json &input);

    /// Compute the element-group stiffness, triggering construction of the
    /// SBFEM eigenvalues and boundary eigenvectors used for enrichments.
    void ComputeSBFemApproximationSpace();

    /// Select modes in the configured radial-exponent ranges (lambda = -s).
    /// Store eigenvector indices, real/imaginary tags and conjugate partners;
    /// this method updates selection metadata rather than returning values.
    void IdentifySingularModes(const json &input);

    /// Build two real modal combinations by minimizing reflection-parity error.
    /// They separate Mode I/II only when the local problem has that symmetry.
    /// Energy normalization is optional and disabled by default.
    void BuildSymmetryAdaptedModes(const json &input);

    /// Apply the two-mode symmetry construction independently to every
    /// requested eigenvalue-exponent range and assemble a block transform.
    void BuildBlockwiseSymmetryAdaptedModes(const json &input);

    /// Build an energy-orthonormal basis from exactly two real enrichments
    /// using a sampled elastic-energy Gram matrix, without requiring
    /// crack-plane reflection symmetry or asserting improved conditioning.
    void BuildEnergyOrthonormalModes(const json &input);

    /// Apply fModeTransform in place to mode values and their derivatives.
    /// Does nothing when disabled; rejects inconsistent mode dimensions.
    void ApplyModeTransform(std::vector<std::array<REAL,2>> &phi,
                            std::vector<std::array<REAL,4>> &dphi) const;

    /// Load an eigenvector already stored by fElementGroup as the current
    /// computational-mesh solution; modeindex indexes the full spectrum.
    void LoadEigenMode(int64_t modeindex);

    /// Export the current mesh solution and geometry using JSON output settings.
    /// rootname sets the filename prefix; refstep selects the VTK output step.
    void PlotSolution(const std::string &rootname, TPZCompMesh &cmesh, int refstep, json &input);

    /// Extend a requested crack-side field across the angular cut using
    /// face evaluation, a first-order angular correction and rotations.
    /// This helper is not called by the current main evaluation path.
    void ExtendEnrichmentFunction(MSectorFamily family,
                                 const std::vector<REAL> &x,
                                 std::vector<std::array<REAL, 2>> &phi,
                                 std::vector<std::array<REAL ,4>> &dphixy);

};

#endif // TPZGFEMENRICHMENT_H
