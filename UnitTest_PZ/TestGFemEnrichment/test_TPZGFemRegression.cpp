#include <catch2/catch_all.hpp>
#include "TPZGFemEnrichment.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {
using Values = std::vector<std::array<REAL,2>>;
using Gradients = std::vector<std::array<REAL,4>>;
using Transform = std::array<std::array<REAL,2>,2>;

json Config(bool bimaterial = false, REAL crackAngle = 0.)
{
    json lower = {{"matid",1}, {"tension state","plane stress"},
                  {"isotropic",{{"young",1.},{"poisson",0.3}}}};
    json config = {{"radius",1.}, {"skeleton porder",6}, {"nref_skeleton",1},
        {"angles",{crackAngle-M_PI,crackAngle-M_PI/2.,crackAngle,
                   crackAngle+M_PI/2.,crackAngle+M_PI}},
        {"sector material ids",{1,1,1,1}}, {"problem_type","Elastic"},
        {"materials",json::array({lower})}};
    if(bimaterial) {
        auto upper = lower;
        upper["matid"] = 2;
        upper["isotropic"]["young"] = 10.;
        config["materials"].push_back(upper);
        config["sector material ids"] = {1,1,2,2};
    }
    return config;
}

struct Sample {
    Values value;
    Gradients gradient;
};

Sample Evaluate(TPZGFemEnrichment &enrichment, const std::vector<REAL> &point,
                TPZGFemEnrichment::MSectorFamily side = TPZGFemEnrichment::EAny)
{
    Sample result;
    enrichment.ComputeEnrichmentWithDerivatives(side,point,result.value,result.gradient);
    REQUIRE(result.value.size() == result.gradient.size());
    REQUIRE(result.value.size() == static_cast<size_t>(enrichment.getNEigenvalues()));
    REQUIRE_FALSE(result.value.empty());
    for(const auto &mode : result.value)
        for(REAL value : mode) REQUIRE(std::isfinite(value));
    for(const auto &mode : result.gradient)
        for(REAL value : mode) REQUIRE(std::isfinite(value));
    return result;
}

REAL SquaredMagnitude(const Sample &sample)
{
    REAL result=0.;
    for(const auto &mode : sample.value) for(REAL entry : mode) result += entry*entry;
    for(const auto &mode : sample.gradient) for(REAL entry : mode) result += entry*entry;
    return result;
}

void CheckSameFace(const Sample &reference, const Sample &forced)
{
    REQUIRE(reference.value.size() == forced.value.size());
    REAL error=0.;
    for(size_t mode=0;mode<reference.value.size();mode++) {
        for(int component=0;component<2;component++)
            error += std::pow(reference.value[mode][component]-forced.value[mode][component],2);
        for(int component=0;component<4;component++)
            error += std::pow(reference.gradient[mode][component]-forced.gradient[mode][component],2);
    }
    const REAL scale=SquaredMagnitude(reference);
    REQUIRE(scale > 0.);
    REQUIRE(std::sqrt(error/scale) < 1.e-7);
}

// Compare gradients with centered differences without crossing sectors or faces.
void CheckDerivatives(TPZGFemEnrichment &enrichment, const std::vector<REAL> &point)
{
    constexpr REAL step = 1.e-6;
    const auto center = Evaluate(enrichment,point);
    REAL squaredError = 0., squaredScale = 0.;
    for(int direction=0; direction<2; direction++) {
        auto plusPoint=point, minusPoint=point;
        plusPoint[direction] += step;
        minusPoint[direction] -= step;
        const auto plus=Evaluate(enrichment,plusPoint);
        const auto minus=Evaluate(enrichment,minusPoint);
        for(size_t mode=0; mode<center.value.size(); mode++) {
            for(int component=0; component<2; component++) {
                const REAL difference=(plus.value[mode][component]-minus.value[mode][component])/(2.*step);
                const REAL analytic=center.gradient[mode][2*component+direction];
                squaredError += (difference-analytic)*(difference-analytic);
                squaredScale += analytic*analytic;
            }
        }
    }
    REQUIRE(squaredScale > 0.);
    INFO("Relative finite-difference gradient error: " << std::sqrt(squaredError/squaredScale));
    REQUIRE(std::sqrt(squaredError/squaredScale) < 2.e-6);
}

// Fit one constant two-mode transform using values only, not gradients.
Transform FitTransform(TPZGFemEnrichment &original, TPZGFemEnrichment &changed)
{
    REAL gram[2][2] = {{0.,0.},{0.,0.}};
    REAL rhs[2][2] = {{0.,0.},{0.,0.}};
    for(const auto &point : std::vector<std::vector<REAL>>{{.31,.22},{-.32,.21},{.18,-.37},{-.43,-.18}}) {
        const auto a=Evaluate(original,point), b=Evaluate(changed,point);
        REQUIRE(a.value.size() == 2);
        REQUIRE(b.value.size() == 2);
        for(int component=0; component<2; component++)
            for(int i=0; i<2; i++) for(int j=0; j<2; j++) {
                gram[i][j] += a.value[i][component]*a.value[j][component];
                rhs[i][j] += a.value[i][component]*b.value[j][component];
            }
    }
    const REAL determinant=gram[0][0]*gram[1][1]-gram[0][1]*gram[1][0];
    REQUIRE(determinant > 1.e-10*gram[0][0]*gram[1][1]);
    Transform transform{};
    for(int output=0; output<2; output++) {
        transform[output][0]=(gram[1][1]*rhs[0][output]-gram[0][1]*rhs[1][output])/determinant;
        transform[output][1]=(gram[0][0]*rhs[1][output]-gram[1][0]*rhs[0][output])/determinant;
    }
    const REAL transformDet=transform[0][0]*transform[1][1]-transform[0][1]*transform[1][0];
    REAL transformScale=0.;
    for(const auto &row : transform) for(REAL entry : row) transformScale += entry*entry;
    REQUIRE(std::abs(transformDet) > 1.e-10*transformScale);
    return transform;
}

// Validate the fitted transform on independent points and all four derivatives.
void CheckSameModalSpace(TPZGFemEnrichment &original, TPZGFemEnrichment &changed)
{
    const auto transform=FitTransform(original,changed);
    REAL valueError=0., valueScale=0., gradientError=0., gradientScale=0.;
    for(const auto &point : std::vector<std::vector<REAL>>{{.27,.41},{-.51,.13},{.23,-.29},{-.24,-.46}}) {
        const auto a=Evaluate(original,point), b=Evaluate(changed,point);
        for(int output=0; output<2; output++) {
            for(int component=0; component<2; component++) {
                const REAL prediction=transform[output][0]*a.value[0][component]+transform[output][1]*a.value[1][component];
                valueError += std::pow(prediction-b.value[output][component],2);
                valueScale += std::pow(b.value[output][component],2);
            }
            for(int derivative=0; derivative<4; derivative++) {
                const REAL prediction=transform[output][0]*a.gradient[0][derivative]+transform[output][1]*a.gradient[1][derivative];
                gradientError += std::pow(prediction-b.gradient[output][derivative],2);
                gradientScale += std::pow(b.gradient[output][derivative],2);
            }
        }
    }
    REQUIRE(valueScale > 0.);
    REQUIRE(gradientScale > 0.);
    REQUIRE(std::sqrt(valueError/valueScale) < 1.e-8);
    REQUIRE(std::sqrt(gradientError/gradientScale) < 1.e-8);
}

json OrthotropicConfig(REAL angle)
{
    auto config=Config();
    config["materials"][0].erase("isotropic");
    config["materials"][0]["orthotropic"]={{"Ex",4.},{"Ey",1.},{"nuxy",.25},
        {"Gxy",.4},{"direction",{std::cos(angle),std::sin(angle)}}};
    return config;
}
}

TEST_CASE("Complex reconstructions span the same real fields and gradients",
          "[gfem][regression][complex]")
{
    auto representativeConfig=Config(true), conjugateConfig=representativeConfig;
    representativeConfig["complex mode reconstruction"]="representative_parts";
    conjugateConfig["complex mode reconstruction"]="conjugate_sum";
    TPZGFemEnrichment representative,conjugate;
    representative.LoadFromJSON(representativeConfig);
    conjugate.LoadFromJSON(conjugateConfig);
    REQUIRE(representative.getNEigenvalues() == 2);
    REQUIRE(conjugate.getNEigenvalues() == 2);
    REQUIRE(std::abs(representative.SelectedEigenvalues()[0].imag()) > 1.e-4);
    CheckSameModalSpace(representative,conjugate);
    CheckDerivatives(representative,{-.32,.21});
    CheckDerivatives(conjugate,{-.43,-.18});
}

TEST_CASE("Optional modal transformations preserve fields and derivative consistency",
          "[gfem][regression][transform]")
{
    for(int variant=0; variant<3; variant++) {
        CAPTURE(variant);
        auto originalConfig=Config(variant==2), transformedConfig=originalConfig;
        if(variant==2) transformedConfig["energy orthonormal modes"]=true;
        else {
            transformedConfig["symmetry adapted modes"]=true;
            transformedConfig["energy normalize adapted modes"]=(variant==1);
        }
        TPZGFemEnrichment original,transformed;
        original.LoadFromJSON(originalConfig);
        transformed.LoadFromJSON(transformedConfig);
        CheckSameModalSpace(original,transformed);
        CheckDerivatives(transformed,{.31,.22});
        CheckDerivatives(transformed,{-.43,-.18});
    }
}

TEST_CASE("Oblique orthotropic modes have consistent Cartesian derivatives",
          "[gfem][regression][orthotropic][derivatives]")
{
    for(REAL angle : {M_PI/6.,M_PI/4.}) {
        CAPTURE(angle);
        auto config=OrthotropicConfig(angle);
        TPZGFemEnrichment enrichment;
        enrichment.LoadFromJSON(config);
        REQUIRE(enrichment.getNEigenvalues() == 2);
        CheckDerivatives(enrichment,{.31,.22});
        CheckDerivatives(enrichment,{.18,-.37});
    }
}

TEST_CASE("Circular geometry and crack-face selection support rotated angular cuts",
          "[gfem][regression][geometry]")
{
    for(REAL angle : {0.,M_PI/6.}) {
        CAPTURE(angle);
        auto config=Config(false,angle);
        TPZGFemEnrichment enrichment;
        enrichment.LoadFromJSON(config);
        REQUIRE(enrichment.CrackAxisAngle() == Catch::Approx(angle).margin(1.e-12));
        auto *mesh=enrichment.GetGeoMesh();
        REQUIRE(mesh != nullptr);
        int crackFaceEndpoints=0;
        for(int64_t node=0; node<mesh->NodeVec().NElements(); node++) {
            const auto &coordinate=mesh->NodeVec()[node];
            const REAL radius=std::hypot(coordinate.Coord(0),coordinate.Coord(1));
            REQUIRE((radius < 1.e-12 || std::abs(radius-1.) < 1.e-10));
            if(std::hypot(coordinate.Coord(0)+std::cos(angle),coordinate.Coord(1)+std::sin(angle)) < 1.e-10)
                crackFaceEndpoints++;
        }
        // The two faces meet geometrically but retain distinct endpoint nodes.
        REQUIRE(crackFaceEndpoints == 2);
        const REAL firstAngle=angle-M_PI+1.e-10;
        const REAL lastAngle=angle+M_PI-1.e-10;
        const std::vector<REAL> firstPoint={.43*std::cos(firstAngle),.43*std::sin(firstAngle)};
        const std::vector<REAL> lastPoint={.43*std::cos(lastAngle),.43*std::sin(lastAngle)};
        // Force the opposite geometric side and verify selection of the
        // requested angular endpoint, including its Cartesian derivatives.
        CheckSameFace(Evaluate(enrichment,firstPoint),
                      Evaluate(enrichment,lastPoint,TPZGFemEnrichment::EUpper));
        CheckSameFace(Evaluate(enrichment,lastPoint),
                      Evaluate(enrichment,firstPoint,TPZGFemEnrichment::ELower));
        for(REAL offset : {-1.e-10,0.,1.e-10}) {
            const REAL theta=angle+M_PI+offset;
            const std::vector<REAL> point={.43*std::cos(theta),.43*std::sin(theta)};
            Evaluate(enrichment,point);
            Evaluate(enrichment,point,TPZGFemEnrichment::EUpper);
            Evaluate(enrichment,point,TPZGFemEnrichment::ELower);
        }
        if(angle==0.) {
            Evaluate(enrichment,{-.43,0.});
            Evaluate(enrichment,{-.43,-0.});
        }
        CheckDerivatives(enrichment,{.31,.22});
    }
}

TEST_CASE("Exponent ranges select both independent modes of each Williams order",
          "[gfem][regression][selection]")
{
    for(int variant=0; variant<4; variant++) {
        CAPTURE(variant);
        auto config=Config();
        config["skeleton porder"]=8;
        config["eigenvalue exponent ranges"]=variant==0 ? json::array({{.4,.6}})
            : variant==1 ? json::array({{1.4,1.6}}) : json::array({{.4,.6},{1.4,1.6}});
        if(variant==3) config["blockwise symmetry adapted modes"]=true;
        TPZGFemEnrichment enrichment;
        enrichment.LoadFromJSON(config);
        const auto eigenvalues=enrichment.SelectedEigenvalues();
        REQUIRE(eigenvalues.size() == (variant<2 ? 2 : 4));
        int leading=0,next=0;
        for(const auto &eigenvalue : eigenvalues) {
            REQUIRE(std::abs(eigenvalue.imag()) < 1.e-5);
            if(std::abs(eigenvalue.real()+.5) < 1.e-5) leading++;
            else if(std::abs(eigenvalue.real()+1.5) < 1.e-4) next++;
            else FAIL("Unexpected selected Williams exponent: " << eigenvalue);
        }
        REQUIRE(leading == (variant==1 ? 0 : 2));
        REQUIRE(next == (variant==0 ? 0 : 2));
        CheckDerivatives(enrichment,{.31,.22});
    }
}

TEST_CASE("Omitted reconstruction and normalization options match explicit defaults",
          "[gfem][regression][defaults]")
{
    auto implicitConfig=Config(true), explicitConfig=implicitConfig;
    explicitConfig["complex mode reconstruction"]="representative_parts";
    explicitConfig["symmetry adapted modes"]=false;
    explicitConfig["blockwise symmetry adapted modes"]=false;
    explicitConfig["energy orthonormal modes"]=false;
    explicitConfig["energy normalize adapted modes"]=false;
    TPZGFemEnrichment implicit,explicitModes;
    implicit.LoadFromJSON(implicitConfig);
    explicitModes.LoadFromJSON(explicitConfig);
    CheckSameModalSpace(implicit,explicitModes);
    // A span check alone would miss a changed default scale, such as the
    // factor introduced by conjugate sums. Check sampled magnitudes too.
    for(const auto &point : std::vector<std::vector<REAL>>{{.31,.22},{-.43,-.18}}) {
        const REAL implicitMagnitude=SquaredMagnitude(Evaluate(implicit,point));
        const REAL explicitMagnitude=SquaredMagnitude(Evaluate(explicitModes,point));
        REQUIRE(implicitMagnitude == Catch::Approx(explicitMagnitude).epsilon(1.e-9));
    }
    auto adaptedImplicit=Config(), adaptedExplicit=adaptedImplicit;
    adaptedImplicit["symmetry adapted modes"]=true;
    adaptedExplicit["symmetry adapted modes"]=true;
    adaptedExplicit["energy normalize adapted modes"]=false;
    TPZGFemEnrichment a,b;
    a.LoadFromJSON(adaptedImplicit);
    b.LoadFromJSON(adaptedExplicit);
    // Unlike a span test, compare sampled magnitudes to detect normalization.
    REAL normA=0.,normB=0.;
    for(const auto &point : std::vector<std::vector<REAL>>{{.31,.22},{-.43,-.18}}) {
        const auto sa=Evaluate(a,point),sb=Evaluate(b,point);
        for(size_t mode=0; mode<sa.value.size(); mode++) for(int k=0;k<2;k++) {
            normA += sa.value[mode][k]*sa.value[mode][k];
            normB += sb.value[mode][k]*sb.value[mode][k];
        }
    }
    REQUIRE(normA == Catch::Approx(normB).epsilon(1.e-9));
}

TEST_CASE("Malformed optional configurations fail with explicit exceptions",
          "[gfem][regression][validation]")
{
    std::vector<json> invalid;
    auto config=Config();config["complex mode reconstruction"]="unknown";invalid.push_back(config);
    config=Config();config["eigenvalue exponent ranges"]=json::array();invalid.push_back(config);
    config=Config();config["eigenvalue exponent ranges"]={{.1}};invalid.push_back(config);
    config=Config();config["eigenvalue exponent ranges"]={{-.1,.6}};invalid.push_back(config);
    config=Config();config["eigenvalue exponent ranges"]={{.6,.4}};invalid.push_back(config);
    config=OrthotropicConfig(0.);config["materials"][0]["orthotropic"]["direction"]={0.,0.};invalid.push_back(config);
    config=OrthotropicConfig(0.);config["materials"][0]["orthotropic"]["direction"]={1.};invalid.push_back(config);
    config=Config();config["symmetry adapted modes"]=true;config["energy orthonormal modes"]=true;invalid.push_back(config);
    config=Config();config["symmetry adapted modes"]=true;config["blockwise symmetry adapted modes"]=true;invalid.push_back(config);
    config=Config();config["blockwise symmetry adapted modes"]=true;invalid.push_back(config);
    for(size_t index=0;index<invalid.size();index++) {
        CAPTURE(index);
        TPZGFemEnrichment enrichment;
        REQUIRE_THROWS_AS(enrichment.LoadFromJSON(invalid[index]),std::invalid_argument);
    }
}
