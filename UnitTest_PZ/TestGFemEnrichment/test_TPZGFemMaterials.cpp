#include <catch2/catch_all.hpp>

#include "TPZGFemEnrichment.h"

#include <algorithm>
#include <cmath>
#include <nlohmann/json.hpp>

using nlohmann::json;

namespace {

json HomogeneousConfig(const json &material)
{
    return {
        {"radius", 1.0},
        {"skeleton porder", 6},
        {"nref_skeleton", 1},
        {"angles", {-M_PI, -M_PI/2., 0., M_PI/2., M_PI}},
        {"sector material ids", {1, 1, 1, 1}},
        {"problem_type", "Elastic"},
        {"materials", {material}}
    };
}

void RequireFiniteEnrichment(TPZGFemEnrichment &enrichment,
                             const std::vector<REAL> &point)
{
    std::vector<std::array<REAL,2>> phi;
    std::vector<std::array<REAL,4>> dphi;
    enrichment.ComputeEnrichmentWithDerivatives(
        TPZGFemEnrichment::EAny, point, phi, dphi);

    REQUIRE(phi.size() == dphi.size());
    REQUIRE(phi.size() == enrichment.getNEigenvalues());
    REQUIRE_FALSE(phi.empty());
    for(const auto &mode : phi) {
        REQUIRE(std::isfinite(mode[0]));
        REQUIRE(std::isfinite(mode[1]));
    }
    for(const auto &gradient : dphi) {
        for(const auto value : gradient) {
            REQUIRE(std::isfinite(value));
        }
    }
}

}

TEST_CASE("Orthotropic material recovers the isotropic singular exponents",
          "[gfem][sbfem][orthotropic]")
{
    constexpr REAL young = 2.0;
    constexpr REAL poisson = 0.25;
    constexpr REAL shear = young/(2.*(1.+poisson));

    json isotropicMaterial = {
        {"matid", 1},
        {"tension state", "plane stress"},
        {"isotropic", {{"young", young}, {"poisson", poisson}}}
    };
    json orthotropicMaterial = {
        {"matid", 1},
        {"tension state", "plane stress"},
        {"orthotropic", {
            {"Ex", young}, {"Ey", young}, {"nuxy", poisson},
            {"Gxy", shear}, {"direction", {1.0, 0.0}}
        }}
    };

    auto isotropicConfig = HomogeneousConfig(isotropicMaterial);
    auto orthotropicConfig = HomogeneousConfig(orthotropicMaterial);
    TPZGFemEnrichment isotropic;
    TPZGFemEnrichment orthotropic;
    isotropic.LoadFromJSON(isotropicConfig);
    orthotropic.LoadFromJSON(orthotropicConfig);

    const auto isotropicEigenvalues = isotropic.SelectedEigenvalues();
    const auto orthotropicEigenvalues = orthotropic.SelectedEigenvalues();
    REQUIRE(isotropicEigenvalues.size() == orthotropicEigenvalues.size());
    REQUIRE(isotropicEigenvalues.size() == 2);
    for(int i = 0; i < isotropicEigenvalues.size(); i++) {
        REQUIRE(orthotropicEigenvalues[i].real() ==
                Catch::Approx(isotropicEigenvalues[i].real()).margin(1.e-7));
        REQUIRE(orthotropicEigenvalues[i].imag() ==
                Catch::Approx(isotropicEigenvalues[i].imag()).margin(1.e-7));
    }
    RequireFiniteEnrichment(orthotropic, {0.35, 0.2});
}

TEST_CASE("Rotated orthotropic material produces finite SBFEM enrichments",
          "[gfem][sbfem][orthotropic]")
{
    json material = {
        {"matid", 1},
        {"tension state", "plane stress"},
        {"orthotropic", {
            {"Ex", 4.0}, {"Ey", 1.0}, {"nuxy", 0.2},
            {"Gxy", 0.5}, {"direction", {0.0, 1.0}}
        }}
    };
    auto config = HomogeneousConfig(material);
    TPZGFemEnrichment enrichment;
    enrichment.LoadFromJSON(config);

    REQUIRE(enrichment.getNEigenvalues() == 2);
    RequireFiniteEnrichment(enrichment, {0.35, 0.2});
    RequireFiniteEnrichment(enrichment, {-0.35, 0.0});
}

TEST_CASE("Repeated aligned orthotropic modes are separated by crack-plane symmetry",
          "[gfem][sbfem][orthotropic][symmetry]")
{
    json material = {
        {"matid", 1},
        {"tension state", "plane stress"},
        {"orthotropic", {
            {"Ex", 4.0}, {"Ey", 1.0}, {"nuxy", 0.25},
            {"Gxy", 0.4}, {"direction", {1.0, 0.0}}
        }}
    };
    auto config = HomogeneousConfig(material);
    config["symmetry adapted modes"] = true;
    TPZGFemEnrichment enrichment;
    enrichment.LoadFromJSON(config);

    const REAL theta = 0.73;
    const REAL radius = 0.47;
    std::vector<std::array<REAL,2>> plus, minus;
    std::vector<std::array<REAL,4>> dplus, dminus;
    enrichment.ComputeEnrichmentWithDerivatives(
        TPZGFemEnrichment::EAny,
        {radius*std::cos(theta),radius*std::sin(theta)}, plus, dplus);
    enrichment.ComputeEnrichmentWithDerivatives(
        TPZGFemEnrichment::EAny,
        {radius*std::cos(theta),-radius*std::sin(theta)}, minus, dminus);

    REQUIRE(plus.size() == 2);
    const REAL scaleI = std::max({std::abs(plus[0][0]),std::abs(plus[0][1]),1.e-12});
    const REAL scaleII = std::max({std::abs(plus[1][0]),std::abs(plus[1][1]),1.e-12});
    // Mode I: u_x even and u_y odd across the crack plane.
    REQUIRE(std::abs(plus[0][0]-minus[0][0])/scaleI < 1.e-6);
    REQUIRE(std::abs(plus[0][1]+minus[0][1])/scaleI < 1.e-6);
    // Mode II: u_x odd and u_y even across the crack plane.
    REQUIRE(std::abs(plus[1][0]+minus[1][0])/scaleII < 1.e-6);
    REQUIRE(std::abs(plus[1][1]-minus[1][1])/scaleII < 1.e-6);
}

TEST_CASE("Bimaterial oscillatory mode is represented by two real enrichments",
          "[gfem][sbfem][bimaterial]")
{
    json config = {
        {"radius", 1.0},
        {"skeleton porder", 8},
        {"nref_skeleton", 2},
        {"angles", {-M_PI, 0., M_PI}},
        {"sector material ids", {1, 2}},
        {"problem_type", "Elastic"},
        {"materials", {
            {
                {"matid", 1}, {"tension state", "plane strain"},
                {"isotropic", {{"young", 1.0}, {"poisson", 0.3}}}
            },
            {
                {"matid", 2}, {"tension state", "plane strain"},
                {"isotropic", {{"young", 10.0}, {"poisson", 0.3}}}
            }
        }}
    };

    TPZGFemEnrichment enrichment;
    enrichment.LoadFromJSON(config);
    const auto eigenvalues = enrichment.SelectedEigenvalues();

    REQUIRE(eigenvalues.size() == 2);
    REQUIRE(eigenvalues[0].real() == Catch::Approx(-0.5).margin(1.e-6));
    REQUIRE(std::abs(eigenvalues[0].imag()) > 1.e-4);
    REQUIRE(eigenvalues[1].real() == Catch::Approx(eigenvalues[0].real()).margin(1.e-12));
    REQUIRE(eigenvalues[1].imag() == Catch::Approx(eigenvalues[0].imag()).margin(1.e-12));
    RequireFiniteEnrichment(enrichment, {0.5, 0.5});
}
