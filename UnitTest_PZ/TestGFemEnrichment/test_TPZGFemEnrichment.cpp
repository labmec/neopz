#include <catch2/catch_all.hpp>

#include "TPZGFemEnrichment.h"
#include <nlohmann/json.hpp>
using nlohmann::json;
#include <random>
#include <cmath>
#include <vector>
#include <array>
#include <iostream>


// Basic test scaffold for TPZGFemEnrichment.
// The actual tests will be added one-by-one.

TEST_CASE("TPZGFemEnrichment - basic construction and API", "[gfem][enrichment]") {
        // Minimal CrackDef-like JSON to be used as input for construction tests
        const char *json_text = R"({
            "radius": 1.0,
            "angles": [0.0, 3.141592653589793, 6.283185307179586],
            "sector material ids": [1, 1],
            "problem_type" : "Elastic",

            "materials": [
                {
                    "matid": 1,
                    "tension state": "plane stress",
                    "isotropic": { "young": 210000.0, "poisson": 0.3 }
                }
            ]
        })";
    
        std::cout << "Parsing JSON for TPZGFemEnrichment test..." << std::endl;
        json j = json::parse(json_text);
            // Instantiate TPZGFemEnrichment object
        TPZGFemEnrichment enrichment;
        enrichment.LoadFromJSON(j);

    SECTION("test two ways of computing a solution") {
        REQUIRE(enrichment.Radius() == Catch::Approx(1.0));
        // TODO: create minimal objects needed to construct TPZGFemEnrichment
        // e.g., call LoadFromJSON(j) when mesh/NeoPZ context is available.

        // generate 10 random points uniformly inside circle of radius enrichment.Radius()
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> ang_dist(0.0, 2.0 * M_PI);
        std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

        const double R = enrichment.Radius();
        std::vector<std::array<double,2>> points;
        points.reserve(10);
        for(int i = 0; i < 10; ++i) {
            double theta = ang_dist(gen);
            double rho = std::sqrt(unit_dist(gen)) * R; // sqrt for uniform sampling in disk
            double x = rho * std::cos(theta);
            double y = rho * std::sin(theta);
            points.push_back({x, y});
            REQUIRE(std::hypot(x, y) <= R + 1e-12);
        }

        // optional: print generated points for debug
        for(size_t i = 0; i < points.size(); ++i) {
            std::cout << "point " << i << ": (" << points[i][0] << ", " << points[i][1] << ")\n";
        }

        // For each generated point compare the three enrichment APIs
        for(size_t ip=0; ip < points.size(); ++ip) {
            std::vector<REAL> xp = {points[ip][0], points[ip][1]};
            TPZGFemEnrichment::MSectorFamily family = TPZGFemEnrichment::ENone;

            std::vector<std::array<REAL,2>> phiA, phiB, phiC;
            std::vector<std::array<REAL,4>> dA, dB, dC;

            // 1) Direct call
            enrichment.ComputeEnrichmentWithDerivatives(family, xp, phiA, dA);

            // 2) Via GetEnrichmentFunction()
            auto f = enrichment.GetEnrichmentFunction();
            f(family, xp, phiB, dB);

            // 3) Direct EnrichmentFunction
            enrichment.EnrichmentFunction(family, xp, phiC, dC);

            REQUIRE(phiA.size() == phiB.size());
            REQUIRE(phiA.size() == phiC.size());
            REQUIRE(dA.size() == dB.size());
            REQUIRE(dA.size() == dC.size());
            REQUIRE(phiA.size() > 0);

            for(size_t im=0; im<phiA.size(); ++im) {
                REQUIRE(phiA[im][0] == Catch::Approx(phiB[im][0]).margin(1e-8));
                REQUIRE(phiA[im][1] == Catch::Approx(phiB[im][1]).margin(1e-8));
                REQUIRE(phiA[im][0] == Catch::Approx(phiC[im][0]).margin(1e-8));
                REQUIRE(phiA[im][1] == Catch::Approx(phiC[im][1]).margin(1e-8));
                for(int id=0; id<4; ++id) {
                    REQUIRE(dA[im][id] == Catch::Approx(dB[im][id]).margin(1e-8));
                    REQUIRE(dA[im][id] == Catch::Approx(dC[im][id]).margin(1e-8));
                }
            }
        }
    }

    SECTION("verify the derivatives") {
        std::array<REAL,2> point = {0.25, 0.5};
        double del=0.1;
        for(int pindex = 0; pindex < 18; ++pindex) {
            double pangle = -M_PI/2.0 + pindex * (M_PI/6.0);
            std::array<REAL,2> point = {0.5*std::cos(pangle), 0.5*std::sin(pangle)};
            std::cout << "Testing pindex " << pindex << " pangle " << pangle << " point: (" << point[0] << ", " << point[1] << ")"  << std::endl;
            // compute base value and derivatives at point
            std::vector<REAL> xp = {point[0], point[1]};
            TPZGFemEnrichment::MSectorFamily family = TPZGFemEnrichment::ENone;
            if(pangle < M_PI/2.0) {
                family = TPZGFemEnrichment::EUpper;
                std::cout << "Upper sector: " << pangle << std::endl;
            } else if (pangle >= 1.5*M_PI) {
                family = TPZGFemEnrichment::ELower;
                std::cout << "Lower sector: " << pangle << std::endl;
            } else {
                std::cout << "No sector: " << pangle << std::endl;
            }
            
            std::vector<std::array<REAL,2>> phi0;
            std::vector<std::array<REAL,4>> d0;
            enrichment.ComputeEnrichmentWithDerivatives(family, xp, phi0, d0);
            REQUIRE(phi0.size() > 0);
            REQUIRE(phi0.size() == d0.size());

            // directions: 10 equally spaced between 0 and 2*pi
            const int ndirs = 10;
            std::vector<std::array<REAL,2>> dirs;
            dirs.reserve(ndirs);
            for(int i=0;i<ndirs;++i){
                double theta = (2.0*M_PI*double(i))/double(ndirs);
                dirs.push_back({REAL(std::cos(theta)), REAL(std::sin(theta))});
            }
            // alphas to test (halving sequence -> expect quadratic scaling)
            std::vector<double> alphas = {1.0, 0.5, 0.25, 0.125};

            // for each direction, compute residual norms for each alpha and verify scaling ~ alpha^2
            for(int id=0; id<ndirs; ++id) {
                auto dir = dirs[id];
                std::vector<double> norms;
                norms.reserve(alphas.size());
                for(double alpha : alphas) {
                    // displacement dx = alpha * del * dir
                    const double sx = alpha * del * double(dir[0]);
                    const double sy = alpha * del * double(dir[1]);
                    std::vector<REAL> xq = { REAL(point[0] + sx), REAL(point[1] + sy) };

                    std::vector<std::array<REAL,2>> phi_q;
                    std::vector<std::array<REAL,4>> d_q;
                    enrichment.ComputeEnrichmentWithDerivatives(family, xq, phi_q, d_q);
                    REQUIRE(phi_q.size() == phi0.size());
                    REQUIRE(d_q.size() == d0.size());

                    // compute residual r = phi(x+dx) - phi(x) - J(x)*dx
                    double accum = 0.0;
                    for(size_t im=0; im<phi0.size(); ++im) {
                        // J*dx using derivatives in d0[im]: [du/dx, du/dy, dv/dx, dv/dy]
                        double j0 = double(d0[im][0]) * sx + double(d0[im][1]) * sy;
                        double j1 = double(d0[im][2]) * sx + double(d0[im][3]) * sy;
                        std::cout << "im " << im << " phi_q: (" << double(phi_q[im][0]) << ", " << double(phi_q[im][1]) << "), "
                                  << "phi_0: (" << double(phi0[im][0]) << ", " << double(phi0[im][1]) << "), "
                                  << "J*dx: (" << j0 << ", " << j1 << ")\n";
                        double r0 = double(phi_q[im][0]) - double(phi0[im][0]) - j0;
                        double r1 = double(phi_q[im][1]) - double(phi0[im][1]) - j1;
                        accum += r0*r0 + r1*r1;
                    }
                    norms.push_back(std::sqrt(accum));
                }
                // print norms for debug
                std::cout << "pangle " << pangle << " Direction " << id << " norms: ";
                for(auto n : norms) {
                    std::cout << n << " ";
                }
                std::cout << std::endl;

                // verify quadratic scaling: for halving alpha residual should reduce ~4x
                for(size_t k=0; k+1<norms.size(); ++k) {
                    double n1 = norms[k];
                    double n2 = norms[k+1];
                    // if residuals are extremely small, skip strict ratio check
                    if(n2 < 1e-14 && n1 < 1e-14) {
                        continue;
                    }
                    // avoid division by zero
                    REQUIRE(n2 > 0.0);
                    double ratio = n1 / n2;
                    // print ratio for debug
                    std::cout << " Direction " << id << " ratio[" << k << "] = " << ratio << std::endl;
                    // expect ratio around 4 for quadratic behaviour; allow loose bounds
                    REQUIRE(ratio > 2.0);
                    REQUIRE(ratio < 8.0);
                }
            }
        }
    }
}


TEST_CASE("TPZGFemEnrichment - derivatives and rotation", "[gfem][enrichment][derivatives]") {
    SECTION("compute rotation matrix") {
        // Minimal CrackDef-like JSON used as input for rotation helper tests
        const char *json_text = R"({
            "radius": 1.0,
            "angles": [0.0, 3.141592653589793, 6.283185307179586],
            "sector material ids": [1, 2],
            "materials": [
                {
                    "matid": 1,
                    "tension state": "plane stress",
                    "isotropic": { "young": 210000.0, "poisson": 0.3 }
                }
            ],
            "boundary_conditions": [
                { "matid": 2, "name": "fixed", "value": [0.0, 0.0], "type": 0 }
            ]
        })";

        json j = json::parse(json_text);
        // Verify rotated derivatives behave consistently with rotated displacements
        // Check that R*(J*dx) == J'*(R*dx) where J' = R*J*R^T
        const std::vector<double> angles = {0.0, 0.3, 0.7, 1.2, -0.5};
        std::mt19937 rng(12345);
        std::uniform_real_distribution<double> val(-10.0, 10.0);
        std::uniform_real_distribution<double> dxy(-1.0, 1.0);

        for(double ang : angles) {
            for(int trial=0; trial<50; ++trial) {
                // random tensor J in form [t0,t1,t2,t3]
                std::array<REAL,4> J = {REAL(val(rng)), REAL(val(rng)), REAL(val(rng)), REAL(val(rng))};
                // random displacement dx
                std::array<REAL,2> dx = {REAL(dxy(rng)), REAL(dxy(rng))};

                // compute J*dx
                double j0 = double(J[0])*double(dx[0]) + double(J[1])*double(dx[1]);
                double j1 = double(J[2])*double(dx[0]) + double(J[3])*double(dx[1]);
                std::array<REAL,2> jvec = {REAL(j0), REAL(j1)};

                // rotate jvec: R * (J*dx)
                auto jrot = jvec;
                TPZGFemEnrichment::RotateVector(ang, jrot);

                // rotate dx and J to compute J'*(R*dx)
                auto dxrot = dx;
                TPZGFemEnrichment::RotateVector(ang, dxrot);
                auto Jprime = J;
                TPZGFemEnrichment::RotateTensor(ang, Jprime);

                // compute J' * dxrot
                double k0 = double(Jprime[0])*double(dxrot[0]) + double(Jprime[1])*double(dxrot[1]);
                double k1 = double(Jprime[2])*double(dxrot[0]) + double(Jprime[3])*double(dxrot[1]);

                // compare
                REQUIRE(double(jrot[0]) == Catch::Approx(k0).margin(1e-10));
                REQUIRE(double(jrot[1]) == Catch::Approx(k1).margin(1e-10));
            }
        }
    }
}
