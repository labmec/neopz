#include "ExampleConfig.h"
#include <algorithm>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv)
{
    if(argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: sbfem_evaluate [config.json [x y]] [--verify-gradients]\n";
        return 0;
    }
    try {
        std::vector<std::string> args(argv+1,argv+argc);
        const bool verify = !args.empty() && args.back() == "--verify-gradients";
        if(verify) args.pop_back();
        if(args.size() != 0 && args.size() != 1 && args.size() != 3)
            throw std::invalid_argument("Usage: sbfem_evaluate [config.json [x y]] [--verify-gradients]");
        auto config = ReadConfiguration(args.empty() ? "isotropic.json" : args[0]);
        std::vector<REAL> point = {.35,.20};
        if(args.size() == 3) point = {ReadCoordinate(args[1]),ReadCoordinate(args[2])};
        const REAL r = std::hypot(point[0],point[1]);
        if(!(r > 0.) || r > config.at("radius").get<REAL>())
            throw std::invalid_argument("Choose a point inside the circle, excluding the singular tip (0,0)");

        TPZGFemEnrichment enrichment;
        enrichment.LoadFromJSON(config); // Generate once; reuse for many points.
        std::vector<std::array<REAL,2>> phi;
        std::vector<std::array<REAL,4>> dphi;
        enrichment.ComputeEnrichmentWithDerivatives(
            TPZGFemEnrichment::EAny,point,phi,dphi);
        if(phi.empty() || phi.size() != dphi.size())
            throw std::runtime_error("Empty or inconsistent mode arrays");
        std::cout << "\nPoint = (" << point[0] << ", " << point[1] << ")\n"
                  << "index  phi_x  phi_y  dphi_x/dx  dphi_x/dy  dphi_y/dx  dphi_y/dy\n"
                  << std::scientific << std::setprecision(10);
        for(size_t j=0; j<phi.size(); ++j) {
            std::cout << j;
            for(auto value : phi[j]) {
                if(!std::isfinite(value)) throw std::runtime_error("Non-finite mode value");
                std::cout << "  " << value;
            }
            for(auto value : dphi[j]) {
                if(!std::isfinite(value)) throw std::runtime_error("Non-finite derivative");
                std::cout << "  " << value;
            }
            std::cout << '\n';
        }
        if(verify) {
            // The supplied verification points are away from sector/material
            // boundaries. Do not use centered differences across a crack.
            const REAL h = std::min<REAL>(1.e-6*enrichment.Radius(),.01*r);
            REAL error2=0., norm2=0.;
            for(int axis=0; axis<2; ++axis) {
                auto plus=point, minus=point;
                plus[axis]+=h; minus[axis]-=h;
                std::vector<std::array<REAL,2>> vp,vm;
                std::vector<std::array<REAL,4>> unused;
                enrichment.ComputeEnrichmentWithDerivatives(TPZGFemEnrichment::EAny,plus,vp,unused);
                enrichment.ComputeEnrichmentWithDerivatives(TPZGFemEnrichment::EAny,minus,vm,unused);
                for(size_t j=0; j<phi.size(); ++j) for(int component=0; component<2; ++component) {
                    const REAL difference = (vp[j][component]-vm[j][component])/(2.*h);
                    const REAL exact = dphi[j][2*component+axis];
                    error2 += (difference-exact)*(difference-exact);
                    norm2 += exact*exact;
                }
            }
            const REAL relative = std::sqrt(error2/std::max<REAL>(norm2,1.e-30));
            std::cout << "Finite-difference relative gradient discrepancy = " << relative << '\n';
            if(!std::isfinite(relative) || relative > 2.e-6)
                throw std::runtime_error("Gradient verification failed");
        }
        std::cout << "sbfem_evaluate: OK\n";
        return 0;
    } catch(const std::exception &error) {
        std::cerr << "sbfem_evaluate: " << error.what() << '\n';
        return 1;
    }
}
