#include "ExampleConfig.h"
#include <iomanip>
#include <iostream>

int main(int argc, char **argv)
{
    if((argc == 2 && std::string(argv[1]) == "--help") || argc > 3) {
        std::cout << "Usage: sbfem_eigenvalues [config.json] [--verify-isotropic]\n";
        return argc > 3 ? 2 : 0;
    }
    try {
        if(argc == 3 && std::string(argv[2]) != "--verify-isotropic")
            throw std::invalid_argument("Unknown option: " + std::string(argv[2]));
        auto config = ReadConfiguration(argc > 1 ? argv[1] : "isotropic.json");
        TPZGFemEnrichment enrichment;
        enrichment.LoadFromJSON(config);
        const auto eigenvalues = enrichment.SelectedEigenvalues();
        if(eigenvalues.size() == 0) throw std::runtime_error("No modes selected");
        std::cout << "\nSelected crack-tip modes (zero-based index)\n"
                  << "index       Re(s)             Im(s)             Re(lambda)        Im(lambda)\n"
                  << std::scientific << std::setprecision(10);
        for(int j=0; j<eigenvalues.size(); ++j) {
            const auto s = eigenvalues[j];
            const auto lambda = -s; // Evaluator uses xi^lambda, with xi = r/R.
            if(!std::isfinite(s.real()) || !std::isfinite(s.imag()))
                throw std::runtime_error("Non-finite eigenvalue");
            std::cout << j << "  " << s.real() << "  " << s.imag()
                      << "  " << lambda.real() << "  " << lambda.imag() << '\n';
            if(argc == 3 && (std::abs(s.real()+.5) > 1.e-4 || std::abs(s.imag()) > 1.e-4))
                throw std::runtime_error("Expected a leading isotropic eigenvalue near -0.5");
        }
        if(argc == 3 && eigenvalues.size() != 2)
            throw std::runtime_error("Expected two leading isotropic modes");
        std::cout << "sbfem_eigenvalues: OK\n";
        return 0;
    } catch(const std::exception &error) {
        std::cerr << "sbfem_eigenvalues: " << error.what() << '\n';
        return 1;
    }
}
