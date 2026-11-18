#ifndef SBFEM_ENRICHMENT_EXAMPLE_CONFIG_H
#define SBFEM_ENRICHMENT_EXAMPLE_CONFIG_H

#include "TPZGFemEnrichment.h"
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

// Basic guards for these examples; the generator does not validate all JSON.
inline nlohmann::json ReadConfiguration(const std::string &filename)
{
    std::ifstream input(filename);
    if(!input) throw std::runtime_error("Cannot open configuration: " + filename);
    nlohmann::json config;
    input >> config;
    if(config.at("problem_type") != "Elastic")
        throw std::invalid_argument("These examples support problem_type Elastic only");
    const double radius = config.at("radius").get<double>();
    if(!std::isfinite(radius) || radius <= 0.)
        throw std::invalid_argument("radius must be finite and positive");
    const auto angles = config.at("angles").get<std::vector<double>>();
    const auto ids = config.at("sector material ids").get<std::vector<int>>();
    if(angles.size() < 3 || ids.size() + 1 != angles.size())
        throw std::invalid_argument("Provide increasing sector endpoints and one material ID per interval");
    for(size_t i=0; i<angles.size(); ++i)
        if(!std::isfinite(angles[i]) || (i && angles[i] <= angles[i-1]))
            throw std::invalid_argument("angles must be finite and strictly increasing");
    constexpr double pi = 3.14159265358979323846;
    if(std::abs(angles.back()-angles.front()-2.*pi) > 1.e-6)
        throw std::invalid_argument("angles must span exactly 2*pi radians");
    return config;
}

inline double ReadCoordinate(const std::string &text)
{
    size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if(consumed != text.size() || !std::isfinite(value))
        throw std::invalid_argument("Coordinates must be finite numbers");
    return value;
}

#endif
