// ============================================================================
//  params.cpp -- input parsing and validation
// ============================================================================
#include "params.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace lbm {

namespace {
std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}
std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}
bool toBool(const std::string& v) {
    const std::string l = lower(v);
    if (l == "1" || l == "true" || l == "yes" || l == "on") return true;
    if (l == "0" || l == "false" || l == "no" || l == "off") return false;
    throw std::runtime_error("cannot interpret '" + v + "' as a boolean");
}
}  // namespace


//===========================================//    
// Reading input and assigning to parameters //
//==========================================//
void Params::set(const std::string& keyIn, const std::string& v) {
    const std::string key = lower(trim(keyIn));
    try {
        if      (key == "nx") NX = std::stoi(v);
        else if (key == "ny") NY = std::stoi(v);
        else if (key == "nz") NZ = std::stoi(v);
        else if (key == "tau") tau = std::stod(v);
        else if (key == "lambda") Lambda = std::stod(v);
        else if (key == "fx") Fx = std::stod(v);
        else if (key == "fy") Fy = std::stod(v);
        else if (key == "fz") Fz = std::stod(v);
        else if (key == "walls_y" || key == "wallsy") wallsY = toBool(v);
        else if (key == "walls_z" || key == "wallsz") wallsZ = toBool(v);
        else if (key == "steps") steps = std::stoi(v);
        else if (key == "print_every" || key == "print") printEvery = std::stoi(v);
        else if (key == "vtk_every" || key == "vtk") vtkEvery = std::stoi(v);
        else if (key == "traj_every" || key == "traj") trajEvery = std::stoi(v);
        else if (key == "tol") tol = std::stod(v);
        else if (key == "mass_correction_every") massCorrectionEvery = std::stoi(v);
        else if (key == "output_dir") outputDir = v;
        else if (key == "prefix") prefix = v;
        else if (key == "squirmer" || key == "sq") squirmer = toBool(v);
        else if (key == "radius" || key == "r") R = std::stod(v);
        else if (key == "b1") B1 = std::stod(v);
        else if (key == "beta") beta = std::stod(v);
        else if (key == "rho_p" || key == "rhop") rhoP = std::stod(v);
        else if (key == "x0") x0 = std::stod(v);
        else if (key == "y0") y0 = std::stod(v);
        else if (key == "z0") z0 = std::stod(v);
        else if (key == "angle") angle = std::stod(v);
        else if (key == "planar") planar = toBool(v);
        else if (key == "eps") eps = std::stod(v);
        else if (key == "hc") hc = std::stod(v);
        else if (key == "hmin") hmin = std::stod(v);
        else throw std::runtime_error("unknown parameter '" + key + "'");
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("invalid value '" + v + "' for parameter '" + key + "'");
    }
}


//=====================//    
// Loading  input file //
//=====================//
void Params::loadFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open input file '" + path + "'");
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        const auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = trim(line);
        if (line.empty()) continue;
        const auto eq = line.find('=');
        if (eq == std::string::npos)
            throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": expected 'key = value'");
        try {
            set(line.substr(0, eq), trim(line.substr(eq + 1)));
        } catch (const std::runtime_error& e) {
            throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": " + e.what());
        }
    }
}


//=========================//    
// Overiding input params. //
//=========================//
void Params::parseCommandLine(int argc, char** argv) {
    for (int k = 1; k < argc; ++k) {
        const std::string a = argv[k];
        const auto eq = a.find('=');
        if (eq == std::string::npos) loadFile(a);            // an input file
        else set(a.substr(0, eq), a.substr(eq + 1));         // an override
    }
}


//=======================//    
// Validating constants  //
//=======================//
void Params::validate() const {
    auto fail = [](const std::string& m) { throw std::runtime_error(m); };
    if (NX < 1 || NY < 1 || NZ < 1) fail("grid dimensions must be positive");
    if (wallsY && NY < 4) fail("ny must be >= 4 with walls in y");
    if (wallsZ && NZ < 4) fail("nz must be >= 4 with walls in z");
    if (tau <= 0.5) fail("tau must be > 0.5");
    if (Lambda <= 0.0) fail("lambda must be > 0");
    if (steps < 0) fail("steps must be >= 0");
    if (printEvery < 1 || trajEvery < 1) fail("print_every and traj_every must be >= 1");
    if (!squirmer) return;
    if (R <= 1.0) fail("squirmer radius must be > 1 lattice unit");
    if (std::fabs(angle) >= 45.0) fail("initial angle must satisfy |angle| < 45 deg from the x-axis");
    if (rhoP <= 0.0) fail("rho_p must be > 0");
    if (eps < 0.0 || hc <= 0.0 || hmin < 0.0) fail("eps >= 0, hc > 0, hmin >= 0 required");
    const double y = y0 >= 0 ? y0 : 0.5 * (NY - 1);
    if (wallsY && (y - R - 0.5 < hmin || (NY - 1.5) - y - R < hmin))
        fail("squirmer initially overlaps (or is too close to) a y-wall");
    if (2.0 * R + 2.0 > NX || (!wallsY && 2.0 * R + 2.0 > NY) || 2.0 * R + 2.0 > NZ)
        fail("domain too small for the squirmer");
    if (2.0 * B1 / 3.0 > 0.05) fail("swimming speed U0 = 2 B1 / 3 too large (Mach number); keep U0 < 0.05");
}

//=============================//    
// Printing system parameters  //
//=============================//
void Params::print(std::FILE* fp) const {
    std::fprintf(fp, "grid        : %d x %d x %d\n", NX, NY, NZ);
    std::fprintf(fp, "fluid       : tau = %.4f (nu = %.5f), TRT Lambda = %.4f, F = (%.3e, %.3e, %.3e)\n",
                 tau, (tau - 0.5) / 3.0, Lambda, Fx, Fy, Fz);
    std::fprintf(fp, "walls       : y %s, z %s (x periodic)\n", wallsY ? "on" : "off", wallsZ ? "on" : "off");
    std::fprintf(fp, "run         : %d steps, output '%s/%s_*'\n", steps, outputDir.c_str(), prefix.c_str());
    if (!squirmer) return;
    const double U0 = 2.0 * B1 / 3.0, nu = (tau - 0.5) / 3.0;
    std::fprintf(fp, "squirmer    : R = %.2f, B1 = %.4f, beta = %+.2f (%s), U0 = %.5f, Re = 2RU0/nu = %.3f\n",
                 R, B1, beta, beta < 0 ? "pusher" : (beta > 0 ? "puller" : "neutral"), U0, 2 * R * U0 / nu);
    std::fprintf(fp, "              angle = %.1f deg, rho_p = %.2f, planar = %d\n", angle, rhoP, (int)planar);
    std::fprintf(fp, "wall pot.   : eps = %.3f (F_max = eps F_S), hc = %.2f, hmin = %.2f\n", eps, hc, hmin);
}

}  // namespace lbm
