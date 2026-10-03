// ============================================================================
//  analytic.hpp -- exact Poiseuille solutions used for validation
// ============================================================================
#pragma once
#include <string>

#include "lbm.hpp"

namespace lbm {

// Rectangular duct -a < y < a, -b < z < b, driven by force density G along x
// (series solution, e.g. F. M. White, "Viscous Fluid Flow").
double ductVelocity(double y, double z, double a, double b, double G, double mu);

// Plane channel -a < y < a:  u = G (a^2 - y^2) / (2 mu).
double planeVelocity(double y, double a, double G, double mu);

// Relative L2 error of u_x over the cross-section x = 0 against the exact
// solution (duct if walls in y and z, plane channel if walls in one direction).
// Writes centre-line profiles to <base>_profile_y.dat and <base>_profile_z.dat.
// Returns a negative value if there are no walls.
double poiseuilleError(const LBM& L, const std::string& base);

}  // namespace lbm
