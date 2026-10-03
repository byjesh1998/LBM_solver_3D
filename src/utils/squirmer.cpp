// ============================================================================
//  squirmer.cpp
// 
// Initialize squirer and define properties like surface velocity etc...
// ============================================================================

#include "squirmer.hpp"

#include <algorithm>
#include <cmath>

#include "lattice.hpp"
#include "params.hpp"

namespace lbm {

//=======================//    
// Initializing squirmer //
//=======================//  
void Squirmer::init(
    const Params& p, double nu, double xc, double yc, double zc) {
    R = p.R; B1 = p.B1; beta = p.beta; planar = p.planar; //reading contstants
    M = p.rhoP * 4.0 / 3.0 * M_PI * R * R * R; //mass
    I = 0.4 * M * R * R; //moment of inertia
    X[0] = xc; X[1] = yc; X[2] = zc; //initial positions
    const double a = p.angle * M_PI / 180.0; //orientation angle measured from x-axis
    e[0] = std::cos(a); e[1] = std::sin(a); e[2] = 0.0; //orientation vector (change it to general 3D case)
    eps = p.eps; hc = p.hc; hmin = p.hmin; //constant for soft potential 
    FS = 6.0 * M_PI * nu * R * U0();   // rho0 = 1
}

//=======================//    
// Define surface slip   //
//=======================//
void Squirmer::surfaceVelocity(const double r[3], double ub[3], bool withSlip) const {
    double wr[3];
    cross(W, r, wr);
    for (int d = 0; d < 3; ++d) ub[d] = U[d] + wr[d];
    if (!withSlip) return;
    const double rn = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]); //magnitude of radial vector
    if (rn < 1e-12) return;
    const double rh[3] = {r[0] / rn, r[1] / rn, r[2] / rn}; //radial unit vector
    const double ct = e[0] * rh[0] + e[1] * rh[1] + e[2] * rh[2]; // cos \theta
    const double amp = B1 * (1.0 + beta * ct); 
    for (int d = 0; d < 3; ++d) ub[d] += amp * (ct * rh[d] - e[d]); //surface slip velocity
}


//======================================//    
// Updating position and orientations   //
//======================================//
void Squirmer::integrate() {
    for (int d = 0; d < 3; ++d) {
        U[d] += (Fh[d] + Fc[d] + Fw[d]) / M; //Velocity updates
        W[d] += (Th[d] + Tc[d]) / I; //angular velocity updates
        Fc[d] = Tc[d] = 0.0;
    }
    if (planar) { U[2] = 0.0; W[0] = 0.0; W[1] = 0.0; } //for 2D case

    for (int d = 0; d < 3; ++d) X[d] += U[d]; //position updates

    double we[3];
    cross(W, e, we);  //cross for orientation rate
    for (int d = 0; d < 3; ++d) e[d] += we[d]; //orientation update 
    if (planar) e[2] = 0.0;
    const double n = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
    for (int d = 0; d < 3; ++d) e[d] /= n; //normalizing updated orientation vector
}

//=======================//    
// Soft wall potential   //
//=======================//
double softWallForce(double h, double eps, double FS, double hc) {
    if (h >= hc) return 0.0;
    const double s = 1.0 - std::max(h, 0.0) / hc;
    return eps * FS * s * s;
}

}  // namespace lbm
