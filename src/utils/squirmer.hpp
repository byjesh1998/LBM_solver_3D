// ============================================================================
//  squirmer.hpp -- spherical squirmer (Lighthill 1952, Blake 1971)
//
//  Surface slip:  u_s = B1 (1 + beta cos th) (cos th r^ - e),   cos th = e . r^
//  Swimming speed in free space: U0 = 2 B1 / 3.
//  beta < 0: pusher,  beta = 0: neutral,  beta > 0: puller.
// ============================================================================
#pragma once

namespace lbm {

struct Params;

struct Squirmer {
    // model
    double R = 4.0, B1 = 0.015, beta = 0.0; //default values
    double M = 0.0, I = 0.0;          // mass and moment of inertia
    bool   planar = true;


    // state (X[0] is kept unwrapped so trajectories are continuous)
    double X[3] = {0, 0, 0};  //positions
    double e[3] = {1, 0, 0};  //orientations
    double U[3] = {0, 0, 0};  //velocities
    double W[3] = {0, 0, 0};  //angular velocities


    // forces & torques acting during the current step
    double Fh[3] = {0, 0, 0}, Th[3] = {0, 0, 0};  // hydrodynamic (momentum exchange)
    double Fc[3] = {0, 0, 0}, Tc[3] = {0, 0, 0};  // node covering / uncovering
    double Fw[3] = {0, 0, 0};                     // soft wall potential


    // soft wall potential
    double eps = 0.2, hc = 2.0, hmin = 0.5;
    double FS = 0.0;                  // Stokes drag scale 6 pi mu R U0
    long   contacts = 0;              // steps in which the hard gap had to act


    void init(const Params& p, double nu, double xc, double yc, double zc);

    double U0() const { return 2.0 * B1 / 3.0; }

    // Rigid-body velocity U + W x r at offset r from the centre, plus the
    // squirming slip evaluated at the direction r^ if withSlip is true.
    void surfaceVelocity(const double r[3], double ub[3], bool withSlip) const;

    // Newton-Euler update using Fh + Fc + Fw and Th + Tc; resets Fc, Tc.
    void integrate();
};


// Magnitude of the soft wall force at surface-wall gap h:  eps F_S (1 - h/hc)^2 for h < hc.
double softWallForce(double h, double eps, double FS, double hc);

}  // namespace lbm
