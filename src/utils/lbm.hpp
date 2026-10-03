// ============================================================================
//  lbm.hpp -- D3Q19 TRT lattice Boltzmann solver with walls and a squirmer
//
//  Layout: populations are stored structure-of-arrays, f[i * N + n], with the
//  node index n = x + NX * (y + NY * z). x is always periodic; y and z are
//  periodic unless walls are switched on, in which case the outermost node
//  layer is solid and the no-slip plane lies half-way to the first fluid node.
// ============================================================================
#pragma once
#include <cstdint>
#include <vector>

#include "lattice.hpp"
#include "params.hpp"
#include "squirmer.hpp"

namespace lbm {

class LBM {
public:
    explicit LBM(const Params& p);

    // One LBM time step: TRT collision (+ Guo forcing) fused with push
    // streaming and halfway bounce-back on WALL and PARTICLE nodes. Also
    // accumulates the momentum-exchange force on the walls and the squirmer.

    void step();

    // Squirmer update: soft wall force, Newton-Euler integration, hard gap,
    // and re-mapping of the sphere onto the lattice.
    void advanceSquirmer();

    // Restore the mean fluid density to 1 (moving boundaries leak a little mass).
    void correctMass();


    // ---- diagnostics ----
    long   countFluid() const;
    double meanDensity() const;
    double maxSpeed() const;
    double sumUx() const;




//change 3d vector into 1d vector
// (x,y,z)
//    │
//    ↓
// idx(x,y,z)
//    │
//    ↓
// single number n, n=x + N_x (y+ N_y z)
//    │
//    ↓
// vector[n]

    inline size_t idx(int x, int y, int z) const { return x + (size_t)NX * (y + (size_t)NY * z); }
    inline size_t fi(int i, size_t n) const { return (size_t)i * N + n; } //Get distribution in direction i at lattice node n

    // Minimum-image offset of a point from the squirmer centre.
    void offsetFromSquirmer(double px, double py, double pz, double r[3]) const; // minimum-image calculation finds the shortest periodic displacement rather than treating them as far apart.

    // ---- data (public for I/O and tests) ----
    int NX, NY, NZ;
    size_t N;
    double tau, omPlus, omMinus, nu;
    double F[3]; //external force
    bool wallsY, wallsZ;

    std::vector<double>  f, fnew;
    std::vector<double>  rho, ux, uy, uz;
    std::vector<uint8_t> flag;

    double wallForce[3] = {0, 0, 0};

    bool hasSquirmer = false;
    Squirmer sq;

private:
    void updateParticleFlags(bool initial);
    void applyWallPotential();
    void enforceMinimumGap();
};

}  // namespace lbm
