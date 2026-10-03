// ============================================================================
//  lbm.cpp -- collision, streaming, bounce-back and squirmer coupling
// ============================================================================


// This code simulates a 3D fluid using TRT-LBM, applies external forcing, 
// represents a spherical squirmer as moving solid lattice nodes, 
// uses bounce-back to couple the particle to the fluid, calculates hydrodynamic 
// force/torque on the squirmer, and then moves the squirmer accordingly.


//  ┌────────────────────────────┐
//  │          LBM FLUID         │
//  │                            │
//  │  rho, ux, uy, uz           │
//  │       ↓                    │
//  │  collision + forcing       │
//  │       ↓                    │
//  │  streaming                 │
//  │       ↓                    │
//  │  bounce-back               │
//  └──────────────┬─────────────┘
//                 │
//            force + torque
//                 │
//                 |
//  ┌─────────────────────────────┐
//  │          SQUIRMER           │
//  │                             │
//  │  hydrodynamic force Fh      │
//  │  hydrodynamic torque Th     │
//  │  wall force Fw              │
//  │       ↓                     │
//  │  integrate()                │
//  │       ↓                     │
//  │  new position/orientation   │
//  └──────────────┬──────────────┘
//                 │
//            new sphere
//                 │
//                 ↓
//  ┌─────────────────────────────┐
//  │  updateParticleFlags()      │
//  │                             │
//  │  fluid ↔ particle lattice   │
//  │  nodes are updated          │
//  └─────────────────────────────┘



#include "lbm.hpp"

#include <algorithm>
#include <cmath>

namespace lbm {

using namespace d3q19;

LBM::LBM(const Params& p)
    : NX(p.NX), NY(p.NY), NZ(p.NZ), N((size_t)p.NX * p.NY * p.NZ), tau(p.tau),
      F{p.Fx, p.Fy, p.Fz}, wallsY(p.wallsY), wallsZ(p.wallsZ),
      f(Q * N), fnew(Q * N), rho(N, 1.0), ux(N, 0.0), uy(N, 0.0), uz(N, 0.0), flag(N, FLUID)
{
    // TRT: omega+ sets the viscosity, omega- follows from the magic parameter
    // Lambda = (1/omega+ - 1/2)(1/omega- - 1/2).
    omPlus = 1.0 / tau;
    omMinus = 1.0 / (0.5 + p.Lambda / (tau - 0.5));
    nu = (tau - 0.5) / 3.0;

    //lattice cell is classified as either a wall or fluid
    for (int z = 0; z < NZ; ++z)
        for (int y = 0; y < NY; ++y)
            for (int x = 0; x < NX; ++x) {
                const bool wall = (wallsY && (y == 0 || y == NY - 1)) ||
                                  (wallsZ && (z == 0 || z == NZ - 1));
                flag[idx(x, y, z)] = wall ? WALL : FLUID;
            }


    // fluid at rest, rho = 1
    for (size_t n = 0; n < N; ++n)
        for (int i = 0; i < Q; ++i) f[fi(i, n)] = w[i];
    fnew = f;

    if (p.squirmer) {
        hasSquirmer = true;
        sq.init(p, nu,
                p.x0 >= 0 ? p.x0 : 0.5 * (NX - 1),
                p.y0 >= 0 ? p.y0 : 0.5 * (NY - 1),
                p.z0 >= 0 ? p.z0 : 0.5 * (NZ - 1));
        updateParticleFlags(true);
    }
}


//=======================//    
// Initializing squirmer //
//=======================//  
void LBM::offsetFromSquirmer(double px, double py, double pz, double r[3]) const {
    r[0] = px - sq.X[0]; r[0] -= NX * std::round(r[0] / NX);
    r[1] = py - sq.X[1]; if (!wallsY) r[1] -= NY * std::round(r[1] / NY);
    r[2] = pz - sq.X[2]; if (!wallsZ) r[2] -= NZ * std::round(r[2] / NZ);
}


//================================//    
// Time steps                     //
//================================// 
//          ┌────────────────┐
//          │ f distributions│
//          └───────┬────────┘
//                  ↓
//          calculate rho, u
//                  ↓
//              collision
//                  ↓
//              streaming
//                  ↓
//        ┌─────────┴─────────┐
//        ↓                   ↓
//     fluid cell          solid cell
//        ↓                   ↓
//   move normally       bounce back
//                            ↓
//                      force/torque
// ----------------------------------------------------------------------------
void LBM::step() {
    const double gP = 1.0 - 0.5 * omPlus, gM = 1.0 - 0.5 * omMinus;
    const double Fx = F[0], Fy = F[1], Fz = F[2];
    double FWx = 0, FWy = 0, FWz = 0;
    double FPx = 0, FPy = 0, FPz = 0, TPx = 0, TPy = 0, TPz = 0;

    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:FWx,FWy,FWz,FPx,FPy,FPz,TPx,TPy,TPz)
    for (int z = 0; z < NZ; ++z)
        for (int y = 0; y < NY; ++y)
            for (int x = 0; x < NX; ++x) {
                const size_t n = idx(x, y, z);
                if (flag[n] != FLUID) continue;

                // ---- moments (Guo: velocity includes half the force) ----
                double fl[Q];
                double r = 0, mx = 0, my = 0, mz = 0;
                for (int i = 0; i < Q; ++i) {
                    fl[i] = f[fi(i, n)]; //density
                    r += fl[i]; mx += cx[i] * fl[i]; my += cy[i] * fl[i]; mz += cz[i] * fl[i]; //velocity --> rho * u
                }


                const double u = (mx + 0.5 * Fx) / r, v = (my + 0.5 * Fy) / r, wz = (mz + 0.5 * Fz) / r; //velocity from force
                rho[n] = r; ux[n] = u; uy[n] = v; uz[n] = wz; //u=(\sum f_i * c_i +.5 F)/rho
                const double usq = u * u + v * v + wz * wz, uF = u * Fx + v * Fy + wz * Fz;



                //==========================================================================//    
                // Collision                                                                //
                //                                                                          //
                // f_i(post)=f_i - w_+ (f_i^+ - f_i^+(eq)- w_- (f_i^- - f_i^-(eq)) +forcing //
                //==========================================================================// 
                for (int i = 0; i < Q; ++i) {
                    // ---- TRT collision with Guo forcing ----
                    const int io = opp[i];
                    const double cu = cx[i] * u + cy[i] * v + cz[i] * wz;
                    const double cF = cx[i] * Fx + cy[i] * Fy + cz[i] * Fz;
                    const double feqP = w[i] * r * (1.0 + 4.5 * cu * cu - 1.5 * usq); //split into symmetric and antisymmetric
                    const double feqM = w[i] * r * 3.0 * cu;
                    const double fP = 0.5 * (fl[i] + fl[io]), fM = 0.5 * (fl[i] - fl[io]);
                    const double SP = w[i] * (9.0 * cu * cF - 3.0 * uF), SM = w[i] * 3.0 * cF;
                    const double fpost = fl[i] - omPlus * (fP - feqP) - omMinus * (fM - feqM)
                                         + gP * SP + gM * SM; //forcing correction



                    //=======================//    
                    // Streaming             //
                    //=======================// 
                    // ---- streaming (periodic wrap; walls are solid nodes) ----
                    int xn = x + cx[i], yn = y + cy[i], zn = z + cz[i];
                    if (xn < 0) xn += NX; else if (xn >= NX) xn -= NX; //periodic boundary condition
                    if (yn < 0) yn += NY; else if (yn >= NY) yn -= NY;
                    if (zn < 0) zn += NZ; else if (zn >= NZ) zn -= NZ;
                    const size_t nb = idx(xn, yn, zn);
                    const uint8_t t = flag[nb];
                    if (t == FLUID) { fnew[fi(i, nb)] = fpost; continue; }


                    // ---- halfway bounce-back with moving-boundary correction ----
                    double ub[3] = {0, 0, 0}, rb[3] = {0, 0, 0};
                    if (t == PARTICLE) {
                        offsetFromSquirmer(x + 0.5 * cx[i], y + 0.5 * cy[i], z + 0.5 * cz[i], rb);
                        sq.surfaceVelocity(rb, ub, true);// moving-boundary halfway bounce-back correction.
                    } 


                    const double cub = cx[i] * ub[0] + cy[i] * ub[1] + cz[i] * ub[2];
                    const double fb = fpost - 6.0 * w[i] * r * cub;
                    fnew[fi(io, n)] = fb;


                    // ---- Galilean-invariant momentum exchange (Wen et al. 2014) ----
                    const double s = fpost + fb, dl = fpost - fb;
                    const double dpx = s * cx[i] - dl * ub[0];
                    const double dpy = s * cy[i] - dl * ub[1];
                    const double dpz = s * cz[i] - dl * ub[2];


                    if (t == WALL) { FWx += dpx; FWy += dpy; FWz += dpz; } //momentum transferred between the fluid and boundary.
                    else {
                        FPx += dpx; FPy += dpy; FPz += dpz;
                        TPx += rb[1] * dpz - rb[2] * dpy; //hydrodynamic torque
                        TPy += rb[2] * dpx - rb[0] * dpz;
                        TPz += rb[0] * dpy - rb[1] * dpx;
                    }
                }
            }
    f.swap(fnew);
    wallForce[0] = FWx; wallForce[1] = FWy; wallForce[2] = FWz; //if its squirmer, both force and torque handed to swimmer
    if (hasSquirmer) {
        sq.Fh[0] = FPx; sq.Fh[1] = FPy; sq.Fh[2] = FPz;
        sq.Th[0] = TPx; sq.Th[1] = TPy; sq.Th[2] = TPz;
    }
}





// ----------------------------------------------------------------------------
void LBM::applyWallPotential() {
    sq.Fw[0] = sq.Fw[1] = sq.Fw[2] = 0.0;
    if (wallsY) {
        sq.Fw[1] += softWallForce(sq.X[1] - sq.R - 0.5, sq.eps, sq.FS, sq.hc);
        sq.Fw[1] -= softWallForce((NY - 1.5) - sq.X[1] - sq.R, sq.eps, sq.FS, sq.hc);
    }
    if (wallsZ) {
        sq.Fw[2] += softWallForce(sq.X[2] - sq.R - 0.5, sq.eps, sq.FS, sq.hc);
        sq.Fw[2] -= softWallForce((NZ - 1.5) - sq.X[2] - sq.R, sq.eps, sq.FS, sq.hc);
    }
}


//Preventing the squirmer from getting too close
void LBM::enforceMinimumGap() {
    auto clampAxis = [&](int d, int Nd) {
        const double lo = 0.5 + sq.R + sq.hmin, hi = (Nd - 1.5) - sq.R - sq.hmin;
        if (sq.X[d] < lo) { sq.X[d] = lo; if (sq.U[d] < 0) sq.U[d] = 0; ++sq.contacts; }
        if (sq.X[d] > hi) { sq.X[d] = hi; if (sq.U[d] > 0) sq.U[d] = 0; ++sq.contacts; }
    };
    if (wallsY) clampAxis(1, NY);
    if (wallsZ) clampAxis(2, NZ);
}



void LBM::advanceSquirmer() {
    if (!hasSquirmer) return;
    applyWallPotential();
    sq.integrate();
    enforceMinimumGap();
    updateParticleFlags(false);
}

//=======================//    
// Remapping             //
//=======================//  
// Re-map the sphere onto the lattice. Nodes swallowed by the sphere give their
// momentum to it; released nodes are filled with the equilibrium at the local
// rigid-body velocity, whose momentum is taken from the sphere (Aidun et al. 1998).
void LBM::updateParticleFlags(bool initial) {
    const double R = sq.R;
    const int x0 = (int)std::floor(sq.X[0] - R) - 2, x1 = (int)std::ceil(sq.X[0] + R) + 2;
    int y0 = (int)std::floor(sq.X[1] - R) - 2, y1 = (int)std::ceil(sq.X[1] + R) + 2;
    const int z0 = (int)std::floor(sq.X[2] - R) - 2, z1 = (int)std::ceil(sq.X[2] + R) + 2;
    if (wallsY) { y0 = std::max(0, y0); y1 = std::min(NY - 1, y1); }

    for (int zz = z0; zz <= z1; ++zz)
        for (int yy = y0; yy <= y1; ++yy)
            for (int xx = x0; xx <= x1; ++xx) {
                const int x = ((xx % NX) + NX) % NX;
                const int y = ((yy % NY) + NY) % NY;
                const int z = ((zz % NZ) + NZ) % NZ;
                const size_t n = idx(x, y, z);
                if (flag[n] == WALL) continue;
                double r[3];
                offsetFromSquirmer(x, y, z, r);
                const bool inside = r[0] * r[0] + r[1] * r[1] + r[2] * r[2] < R * R;
                const bool was = flag[n] == PARTICLE;
                if (inside == was) continue;
                if (initial) { flag[n] = inside ? PARTICLE : FLUID; continue; }

                if (inside) {
                    double p[3] = {0, 0, 0}, t[3]; //For a node swallowed by the particle, the code extracts its fluid momentum
                    for (int i = 0; i < Q; ++i) {
                        const double fv = f[fi(i, n)];
                        p[0] += fv * cx[i]; p[1] += fv * cy[i]; p[2] += fv * cz[i];
                    }
                    cross(r, p, t);
                    for (int d = 0; d < 3; ++d) { sq.Fc[d] += p[d]; sq.Tc[d] += t[d]; }
                    flag[n] = PARTICLE;
                } else {
                    double rs = 0; int cnt = 0;
                    for (int i = 1; i < Q; ++i) {
                        const int xn = (x + cx[i] + NX) % NX, yn = (y + cy[i] + NY) % NY, zn = (z + cz[i] + NZ) % NZ;
                        const size_t nb = idx(xn, yn, zn);
                        if (flag[nb] == FLUID) { rs += rho[nb]; ++cnt; }
                    }
                    const double rr = cnt ? rs / cnt : 1.0;
                    double ub[3];
                    sq.surfaceVelocity(r, ub, false);
                    for (int i = 0; i < Q; ++i) f[fi(i, n)] = equilibrium(i, rr, ub[0], ub[1], ub[2]);
                    rho[n] = rr; ux[n] = ub[0]; uy[n] = ub[1]; uz[n] = ub[2];
                    double p[3] = {rr * ub[0], rr * ub[1], rr * ub[2]}, t[3];
                    cross(r, p, t);
                    for (int d = 0; d < 3; ++d) { sq.Fc[d] -= p[d]; sq.Tc[d] -= t[d]; }
                    flag[n] = FLUID;
                }
            }
}


// ----------------------------------------------------------------------------
void LBM::correctMass() {
    double m = 0; long nf = 0;
    #pragma omp parallel for reduction(+:m,nf)
    for (long n = 0; n < (long)N; ++n) {
        if (flag[n] != FLUID) continue;
        for (int i = 0; i < Q; ++i) m += f[fi(i, n)];
        ++nf;
    }
    if (nf == 0) return;
    const double d = 1.0 - m / nf;
    #pragma omp parallel for
    for (long n = 0; n < (long)N; ++n) {
        if (flag[n] != FLUID) continue;
        for (int i = 0; i < Q; ++i) f[fi(i, n)] += w[i] * d;// bring the average fluid density back toward rho=1
    }
}

// cout fluid nodes 
long LBM::countFluid() const { return (long)std::count(flag.begin(), flag.end(), (uint8_t)FLUID); }

// calculate mean density
double LBM::meanDensity() const {
    double m = 0; long nf = 0;
    for (size_t n = 0; n < N; ++n) if (flag[n] == FLUID) { m += rho[n]; ++nf; }
    return nf ? m / nf : 0.0;
}

// calculate maximum speed
double LBM::maxSpeed() const {
    double m = 0;
    for (size_t n = 0; n < N; ++n)
        if (flag[n] == FLUID) m = std::max(m, std::sqrt(ux[n] * ux[n] + uy[n] * uy[n] + uz[n] * uz[n]));
    return m;
}

// calculate maximum speed in x direction
double LBM::sumUx() const {
    double s = 0;
    for (size_t n = 0; n < N; ++n) if (flag[n] == FLUID) s += ux[n];
    return s;
}

}  // namespace lbm
