// ============================================================================
//  lbm-squirmer: 3D lattice Boltzmann (D3Q19, TRT) solver for channel flow
//  with a self-propelled spherical squirmer.
//
//  Usage:  ./bin/lbm_squirmer <input_file> [key=value ...]
//
//  Outputs (in output_dir):
//    <prefix>_final.vtk          final flow field
//    <prefix>_<step>.vtk         snapshots (vtk_every > 0)
//    <prefix>_traj.dat           squirmer trajectory (squirmer = 1)
//    <prefix>_squirmer.dat       squirmer state at the end of the run
//    <prefix>_profile_{y,z}.dat  Poiseuille profiles vs exact (no squirmer)
//    <prefix>_summary.txt        key diagnostics, machine readable
// ============================================================================
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

#include "utils/analytic.hpp"
#include "utils/io.hpp"
#include "utils/lbm.hpp"
#include "utils/params.hpp"

using namespace lbm;

//=========================//    
// Running the simulation  //
//=========================//
int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <input_file> [key=value ...]\n", argv[0]);
        return 1;
    }
    Params p;
    try {
        p.parseCommandLine(argc, argv);
        p.validate();
        ensureDirectory(p.outputDir);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    const std::string base = p.outputDir + "/" + p.prefix;

    std::printf("================ lbm-squirmer ================\n");
    p.print(stdout);
    std::printf("==============================================\n");

    LBM L(p);
    std::unique_ptr<TrajectoryWriter> traj;
    if (L.hasSquirmer) traj = std::make_unique<TrajectoryWriter>(base + "_traj.dat");

    const auto t0 = std::chrono::steady_clock::now();
    double prevSum = 0.0;
    int t = 0;
    bool converged = false;
    for (t = 1; t <= p.steps; ++t) {
        L.step();
        if (L.hasSquirmer) {
            L.advanceSquirmer();
            if (p.massCorrectionEvery > 0 && t % p.massCorrectionEvery == 0) L.correctMass();
            if (t % p.trajEvery == 0) traj->write(t, L.sq);
        }
        if (p.vtkEvery > 0 && t % p.vtkEvery == 0) writeVTK(L, base + "_" + std::to_string(t) + ".vtk");

        if (t % p.printEvery == 0) {
            if (L.hasSquirmer) {
                const Squirmer& s = L.sq;
                std::printf("t=%8d  X=(%9.2f,%7.2f,%7.2f)  angle=%7.2f deg  |U|/U0=%.3f  <rho>=%.6f\n", t,
                            s.X[0], s.X[1], s.X[2], std::atan2(s.e[1], s.e[0]) * 180.0 / M_PI,
                            std::sqrt(s.U[0] * s.U[0] + s.U[1] * s.U[1] + s.U[2] * s.U[2]) / s.U0(),
                            L.meanDensity());
                traj->flush();
            } else {
                const double sx = L.sumUx();
                const double change = std::fabs(sx - prevSum) / (std::fabs(sx) + 1e-30);
                prevSum = sx;
                std::printf("t=%8d  u_max=%.6e  F_wall_x=%.6e  rel.change=%.3e\n", t, L.maxSpeed(),
                            L.wallForce[0], change);
                if (p.tol > 0 && change < p.tol) { converged = true; std::fflush(stdout); break; }
            }
            std::fflush(stdout);
        }
    }
    const int done = converged ? t : p.steps;
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("%d steps in %.1f s (%.2f MLUPS)\n", done, sec, sec > 0 ? 1e-6 * (double)L.N * done / sec : 0.0);

    // ---- final output ----
    writeVTK(L, base + "_final.vtk");
    std::FILE* sum = std::fopen((base + "_summary.txt").c_str(), "w");
    std::fprintf(sum, "steps %d\nconverged %d\nnu %.10g\nmean_density %.10g\nmax_speed %.10g\n",
                 done, (int)converged, L.nu, L.meanDensity(), L.maxSpeed());
    std::fprintf(sum, "wall_force %.10g %.10g %.10g\n", L.wallForce[0], L.wallForce[1], L.wallForce[2]);
    std::fprintf(sum, "body_force_total %.10g %.10g %.10g\n",
                 p.Fx * L.countFluid(), p.Fy * L.countFluid(), p.Fz * L.countFluid());
    if (L.hasSquirmer) {
        writeSquirmerState(L, done, base + "_squirmer.dat");
        std::printf("hard-gap safety activations: %ld steps\n", L.sq.contacts);
        std::fprintf(sum, "contacts %ld\n", L.sq.contacts);
    } else {
        const double err = poiseuilleError(L, base);
        if (err >= 0) {
            std::printf("relative L2 error vs exact Poiseuille solution: %.4e\n", err);
            std::fprintf(sum, "l2_error %.10g\n", err);
        }
    }
    std::fclose(sum);
    std::printf("output written to %s_*\n", base.c_str());
    return 0;
}
