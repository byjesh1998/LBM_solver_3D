// ============================================================================
//  params.hpp -- simulation parameters, input-file parser, validation
//
//  Input files are plain text, one "key = value" per line, '#' starts a comment.
//  Any key can be overridden on the command line as key=value, e.g.
//      ./bin/lbm_squirmer inputs/channel_neutral.in beta=1 steps=90000
// ============================================================================
#pragma once
#include <cstdio>
#include <string>

namespace lbm {

struct Params {
    // ---- grid (lattice units; walls occupy the outermost layer when enabled) ----
    int NX = 64, NY = 34, NZ = 32;


    // ---- fluid ----
    double tau = 1.0;                 // relaxation time, nu = (tau - 1/2)/3
    double Lambda = 3.0 / 16.0;       // TRT magic parameter
    double Fx = 0.0, Fy = 0.0, Fz = 0.0;  // body-force density (drives background flow)


    // ---- boundaries (x is always periodic) ----
    bool wallsY = true;               // no-slip walls at y = 0.5 and y = NY-1.5
    bool wallsZ = false;              // no-slip walls at z = 0.5 and z = NZ-1.5


    // ---- run control ----
    int    steps = 10000;
    int    printEvery = 1000;
    int    vtkEvery = 0;              // 0: only the final field
    int    trajEvery = 50;
    double tol = 0.0;                 // steady-state tolerance (no squirmer); 0 = run all steps
    int    massCorrectionEvery = 10;  // with a squirmer; 0 = off
    std::string outputDir = "outputs";
    std::string prefix = "run";


    // ---- squirmer ----
    bool   squirmer = false;
    double R = 4.0;                   // radius
    double B1 = 0.015;                // first squirming mode, U0 = 2 B1 / 3
    double beta = 0.0;                // B2/B1: < 0 pusher, > 0 puller
    double rhoP = 1.0;                // particle / fluid density ratio
    double x0 = -1, y0 = -1, z0 = -1; // initial centre (-1: domain centre)
    double angle = 30.0;              // initial orientation in x-y plane [deg]
    //    double angle = 30.0;              // initial orientation in x-z plane [deg]
    bool   planar = true;             // restrict motion to the x-y plane


    // ---- soft wall potential: F = eps F_S (1 - h/hc)^2 for gap h < hc ----
    double eps = 0.2;                 // max force in units of F_S = 6 pi mu R U0
    double hc = 2.0;                  // range (lattice units)
    double hmin = 0.5;                // hard safety gap

    void set(const std::string& key, const std::string& value);   // throws on unknown key
    void loadFile(const std::string& path);
    void parseCommandLine(int argc, char** argv);                  // [input_file] [key=value ...]
    void validate() const;                                         // throws std::runtime_error
    void print(std::FILE* fp) const;
};

}  // namespace lbm
