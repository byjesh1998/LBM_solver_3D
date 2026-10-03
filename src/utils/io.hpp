// ============================================================================
//  io.hpp -- output: VTK fields, squirmer trajectory, final squirmer state
// ============================================================================
#pragma once
#include <cstdio>
#include <string>

#include "lbm.hpp"


namespace lbm {

// Legacy VTK (STRUCTURED_POINTS, binary big-endian): flag (uint8),
// density (float64), velocity (3 x float64). Opens in ParaView/VisIt.
void writeVTK(const LBM& L, const std::string& path);

// Trajectory file: one row every traj_every steps, columns
//   t  X Y Z  ex ey ez  Ux Uy Uz  Wx Wy Wz  Fwall_y Fwall_z
class TrajectoryWriter {
public:
    explicit TrajectoryWriter(const std::string& path);
    ~TrajectoryWriter();
    void write(int t, const Squirmer& s);
    void flush();
private:
    std::FILE* fp_ = nullptr;
};

// "key value" text file with the squirmer state at the end of the run
// (used by the post-processing to locate the sphere in the final VTK field).
void writeSquirmerState(const LBM& L, int t, const std::string& path);

// Make sure the output directory exists.
void ensureDirectory(const std::string& dir);

}  // namespace lbm
