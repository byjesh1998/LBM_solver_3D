// ============================================================================
//  io.cpp 
//
// defining how to write output files
// ============================================================================
#include "io.hpp"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace lbm {

namespace {
// VTK legacy binary is big-endian.
bool hostIsLittleEndian() {
    const uint16_t one = 1;
    return *reinterpret_cast<const uint8_t*>(&one) == 1;
}
void writeBigEndianDoubles(std::FILE* fp, const std::vector<double>& v) {
    std::vector<uint64_t> buf(v.size());
    const bool swap = hostIsLittleEndian();
    for (size_t k = 0; k < v.size(); ++k) {
        uint64_t u;
        std::memcpy(&u, &v[k], 8);
        buf[k] = swap ? __builtin_bswap64(u) : u;
    }
    std::fwrite(buf.data(), 8, buf.size(), fp);
}
}  // namespace

void ensureDirectory(const std::string& dir) {
    if (dir.empty()) return;
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) throw std::runtime_error("cannot create output directory '" + dir + "': " + ec.message());
}

//=======================//    
// vtk files.            //
//=======================//
void writeVTK(const LBM& L, const std::string& path) {
    std::FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) throw std::runtime_error("cannot write '" + path + "'");
    std::fprintf(fp, "# vtk DataFile Version 3.0\nlbm-squirmer\nBINARY\nDATASET STRUCTURED_POINTS\n");
    std::fprintf(fp, "DIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA %zu\n", L.NX, L.NY, L.NZ, L.N);

    std::fprintf(fp, "SCALARS flag unsigned_char 1\nLOOKUP_TABLE default\n");
    std::fwrite(L.flag.data(), 1, L.N, fp);
    std::fprintf(fp, "\n");

    std::vector<double> buf(L.N);
    for (size_t n = 0; n < L.N; ++n) buf[n] = L.flag[n] == FLUID ? L.rho[n] : 1.0;
    std::fprintf(fp, "SCALARS density double 1\nLOOKUP_TABLE default\n");
    writeBigEndianDoubles(fp, buf);
    std::fprintf(fp, "\n");

    buf.assign(3 * L.N, 0.0);
    for (size_t n = 0; n < L.N; ++n)
        if (L.flag[n] == FLUID) { buf[3 * n] = L.ux[n]; buf[3 * n + 1] = L.uy[n]; buf[3 * n + 2] = L.uz[n]; }
    std::fprintf(fp, "VECTORS velocity double\n");
    writeBigEndianDoubles(fp, buf);
    std::fprintf(fp, "\n");
    std::fclose(fp);
}

//==========================//    
// trajectory of squirmer   //
//==========================//
TrajectoryWriter::TrajectoryWriter(const std::string& path) {
    fp_ = std::fopen(path.c_str(), "w");
    if (!fp_) throw std::runtime_error("cannot write '" + path + "'");
    std::fprintf(fp_, "# t X Y Z ex ey ez Ux Uy Uz Wx Wy Wz Fwall_y Fwall_z\n");
}
TrajectoryWriter::~TrajectoryWriter() { if (fp_) std::fclose(fp_); }
void TrajectoryWriter::flush() { if (fp_) std::fflush(fp_); }
void TrajectoryWriter::write(int t, const Squirmer& s) {
    std::fprintf(fp_, "%d %.6f %.6f %.6f %.6f %.6f %.6f %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                 t, s.X[0], s.X[1], s.X[2], s.e[0], s.e[1], s.e[2],
                 s.U[0], s.U[1], s.U[2], s.W[0], s.W[1], s.W[2], s.Fw[1], s.Fw[2]);
}



void writeSquirmerState(const LBM& L, int t, const std::string& path) {
    std::FILE* fp = std::fopen(path.c_str(), "w");
    if (!fp) throw std::runtime_error("cannot write '" + path + "'");
    const Squirmer& s = L.sq;
    std::fprintf(fp, "t %d\nR %.10g\nB1 %.10g\nbeta %.10g\nU0 %.10g\nnu %.10g\nFS %.10g\n",
                 t, s.R, s.B1, s.beta, s.U0(), L.nu, s.FS);
    std::fprintf(fp, "X %.10g %.10g %.10g\ne %.10g %.10g %.10g\n", s.X[0], s.X[1], s.X[2], s.e[0], s.e[1], s.e[2]);
    std::fprintf(fp, "U %.10g %.10g %.10g\nW %.10g %.10g %.10g\ncontacts %ld\n",
                 s.U[0], s.U[1], s.U[2], s.W[0], s.W[1], s.W[2], s.contacts);
    std::fclose(fp);
}

}  // namespace lbm
