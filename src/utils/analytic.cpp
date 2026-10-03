// ============================================================================
//  analytic.cpp
// ============================================================================
#include "analytic.hpp"

#include <cmath>
#include <cstdio>

namespace lbm {

double ductVelocity(double y, double z, double a, double b, double G, double mu) {
    double s = 0.0;
    for (int n = 1; n < 400; n += 2) {
        const double k = n * M_PI / (2.0 * a), A = k * std::fabs(z), B = k * b;
        // cosh(A)/cosh(B) without overflow
        const double ratio = std::exp(A - B) * (1.0 + std::exp(-2.0 * A)) / (1.0 + std::exp(-2.0 * B));
        const double sign = ((n - 1) / 2) % 2 == 0 ? 1.0 : -1.0;
        s += sign * (1.0 - ratio) * std::cos(k * y) / ((double)n * n * n);
    }
    return 16.0 * a * a * G / (mu * M_PI * M_PI * M_PI) * s;
}

double planeVelocity(double y, double a, double G, double mu) {
    return G / (2.0 * mu) * (a * a - y * y);
}

double poiseuilleError(const LBM& L, const std::string& base) {
    if (!L.wallsY && !L.wallsZ) return -1.0;
    const double mu = L.nu, G = L.F[0];
    const double yc = 0.5 * (L.NY - 1), zc = 0.5 * (L.NZ - 1);
    const double a = 0.5 * (L.NY - 2), b = 0.5 * (L.NZ - 2);   // walls half-way to the first fluid node
    auto exact = [&](int y, int z) {
        if (L.wallsY && L.wallsZ) return ductVelocity(y - yc, z - zc, a, b, G, mu);
        if (L.wallsY) return planeVelocity(y - yc, a, G, mu);
        return planeVelocity(z - zc, b, G, mu);
    };

    double num = 0, den = 0;
    for (int z = 0; z < L.NZ; ++z)
        for (int y = 0; y < L.NY; ++y) {
            const size_t n = L.idx(0, y, z);
            if (L.flag[n] != FLUID) continue;
            const double ue = exact(y, z);
            num += (L.ux[n] - ue) * (L.ux[n] - ue);
            den += ue * ue;
        }

    const int zm = L.NZ / 2, ym = L.NY / 2;
    if (std::FILE* fp = std::fopen((base + "_profile_y.dat").c_str(), "w")) {
        std::fprintf(fp, "# y-yc  ux_lbm  ux_exact   (z = %d)\n", zm);
        for (int y = 0; y < L.NY; ++y)
            if (L.flag[L.idx(0, y, zm)] == FLUID)
                std::fprintf(fp, "%9.3f %.10e %.10e\n", y - yc, L.ux[L.idx(0, y, zm)], exact(y, zm));
        std::fclose(fp);
    }
    if (std::FILE* fp = std::fopen((base + "_profile_z.dat").c_str(), "w")) {
        std::fprintf(fp, "# z-zc  ux_lbm  ux_exact   (y = %d)\n", ym);
        for (int z = 0; z < L.NZ; ++z)
            if (L.flag[L.idx(0, ym, z)] == FLUID)
                std::fprintf(fp, "%9.3f %.10e %.10e\n", z - zc, L.ux[L.idx(0, ym, z)], exact(ym, z));
        std::fclose(fp);
    }
    return den > 0 ? std::sqrt(num / den) : -1.0;
}

}  // namespace lbm
