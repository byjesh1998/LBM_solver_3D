// ============================================================================
//  lattice.hpp -- D3Q19 velocity set, weights and node types
//
// 3 spatial dimensions and 19 decrete velocity points
// ============================================================================



#pragma once
#include <cstdint>


namespace lbm {
//=============================//    
// Define the weights in d3q19 //
//=============================//  
// Discrete velocities are stored in opposite pairs: opp[i] = i + 1 (odd i), i - 1 (even i).
namespace d3q19 {
inline constexpr int Q = 19;
inline constexpr int cx[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0}; // (cx,cy,cz) gives the direction vector
inline constexpr int cy[Q] = {0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1};
inline constexpr int cz[Q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1};
inline constexpr int opp[Q] = {0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17};
inline constexpr double w[Q] = {1.0/3.0, //gives the weights for each direction, their sum should be unity
    1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18,
    1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36,
    1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36};

inline constexpr double cs2 = 1.0 / 3.0;   // lattice speed of sound squared
}  // namespace d3q19


//=======================//    
// Three node types.     //
//=======================//  
enum NodeType : uint8_t { FLUID = 0, WALL = 1, PARTICLE = 2 };

//======================================================================================================//    
// Equilibrium distribution  f_i(eq)=w_i rho [ 1 + (c_i . u)/cs2 + (c_i . u)^2/2cs4 - (u . u)/2cs2]. //
//======================================================================================================//  
// Second-order equilibrium distribution for direction i.
inline double equilibrium(int i, double rho, double ux, double uy, double uz) {
    using namespace d3q19;
    const double cu = cx[i] * ux + cy[i] * uy + cz[i] * uz;
    const double usq = ux * ux + uy * uy + uz * uz;
    return w[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usq);
}

//=================================//    
// Definition of cross product     //
//=================================//
inline void cross(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

}  // namespace lbm
