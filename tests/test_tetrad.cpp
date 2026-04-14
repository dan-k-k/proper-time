// tests/test_tetrad.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "RelSpacetime/LocalFrame.h"
#include "RelSpacetime/MetricTensor.h"

using namespace RelSpacetime;

bool isClose(double a, double b, double epsilon = 1e-6) {
    return std::abs(a - b) < epsilon;
}

// Helper <A, B> = A^T * g * B
double innerProduct(const Eigen::Vector4d& a, const Eigen::Vector4d& b, const Eigen::Matrix4d& g) {
    return a.transpose() * g * b;
}

int main() {
    std::cout << "Running LocalFrame & Tetrad Tests...\n";

    MinkowskiMetric flatSpace;
    Eigen::Matrix4d g_flat = flatSpace.getMetricAt(SpacetimeVector(0,0,0,0));

    // Test 1: Orthonormality at High Velocities (Relativistic Gram-Schmidt)
    // Ship moving at 0.6c along the X axis: gamma = 1.25.
    Eigen::Vector4d highSpeedU(1.25, 0.75, 0.0, 0.0);
    Eigen::Vector3d fwd(1.0, 0.0, 0.0);
    Eigen::Vector3d rgt(0.0, 1.0, 0.0);
    Eigen::Vector3d up(0.0, 0.0, 1.0);

    LocalFrame speedFrame(g_flat, highSpeedU, fwd, rgt, up);

    Eigen::Vector4d e0 = speedFrame.getE0();
    Eigen::Vector4d e1 = speedFrame.getE1();
    Eigen::Vector4d e2 = speedFrame.getE2();
    Eigen::Vector4d e3 = speedFrame.getE3();

    assert(isClose(innerProduct(e0, e0, g_flat), -1.0) && "FAIL: e0 must have magnitude -1 (Timelike)");
    assert(isClose(innerProduct(e1, e1, g_flat), 1.0) && "FAIL: e1 must have magnitude +1 (Spacelike)");
    assert(isClose(innerProduct(e2, e2, g_flat), 1.0) && "FAIL: e2 must have magnitude +1 (Spacelike)");
    assert(isClose(innerProduct(e3, e3, g_flat), 1.0) && "FAIL: e3 must have magnitude +1 (Spacelike)");

    assert(isClose(innerProduct(e0, e1, g_flat), 0.0) && "FAIL: e1 is not orthogonal to 4-velocity!");
    assert(isClose(innerProduct(e1, e2, g_flat), 0.0) && "FAIL: e2 is not orthogonal to e1!");
    assert(isClose(innerProduct(e0, e3, g_flat), 0.0) && "FAIL: e3 is not orthogonal to 4-velocity!");

    assert(isClose(e1[0], 0.75) && isClose(e1[1], 1.25) && "FAIL: e1 (forward) does not match exact analytical Lorentz boost!");
    assert(isClose(e1[2], 0.0) && isClose(e1[3], 0.0) && "FAIL: e1 has ghost components in y/z!");

    assert(isClose(e2[0], 0.0) && isClose(e2[1], 0.0) && isClose(e2[2], 1.0) && "FAIL: e2 (right) was incorrectly boosted!");
    assert(isClose(e3[0], 0.0) && isClose(e3[3], 1.0) && "FAIL: e3 (up) was incorrectly boosted!");

    // Test 2: 4-Acceleration Orthogonality
    Eigen::Vector3d localThrust(10.0, -5.5, 3.2);
    Eigen::Vector4d globalAccel = speedFrame.getGlobalFourAcceleration(localThrust);
    // A_mu * U^mu = 0
    double accelVelocityDot = innerProduct(globalAccel, e0, g_flat);
    
    std::cout << "4-Acceleration * 4-Velocity Dot Product: " << accelVelocityDot << "\n";
    assert(isClose(accelVelocityDot, 0.0) && "CRITICAL FAIL: Acceleration added energy to the rest mass! A_mu * U^mu != 0");

    double expectedAccelSq = localThrust.squaredNorm(); // 10^2 + (-5.5)^2 + 3.2^2 = 140.49
    double actualAccelSq = innerProduct(globalAccel, globalAccel, g_flat);
    
    std::cout << "Expected Proper Accel Sq: " << expectedAccelSq << " | Computed: " << actualAccelSq << "\n";
    assert(isClose(actualAccelSq, expectedAccelSq, 1e-5) && "FAIL: The 4-Acceleration magnitude drifted! Tetrad mapping is stretching the vector.");

    Eigen::Vector3d degenFwd(1.0, 0.0, 0.0);
    Eigen::Vector3d degenRgt(1.0, 0.0, 0.0); 
    
    LocalFrame degenFrame(g_flat, highSpeedU, degenFwd, degenRgt, up);
    
    Eigen::Vector4d safeE2 = degenFrame.getE2();
    
    assert(!std::isnan(safeE2[0]) && "FAIL: Degenerate Gram-Schmidt projection generated NaN!");
    assert(isClose(safeE2.norm(), 0.0) && "FAIL: Degenerate vector should be cleanly zeroed out.");

    std::cout << "All Tetrad & Local Frame tests passed! Frame mapping is mathematically bulletproof.\n";
    return 0;
}

