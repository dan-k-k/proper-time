// tests/test_kerr.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "RelSpacetime/SpacetimeVector.h"
#include "RelSpacetime/MetricTensor.h"
#include "RelSpacetime/WorldlineIntegrator.h" // UPDATED HEADER

using namespace RelSpacetime;

bool isClose(double a, double b, double epsilon = 1e-5) {
    return std::abs(a - b) < epsilon;
}

EntityState createRestState(double x, double y, double z, const MetricTensor& metric) {
    EntityState state;
    state.position = SpacetimeVector(0.0, x, y, z);
    
    // to stay on the mass shell, u_mu u^mu = -1
    Eigen::Matrix4d g = metric.getMetricAt(state.position);
    double u0 = 1.0 / std::sqrt(-g(0, 0));
    
    state.fourVelocity = Eigen::Vector4d(u0, 0.0, 0.0, 0.0);
    state.properTime = 0.0;
    state.coordinateTime = 0.0;
    state.timeScale = 1.0 / u0; 
    return state;
}

int main() {
    std::cout << "Running Kerr Metric & Geodesic Integrator Tests...\n";

    double server_dt = 0.1; 

    // zero-thrust inputs
    EngineInputs zeroInputs;
    zeroInputs.localProperThrust = Eigen::Vector3d::Zero();
    zeroInputs.forward = Eigen::Vector3d(1.0, 0.0, 0.0);
    zeroInputs.right   = Eigen::Vector3d(0.0, 1.0, 0.0);
    zeroInputs.up      = Eigen::Vector3d(0.0, 0.0, 1.0);

    // Test 1: Radial Free-Fall (Schwarzschild-like, a = 0)
    KerrMetric staticBlackHole(1.0, 0.0); // Mass = 1, Spin = 0
    EntityState fallingObj = createRestState(10.0, 0.0, 0.0, staticBlackHole);

    fallingObj = WorldlineIntegrator::computeNextState(fallingObj, zeroInputs, server_dt, staticBlackHole);
    
    std::cout << "Test 1 - Dropped from x=10.0, new x=" << fallingObj.position.x() << "\n";
    
    assert(fallingObj.position.x() < 10.0 && "FAIL: Object did not fall towards the black hole.");
    assert(isClose(fallingObj.position.y(), 0.0) && "FAIL: Object moved sideways without black hole spin.");
    assert(isClose(fallingObj.position.z(), 0.0) && "FAIL: Object moved in z axis.");

    // Test 2: Exact Gravitational Time Dilation (Stationary Approximation)
    double test2_dt = 1e-6; 
    EntityState deepSpaceObj = createRestState(1000.0, 0.0, 0.0, staticBlackHole);
    EntityState deepGravityObj = createRestState(2.5, 0.0, 0.0, staticBlackHole);

    deepSpaceObj = WorldlineIntegrator::computeNextState(deepSpaceObj, zeroInputs, test2_dt, staticBlackHole);
    deepGravityObj = WorldlineIntegrator::computeNextState(deepGravityObj, zeroInputs, test2_dt, staticBlackHole);

    // Analytical truth
    Eigen::Matrix4d g_deepSpace = staticBlackHole.getMetricAt(SpacetimeVector(0.0, 1000.0, 0.0, 0.0));
    Eigen::Matrix4d g_deepGravity = staticBlackHole.getMetricAt(SpacetimeVector(0.0, 2.5, 0.0, 0.0));

    double expected_space_tau = test2_dt * std::sqrt(-g_deepSpace(0,0));
    double expected_gravity_tau = test2_dt * std::sqrt(-g_deepGravity(0,0));

    std::cout << "Test 2 - Deep Space dTau: " << deepSpaceObj.properTime 
            << " | Expected: " << expected_space_tau << "\n";
    std::cout << "         Near Horizon dTau: " << deepGravityObj.properTime 
            << " | Expected: " << expected_gravity_tau << "\n";

    assert(isClose(deepSpaceObj.properTime, expected_space_tau, 1e-8) && "FAIL: Deep space proper time does not match.");
    assert(isClose(deepGravityObj.properTime, expected_gravity_tau, 1e-8) && "FAIL: Near-horizon proper time does not match.");

    // Test 3: Frame-Dragging & Energy Conservation (Spinning Black Hole, a > 0)
    KerrMetric spinningBlackHole(1.0, 0.9); // Mass = 1, Spin = 0.9
    EntityState draggedObj = createRestState(5.0, 0.0, 0.0, spinningBlackHole);

    // Initial energy
    Eigen::Matrix4d g_initial = spinningBlackHole.getMetricAt(draggedObj.position);
    Eigen::Vector4d u_cov_initial = g_initial * draggedObj.fourVelocity;
    double E_initial = -u_cov_initial[0]; // Specific energy E = -u_t

    // Update over a few frames
    for(int i = 0; i < 5; ++i) {
        draggedObj = WorldlineIntegrator::computeNextState(draggedObj, zeroInputs, server_dt, spinningBlackHole);
    }

    // Final Energy & Mass Shell 
    Eigen::Matrix4d g_final = spinningBlackHole.getMetricAt(draggedObj.position);
    Eigen::Vector4d u_cov_final = g_final * draggedObj.fourVelocity;
    double E_final = -u_cov_final[0];

    double norm_sq = draggedObj.fourVelocity.transpose() * g_final * draggedObj.fourVelocity;

    std::cout << "Test 3 - Dropped from x=5.0 near spinning BH.\n";
    std::cout << "         New Position: x=" << draggedObj.position.x() << ", y=" << draggedObj.position.y() << "\n";
    std::cout << "         Initial E: " << E_initial << " | Final E: " << E_final << "\n";

    assert(std::abs(draggedObj.position.y()) > 1e-7 && "FAIL: Frame-dragging did not alter the Y position.");
    assert(isClose(E_initial, E_final, 1e-5) && "FAIL: Specific energy E was not conserved during frame-dragging!");
    assert(isClose(norm_sq, -1.0, 1e-5) && "FAIL: 4-Velocity drifted off the mass shell during frame-dragging!");

    std::cout << "\nAll Kerr Metric tests passed! General Relativity is fully operational.\n";
    return 0;
}

