// tests/test_integrator.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "RelSpacetime/WorldlineIntegrator.h"
#include "RelSpacetime/MetricTensor.h"

using namespace RelSpacetime;

bool isClose(double a, double b, double epsilon = 1e-6) {
    return std::abs(a - b) < epsilon;
}

int main() {
    std::cout << "Running Stateless Integrator Tests...\n";
    
    MinkowskiMetric flatSpace;

    // Old truth
    EntityState startState;
    startState.position = SpacetimeVector(0.0, 0.0, 0.0, 0.0); // Origin at t=0
    startState.fourVelocity = Eigen::Vector4d(1.0, 0.0, 0.0, 0.0); // At rest (gamma = 1)
    startState.properTime = 0.0;
    startState.coordinateTime = 0.0;
    startState.timeScale = 1.0;

    EngineInputs inputs;
    inputs.localProperThrust = Eigen::Vector3d(0.0, 0.0, 0.0);
    inputs.forward = Eigen::Vector3d(1.0, 0.0, 0.0);
    inputs.right   = Eigen::Vector3d(0.0, 1.0, 0.0);
    inputs.up      = Eigen::Vector3d(0.0, 0.0, 1.0);

    double server_dt = 1.0; // 1 second server tick

    EntityState nextState = WorldlineIntegrator::computeNextState(startState, inputs, server_dt, flatSpace);

    // In flat space with zero thrust, proper time should match coordinate time exactly.
    assert(isClose(nextState.coordinateTime, 1.0) && "Coordinate time failed to advance");
    assert(isClose(nextState.properTime, 1.0) && "Proper time should equal coordinate time at rest");
    assert(isClose(nextState.position.x(), 0.0) && "Entity should not have moved spatially");
    assert(isClose(nextState.timeScale, 1.0) && "Time scale should remain 1.0 at rest");

    std::cout << "\nRunning Analytical Hyperbolic Motion Test...\n";
    
    double a = 1.0; // Constant proper acceleration. c=1 units
    inputs.localProperThrust = Eigen::Vector3d(a, 0.0, 0.0); 
    
    double tick_dt = 0.016667;
    int num_ticks = 300;
    double sim_time = tick_dt * num_ticks; 
    
    EntityState thrustState = startState;

    for(int i = 0; i < num_ticks; ++i) {
        thrustState = WorldlineIntegrator::computeNextState(thrustState, inputs, tick_dt, flatSpace);
    }
    
    // The Analytical Ground Truth
    double expected_U0 = std::sqrt(1.0 + a * a * sim_time * sim_time);
    double expected_U1 = a * sim_time;
    double expected_x = (std::sqrt(1.0 + a * a * sim_time * sim_time) - 1.0) / a;
    double expected_tau = std::asinh(a * sim_time) / a;
    double expected_timeScale = 1.0 / expected_U0;

    double epsilon = 1e-6;

    assert(isClose(thrustState.fourVelocity[0], expected_U0, epsilon) && "FAIL: U^0 (Gamma) drift detected!");
    assert(isClose(thrustState.fourVelocity[1], expected_U1, epsilon) && "FAIL: U^1 drift detected!");
    assert(isClose(thrustState.position.x(), expected_x, epsilon) && "FAIL: Spatial position diverges from hyperbolic trajectory!");
    assert(isClose(thrustState.properTime, expected_tau, epsilon) && "FAIL: Proper time integration failed!");
    assert(isClose(thrustState.timeScale, expected_timeScale, epsilon) && "FAIL: Time scale calculation is incorrect!");

    std::cout << "Test 5 - Hyperbolic Motion passed against exact analytical solutions!\n";

    // Test Speed of Light Limit and Extreme Time Dilation Drift
    std::cout << "\nRunning Extreme Relativistic Limit (Test 6)...\n";
    
    EntityState speedState = startState;
    double extreme_a = 5.0; // 5g acceleration
    inputs.localProperThrust = Eigen::Vector3d(extreme_a, 0.0, 0.0); 
    
    double extreme_dt = 0.016667;
    int extreme_steps = 6000;
    double total_sim_time = extreme_dt * extreme_steps; // 100 coordinate seconds
    
    for (int i = 0; i < extreme_steps; ++i) {
        speedState = WorldlineIntegrator::computeNextState(speedState, inputs, extreme_dt, flatSpace);
    }

    // Analytical Truth
    double expected_extreme_U0 = std::sqrt(1.0 + extreme_a * extreme_a * total_sim_time * total_sim_time);
    // v = U^1 / U^0. As t -> infinity, v -> 1.0 (c)
    double expected_spatial_vel = (extreme_a * total_sim_time) / expected_extreme_U0;
    
    // tau = (1/a) * arcsinh(a * t)
    // At a=5 and t=100, tau should be roughly 1.38 seconds, despite 100 seconds passing for the server!
    double expected_extreme_tau = std::asinh(extreme_a * total_sim_time) / extreme_a;

    double extreme_epsilon = 1e-6;

    double actual_spatial_vel = speedState.fourVelocity[1] / speedState.fourVelocity[0];

    assert(actual_spatial_vel < 1.0 && "FAIL: Entity broke the speed of light!");
    assert(isClose(actual_spatial_vel, expected_spatial_vel, extreme_epsilon) && "FAIL: Spatial velocity curve is incorrect near c.");
    
    assert(isClose(speedState.properTime, expected_extreme_tau, extreme_epsilon) && "FAIL: Proper time integration broke down at extreme relativistic speeds!");
    
    double extremeNormSq = -std::pow(speedState.fourVelocity[0], 2) + 
                            std::pow(speedState.fourVelocity[1], 2) + 
                            std::pow(speedState.fourVelocity[2], 2) + 
                            std::pow(speedState.fourVelocity[3], 2);
    assert(isClose(extremeNormSq, -1.0, extreme_epsilon) && "FAIL: 4-Velocity drifted off the mass shell at extreme speeds!");

    std::cout << "Test 6 - Extreme limit passed! Server Time: " << total_sim_time 
              << "s | Ship Proper Time: " << speedState.properTime 
              << "s | Velocity: " << actual_spatial_vel << "c\n";

    // 7. Test General Relativity: Kerr Metric Conserved Quantities
    std::cout << "\nRunning Kerr Metric Gravity Drop (Conserved Quantities) test...\n";

    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.5; // Rotating black hole
    KerrMetric kerrSpace(blackHoleMass, blackHoleSpin);

    EntityState dropState;
    // Place entity at x = 10.0 in the equatorial plane (z=0, y=0)
    dropState.position = SpacetimeVector(0.0, 10.0, 0.0, 0.0);
    // Initial u^0 for a stationary observer in curved space
    Eigen::Matrix4d initialMetric = kerrSpace.getMetricAt(dropState.position);
    double u0_initial = 1.0 / std::sqrt(-initialMetric(0, 0));
    
    dropState.fourVelocity = Eigen::Vector4d(u0_initial, 0.0, 0.0, 0.0);
    dropState.properTime = 0.0;
    dropState.coordinateTime = 0.0;
    dropState.timeScale = 1.0 / u0_initial; 

    inputs.localProperThrust = Eigen::Vector3d::Zero(); // free-fall

    // 4-velocity: u_\mu = g_{\mu\nu} u^\nu
    Eigen::Vector4d u_cov_initial = initialMetric * dropState.fourVelocity;
    
    double E_initial = -u_cov_initial[0]; // Specific Energy (Time translation symmetry)
    double Lz_initial = u_cov_initial[3]; // Specific Axial Angular Momentum (Axisymmetry)

    double kerr_dt = 1.0; 
    EntityState fallenState = WorldlineIntegrator::computeNextState(dropState, inputs, kerr_dt, kerrSpace);

    Eigen::Matrix4d fallenMetric = kerrSpace.getMetricAt(fallenState.position);
    Eigen::Vector4d u_cov_final = fallenMetric * fallenState.fourVelocity;
    
    double E_final = -u_cov_final[0];
    double Lz_final = u_cov_final[3];

    double kerr_epsilon = 1e-4; // Tolerance for RK4 numerical drift

    assert(fallenState.position.x() < 10.0 && "FAIL: Entity did not fall towards the black hole.");
    assert(fallenState.fourVelocity[1] < 0.0 && "FAIL: Entity did not accelerate inward due to gravity.");

    double fallenNormSq = fallenState.fourVelocity.transpose() * fallenMetric * fallenState.fourVelocity;
    assert(isClose(fallenNormSq, -1.0, kerr_epsilon) && "FAIL: 4-Velocity drifted off the mass shell!");

    assert(isClose(E_initial, E_final, kerr_epsilon) && "FAIL: Specific Energy (E) is not conserved! Christoffel/RK4 math is flawed.");
    assert(isClose(Lz_initial, Lz_final, kerr_epsilon) && "FAIL: Specific Axial Angular Momentum (L_z) is not conserved!");

    std::cout << "Test 7 - Conserved Quantities verified! E=" << E_final << ", Lz=" << Lz_final << " remained perfectly stable.\n";

    std::cout << "\nAll Integrator tests passed! Stateless boundary is solid.\n";
    return 0;
}

