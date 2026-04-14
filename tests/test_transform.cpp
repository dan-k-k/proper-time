// tests/test_transform.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "RelSpacetime/SpacetimeVector.h"
#include "RelSpacetime/MetricTensor.h"
#include "RelSpacetime/Transform.h"

using namespace RelSpacetime;

bool isClose(double a, double b, double epsilon = 1e-9) {
    return std::abs(a - b) < epsilon;
}

int main() {
    std::cout << "Running RelSpacetime Transform Tests...\n";

    MinkowskiMetric flatSpace;

    // Test 1: Simple 1D Boost (X-axis)
    // v = 0.6c. gamma should be 1.25.
    Eigen::Vector3d xVelocity(0.6, 0.0, 0.0); 
    Eigen::Matrix4d boostX = Transform::createBoost(xVelocity);
    SpacetimeVector eventA(1.0, 1.0, 0.0, 0.0);
    SpacetimeVector boostedA = Transform::apply(boostX, eventA);
    
    // t' = gamma * (t - vx) = 1.25 * (1.0 - 0.6(1.0)) = 0.5
    // x' = gamma * (-vt + x) = 1.25 * (-0.6(1.0) + 1.0) = 0.5
    assert(isClose(boostedA.t(), 0.5) && "FAIL: Boost X Time component incorrect");
    assert(isClose(boostedA.x(), 0.5) && "FAIL: Boost X Spatial component incorrect");

    // Test 2: Lorentz Invariance
    SpacetimeVector event1(0.0, 0.0, 0.0, 0.0);
    SpacetimeVector event2(5.0, 3.0, -2.0, 4.0); // Arbitrary event
    double original_ds2 = flatSpace.calculateIntervalSq(event1, event2);

    // Create a 3D boost (vx=0.2, vy=-0.5, vz=0.4)
    Eigen::Vector3d wildVelocity(0.2, -0.5, 0.4); 
    Eigen::Matrix4d wildBoost = Transform::createBoost(wildVelocity);

    SpacetimeVector boosted1 = Transform::apply(wildBoost, event1);
    SpacetimeVector boosted2 = Transform::apply(wildBoost, event2);
    double boosted_ds2 = flatSpace.calculateIntervalSq(boosted1, boosted2);

    assert(isClose(original_ds2, boosted_ds2) && "FAIL: Lorentz Invariance violated! ds^2 changed after 3D boost.");

    // Test 3: Speed of Light Exceptions
    bool caught1D = false;
    try { Transform::createBoost(Eigen::Vector3d(1.0, 0.0, 0.0)); } 
    catch (const std::invalid_argument&) { caught1D = true; }
    assert(caught1D && "FAIL: createBoost failed to throw on v=1.0");

    bool caught3D = false;
    try { 
        // Magnitude of this vector is sqrt(0.8^2 + 0.8^2) = sqrt(1.28) > 1.0
        Transform::createBoost(Eigen::Vector3d(0.8, 0.8, 0.0)); 
    } 
    catch (const std::invalid_argument&) { caught3D = true; }
    assert(caught3D && "FAIL: createBoost failed to throw on total v > 1.0");

    // Test 4: 1D Velocity Addition (Firing Forward) 
    // Ship moving 0.9c X, fires projectile 0.3c X. 
    // u = (0.9 + 0.3) / (1 + 0.9 * 0.3) = 1.2 / 1.27 = 0.944881889c
    Eigen::Vector3d ship_v(0.9, 0.0, 0.0);
    Eigen::Vector3d proj_u_prime(0.3, 0.0, 0.0);
    Eigen::Vector3d combined = Transform::addVelocities(ship_v, proj_u_prime);
    
    assert(combined.x() < 1.0 && "FAIL: Projectile broke the speed of light!");
    assert(isClose(combined.x(), 0.944881889) && "FAIL: 1D Velocity addition math is incorrect");
    // assert(isClose(combined.x(), 0.3) && "FAIL: 1D Velocity addition math is incorrect");

    std::cout << "All Transform tests passed! Lorentz Invariance is confirmed.\n";
    return 0;
}

