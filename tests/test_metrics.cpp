// tests/test_metrics.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "RelSpacetime/SpacetimeVector.h"
#include "RelSpacetime/MetricTensor.h" 
#include "RelSpacetime/RelativityTypes.h" 
#include "RelSpacetime/AnalyticalRaytracer.h"

using namespace RelSpacetime;

bool isClose(double a, double b, double epsilon = 1e-9) {
    return std::abs(a - b) < epsilon;
}

int main() {
    std::cout << "Running RelSpacetime Metrics Tests...\n";
    MinkowskiMetric flatSpace;

    // Test 1: Lightlike (Null) Interval
    SpacetimeVector emission(0.0, 0.0, 0.0, 0.0);
    SpacetimeVector absorption(1.0, 1.0, 0.0, 0.0);
    double ds2_light = flatSpace.calculateIntervalSq(emission, absorption);
    assert(isClose(ds2_light, 0.0) && "FAIL: Lightlike interval must be 0");

    // Test 2: Timelike Interval
    SpacetimeVector start_wait(0.0, 5.0, 0.0, 0.0);
    SpacetimeVector end_wait(2.0, 5.0, 0.0, 0.0);
    double ds2_time = flatSpace.calculateIntervalSq(start_wait, end_wait);
    assert(isClose(ds2_time, -4.0) && "FAIL: Timelike interval must be negative proper time squared");

    // Test 3: Spacelike Interval
    SpacetimeVector event_A(0.0, 0.0, 0.0, 0.0);
    SpacetimeVector event_B(0.0, 3.0, 0.0, 0.0);
    double ds2_space = flatSpace.calculateIntervalSq(event_A, event_B);
    assert(isClose(ds2_space, 9.0) && "FAIL: Spacelike interval must be positive proper length squared");

    std::cout << "All Metric Phase 1 tests passed! Math is sound.\n";
    
    // Test 4: Flat Space Light Travel Time (Engine Lookback)
    SpacetimeVector observer(10.0, 10.0, 0.0, 0.0); // Observer at x=10
    SpacetimeVector target(10.0, 0.0, 0.0, 0.0);    // Target at x=0
    
    NullGeodesicResult flatLookback = AnalyticalRaytracer::traceRays(observer, target, flatSpace);

    assert(flatLookback.primary.isValid && "FAIL: Flat space lookback must be valid");
    assert(isClose(flatLookback.primary.travelTime, 10.0) && "FAIL: Light travel time in flat space should equal spatial distance (c=1)");
    
    assert(isClose(flatLookback.primary.arrivalK[0], -1.0) && "FAIL: Arriving photon must point backward in time");
    assert(isClose(flatLookback.primary.arrivalK[1], -1.0) && "FAIL: Arriving photon spatial vector must point toward target");
    assert(!flatLookback.secondary.isValid && "FAIL: Flat space should not produce a secondary image");

    // Test 5: Exact Kerr Metric Tensor Evaluation (Equatorial Plane)
    std::cout << "\nRunning Exact Kerr Metric Evaluation...\n";
    
    double M = 1.0;
    double a = 0.5;
    KerrMetric kerr(M, a); 
    
    SpacetimeVector deepSpace(0.0, 100000.0, 0.0, 0.0); // g_00 will be ~ -0.99998
    SpacetimeVector nearBH(0.0, 2.5, 0.0, 0.0); // x = 2.5

    Eigen::Matrix4d g_deep = kerr.getMetricAt(deepSpace);
    Eigen::Matrix4d g_near = kerr.getMetricAt(nearBH);

    assert(isClose(g_deep(0,0), -1.0, 1e-4) && "FAIL: Kerr metric does not flatten out at infinite distance");

    double r = std::sqrt(2.5 * 2.5 - a * a); 
    double expected_g00 = -1.0 + (2.0 * M) / r; 

    assert(isClose(g_near(0,0), expected_g00, 1e-6) && "FAIL: Kerr g_00 does not match exact analytical formula!");

    // Frame-dragging check. Dragging in the XY plane creates a g_0y (or g_01/g_02) term.
    bool hasFrameDragging = std::abs(g_near(0,1)) > 1e-6 || std::abs(g_near(0,2)) > 1e-6;
    assert(hasFrameDragging && "FAIL: Spinning black hole must have non-zero shift vector (frame dragging)");

    std::cout << "Test 5 - Exact Kerr g_00 verified: " << g_near(0,0) << " == " << expected_g00 << "\n";

    // Test 6: Kerr Metric Secondary Image Travel Time
    std::cout << "\nRunning Kerr Secondary Image Tests...\n";
    
    // Observer and Target are on opposite sides of the Black Hole
    SpacetimeVector obsOpposite(0.0, 10.0, 0.0, 0.0);
    SpacetimeVector tgtOpposite(0.0, -10.0, 0.0, 0.0);
    
    NullGeodesicResult oppositeGeodesic = AnalyticalRaytracer::traceRays(obsOpposite, tgtOpposite, kerr);

    assert(oppositeGeodesic.primary.isValid && "FAIL: Primary image should be valid.");
    assert(oppositeGeodesic.secondary.isValid && "FAIL: Secondary image should survive at this distance.");
    
    std::cout << "Primary Travel Time: " << oppositeGeodesic.primary.travelTime << "\n";
    std::cout << "Secondary Travel Time: " << oppositeGeodesic.secondary.travelTime << "\n";

    assert(oppositeGeodesic.secondary.travelTime > oppositeGeodesic.primary.travelTime && 
           "FAIL: Secondary image must take longer because it travels around the back of the BH.");

    std::cout << "Test 6 - Secondary image routing and time delay verified.\n";

    std::cout << "All Metric Phase 2 tests passed! Matrices are perfectly precise.\n";
    return 0;
}

