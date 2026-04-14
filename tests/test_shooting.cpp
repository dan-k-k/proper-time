// tests/test_shooting.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "RelSpacetime/MetricTensor.h"
#include "RelSpacetime/RelativityTypes.h"
#include "RelSpacetime/AnalyticalRaytracer.h"

using namespace RelSpacetime;

bool isClose(double a, double b, double epsilon = 1e-4) {
    return std::abs(a - b) < epsilon;
}

int main() {
    std::cout << "Running Lookback Key (Shapiro Delay) & Spatial Bending Tests...\n";

    MinkowskiMetric flatSpace;
    KerrMetric kerrSpace(1.0, 0.9); // Mass = 1, Spin = 0.9

    SpacetimeVector observerSafe(0.0, 4.0, 0.0, 0.0);
    SpacetimeVector target(0.0, 0.0, 4.0, 0.0); 

    // Test 1: Flat Space Baseline
    NullGeodesicResult flatRes = AnalyticalRaytracer::traceRays(observerSafe, target, flatSpace);
    
    double exactDist = std::sqrt(4.0*4.0 + 4.0*4.0); // ~5.65685
    std::cout << "Flat Space Travel Time: " << flatRes.primary.travelTime << " (Expected: " << exactDist << ")\n";
    
    assert(flatRes.primary.isValid && "FAIL: Flat space should always be valid.");
    assert(std::abs(flatRes.primary.travelTime - exactDist) < 1e-4 && "FAIL: Flat space travel time must equal Euclidean distance.");

    // Test 2: Analytical Kerr Space (Shapiro Delay)
    NullGeodesicResult kerrRes = AnalyticalRaytracer::traceRays(observerSafe, target, kerrSpace);

    assert(kerrRes.primary.isValid && "FAIL: Kerr space Lookback should be valid.");

    double r_obs = std::sqrt(4.0*4.0); // 4.0
    double r_tgt = std::sqrt(4.0*4.0); // 4.0
    double d = exactDist; // ~5.65685
    double M = 1.0;

    double expected_shapiro_time = d + 2.0 * M * std::log((r_obs + r_tgt + d) / (r_obs + r_tgt - d));

    std::cout << "Flat Space Time: " << d << "\n";
    std::cout << "Theoretical GR Shapiro Time: " << expected_shapiro_time << "\n";
    std::cout << "Engine Computed Time: " << kerrRes.primary.travelTime << "\n";

    assert(kerrRes.primary.travelTime > d && "FAIL: Shapiro Time Delay did not slow down the light!");
    assert(isClose(kerrRes.primary.travelTime, expected_shapiro_time, 0.1) && "FAIL: Shapiro delay diverges heavily");
    
    // Test 3: Collinear Test (Einstein Ring / Black Hole Transit)    
    std::cout << "\nRunning Collinear Test...\n";
    SpacetimeVector obsCollinear(0.0, 4.0, 0.0, 0.0);
    SpacetimeVector tgtCollinear(0.0, -4.0, 0.0, 0.0);
    
    NullGeodesicResult collinearRes = AnalyticalRaytracer::traceRays(obsCollinear, tgtCollinear, kerrSpace);

    assert(collinearRes.primary.isValid && "FAIL: Collinear transit incorrectly flagged as invalid! Light should survive to form an Einstein Ring.");
    std::cout << "Collinear Tangential Stretch: " << collinearRes.primary.tangentialStretch << "\n";
    std::cout << "Collinear Radial Squish: " << collinearRes.primary.radialSquish << "\n";
    assert(isClose(collinearRes.primary.tangentialStretch, 15.0, 1e-5) && "FAIL: Collinear transit did not clamp stretch to 15.0.");
    assert(isClose(collinearRes.primary.radialSquish, 0.05, 1e-5) && "FAIL: Perfect ring squish should be 0.05.");
    assert(isClose(collinearRes.primary.ringFactor, 1.0, 1e-5) && "FAIL: Perfect ring factor should be 1.0.");

    // Test 5: Weak-Field Lensing Magnification (Jacobian Derivation)
    std::cout << "\nRunning Thin-Lens Jacobian (Magnification) Test...\n";

    SpacetimeVector obsLens(0.0, 25.0, -100.0, 0.0);
    SpacetimeVector tgtLens(0.0, 25.0, 100.0, 0.0);

    NullGeodesicResult lensRes = AnalyticalRaytracer::traceRays(obsLens, tgtLens, kerrSpace);

    Eigen::Vector3d obs_vec(obsLens.x(), obsLens.y(), obsLens.z());
    Eigen::Vector3d tgt_vec(tgtLens.x(), tgtLens.y(), tgtLens.z());

    double D_l = std::sqrt(25.0*25.0 + 100.0*100.0); 
    double D_ls = std::sqrt(25.0*25.0 + 100.0*100.0); 
    double D_s = D_l + D_ls;
    
    Eigen::Vector3d dirToBH = (-obs_vec).normalized();
    Eigen::Vector3d viewDir = (tgt_vec - obs_vec).normalized();
    double beta = std::acos(dirToBH.dot(viewDir)); 
    
    double theta_E_sq = (4.0 * kerrSpace.getMass() * D_ls) / (D_l * D_s);
    double expected_apparent_angle = (beta + std::sqrt(beta * beta + 4.0 * theta_E_sq)) / 2.0;
    
    double expected_stretch = expected_apparent_angle / beta;
    double expected_squish = 1.0 / std::sqrt(1.0 + (theta_E_sq / (beta * beta)));
    double expected_ring = std::sqrt(theta_E_sq) / std::sqrt(beta * beta + theta_E_sq);

    std::cout << "Analytical Flat Angle (beta): " << beta << " rad\n";
    std::cout << "Expected Stretch: " << expected_stretch << " | Engine: " << lensRes.primary.tangentialStretch << "\n";
    std::cout << "Expected Squish:  " << expected_squish  << " | Engine: " << lensRes.primary.radialSquish << "\n";
    std::cout << "Expected Ring:    " << expected_ring    << " | Engine: " << lensRes.primary.ringFactor << "\n";

    assert(lensRes.primary.isValid && "FAIL: Light should easily reach the observer.");
    assert(isClose(lensRes.primary.tangentialStretch, expected_stretch, 0.01) && "FAIL: Tangential stretch deviates from analytical baseline!");
    assert(isClose(lensRes.primary.radialSquish, expected_squish, 0.01) && "FAIL: Radial squish deviates from analytical baseline!");
    assert(isClose(lensRes.primary.ringFactor, expected_ring, 0.01) && "FAIL: Ring factor deviates from analytical baseline!");
    
    Eigen::Vector3d expectedAxis = dirToBH.cross(viewDir).normalized();
    assert(isClose(std::abs(lensRes.primary.arcTangentAxis.dot(expectedAxis)), 1.0, 1e-4) && "FAIL: Arc tangent rotation axis is mathematically incorrect!");

    assert(lensRes.primary.arrivalK[1] > 0.0 && "FAIL: The light bent the wrong way!");
    
    // Test 6: Strong Field Beloborodov
    std::cout << "\nRunning Regression Test: Strong-Field Drop-In...\n";
    
    SpacetimeVector obsDropIn(0.0, 50.0, 0.0, 0.0);
    SpacetimeVector tgtDropIn(0.0, 0.0, 2.1, 0.0);
    
    NullGeodesicResult dropInRes = AnalyticalRaytracer::traceRays(obsDropIn, tgtDropIn, kerrSpace);

    std::cout << "Drop-In Stretch: " << dropInRes.primary.tangentialStretch << "\n";
    
    assert(dropInRes.primary.tangentialStretch > 1.05 && "FAIL: Tangential stretch collapsed! Beloborodov is failing on deep-well targets.");
    assert(dropInRes.primary.radialSquish < 0.95 && "FAIL: Deep well target failed to radially squish.");

    // Test 7: Continuous Curvature
    std::cout << "\nRunning Regression Test: Non-Sandwiched Curvature...\n";
    
    SpacetimeVector obsSameSide(0.0, 10.0, 5.0, 0.0);
    SpacetimeVector tgtSameSide(0.0, 20.0, 5.0, 0.0);
    
    NullGeodesicResult sameSideRes = AnalyticalRaytracer::traceRays(obsSameSide, tgtSameSide, kerrSpace);

    assert(sameSideRes.primary.tangentialStretch > 1.0 && "FAIL: Light bending was prematurely culled! The '7-second jump' bug has returned.");
    assert(sameSideRes.primary.radialSquish < 1.0 && "FAIL: Same-side target should still experience slight squish.");

    // Test 8: Secondary Image Arrival Angle (Spatial Inversion)
    std::cout << "\nRunning Secondary Image Optics & Arrival Test...\n";
    
    SpacetimeVector obsSec(0.0, 0.0, 0.0, -100.0);
    SpacetimeVector tgtSec(0.0, 5.0, 0.0, 100.0);
    
    NullGeodesicResult secBendRes = AnalyticalRaytracer::traceRays(obsSec, tgtSec, kerrSpace);

    assert(secBendRes.secondary.tangentialStretch > 1.0 && "FAIL: Secondary image should undergo tangential stretch.");
    assert(secBendRes.secondary.radialSquish < 1.0 && "FAIL: Secondary image should undergo severe radial squish.");
    
    std::cout << "Primary Image Apparent X:   " << secBendRes.primary.arrivalK[1] << "\n";
    std::cout << "Secondary Image Apparent X: " << secBendRes.secondary.arrivalK[1] << "\n";

    assert(secBendRes.primary.arrivalK[1] > 0.0 && "FAIL: Primary image should appear on the same side as the target.");
    assert(secBendRes.secondary.arrivalK[1] < 0.0 && "FAIL: Secondary image failed to invert! It should arrive from the opposite side of the Black Hole.");

    std::cout << "\nALL TESTS PASSED! Lookback Key & Visuals API are ready.\n";
    return 0;
}

