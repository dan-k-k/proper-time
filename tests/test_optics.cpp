// tests/test_optics.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "RelSpacetime/VisualFX.h"
#include "RelSpacetime/MetricTensor.h"
#include "RelSpacetime/AnalyticalRaytracer.h"

using namespace RelSpacetime;

bool isClose(double a, double b, double epsilon = 1e-4) {
    return std::abs(a - b) < epsilon;
}

int main() {
    std::cout << "Running Optics & Visuals Tests...\n";
    MinkowskiMetric flatSpace;
    VisualFX optics(1.0); // spectral index = 1.0

    // Observer at rest at the origin
    EntityState obsState;
    obsState.position = SpacetimeVector(0.0, 0.0, 0.0, 0.0);
    obsState.fourVelocity = Eigen::Vector4d(1.0, 0.0, 0.0, 0.0);
    
    LocalFrame obsFrame(flatSpace.getMetricAt(obsState.position), obsState.fourVelocity, 
                        Eigen::Vector3d(1,0,0), Eigen::Vector3d(0,1,0), Eigen::Vector3d(0,0,1));

    // Target directly in front (+X axis) moving AWAY at 0.6c. gamma = 1.25. U^mu = (1.25, 0.75, 0, 0)
    EntityState tgtAway;
    tgtAway.position = SpacetimeVector(-10.0, 10.0, 0.0, 0.0); // Emitted 10 seconds ago at x=10
    tgtAway.fourVelocity = Eigen::Vector4d(1.25, 0.75, 0.0, 0.0);

    NullGeodesicResult mockGeodesic;
    mockGeodesic.primary.isValid = true;
    mockGeodesic.primary.travelTime = 10.0; 
    mockGeodesic.primary.arrivalK = Eigen::Vector4d(-1.0, 1.0, 0.0, 0.0);
    mockGeodesic.primary.tangentialStretch = 1.0;
    mockGeodesic.primary.radialSquish = 1.0;
    mockGeodesic.primary.ringFactor = 0.0;
    mockGeodesic.primary.arcTangentAxis = Eigen::Vector3d::Zero();
    mockGeodesic.secondary.isValid = false; 

    // Test 1: Redshift (Moving Away) & Flat Space Optics
    VisualRenderState visAway = optics.computeVisuals(obsState, obsFrame, tgtAway, tgtAway, mockGeodesic, flatSpace);
    
    std::cout << "Target moving away Doppler: " << visAway.primary.dopplerShift << "\n";
    assert(visAway.primary.dopplerShift < 1.0 && "FAIL: Object moving away should be Redshifted (< 1.0)");
    assert(isClose(visAway.primary.dopplerShift, 0.5) && "FAIL: Exact Redshift math is wrong.");
    
    // Test 2: Blueshift (Moving Towards)
    EntityState tgtTowards = tgtAway;
    tgtTowards.fourVelocity = Eigen::Vector4d(1.25, -0.75, 0.0, 0.0);
    VisualRenderState visTowards = optics.computeVisuals(obsState, obsFrame, tgtTowards, tgtTowards, mockGeodesic, flatSpace);

    std::cout << "Target moving towards Doppler: " << visTowards.primary.dopplerShift << "\n";
    assert(visTowards.primary.dopplerShift > 1.0 && "FAIL: Object moving towards should be Blueshifted (> 1.0)");
    assert(isClose(visTowards.primary.dopplerShift, 2.0) && "FAIL: Exact Blueshift math is wrong.");

    // Test 3: Exact Gravitational Redshift in Curved Space
    KerrMetric staticBH(1.0, 0.0); // Mass = 1, Spin = 0 (Schwarzschild geometry)
    
    EntityState obsDeep;
    obsDeep.position = SpacetimeVector(0.0, 100.0, 0.0, 0.0); // r = 100
    obsDeep.fourVelocity = Eigen::Vector4d(1.0, 0.0, 0.0, 0.0); 
    
    EntityState tgtGravity;
    tgtGravity.position = SpacetimeVector(-50.0, 3.0, 0.0, 0.0); // r = 3
    tgtGravity.fourVelocity = Eigen::Vector4d(1.0, 0.0, 0.0, 0.0); 

    Eigen::Matrix4d g_obs = staticBH.getMetricAt(obsDeep.position);
    Eigen::Matrix4d g_emit = staticBH.getMetricAt(tgtGravity.position);
    
    // Lock onto the mass shell
    obsDeep.fourVelocity[0] = 1.0 / std::sqrt(-g_obs(0,0));
    tgtGravity.fourVelocity[0] = 1.0 / std::sqrt(-g_emit(0,0));

    LocalFrame obsDeepFrame(g_obs, obsDeep.fourVelocity, 
                            Eigen::Vector3d(1,0,0), Eigen::Vector3d(0,1,0), Eigen::Vector3d(0,0,1));

    NullGeodesicResult mockGravityGeodesic;
    mockGravityGeodesic.primary.isValid = true;
    mockGravityGeodesic.primary.travelTime = 97.0; // Arbitrary
    mockGravityGeodesic.primary.arrivalK = Eigen::Vector4d(-1.0, 1.0, 0.0, 0.0);
    mockGravityGeodesic.primary.tangentialStretch = 1.0;
    mockGravityGeodesic.primary.radialSquish = 1.0;
    mockGravityGeodesic.primary.ringFactor = 0.0;
    mockGravityGeodesic.primary.arcTangentAxis = Eigen::Vector3d::Zero();
    mockGravityGeodesic.secondary.isValid = false; 

    VisualRenderState visGravity = optics.computeVisuals(obsDeep, obsDeepFrame, tgtGravity, tgtGravity, mockGravityGeodesic, staticBH);

    // Analytical truth. g_00 = -1 + 2M/r
    double expected_g00_emit = -1.0 + (2.0 * 1.0) / 3.0;   // -0.333333...
    double expected_g00_obs  = -1.0 + (2.0 * 1.0) / 100.0; // -0.98
    double expected_grav_doppler = std::sqrt(expected_g00_emit / expected_g00_obs);

    std::cout << "Gravitational Doppler (Target at r=3, Obs at r=100): " << visGravity.primary.dopplerShift << "\n";
    std::cout << "Expected Analytical Doppler: " << expected_grav_doppler << "\n";

    assert(visGravity.primary.dopplerShift < 1.0 && "GOTCHA C FAIL: Gravitational redshift ignored.");
    assert(isClose(visGravity.primary.dopplerShift, expected_grav_doppler, 1e-6) && "FAIL: Exact gravitational redshift calculation is incorrect!");

    // Test 4: Secondary Image Independent Visuals
    std::cout << "\nRunning Secondary Image Optics Tests...\n";
    
    // 2nd image took 5 seconds longer to arrive
    NullGeodesicResult dualGeodesic;
    dualGeodesic.primary.isValid = true;
    dualGeodesic.primary.travelTime = 10.0;
    dualGeodesic.primary.arrivalK = Eigen::Vector4d(-1.0, 1.0, 0.0, 0.0);
    dualGeodesic.primary.tangentialStretch = 2.0;
    dualGeodesic.primary.radialSquish = 0.8;
    dualGeodesic.primary.ringFactor = 0.1;
    dualGeodesic.primary.arcTangentAxis = Eigen::Vector3d(0, 1, 0);
    
    dualGeodesic.secondary.isValid = true;
    dualGeodesic.secondary.travelTime = 15.0; // Took longer!
    dualGeodesic.secondary.arrivalK = Eigen::Vector4d(-1.0, -1.0, 0.0, 0.0); // Arriving from the other side
    dualGeodesic.secondary.tangentialStretch = 0.5;
    dualGeodesic.secondary.radialSquish = 0.2;
    dualGeodesic.secondary.ringFactor = 0.9;
    dualGeodesic.secondary.arcTangentAxis = Eigen::Vector3d(0, 1, 0);

    // The primary past state (10 seconds ago, moving away slowly)
    EntityState primaryPast;
    primaryPast.fourVelocity = Eigen::Vector4d(1.05, 0.3, 0.0, 0.0);
    // The secondary past state (15 seconds ago, before it accelerated, moving away fast)
    EntityState secondaryPast;
    secondaryPast.fourVelocity = Eigen::Vector4d(1.25, 0.75, 0.0, 0.0);

    VisualRenderState dualVis = optics.computeVisuals(obsState, obsFrame, primaryPast, secondaryPast, dualGeodesic, flatSpace);

    assert(dualVis.secondary.isVisible && "FAIL: Secondary image should be marked visible.");
    assert(isClose(dualVis.secondary.tangentialStretch, 0.5) && "FAIL: Secondary tangential stretch not mapped.");
    assert(isClose(dualVis.secondary.radialSquish, 0.2) && "FAIL: Secondary radial squish not mapped.");
    assert(isClose(dualVis.secondary.ringFactor, 0.9) && "FAIL: Secondary ring factor not mapped.");
    
    std::cout << "Primary Doppler (recent past): " << dualVis.primary.dopplerShift << "\n";
    std::cout << "Secondary Doppler (deeper past): " << dualVis.secondary.dopplerShift << "\n";

    assert(dualVis.primary.dopplerShift != dualVis.secondary.dopplerShift && 
           "FAIL: Primary and Secondary images must process their distinct historical velocities independently!");

    std::cout << "Test 4 - Secondary image mapping verified.\n";

    // Test 5: Analytical Light Bending - Quantitative Verification
    std::cout << "\nRunning Analytical Light Bending Tests...\n";

    KerrMetric testBH(1.0, 0.0);

    // 5A: Einstein ring
    SpacetimeVector ringObs(0.0, 0.0, 0.0, 100.0);
    SpacetimeVector ringTgt(0.0, 0.0, 0.0, -100.0);
    
    NullGeodesicResult ringResult = AnalyticalRaytracer::traceRays(ringObs, ringTgt, testBH);
    
    assert(isClose(ringResult.primary.tangentialStretch, 15.0) && "FAIL: Perfect ring primary should clamp stretch to 15.0");
    assert(isClose(ringResult.primary.radialSquish, 0.05) && "FAIL: Perfect ring primary squish should be razor thin (0.05)");
    assert(isClose(ringResult.primary.ringFactor, 1.0) && "FAIL: Perfect ring primary should set ringFactor to 1.0");
    
    assert(isClose(ringResult.secondary.tangentialStretch, 15.0) && "FAIL: Perfect ring secondary should clamp stretch to 15.0");
    assert(isClose(ringResult.secondary.radialSquish, 0.05) && "FAIL: Perfect ring secondary squish should be razor thin (0.05)");
    assert(isClose(ringResult.secondary.ringFactor, 1.0) && "FAIL: Perfect ring secondary should set ringFactor to 1.0");
    
    std::cout << "Test 5A - Perfect Einstein Ring (Both Images) verified.\n";

    // 5B: Weak field deflection
    SpacetimeVector weakObs(0.0, 0.0, 0.0, 100.0);
    SpacetimeVector weakTgt(0.0, 10.0, 0.0, -100.0);

    NullGeodesicResult weakResult = AnalyticalRaytracer::traceRays(weakObs, weakTgt, testBH);

    assert(weakResult.primary.tangentialStretch > 1.0 && "FAIL: Lensed target should stretch tangentially.");
    assert(isClose(weakResult.primary.radialSquish, 0.333, 0.05) && "FAIL: Radial squish calculation differs from analytical expectation.");
    assert(isClose(weakResult.primary.ringFactor, 0.942, 0.05) && "FAIL: Ring factor calculation differs from analytical expectation.");
    
    assert(weakResult.secondary.tangentialStretch > 1.0 && "FAIL: Secondary image should be stretched.");
    assert(isClose(weakResult.secondary.radialSquish, 0.333, 0.05) && "FAIL: Secondary radial squish calculation differs from analytical expectation.");

    std::cout << "Test 5B - Mild Deflection (Primary Weak Field, Secondary Strong Field) verified.\n";

    // 5C: Looking away, no deflection
    SpacetimeVector awayObs(0.0, 0.0, 0.0, 10.0);
    SpacetimeVector awayTgt(0.0, 0.0, 0.0, 20.0);

    NullGeodesicResult awayResult = AnalyticalRaytracer::traceRays(awayObs, awayTgt, testBH);

    assert(isClose(awayResult.primary.tangentialStretch, 1.0) && "FAIL: Unlensed target should have stretch 1.0");
    assert(isClose(awayResult.primary.radialSquish, 1.0) && "FAIL: Unlensed target should have squish 1.0");
    assert(isClose(awayResult.primary.ringFactor, 0.0) && "FAIL: Unlensed target should have ringFactor 0.0");

    assert(awayResult.secondary.radialSquish <= 0.05 && "FAIL: Secondary image looking away should be extremely clamped.");

    std::cout << "Test 5C - Free space (No Deflection) verified.\n";

    std::cout << "\nOptics tests passed! Visual pipeline is ready for GPU Shaders.\n";
    return 0;
}

