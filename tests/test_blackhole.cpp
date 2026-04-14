// tests/test_blackhole.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "RelSpacetime/VisualFX.h"
#include "RelSpacetime/MetricTensor.h"
#include "RelSpacetime/LocalFrame.h"

using namespace RelSpacetime;

bool isClose(double a, double b, double epsilon = 1e-4) {
    return std::abs(a - b) < epsilon;
}

// Helper to construct a static observer locked to the mass shell
EntityState makeStaticObserver(double r, const MetricTensor& metric) {
    EntityState state;
    state.position = SpacetimeVector(0.0, r, 0.0, 0.0); // Placed on the X-axis
    Eigen::Matrix4d g = metric.getMetricAt(state.position);
    
    // U^t = 1 / sqrt(-g_00) for a static observer
    double Ut = 1.0 / std::sqrt(-g(0, 0));
    state.fourVelocity = Eigen::Vector4d(Ut, 0.0, 0.0, 0.0);
    return state;
}

// Synge's analytical formula for the Schwarzschild shadow angular radius
double getAnalyticalShadowRadius(double r, double M) {
    if (r <= 2.0 * M) return M_PI; // Inside horizon fallback
    
    double bc = 3.0 * std::sqrt(3.0) * M; // Critical impact parameter
    double sin_alpha = (bc / r) * std::sqrt(1.0 - (2.0 * M) / r);
    
    if (r < 3.0 * M) {
        // Inside the photon sphere, the shadow covers more than half the sky
        return M_PI - std::asin(sin_alpha);
    }
    return std::asin(sin_alpha);
}

int main() {
    std::cout << "Running Black Hole Visuals Tests (Schwarzschild Limit)...\n";

    KerrMetric staticBH(1.0, 0.0); // M=1, a=0
    VisualFX visFX(1.0);
    double M = staticBH.getMass();

    // Test 1: Weak Field Limit (Far Observer)
    EntityState obsFar = makeStaticObserver(100.0, staticBH);
    // Orient the engine to look directly at the origin (-X direction)
    LocalFrame frameFar(staticBH.getMetricAt(obsFar.position), obsFar.fourVelocity, 
                        Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
    
    BlackHoleRenderState stateFar = visFX.computeBlackHoleVisuals(obsFar, frameFar, staticBH);
    double expectedFar = getAnalyticalShadowRadius(100.0, M);
    
    std::cout << "Test 1 - Far Observer (r=100) Shadow Radius: " << stateFar.shadowAngularRadius 
              << " (Expected: " << expectedFar << ")\n";
              
    assert(isClose(stateFar.shadowAngularRadius, expectedFar) && "FAIL: Weak field shadow radius diverges from analytical truth.");
    assert(isClose(stateFar.apparentDirection.x(), 1.0) && "FAIL: Fictitious principal ray should map directly to the engine's forward axis.");

    // Test 2: Strong Field (Near Photon Sphere)
    EntityState obsStrong = makeStaticObserver(6.0, staticBH);
    LocalFrame frameStrong(staticBH.getMetricAt(obsStrong.position), obsStrong.fourVelocity, 
                           Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
    
    BlackHoleRenderState stateStrong = visFX.computeBlackHoleVisuals(obsStrong, frameStrong, staticBH);
    double expectedStrong = getAnalyticalShadowRadius(6.0, M);
    
    std::cout << "Test 2 - Strong Field (r=6) Shadow Radius: " << stateStrong.shadowAngularRadius 
              << " (Expected: " << expectedStrong << ")\n";
              
    assert(isClose(stateStrong.shadowAngularRadius, expectedStrong) && "FAIL: Strong field shadow radius differs from analytical truth.");

    // Test 3: Exactly on the Photon Sphere
    EntityState obsSphere = makeStaticObserver(3.0, staticBH);
    LocalFrame frameSphere(staticBH.getMetricAt(obsSphere.position), obsSphere.fourVelocity, 
                           Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
    
    BlackHoleRenderState stateSphere = visFX.computeBlackHoleVisuals(obsSphere, frameSphere, staticBH);
    
    std::cout << "Test 3 - Photon Sphere (r=3) Shadow Radius: " << stateSphere.shadowAngularRadius 
              << " (Expected: " << M_PI / 2.0 << ")\n";
              
    assert(isClose(stateSphere.shadowAngularRadius, M_PI / 2.0) && "FAIL: At the photon sphere, the shadow MUST cover exactly half the sky (pi/2).");

    // Test 4: Interior of Photon Sphere (Sky Inversion)
    EntityState obsInside = makeStaticObserver(2.5, staticBH);
    LocalFrame frameInside(staticBH.getMetricAt(obsInside.position), obsInside.fourVelocity, 
                           Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
    
    BlackHoleRenderState stateInside = visFX.computeBlackHoleVisuals(obsInside, frameInside, staticBH);
    double expectedInside = getAnalyticalShadowRadius(2.5, M);
    
    std::cout << "Test 4 - Interior (r=2.5) Shadow Radius: " << stateInside.shadowAngularRadius 
              << " (Expected: " << expectedInside << ")\n";
              
    assert(stateInside.shadowAngularRadius > M_PI / 2.0 && "FAIL: Inside the photon sphere, the shadow should envelope the observer (> pi/2).");
    assert(isClose(stateInside.shadowAngularRadius, expectedInside) && "FAIL: Sigma-flip math failed for interior radial bounds.");

    // Test 5: Relativistic Aberration (Falling Inwards)
    std::cout << "\nRunning Relativistic Aberration Test...\n";
    
    // We take the static observer at r=6, and give them a local inward boost.
    // Boost towards the Black Hole (-X direction, which is +E1 in the frame).
    double v_fall = 0.6; // 0.6c
    double gamma = 1.0 / std::sqrt(1.0 - v_fall * v_fall); // 1.25
    
    EntityState obsFalling = obsStrong; 
    // Construct the boosted 4-velocity using the static observer's tetrad
    obsFalling.fourVelocity = gamma * frameStrong.getE0() + (gamma * v_fall) * frameStrong.getE1();
    
    LocalFrame frameFalling(staticBH.getMetricAt(obsFalling.position), obsFalling.fourVelocity, 
                            Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
                            
    BlackHoleRenderState stateFalling = visFX.computeBlackHoleVisuals(obsFalling, frameFalling, staticBH);
    
    // Relativistic Aberration Formula: cos(alpha') = (cos(alpha) + v) / (1 + v*cos(alpha))
    // Note: Light arrives from the front, so we use +v for moving towards the source.
    double expected_cos_aberrated = (std::cos(expectedStrong) + v_fall) / (1.0 + v_fall * std::cos(expectedStrong));
    double expected_aberrated_rad = std::acos(expected_cos_aberrated);

    std::cout << "Static Angle:    " << expectedStrong << "\n";
    std::cout << "Falling Angle:   " << stateFalling.shadowAngularRadius << "\n";
    std::cout << "Expected Angle:  " << expected_aberrated_rad << "\n";

    // Counter-intuitive GR fact: Falling *towards* a black hole causes its apparent shadow to SHRINK due to aberration.
    assert(stateFalling.shadowAngularRadius < stateStrong.shadowAngularRadius && "FAIL: Moving towards the BH should shrink its apparent size via aberration!");
    assert(isClose(stateFalling.shadowAngularRadius, expected_aberrated_rad) && "FAIL: Exact aberration math mapping failed.");

    // Test 6: Polar vs Equatorial Symmetry
    EntityState obsPolar;
    obsPolar.position = SpacetimeVector(0.0, 0.0, 0.0, 6.0); // Z=6 (Polar)
    Eigen::Matrix4d gPolar = staticBH.getMetricAt(obsPolar.position);
    obsPolar.fourVelocity = Eigen::Vector4d(1.0 / std::sqrt(-gPolar(0, 0)), 0.0, 0.0, 0.0);
    
    // Orient engine looking down (-Z)
    LocalFrame framePolar(gPolar, obsPolar.fourVelocity, Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0)); 
    BlackHoleRenderState statePolar = visFX.computeBlackHoleVisuals(obsPolar, framePolar, staticBH);
    
    assert(isClose(statePolar.shadowAngularRadius, expectedStrong) && "FAIL: Schwarzschild metric is spherically symmetric. Polar observer must see identical shadow size as equatorial.");

    std::cout << "\nALL BLACK HOLE VISUAL TESTS PASSED!\n";
    return 0;
}

