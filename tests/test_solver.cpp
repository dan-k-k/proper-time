// tests/test_solver.cpp
#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include "RelSpacetime/RelativityTypes.h"
#include "RelSpacetime/MetricTensor.h"
#include "RelSpacetime/SpacetimeInterpolator.h"
#include "RelSpacetime/AnalyticalRaytracer.h"

using namespace RelSpacetime;

int main() {
    std::cout << "Running Iterative Solver & Buffer Tests...\n";

    MinkowskiMetric flatSpace;

    EntityState observerState;
    observerState.position = SpacetimeVector(10.0, 0.0, 0.0, 0.0);
    observerState.coordinateTime = 10.0;

    // Target's Ring Buffer. The target is moving along the X-axis at v = +0.5c. x(t) = 0.5 * t.
    std::vector<EntityState> targetBuffer;
    // add states from t=10 going backwards in time to t=5
    for (int i = 0; i <= 5; ++i) {
        double t = 10.0 - i;
        double x = 0.5 * t;
        
        EntityState state;
        state.coordinateTime = t;
        state.position = SpacetimeVector(t, x, 0.0, 0.0);
        state.fourVelocity = Eigen::Vector4d(1.0, 0.5, 0.0, 0.0);
        targetBuffer.push_back(state);
    }
 
    // Light travels at c=1. Distance D = Travel Time T.
    // Observer is at x=0. Light is emitted when target is at x_emit.
    // T = 10.0 - t_emit.  Also, D = x_emit = 0.5 * t_emit.
    // Therefore: 10.0 - t_emit = 0.5 * t_emit  =>  1.5 * t_emit = 10.0
    // t_emit = 6.6666...
    // T = 10.0 - 6.6666... = 3.3333...
    
    double exact_emission_time = 10.0 / 1.5; // ~6.66667
    double exact_travel_time = 10.0 - exact_emission_time; // ~3.33333
    double exact_x_position = 0.5 * exact_emission_time; // ~3.33333

    Eigen::Vector3d obs_pos(observerState.position.x(), observerState.position.y(), observerState.position.z());
    Eigen::Vector3d tgt_recent(targetBuffer[0].position.x(), targetBuffer[0].position.y(), targetBuffer[0].position.z());
    
    double guess_travel_time = (obs_pos - tgt_recent).norm(); 
    EntityState exactPastState;
    NullGeodesicResult lookbackData;

    for (int i = 0; i < 25; ++i) { 
        double emissionTime = observerState.coordinateTime - guess_travel_time;
        exactPastState = SpacetimeInterpolator::interpolateFromBuffer(targetBuffer, emissionTime, flatSpace);
        lookbackData = AnalyticalRaytracer::traceRays(observerState.position, exactPastState.position, flatSpace);
        
        if (!lookbackData.primary.isValid) break; 

        guess_travel_time = lookbackData.primary.travelTime;
    }

    std::cout << "Expected Emission Time: " << exact_emission_time << "\n";
    std::cout << "Solver Emission Time:   " << exactPastState.coordinateTime << "\n";
    std::cout << "Expected Travel Time:   " << exact_travel_time << "\n";
    std::cout << "Solver Travel Time:     " << guess_travel_time << "\n";

    assert(std::abs(exactPastState.coordinateTime - exact_emission_time) < 1e-4 && "FAIL: Solver failed to converge on the correct historical emission time.");
    assert(std::abs(guess_travel_time - exact_travel_time) < 1e-4 && "FAIL: Solver failed to converge on the correct light travel time.");
    assert(std::abs(exactPastState.position.x() - exact_x_position) < 1e-4 && "FAIL: Interpolated position does not match the intersection point.");

    std::cout << "Iterative Solver tests passed! The ghost target was successfully located in the buffer.\n";
    return 0;
}

