// src/DelayedTimeSolver.cpp 
#include "RelSpacetime/DelayedTimeSolver.h"
#include "RelSpacetime/SpacetimeInterpolator.h"
#include "RelSpacetime/AnalyticalRaytracer.h"

namespace RelSpacetime {

    DelayedTimeResult DelayedTimeSolver::solve(
        const EntityState& observerCurrentState,
        const std::vector<EntityState>& targetHistoryBuffer,
        const MetricTensor& metric,
        double previous_travel_time,
        int max_iterations) 
    {
        DelayedTimeResult result;
        result.converged = false;

        if (targetHistoryBuffer.empty()) return result;

        double guess_travel_time;

        // Initial guess
        if (previous_travel_time > 0.0) {
            guess_travel_time = previous_travel_time;
        } else {
            Eigen::Vector3d obs_pos(observerCurrentState.position.x(), observerCurrentState.position.y(), observerCurrentState.position.z());
            Eigen::Vector3d tgt_recent(targetHistoryBuffer[0].position.x(), targetHistoryBuffer[0].position.y(), targetHistoryBuffer[0].position.z());
            guess_travel_time = (obs_pos - tgt_recent).norm();
        }

        // Iterative solving 
        for (int i = 0; i < max_iterations; ++i) {
            double emissionTime = observerCurrentState.coordinateTime - guess_travel_time;
            result.primaryPastState = SpacetimeInterpolator::interpolateFromBuffer(targetHistoryBuffer, emissionTime, metric);

            result.geodesicData = AnalyticalRaytracer::traceRays(observerCurrentState.position, result.primaryPastState.position, metric);
            
            if (!result.geodesicData.primary.isValid) break; 
            
            // Check for convergence 
            if (std::abs(guess_travel_time - result.geodesicData.primary.travelTime) < 1e-5) {
                result.converged = true;
                break;
            }

            guess_travel_time = result.geodesicData.primary.travelTime;
        }

        // Secondary image's past state
        result.secondaryPastState = result.primaryPastState; 
        if (result.geodesicData.secondary.isValid) {
            double secEmissionTime = observerCurrentState.coordinateTime - result.geodesicData.secondary.travelTime;
            result.secondaryPastState = SpacetimeInterpolator::interpolateFromBuffer(targetHistoryBuffer, secEmissionTime, metric);
        }

        return result;
    }

} // namespace RelSpacetime.

