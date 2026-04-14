// include/RelSpacetime/DelayedTimeSolver.h
#pragma once
#include <vector>
#include "RelativityTypes.h"
#include "MetricTensor.h"

namespace RelSpacetime {

    struct DelayedTimeResult {
        EntityState primaryPastState;
        EntityState secondaryPastState;
        NullGeodesicResult geodesicData;
        bool converged; 
    };

    class DelayedTimeSolver {
    public:
        static DelayedTimeResult solve(
            const EntityState& observerCurrentState,
            const std::vector<EntityState>& targetHistoryBuffer,
            const MetricTensor& metric,
            double previous_travel_time = -1.0,
            int max_iterations = 4);
    };

} // namespace RelSpacetime

