// include/RelSpacetime/AnalyticalRaytracer.h
#pragma once 
#include <Eigen/Dense>
#include "RelativityTypes.h"
#include "MetricTensor.h"

namespace RelSpacetime {

    class AnalyticalRaytracer {
    public:
        // Unified function that calculates both travel time and apparent visual bending
        static NullGeodesicResult traceRays(
            const SpacetimeVector& observerPos, 
            const SpacetimeVector& targetPos, 
            const MetricTensor& metric);
    };

} // namespace RelSpacetime

