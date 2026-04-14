// include/RelSpacetime/SpacetimeInterpolator.h
#pragma once
#include "RelativityTypes.h"
#include "MetricTensor.h"
#include <vector> 

namespace RelSpacetime {

    class SpacetimeInterpolator {
    public:
        static EntityState interpolate(
            const EntityState& stateA, 
            const EntityState& stateB, 
            double targetCoordinateTime, 
            const MetricTensor& metric);
        
        static EntityState interpolateFromBuffer(
            const std::vector<EntityState>& historyBuffer, 
            double targetCoordinateTime, 
            const MetricTensor& metric);
    };

} // namespace RelSpacetime

