// include/RelSpacetime/WorldlineIntegrator.h
#pragma once
#include "RelativityTypes.h"
#include "MetricTensor.h"

namespace RelSpacetime {

    class WorldlineIntegrator {
    public:
        // For the next server tick
        static EntityState computeNextState(
            const EntityState& currentState, 
            const EngineInputs& inputs,
            double server_dt, 
            const MetricTensor& metric);
    };

} // namespace RelSpacetime

