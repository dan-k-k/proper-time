// include/RelSpacetime/RelativityVisualsAPI.h
#pragma once
#include <vector>
#include "RelativityTypes.h"
#include "LocalFrame.h"
#include "MetricTensor.h"
#include "VisualFX.h"
#include "DelayedTimeSolver.h"

namespace RelSpacetime {

    class RelativityVisualsAPI {
    private:
        VisualFX visFX;

    public:
        RelativityVisualsAPI(double spectralIndex) : visFX(spectralIndex) {}

        VisualRenderState getVisualStateForOpponent(
            const EntityState& observerCurrentState,
            const LocalFrame& observerFrame,
            const std::vector<EntityState>& targetHistoryBuffer,
            const MetricTensor& metric,
            double previous_travel_time = -1.0) const
        {
            if (targetHistoryBuffer.empty()) return VisualRenderState();

            DelayedTimeResult solvedState = DelayedTimeSolver::solve(
                observerCurrentState, targetHistoryBuffer, metric, previous_travel_time);

            return visFX.computeVisuals(
                observerCurrentState, 
                observerFrame, 
                solvedState.primaryPastState, 
                solvedState.secondaryPastState, 
                solvedState.geodesicData, 
                metric);
        }

        BlackHoleRenderState getBlackHoleVisualState(
            const EntityState& observerCurrentState,
            const LocalFrame& observerFrame,
            const MetricTensor& metric) const
        {
            return visFX.computeBlackHoleVisuals(observerCurrentState, observerFrame, metric);
        }
    };

} // namespace RelSpacetime

