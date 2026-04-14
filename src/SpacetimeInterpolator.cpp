// src/SpacetimeInterpolator.cpp
#include "RelSpacetime/SpacetimeInterpolator.h"
#include <cmath>
#include <algorithm>

namespace RelSpacetime {

    EntityState SpacetimeInterpolator::interpolate(
        const EntityState& stateA, 
        const EntityState& stateB, 
        double targetTime, 
        const MetricTensor& metric) 
    {
        double tA = stateA.coordinateTime;
        double tB = stateB.coordinateTime;
        double dt = tB - tA;

        if (dt <= 1e-9) return stateB;
        if (targetTime <= tA) return stateA;
        if (targetTime >= tB) return stateB;

        double f = (targetTime - tA) / dt;

        // 3-Positions
        Eigen::Vector3d pA(stateA.position.x(), stateA.position.y(), stateA.position.z());
        Eigen::Vector3d pB(stateB.position.x(), stateB.position.y(), stateB.position.z());

        // 3-Velocities (v^i = U^i / U^0)
        Eigen::Vector3d vA(stateA.fourVelocity[1] / stateA.fourVelocity[0],
                           stateA.fourVelocity[2] / stateA.fourVelocity[0],
                           stateA.fourVelocity[3] / stateA.fourVelocity[0]);
                           
        Eigen::Vector3d vB(stateB.fourVelocity[1] / stateB.fourVelocity[0],
                           stateB.fourVelocity[2] / stateB.fourVelocity[0],
                           stateB.fourVelocity[3] / stateB.fourVelocity[0]);

        // Tangents scaled by dt for Hermite Spline (dx/df = dx/dt * dt/df)
        Eigen::Vector3d mA = vA * dt;
        Eigen::Vector3d mB = vB * dt;

        // Cubic Hermite basis functions
        double f2 = f * f;
        double f3 = f2 * f;

        double h00 =  2.0 * f3 - 3.0 * f2 + 1.0;
        double h10 =  f3 - 2.0 * f2 + f;
        double h01 = -2.0 * f3 + 3.0 * f2;
        double h11 =  f3 - f2;

        // Derivatives with respect to f
        double dh00 =  6.0 * f2 - 6.0 * f;
        double dh10 =  3.0 * f2 - 4.0 * f + 1.0;
        double dh01 = -6.0 * f2 + 6.0 * f;
        double dh11 =  3.0 * f2 - 2.0 * f;

        // Interpolation
        Eigen::Vector3d interpPos = h00 * pA + h10 * mA + h01 * pB + h11 * mB;
        // v = dP/dt = (dP/df) * (df/dt) = (dP/df) / dt
        Eigen::Vector3d interpVel = (dh00 * pA + dh10 * mA + dh01 * pB + dh11 * mB) / dt;

        // Linear interpolation for proper time
        double interpProperTime = stateA.properTime + f * (stateB.properTime - stateA.properTime);

        EntityState result;
        result.coordinateTime = targetTime;
        result.properTime = interpProperTime;
        result.position = SpacetimeVector(targetTime, interpPos.x(), interpPos.y(), interpPos.z());

        // Renormalise 4-Velocity using the metric
        // V^mu = (1, v^1, v^2, v^3). 
        // U^mu = U^0 * V^mu, such that g_{\mu\nu} U^\mu U^\nu = -1
        Eigen::Vector4d V_coordinate(1.0, interpVel.x(), interpVel.y(), interpVel.z());
        Eigen::Matrix4d g = metric.getMetricAt(result.position);
        // V_normSq = g_{\mu\nu} V^\mu V^\nu
        double V_normSq = V_coordinate.transpose() * g * V_coordinate;

        if (V_normSq < 0.0) {
            // Valid timelike
            double U0 = std::sqrt(-1.0 / V_normSq);
            result.fourVelocity = V_coordinate * U0;
            result.timeScale = 1.0 / U0;
        } else {
            result.fourVelocity = stateA.fourVelocity * (1.0 - f) + stateB.fourVelocity * f;
            double directNormSq = result.fourVelocity.transpose() * g * result.fourVelocity;
            if (directNormSq < 0.0) {
                result.fourVelocity /= std::sqrt(-directNormSq);
            }
            result.timeScale = 1.0 / result.fourVelocity[0];
        }

        return result;
    }

    EntityState SpacetimeInterpolator::interpolateFromBuffer(
        const std::vector<EntityState>& historyBuffer, 
        double targetCoordinateTime, 
        const MetricTensor& metric) 
    {
        if (historyBuffer.empty()) return EntityState();
        if (historyBuffer.size() == 1) return historyBuffer[0];

        // historyBuffer be sorted newest 0 to oldest N.
        for (size_t i = 0; i < historyBuffer.size() - 1; ++i) {
            const EntityState& newer = historyBuffer[i];
            const EntityState& older = historyBuffer[i+1];

            if (targetCoordinateTime <= newer.coordinateTime && targetCoordinateTime >= older.coordinateTime) {
                return interpolate(older, newer, targetCoordinateTime, metric);
            }
        }
        if (targetCoordinateTime < historyBuffer.back().coordinateTime) return historyBuffer.back();
        return historyBuffer.front();
    }

} // namespace RelSpacetime..

