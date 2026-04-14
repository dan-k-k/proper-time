// src/WorldlineIntegrator.cpp
#include "RelSpacetime/WorldlineIntegrator.h"
#include "RelSpacetime/LocalFrame.h"
#include <cmath>

namespace RelSpacetime {

    struct KinematicState {
        Eigen::Vector3d pos; 
        Eigen::Vector3d vel; 
        double tau;          

        KinematicState operator+(const KinematicState& other) const {
            return { pos + other.pos, vel + other.vel, tau + other.tau };
        }
        KinematicState operator*(double scalar) const {
            return { pos * scalar, vel * scalar, tau * scalar };
        }
    };

    EntityState WorldlineIntegrator::computeNextState(
        const EntityState& currentState, 
        const EngineInputs& inputs,
        double server_dt, 
        const MetricTensor& metric) 
    {
        EntityState nextState = currentState; // Start with a copy to mutate

        Eigen::Matrix4d currentMetric = metric.getMetricAt(currentState.position);
        LocalFrame tetrad(currentMetric, currentState.fourVelocity, inputs.forward, inputs.right, inputs.up);
        Eigen::Vector4d a_thrust_4d = tetrad.getGlobalFourAcceleration(inputs.localProperThrust);

        KinematicState rkState = {
            Eigen::Vector3d(currentState.position.x(), currentState.position.y(), currentState.position.z()),
            Eigen::Vector3d(currentState.fourVelocity[1] / currentState.fourVelocity[0], 
                            currentState.fourVelocity[2] / currentState.fourVelocity[0], 
                            currentState.fourVelocity[3] / currentState.fourVelocity[0]),
            currentState.properTime
        };

        auto computeDerivative = [&](double t_offset, const KinematicState& state) -> KinematicState {
            KinematicState derivative;
            double t_eval = currentState.coordinateTime + t_offset;
            SpacetimeVector eval_pos(t_eval, state.pos.x(), state.pos.y(), state.pos.z());

            auto gamma = metric.getChristoffelSymbols(eval_pos);
            auto g = metric.getMetricAt(eval_pos);

            Eigen::Vector4d v(1.0, state.vel.x(), state.vel.y(), state.vel.z());

            // Proper time derivative. intermediate U0 early
            double ds2 = 0.0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    ds2 += g(mu, nu) * v[mu] * v[nu];
                }
            }
            double dtau_dt = (ds2 < 0.0) ? std::sqrt(-ds2) : 1e-9; 
            double U0 = 1.0 / dtau_dt; 
            
            // Intermediate 4-velocity
            Eigen::Vector4d U_intermediate = v * U0;

            // Thrust 4-vector for RK4 sub-step
            LocalFrame intermediateTetrad(g, U_intermediate, inputs.forward, inputs.right, inputs.up);
            Eigen::Vector4d a_thrust_4d = intermediateTetrad.getGlobalFourAcceleration(inputs.localProperThrust);

            // Base acceleration
            Eigen::Vector4d A = Eigen::Vector4d::Zero();
            for (int mu = 0; mu < 4; ++mu) {
                for (int alpha = 0; alpha < 4; ++alpha) {
                    for (int beta = 0; beta < 4; ++beta) {
                        A[mu] -= gamma[mu](alpha, beta) * v[alpha] * v[beta];
                    }
                }
            }
            
            // Coordinate acceleration
            Eigen::Vector3d accel;
            for (int i = 0; i < 3; ++i) {
                // Now using the properly boosted a_thrust_4d
                double thrustTerm = (a_thrust_4d[i + 1] - state.vel[i] * a_thrust_4d[0]) / (U0 * U0);
                accel[i] = A[i + 1] - state.vel[i] * A[0] + thrustTerm;
            }

            derivative.pos = state.vel;       
            derivative.vel = accel;           
            derivative.tau = dtau_dt;         
            
            return derivative;
        };

        // RK4 Integration
        KinematicState k1 = computeDerivative(0.0, rkState);
        KinematicState k2 = computeDerivative(server_dt * 0.5, rkState + k1 * (server_dt * 0.5));
        KinematicState k3 = computeDerivative(server_dt * 0.5, rkState + k2 * (server_dt * 0.5));
        KinematicState k4 = computeDerivative(server_dt, rkState + k3 * server_dt);

        KinematicState finalRkState = rkState + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (server_dt / 6.0);

        nextState.coordinateTime += server_dt;
        nextState.properTime = finalRkState.tau;
        nextState.position = SpacetimeVector(nextState.coordinateTime, finalRkState.pos.x(), finalRkState.pos.y(), finalRkState.pos.z());

        // Reconstruct U^mu
        double final_dtau_dt = computeDerivative(server_dt, finalRkState).tau;
        if (final_dtau_dt > 1e-9) { 
            double gamma_factor = 1.0 / final_dtau_dt; 
            nextState.fourVelocity[0] = gamma_factor; 
            nextState.fourVelocity[1] = finalRkState.vel.x() * gamma_factor;
            nextState.fourVelocity[2] = finalRkState.vel.y() * gamma_factor;
            nextState.fourVelocity[3] = finalRkState.vel.z() * gamma_factor;
        }

        // Re-normalise U^mu
        Eigen::Matrix4d finalMetric = metric.getMetricAt(nextState.position);
        double normSq = nextState.fourVelocity.transpose() * finalMetric * nextState.fourVelocity;
        if (normSq < 0.0) {
            nextState.fourVelocity /= std::sqrt(-normSq);
        }

        nextState.timeScale = 1.0 / nextState.fourVelocity[0];

        return nextState; // Hand the new truth back to game engine
    }

} // namespace RelSpacetime

