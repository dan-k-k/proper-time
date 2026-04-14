// tests/test_integration.cpp
#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm> // Needed for std::reverse
#include "RelSpacetime/RelativityTypes.h"
#include "RelSpacetime/MetricTensor.h"
#include "RelSpacetime/WorldlineIntegrator.h"
#include "RelSpacetime/LocalFrame.h"
#include "RelSpacetime/RelativityVisualsAPI.h"
#include "RelSpacetime/SpacetimeInterpolator.h"

using namespace RelSpacetime;

// Helper to create valid starting states
EntityState createRestState(double x, double y, double z, const MetricTensor& metric) {
    EntityState state;
    state.position = SpacetimeVector(0.0, x, y, z);
    Eigen::Matrix4d g = metric.getMetricAt(state.position);
    double u0 = 1.0 / std::sqrt(-g(0, 0));
    state.fourVelocity = Eigen::Vector4d(u0, 0.0, 0.0, 0.0);
    state.properTime = 0.0;
    state.coordinateTime = 0.0;
    state.timeScale = 1.0 / u0;
    return state;
}

int main() {
    std::cout << " INIT: RelSpacetime SDK Headless Integration Test \n";

    KerrMetric blackHole(1.0, 0.9); // Mass=1, Spin=0.9
    RelativityVisualsAPI visualsAPI(1.0); // Spectral index = 1.0

    // Player (Observer) is orbiting safely at x=15
    EntityState player = createRestState(15.0, 0.0, 0.0, blackHole);
    // Opponent (Target) is deep in the gravity well at x=3
    EntityState opponent = createRestState(3.0, 0.0, 0.0, blackHole);
    double server_dt = 0.1; // 10 ticks per second
    EngineInputs playerInputs = { Eigen::Vector3d::Zero(), Eigen::Vector3d(1,0,0), Eigen::Vector3d(0,1,0), Eigen::Vector3d(0,0,1) };
    // Opponent is thrusting in their local +Y direction
    EngineInputs opponentInputs = { Eigen::Vector3d(0.0, 5.0, 0.0), Eigen::Vector3d(1,0,0), Eigen::Vector3d(0,1,0), Eigen::Vector3d(0,0,1) };

    std::vector<EntityState> opponentHistory;
    opponentHistory.push_back(opponent);

    std::cout << ">>> Simulating 20 Server Ticks...\n";
    for (int i = 0; i < 20; ++i) {
        player = WorldlineIntegrator::computeNextState(player, playerInputs, server_dt, blackHole);
        opponent = WorldlineIntegrator::computeNextState(opponent, opponentInputs, server_dt, blackHole);
        opponentHistory.push_back(opponent);
    }

    std::cout << "\n>>> EXECUTING SDK VISUAL QUERY...\n";

    std::vector<EntityState> godotRingBuffer = opponentHistory;
    std::reverse(godotRingBuffer.begin(), godotRingBuffer.end());

    // Player's Local Frame for visual rendering
    Eigen::Matrix4d playerMetric = blackHole.getMetricAt(player.position);
    LocalFrame playerFrame(playerMetric, player.fourVelocity, playerInputs.forward, playerInputs.right, playerInputs.up);

    // Pass the whole buffer
    VisualRenderState finalVisuals = visualsAPI.getVisualStateForOpponent(
        player, playerFrame, godotRingBuffer, blackHole
    );

    std::cout << "\n[ MECHANICAL OUTPUTS (The Physics Truth) ]\n";
    std::cout << "1. TimeScale (Gameplay Scaler): " << opponent.timeScale << " (Multiplier for Godot's delta_time)\n";
    std::cout << "2. Proper Time (Lore/Age):      " << opponent.properTime << " seconds\n";
    std::cout << "3. Coordinate State (4D Pos):   (t:" << opponent.coordinateTime 
              << ", x:" << opponent.position.x() 
              << ", y:" << opponent.position.y() 
              << ", z:" << opponent.position.z() << ")\n";
    std::cout << "4. 4-Velocity (Server Truth):   (U^0:" << opponent.fourVelocity[0] 
              << ", U^1:" << opponent.fourVelocity[1] 
              << ", U^2:" << opponent.fourVelocity[2] 
              << ", U^3:" << opponent.fourVelocity[3] << ")\n";

    std::cout << "\n[ PRIMARY VISUAL OUTPUTS (The Direct Illusion) ]\n";
    std::cout << "1. Is Visible:                  " << (finalVisuals.primary.isVisible ? "True" : "False (Occluded!)") << "\n";
    if (finalVisuals.primary.isVisible) {
        std::cout << "2. Light Travel Time:           " << finalVisuals.primary.travelTime << " seconds\n";
        std::cout << "3. Apparent Direction (Local):  (" << finalVisuals.primary.apparentDirection.x() << ", " 
                                                        << finalVisuals.primary.apparentDirection.y() << ", " 
                                                        << finalVisuals.primary.apparentDirection.z() << ")\n";
        std::cout << "4. Doppler Shift:               " << finalVisuals.primary.dopplerShift << "\n";
        std::cout << "5. Beaming Factor:              " << finalVisuals.primary.beamingFactor << "\n";
        std::cout << "6. Tangential Stretch:          " << finalVisuals.primary.tangentialStretch << "\n";
        std::cout << "7. Radial Squish:               " << finalVisuals.primary.radialSquish << "\n";
        std::cout << "8. Ring Factor:                 " << finalVisuals.primary.ringFactor << "\n";
        std::cout << "9. Arc Tangent Axis:            (" << finalVisuals.primary.arcTangentAxis.x() << ", " 
                                                        << finalVisuals.primary.arcTangentAxis.y() << ", " 
                                                        << finalVisuals.primary.arcTangentAxis.z() << ")\n";
    }

    std::cout << "\n[ SECONDARY VISUAL OUTPUTS (The Curved Illusion) ]\n";
    std::cout << "1. Is Visible:                  " << (finalVisuals.secondary.isVisible ? "True" : "False (Occluded!)") << "\n";
    if (finalVisuals.secondary.isVisible) {
        std::cout << "2. Light Travel Time:           " << finalVisuals.secondary.travelTime << " seconds\n";
        std::cout << "3. Apparent Direction (Local):  (" << finalVisuals.secondary.apparentDirection.x() << ", " 
                                                        << finalVisuals.secondary.apparentDirection.y() << ", " 
                                                        << finalVisuals.secondary.apparentDirection.z() << ")\n";
        std::cout << "4. Doppler Shift:               " << finalVisuals.secondary.dopplerShift << "\n";
        std::cout << "5. Beaming Factor:              " << finalVisuals.secondary.beamingFactor << "\n";
        std::cout << "6. Tangential Stretch:          " << finalVisuals.secondary.tangentialStretch << "\n";
        std::cout << "7. Radial Squish:               " << finalVisuals.secondary.radialSquish << "\n";
        std::cout << "8. Ring Factor:                 " << finalVisuals.secondary.ringFactor << "\n";
        std::cout << "9. Arc Tangent Axis:            (" << finalVisuals.secondary.arcTangentAxis.x() << ", " 
                                                        << finalVisuals.secondary.arcTangentAxis.y() << ", " 
                                                        << finalVisuals.secondary.arcTangentAxis.z() << ")\n";
    }

    assert(player.properTime > opponent.properTime && "Integration FAIL: Deep well opponent should age slower!");
    if (finalVisuals.primary.isVisible) {
        assert(finalVisuals.primary.dopplerShift > 0.0 && "Integration FAIL: Doppler shift collapsed.");
        assert(finalVisuals.primary.travelTime > 0.0 && "Integration FAIL: Iterative solver failed to find a valid travel time.");
    }
    if (finalVisuals.primary.isVisible && finalVisuals.secondary.isVisible) {
        assert(finalVisuals.secondary.travelTime > finalVisuals.primary.travelTime && "Integration FAIL: Secondary image must take longer to arrive!");
    }
    
    std::cout << " Success. SDK Pipeline is Fully Operational. \n";

    return 0;
}

