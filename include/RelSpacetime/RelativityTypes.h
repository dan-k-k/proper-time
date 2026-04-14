// include/RelSpacetime/RelativityTypes.h
#pragma once
#include <Eigen/Dense>
#include "SpacetimeVector.h"

namespace RelSpacetime {

    // The single source of truth for an entity's physical state. 
    // Game engine processes ring buffers and passes it to the SDK.
    struct EntityState {
        SpacetimeVector position;
        Eigen::Vector4d fourVelocity;
        double properTime;
        double coordinateTime;
        double timeScale; // d\tau / dt 
    };

    // Inputs from the Game Engine.
    struct EngineInputs {
        Eigen::Vector3d localProperThrust;
        Eigen::Vector3d forward;
        Eigen::Vector3d right;
        Eigen::Vector3d up;
    };

    struct RayData {
        bool isValid;                // Did the photon actually reach the observer
        double travelTime;           // Target-observer time/distance
        Eigen::Vector4d arrivalK;    // The 4-momentum/wave-vector of the photon as it arrives at the observer
        Eigen::Vector3d arcTangentAxis; 
        double tangentialStretch;        
        double radialSquish;             
        double ringFactor;               
    };

    struct NullGeodesicResult {
        RayData primary;
        RayData secondary;
    };

} // namespace RelSpacetime

