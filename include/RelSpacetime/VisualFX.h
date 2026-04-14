// include/RelSpacetime/VisualFX.h
#pragma once
#include <Eigen/Dense>
#include "RelativityTypes.h"
#include "LocalFrame.h"
#include "MetricTensor.h"

namespace RelSpacetime {

    struct RenderImage {
        bool isVisible; 
        Eigen::Vector3d apparentDirection; 
        double dopplerShift; 
        double beamingFactor; 
        double travelTime; 
        Eigen::Vector3d arcTangentAxis; 
        double tangentialStretch;        
        double radialSquish;             
        double ringFactor; 
    };

    struct VisualRenderState {
        RenderImage primary;
        RenderImage secondary;
        Eigen::Matrix4d localTetrad; 
    };

    struct BlackHoleRenderState {
        Eigen::Vector3d apparentDirection; 
        double shadowAngularRadius; 
        double dopplerScale;        
        Eigen::Vector3d squishAxis; 
    };
    
    class VisualFX {
    private:
        double spectralIndex; 

    public:
        VisualFX(double initialSpectralIndex = 1.0) : spectralIndex(initialSpectralIndex) {}
        void setSpectralIndex(double alpha) { spectralIndex = alpha; }
            
        // Visuals using the Observer's current state and the Target's past state (from the ring buffer), and the photon that connects them.
        VisualRenderState computeVisuals(
            const EntityState& observerCurrentState,
            const LocalFrame& observerFrame,
            const EntityState& primaryPastState, 
            const EntityState& secondaryPastState, 
            const NullGeodesicResult& geodesicData, 
            const MetricTensor& metric) const;

        BlackHoleRenderState computeBlackHoleVisuals(
            const EntityState& observerCurrentState,
            const LocalFrame& observerFrame,
            const MetricTensor& metric) const;
    };

} // namespace RelSpacetime

