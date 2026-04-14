// include/RelSpacetime/LocalFrame.h
#pragma once
#include <Eigen/Dense>
#include <vector>

namespace RelSpacetime {

    class LocalFrame {
    private:
        Eigen::Matrix4d metric;
        Eigen::Vector4d e0; // Time-like 
        Eigen::Vector4d e1; // 
        Eigen::Vector4d e2; // 
        Eigen::Vector4d e3; // 

        Eigen::Vector4d projectAndNormalize(const Eigen::Vector4d& v, const std::vector<Eigen::Vector4d>& basis) const;

    public:
        // Constructs the orthonormal tetrad from the metric, current 4-velocity, and engine orientation vectors
        LocalFrame(const Eigen::Matrix4d& g, 
                   const Eigen::Vector4d& fourVel, 
                   const Eigen::Vector3d& engineForward, 
                   const Eigen::Vector3d& engineRight, 
                   const Eigen::Vector3d& engineUp);

        // Maps a local 3D proper thrust vector to a global 4D spacetime acceleration vector
        Eigen::Vector4d getGlobalFourAcceleration(const Eigen::Vector3d& localProperThrust) const;

        Eigen::Vector4d getE0() const { return e0; }
        Eigen::Vector4d getE1() const { return e1; }
        Eigen::Vector4d getE2() const { return e2; }
        Eigen::Vector4d getE3() const { return e3; }

        LocalFrame(const Eigen::Matrix4d& g, const Eigen::Vector4d& fourVel);
    };

} // namespace RelSpacetime

