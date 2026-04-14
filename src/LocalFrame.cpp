// src/LocalFrame.cpp
#include "RelSpacetime/LocalFrame.h"
#include <cmath>

namespace RelSpacetime {

    LocalFrame::LocalFrame(const Eigen::Matrix4d& g, 
                           const Eigen::Vector4d& fourVel, 
                           const Eigen::Vector3d& engineForward, 
                           const Eigen::Vector3d& engineRight, 
                           const Eigen::Vector3d& engineUp) 
        : metric(g), e0(fourVel) 
    {
        Eigen::Vector4d v_forward(0.0, engineForward.x(), engineForward.y(), engineForward.z());
        Eigen::Vector4d v_right(0.0, engineRight.x(), engineRight.y(), engineRight.z());
        Eigen::Vector4d v_up(0.0, engineUp.x(), engineUp.y(), engineUp.z());

        // Gram-Schmidt orthogonalisation
        std::vector<Eigen::Vector4d> basis;
        
        basis.push_back(e0); 
        e1 = projectAndNormalize(v_forward, basis);
        basis.push_back(e1);
        e2 = projectAndNormalize(v_right, basis);
        basis.push_back(e2);
        e3 = projectAndNormalize(v_up, basis);
    }

    Eigen::Vector4d LocalFrame::projectAndNormalize(const Eigen::Vector4d& v, const std::vector<Eigen::Vector4d>& basis) const {
        Eigen::Vector4d result = v;
        
        for (const auto& e : basis) {
            double normSq = e.transpose() * metric * e; 
            if (std::abs(normSq) < 1e-12) {
                continue; 
            }
            double dotProd = v.transpose() * metric * e;
            result -= (dotProd / normSq) * e;
        }
        double magSq = result.transpose() * metric * result;
        if (magSq > 1e-12) {
            result /= std::sqrt(magSq);
        } else {
            result.setZero(); 
        }
        return result;
    }

    Eigen::Vector4d LocalFrame::getGlobalFourAcceleration(const Eigen::Vector3d& localProperThrust) const {
        return e1 * localProperThrust.x() + e2 * localProperThrust.y() + e3 * localProperThrust.z();
    }

    LocalFrame::LocalFrame(const Eigen::Matrix4d& g, const Eigen::Vector4d& fourVel) 
        : metric(g), e0(fourVel) 
    {
        Eigen::Vector4d global_x(0.0, 1.0, 0.0, 0.0);
        Eigen::Vector4d global_y(0.0, 0.0, 1.0, 0.0);
        Eigen::Vector4d global_z(0.0, 0.0, 0.0, 1.0);

        std::vector<Eigen::Vector4d> basis;
        
        basis.push_back(e0); 
        e1 = projectAndNormalize(global_x, basis);
        basis.push_back(e1);
        e2 = projectAndNormalize(global_y, basis);
        basis.push_back(e2);
        e3 = projectAndNormalize(global_z, basis);
    } 

} // namespace RelSpacetime

