// include/RelSpacetime/Transform.h
#pragma once
#include <Eigen/Dense>
#include "SpacetimeVector.h"

namespace RelSpacetime {
    class Transform {
    public:
        static Eigen::Matrix4d createBoost(const Eigen::Vector3d& velocity);
        static SpacetimeVector apply(const Eigen::Matrix4d& transformMatrix, const SpacetimeVector& vec);
        static Eigen::Vector3d addVelocities(const Eigen::Vector3d& v, const Eigen::Vector3d& u_prime);
    };
}

