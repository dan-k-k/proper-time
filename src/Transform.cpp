// src/Transform.cpp
#include "RelSpacetime/Transform.h"

namespace RelSpacetime {

    Eigen::Matrix4d Transform::createBoost(const Eigen::Vector3d& velocity) {

        double v2 = velocity.squaredNorm();
        if (v2 < 1e-14) {
            return Eigen::Matrix4d::Identity();
        }

        if (v2 >= 1.0) {
            throw std::invalid_argument("Total velocity magnitude must be strictly less than the speed of light (v < 1.0).");
        }

        double gamma = 1.0 / std::sqrt(1.0 - v2);
        
        double vx = velocity.x();
        double vy = velocity.y();
        double vz = velocity.z();
        double K = (gamma - 1.0) / v2;

        Eigen::Matrix4d lambda;

        // Time
        lambda(0, 0) = gamma;
        lambda(0, 1) = -gamma * vx;
        lambda(0, 2) = -gamma * vy;
        lambda(0, 3) = -gamma * vz;

        // x
        lambda(1, 0) = -gamma * vx;
        lambda(1, 1) = 1.0 + K * vx * vx;
        lambda(1, 2) = K * vx * vy;
        lambda(1, 3) = K * vx * vz;

        // y
        lambda(2, 0) = -gamma * vy;
        lambda(2, 1) = K * vx * vy;
        lambda(2, 2) = 1.0 + K * vy * vy;
        lambda(2, 3) = K * vy * vz;

        // z
        lambda(3, 0) = -gamma * vz;
        lambda(3, 1) = K * vx * vz;
        lambda(3, 2) = K * vy * vz;
        lambda(3, 3) = 1.0 + K * vz * vz;

        return lambda;
    }
    // ---------
    SpacetimeVector Transform::apply(const Eigen::Matrix4d& transformMatrix, const SpacetimeVector& vec) {

        Eigen::Vector4d transformedData = transformMatrix * vec.getEigenVector();
        
        return SpacetimeVector(transformedData[0], transformedData[1], transformedData[2], transformedData[3]);
    }
    // ---------
    Eigen::Vector3d Transform::addVelocities(const Eigen::Vector3d& v, const Eigen::Vector3d& u_prime) {
        double v2 = v.squaredNorm();

        // For simple calculations 
        if (v2 < 1e-14) {
            return u_prime;
        }

        if (v2 >= 1.0) {
            throw std::invalid_argument("Reference frame velocity (v) must be < 1.0c");
        }

        double gamma = 1.0 / std::sqrt(1.0 - v2);
        double dot_product = v.dot(u_prime);
        double denominator = 1.0 + dot_product; 
        Eigen::Vector3d u = (u_prime / gamma + v + (gamma / (1.0 + gamma)) * dot_product * v) / denominator;

        return u;
    }

} // namespace RelSpacetime

