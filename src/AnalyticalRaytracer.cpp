// src/AnalyticalRaytracer.cpp
#include "RelSpacetime/AnalyticalRaytracer.h"
#include <Eigen/Geometry> 
#include <cmath>
#include <algorithm>

namespace RelSpacetime {

    NullGeodesicResult AnalyticalRaytracer::traceRays(
        const SpacetimeVector& observerPos, 
        const SpacetimeVector& targetPos, 
        const MetricTensor& metric) 
    {
        NullGeodesicResult result;
        double mass = metric.getMass();
        
        Eigen::Vector3d obs(observerPos.x(), observerPos.y(), observerPos.z());
        Eigen::Vector3d tgt(targetPos.x(), targetPos.y(), targetPos.z());
        Eigen::Vector3d vecToTarget = tgt - obs;
        double flatDistance = vecToTarget.norm();

        // Touching / Minkowski space (Mass = 0)
        if (flatDistance < 1e-6 || mass < 1e-6) {
            Eigen::Vector3d dir = (flatDistance < 1e-6) ? Eigen::Vector3d::Zero() : vecToTarget.normalized();
            
            result.primary.isValid = true;
            result.primary.travelTime = flatDistance;
            result.primary.arrivalK = metric.enforceNullCondition(observerPos, dir);
            result.primary.arcTangentAxis = Eigen::Vector3d::Zero();
            result.primary.tangentialStretch = 1.0;
            result.primary.radialSquish = 1.0;
            result.primary.ringFactor = 0.0;
            
            result.secondary.isValid = false;
            return result;
        }

        // Base setup for single black hole
        Eigen::Vector3d bh_pos = Eigen::Vector3d::Zero(); 
        Eigen::Vector3d viewDir = vecToTarget.normalized();
        Eigen::Vector3d obs_rel = obs - bh_pos;
        Eigen::Vector3d tgt_rel = tgt - bh_pos; 
        double r_obs = obs_rel.norm(); 
        double r_tgt = tgt_rel.norm(); 
        Eigen::Vector3d dirToBH = -obs_rel.normalized();

        // Horizon survival check
        if (r_obs <= mass || r_tgt <= mass) {
            result.primary.isValid = false; 
            result.secondary.isValid = false; 
            return result;
        }

        double cos_beta = std::clamp(dirToBH.dot(viewDir), -1.0, 1.0);
        double beta = std::acos(cos_beta);

        // Plane of scattering
        Eigen::Vector3d rotationAxis = dirToBH.cross(viewDir);
        double t_closest = -obs_rel.dot(viewDir);
        bool transitsBH = (t_closest > 0.0 && t_closest < flatDistance);

        if (rotationAxis.norm() < 1e-6) {
            if (!transitsBH) {
                // Looking directly away from BH. Primary unlensed, Secondary forms a perfect ring behind you.
                Eigen::Vector3d ringAxis = Eigen::Vector3d(0, 1, 0).cross(dirToBH);
                if (ringAxis.norm() < 1e-6) ringAxis = Eigen::Vector3d(1, 0, 0).cross(dirToBH);
                ringAxis.normalize();

                result.primary.isValid = true;
                result.primary.travelTime = flatDistance; // Neglecting Shapiro looking straight away
                result.primary.arrivalK = metric.enforceNullCondition(observerPos, viewDir);
                result.primary.arcTangentAxis = Eigen::Vector3d::Zero();
                result.primary.tangentialStretch = 1.0;
                result.primary.radialSquish = 1.0;
                result.primary.ringFactor = 0.0;

                result.secondary.isValid = true;
                result.secondary.travelTime = flatDistance + (M_PI * r_obs); // Approximation for ring detour
                result.secondary.arrivalK = metric.enforceNullCondition(observerPos, dirToBH);
                result.secondary.arcTangentAxis = ringAxis;
                result.secondary.tangentialStretch = 15.0;
                result.secondary.radialSquish = 0.05;
                result.secondary.ringFactor = 1.0;
                return result;
            }
            rotationAxis.normalize();
        } else {
            rotationAxis.normalize();
        }

        // PRIMARY IMAGE
        result.primary.isValid = true;
        
        // Primary Shapiro time delay
        double rawDenominator = r_obs + r_tgt - flatDistance;
        double safeDenominator = std::max(rawDenominator, 1e-5 * mass);
        double shapiroDelay = 2.0 * mass * std::log((r_obs + r_tgt + flatDistance) / safeDenominator);
        result.primary.travelTime = flatDistance + shapiroDelay;

        // Primary spatial bending
        double r_s = 2.0 * mass;
        Eigen::Vector3d closestPoint = obs_rel + t_closest * viewDir;
        double b = closestPoint.norm(); 
        double safeB = std::max(b, 0.1 * r_s);
        double deflectionAngle = (4.0 * mass) / safeB;

        double x_obs = -t_closest;
        double x_tgt = flatDistance - t_closest;
        double shiftSameSide = (deflectionAngle / 2.0) * (
            ((x_tgt / flatDistance) * (x_tgt / std::sqrt(b * b + x_tgt * x_tgt) - x_obs / std::sqrt(b * b + x_obs * x_obs))) +
            (((b * b) / flatDistance) * (1.0 / std::sqrt(b * b + x_tgt * x_tgt) - 1.0 / std::sqrt(b * b + x_obs * x_obs)))
        );

        double D_l = r_obs;
        double D_ls = r_tgt;
        double D_s = D_l + D_ls; 
        double theta_E_sq = (4.0 * mass * D_ls) / (D_l * D_s);
        double theta_primary_raw = (beta + std::sqrt(beta * beta + 4.0 * theta_E_sq)) / 2.0;

        // Analytical static shadow boundary
        double bc = 3.0 * std::sqrt(3.0) * mass;
        double sin_alpha = (bc / r_obs) * std::sqrt(1.0 - (2.0 * mass) / r_obs);

        double shadow_angular_radius;
        if (r_obs < 3.0 * mass) {
            // Inside photon sphere, shadow covers > half the sky
            shadow_angular_radius = M_PI - std::asin(std::clamp(sin_alpha, 0.0, 1.0));
        } else {
            shadow_angular_radius = std::asin(std::clamp(sin_alpha, 0.0, 1.0));
        }
        double blend_margin = shadow_angular_radius * 0.5;
        
        double t_p = std::clamp((std::abs(theta_primary_raw) - shadow_angular_radius) / blend_margin, 0.0, 1.0);
        double weak_field_weight_p = t_p * t_p * (3.0 - 2.0 * t_p); 
        double strong_field_theta_p = shadow_angular_radius * 1.01;
        double theta_primary = strong_field_theta_p + weak_field_weight_p * (theta_primary_raw - strong_field_theta_p);

        double apparentShift_p = theta_primary - beta;
        double margin = 5.0 * mass; 
        double transitWeight = std::min(std::clamp(t_closest / margin, 0.0, 1.0), std::clamp((flatDistance - t_closest) / margin, 0.0, 1.0));
        apparentShift_p = (apparentShift_p * transitWeight) + (shiftSameSide * (1.0 - transitWeight));

        Eigen::AngleAxisd primaryRotation(apparentShift_p, rotationAxis);
        Eigen::Vector3d finalArrivalDir_p = primaryRotation * viewDir;

        // Primary GR deformations
        double cos_alpha_p = obs_rel.normalized().dot(tgt_rel.normalized());
        double backgroundWeight = std::clamp(1.0 - cos_alpha_p, 0.0, 1.0);
        double stretch_p = 1.0, squish_p = 1.0, ringFactor_p = 0.0;

        if (backgroundWeight > 0.001) {
            double raw_stretch_p = 1.0, raw_squish_p = 1.0, raw_ringFactor_p = 0.0;
            if (beta > 1e-5) {
                raw_stretch_p = std::abs((beta + apparentShift_p) / beta); 
                raw_squish_p = 1.0 / std::sqrt(1.0 + (theta_E_sq / (beta * beta)));
                raw_ringFactor_p = std::sqrt(theta_E_sq) / std::sqrt((beta * beta) + theta_E_sq);
            } else if (transitsBH) { 
                raw_stretch_p = 15.0; raw_squish_p = 0.05; raw_ringFactor_p = 1.0;
            }
            raw_squish_p = 0.05 + weak_field_weight_p * (raw_squish_p - 0.05);
            stretch_p = 1.0 + (raw_stretch_p - 1.0) * backgroundWeight;
            squish_p = 1.0 + (raw_squish_p - 1.0) * backgroundWeight;
            ringFactor_p = raw_ringFactor_p * backgroundWeight;
        }

        // SECONDARY IMAGE
        result.secondary.isValid = true;
        
        double theta_secondary_raw = (beta - std::sqrt(beta * beta + 4.0 * theta_E_sq)) / 2.0;
        double t_s = std::clamp((std::abs(theta_secondary_raw) - shadow_angular_radius) / blend_margin, 0.0, 1.0);
        double weak_field_weight_s = t_s * t_s * (3.0 - 2.0 * t_s); 
        double strong_field_theta_s = - (shadow_angular_radius * 1.01);
        double theta_secondary = strong_field_theta_s + weak_field_weight_s * (theta_secondary_raw - strong_field_theta_s);

        // Secondary time delay 
        double b_minus = r_obs * std::abs(theta_secondary); 
        double path_l = std::sqrt(r_obs * r_obs + b_minus * b_minus);
        double path_ls = std::sqrt(r_tgt * r_tgt + b_minus * b_minus);
        double secondaryGeometricDist = path_l + path_ls; 
        double secondaryShapiro = 2.0 * mass * std::log((r_obs + r_tgt + secondaryGeometricDist) / safeDenominator);
        result.secondary.travelTime = secondaryGeometricDist + secondaryShapiro;

        // Secondary spatial bending
        double apparentShift_s = ((theta_secondary - beta) * transitWeight) - (shiftSameSide * (1.0 - transitWeight)); 
        Eigen::AngleAxisd secondaryRotation(apparentShift_s, rotationAxis);
        Eigen::Vector3d finalArrivalDir_s = secondaryRotation * viewDir;

        // Secondary GR deformations
        double stretch_s = 1.0, squish_s = 1.0, ringFactor_s = 0.0;  
        if (beta > 1e-5) {
            stretch_s = std::abs((beta + apparentShift_s) / beta); 
            double weak_squish_s = 1.0 / std::sqrt(1.0 + (theta_E_sq / (beta * beta)));
            squish_s = 0.05 + weak_field_weight_s * (weak_squish_s - 0.05);
            ringFactor_s = std::sqrt(theta_E_sq) / std::sqrt((beta * beta) + theta_E_sq);
        } else if (transitsBH) {
            stretch_s = 15.0; squish_s = 0.05; ringFactor_s = 1.0;
        }
        
        // packaging
        result.primary.arrivalK = metric.enforceNullCondition(observerPos, -finalArrivalDir_p);
        result.primary.arcTangentAxis = rotationAxis;
        result.primary.tangentialStretch = std::clamp(stretch_p, 0.01, 15.0);
        result.primary.radialSquish = std::clamp(squish_p, 0.01, 1.0);
        result.primary.ringFactor = std::clamp(ringFactor_p, 0.0, 1.0);

        result.secondary.arrivalK = metric.enforceNullCondition(observerPos, -finalArrivalDir_s);
        result.secondary.arcTangentAxis = rotationAxis;
        result.secondary.tangentialStretch = std::clamp(stretch_s, 0.01, 15.0);
        result.secondary.radialSquish = std::clamp(squish_s, 0.01, 1.0);
        result.secondary.ringFactor = std::clamp(ringFactor_s, 0.0, 1.0);

        return result;
    }

} // namespace RelSpacetime

