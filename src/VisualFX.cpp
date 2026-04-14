// src/VisualFX.cpp
#include "RelSpacetime/VisualFX.h"
#include <cmath>
#include "KerrShadow.inl"

namespace RelSpacetime {

    VisualRenderState VisualFX::computeVisuals(
        const EntityState& observerCurrentState,
        const LocalFrame& observerFrame,
        const EntityState& primaryPastState,
        const EntityState& secondaryPastState,
        const NullGeodesicResult& geodesicData, 
        const MetricTensor& metric) const 
    {
        VisualRenderState output;
        Eigen::Matrix4d g_obs = metric.getMetricAt(observerCurrentState.position);
        Eigen::Vector4d U_obs = observerFrame.getE0();

        output.localTetrad.col(0) = observerFrame.getE0();
        output.localTetrad.col(1) = observerFrame.getE1();
        output.localTetrad.col(2) = observerFrame.getE2();
        output.localTetrad.col(3) = observerFrame.getE3();

        // PRIMARY IMAGE
        output.primary.isVisible = geodesicData.primary.isValid;
        if (output.primary.isVisible) {
            Eigen::Vector4d k_obs_covar_p = g_obs * geodesicData.primary.arrivalK;
            double energy_obs_p = k_obs_covar_p.dot(U_obs);
            double energy_emit_p = k_obs_covar_p.dot(primaryPastState.fourVelocity); 
            
            output.primary.dopplerShift = energy_obs_p / energy_emit_p;
            output.primary.beamingFactor = std::pow(output.primary.dopplerShift, 3.0 + spectralIndex);

            // Project Apparent Direction
            double E_local_p = std::abs(energy_obs_p); 
            double app_x_p = k_obs_covar_p.dot(observerFrame.getE1());
            double app_y_p = k_obs_covar_p.dot(observerFrame.getE2());
            double app_z_p = k_obs_covar_p.dot(observerFrame.getE3());
            output.primary.apparentDirection = Eigen::Vector3d(app_x_p, app_y_p, app_z_p) / std::max(E_local_p, 1e-12);

            Eigen::Vector4d global_axis_p(0.0, 
                                          geodesicData.primary.arcTangentAxis.x(), 
                                          geodesicData.primary.arcTangentAxis.y(), 
                                          geodesicData.primary.arcTangentAxis.z());
            Eigen::Vector4d covar_axis_p = g_obs * global_axis_p;
            
            double axis_x_p = covar_axis_p.dot(observerFrame.getE1());
            double axis_y_p = covar_axis_p.dot(observerFrame.getE2());
            double axis_z_p = covar_axis_p.dot(observerFrame.getE3());
            
            Eigen::Vector3d local_axis_p(axis_x_p, axis_y_p, axis_z_p);
            output.primary.arcTangentAxis = local_axis_p.norm() > 1e-6 ? local_axis_p.normalized() : Eigen::Vector3d::Zero();

            output.primary.travelTime = geodesicData.primary.travelTime;
            output.primary.tangentialStretch = geodesicData.primary.tangentialStretch;
            output.primary.radialSquish = geodesicData.primary.radialSquish;
            output.primary.ringFactor = geodesicData.primary.ringFactor;
        }

        // SECONDARY IMAGE
        output.secondary.isVisible = geodesicData.secondary.isValid;
        if (output.secondary.isVisible) {
            Eigen::Vector4d k_obs_covar_s = g_obs * geodesicData.secondary.arrivalK;
            double energy_obs_s = k_obs_covar_s.dot(U_obs);
            double energy_emit_s = k_obs_covar_s.dot(secondaryPastState.fourVelocity); 
            
            output.secondary.dopplerShift = energy_obs_s / energy_emit_s;
            output.secondary.beamingFactor = std::pow(output.secondary.dopplerShift, 3.0 + spectralIndex);

            // Project Apparent Direction
            double E_local_s = std::abs(energy_obs_s);
            double app_x_s = k_obs_covar_s.dot(observerFrame.getE1());
            double app_y_s = k_obs_covar_s.dot(observerFrame.getE2());
            double app_z_s = k_obs_covar_s.dot(observerFrame.getE3());
            output.secondary.apparentDirection = Eigen::Vector3d(app_x_s, app_y_s, app_z_s) / std::max(E_local_s, 1e-12);

            Eigen::Vector4d global_axis_s(0.0, 
                                          geodesicData.secondary.arcTangentAxis.x(), 
                                          geodesicData.secondary.arcTangentAxis.y(), 
                                          geodesicData.secondary.arcTangentAxis.z());
            Eigen::Vector4d covar_axis_s = g_obs * global_axis_s;
            
            double axis_x_s = covar_axis_s.dot(observerFrame.getE1());
            double axis_y_s = covar_axis_s.dot(observerFrame.getE2());
            double axis_z_s = covar_axis_s.dot(observerFrame.getE3());
            
            Eigen::Vector3d local_axis_s(axis_x_s, axis_y_s, axis_z_s);
            output.secondary.arcTangentAxis = local_axis_s.norm() > 1e-6 ? local_axis_s.normalized() : Eigen::Vector3d::Zero();

            output.secondary.travelTime = geodesicData.secondary.travelTime;
            output.secondary.tangentialStretch = geodesicData.secondary.tangentialStretch;
            output.secondary.radialSquish = geodesicData.secondary.radialSquish;
            output.secondary.ringFactor = geodesicData.secondary.ringFactor;
        }

        return output;
    }

    BlackHoleRenderState VisualFX::computeBlackHoleVisuals(
        const EntityState& observerCurrentState,
        const LocalFrame& observerFrame,
        const MetricTensor& metric) const 
    {
        BlackHoleRenderState output;
        
        double mass = metric.getMass();
        double spin = 0.0; 
        
        double r_obs, theta_obs;
        metric.getBoyerLindquistCoords(observerCurrentState.position, r_obs, theta_obs);

        // Inside singularity 
        if (r_obs < 1e-6) {
            output.apparentDirection = Eigen::Vector3d(1.0, 0.0, 0.0); 
            output.shadowAngularRadius = M_PI; 
            output.squishAxis = Eigen::Vector3d::Zero();
            output.dopplerScale = 1.0;
            return output;
        }

        double a = spin;
        double a2 = a * a;
        
        double x = observerCurrentState.position.x();
        double y = observerCurrentState.position.y();
        double z = observerCurrentState.position.z();
        
        if (x*x + y*y < 1e-12) {
            x = 1e-5; // for multiplication that fails along polar coords
            
            double R2 = x*x + y*y + z*z;
            double half_term = (R2 - a2) * 0.5;
            double r2_new = half_term + std::sqrt(half_term * half_term + a2 * z * z);
            r_obs = std::sqrt(r2_new);
            theta_obs = std::acos(std::clamp(z / r_obs, -1.0, 1.0));
        }

        double rho_cyl2 = x*x + y*y;
        double r2 = r_obs * r_obs;
        double D = (r_obs * rho_cyl2) / ((r2 + a2) * (r2 + a2)) + (z * z) / (r_obs * r2);
        
        Eigen::Vector3d grad_r(
            x / (D * (r2 + a2)),
            y / (D * (r2 + a2)),
            z / (D * r2)
        );
        
        double sin_theta = std::sin(theta_obs); 
        
        Eigen::Vector3d grad_theta(
            (z * grad_r.x()) / (r2 * sin_theta),
            (z * grad_r.y()) / (r2 * sin_theta),
            -1.0 / (r_obs * sin_theta) + (z * grad_r.z()) / (r2 * sin_theta)
        );
        
        Eigen::Vector3d grad_phi(
            -y / rho_cyl2 - (a * grad_r.x()) / (r2 + a2),
            x / rho_cyl2 - (a * grad_r.y()) / (r2 + a2),
            -(a * grad_r.z()) / (r2 + a2)
        );

        auto projectToSky = [&](double k_t, double k_r, double k_theta, double k_phi) -> Eigen::Vector3d {
            double k_x = k_r * grad_r.x() + k_theta * grad_theta.x() + k_phi * grad_phi.x();
            double k_y = k_r * grad_r.y() + k_theta * grad_theta.y() + k_phi * grad_phi.y();
            double k_z = k_r * grad_r.z() + k_theta * grad_theta.z() + k_phi * grad_phi.z();
            
            Eigen::Vector4d k_covar_cart(k_t, k_x, k_y, k_z);
            
            double E_local = -k_covar_cart.dot(observerFrame.getE0());
            double k_1 = k_covar_cart.dot(observerFrame.getE1());
            double k_2 = k_covar_cart.dot(observerFrame.getE2());
            double k_3 = k_covar_cart.dot(observerFrame.getE3());
            
            return Eigen::Vector3d(k_1, k_2, k_3) / std::max(E_local, 1e-12);
        };

        // The BH centre
        // To find where the BH visually appears, we must track a photon ARRIVING from it.
        // An arriving photon from the BH center is moving OUTWARDS (k_r = +1.0).
        // projectToSky returns the arriving photon's momentum. We look in the exact opposite direction.
        Eigen::Vector3d view_center = projectToSky(-1.0, -1.0, 0.0, 0.0).normalized();
        output.apparentDirection = view_center;

        // The shadow edge params (Polar Ray)
        // Using the "top" vertical edge of the shadow (lambda = 0) avoids the equatorial corotating compression.
        // It provides the true, large ~5.2M impact parameter radius.
        double p_cub = a2 - 3.0 * mass * mass;
        double q_cub = 2.0 * mass * (a2 - mass * mass);
        
        double r_peak;
        if (std::abs(a) < 1e-5) {
            r_peak = 3.0 * mass;
        } else {
            // polar photon sphere cubic equation: r^3 - 3Mr^2 + a^2r + a^2M = 0
            double sqrt_minus_p_3 = std::sqrt(std::max(0.0, -p_cub / 3.0));
            double acos_arg = (3.0 * q_cub / (2.0 * p_cub)) * (1.0 / std::max(sqrt_minus_p_3, 1e-12));
            acos_arg = std::clamp(acos_arg, -1.0, 1.0);
            double x_root = 2.0 * sqrt_minus_p_3 * std::cos(std::acos(acos_arg) / 3.0);
            r_peak = x_root + mass;
        }

        double r_peak2 = r_peak * r_peak;
        double delta_peak = r_peak2 - 2.0 * mass * r_peak + a2;
        
        double eta_c = std::pow(r_peak2 + a2, 2.0) / std::max(delta_peak, 1e-12);
        double lambda_c = 0.0;

        // Stable k_r for the edge ray
        double delta_obs = r2 - 2.0 * mass * r_obs + a2;
        double X = r2 + a2; // lambda_c = 0
        double K_val = eta_c + a2; 
        
        double R_val = X * X - delta_obs * K_val;
        R_val = std::max(0.0, R_val);
        
        // Outside the photon sphere, arriving edge photons are OUTGOING (+1).
        // Inside the photon sphere, arriving edge photons are INGOING (-1).
        // never divides by 0 at the event horizon.
        double sigma = (r_obs >= r_peak) ? 1.0 : -1.0;
        
        double k_r_edge = (K_val / (X - sigma * std::sqrt(R_val))) - 1.0;
        double k_theta_edge = std::sqrt(std::max(0.0, eta_c + a2 * std::pow(std::cos(theta_obs), 2.0)));
        double k_phi_edge = 0.0;
        // The edge ray is an ARRIVING photon, look in the opposite direction of its momentum
        Eigen::Vector3d view_edge = -projectToSky(-1.0, k_r_edge, k_theta_edge, k_phi_edge).normalized();

        double cos_alpha = std::clamp(view_center.dot(view_edge), -1.0, 1.0);
        output.shadowAngularRadius = std::acos(cos_alpha);
        
        output.dopplerScale = 1.0; // legacy variable; aberration resolved inside projectToSky

        // Squish axis
        Eigen::Matrix4d g_obs_matrix = metric.getMetricAt(observerCurrentState.position);
        Eigen::Vector4d global_spin_axis(0.0, 0.0, 0.0, 1.0);
        Eigen::Vector4d covar_spin_axis = g_obs_matrix * global_spin_axis;
        
        Eigen::Vector3d screen_squish(
            covar_spin_axis.dot(observerFrame.getE1()),
            covar_spin_axis.dot(observerFrame.getE2()),
            covar_spin_axis.dot(observerFrame.getE3())
        );
        
        screen_squish -= screen_squish.dot(output.apparentDirection) * output.apparentDirection;
        output.squishAxis = screen_squish.norm() > 1e-6 ? screen_squish.normalized() : Eigen::Vector3d::UnitY();

        return output;
    }
} // namespace RelSpacetime

