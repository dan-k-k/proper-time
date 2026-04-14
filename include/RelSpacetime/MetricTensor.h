// include/RelSpacetime/MetricTensor.h
#pragma once
#include <Eigen/Dense>
#include <vector> 
#include "SpacetimeVector.h"
#include "RelativityTypes.h" 

namespace RelSpacetime {

class MetricTensor {
public:
    virtual ~MetricTensor() = default;
    virtual double getMass() const = 0;
    virtual Eigen::Matrix4d getMetricAt(const SpacetimeVector& pos) const = 0;
    Eigen::Matrix4d getInverseMetricAt(const SpacetimeVector& pos) const { return getMetricAt(pos).inverse(); }
    
    Eigen::Vector4d enforceNullCondition(const SpacetimeVector& pos, const Eigen::Vector3d& k_spatial) const;

    virtual std::vector<Eigen::Matrix4d> getChristoffelSymbols(const SpacetimeVector& pos) const;
    double calculateIntervalSq(const SpacetimeVector& pos, const SpacetimeVector& nextPos) const;
    virtual void getBoyerLindquistCoords(const SpacetimeVector& pos, double& r, double& theta) const {
        r = std::sqrt(pos.x()*pos.x() + pos.y()*pos.y() + pos.z()*pos.z());
        theta = (r > 1e-6) ? std::acos(std::clamp(pos.z() / r, -1.0, 1.0)) : (M_PI / 2.0);
    }
};

class KerrMetric : public MetricTensor {
    double mass;
    double spin; // 'a' parameter
public:
    double getMass() const override { return mass; }
    KerrMetric(double M, double a) : mass(M), spin(a) {}
    Eigen::Matrix4d getMetricAt(const SpacetimeVector& pos) const override; 
    std::vector<Eigen::Matrix4d> getChristoffelSymbols(const SpacetimeVector& pos) const override;
    void getBoyerLindquistCoords(const SpacetimeVector& pos, double& r, double& theta) const override;
};

class MinkowskiMetric : public MetricTensor {
private:
    Eigen::Matrix4d eta;
public:
    double getMass() const override { return 0.0; }
    MinkowskiMetric() {
        eta = Eigen::Matrix4d::Zero();
        eta(0, 0) = -1.0;
        eta(1, 1) = 1.0;
        eta(2, 2) = 1.0;
        eta(3, 3) = 1.0;
    }

    Eigen::Matrix4d getMetricAt(const SpacetimeVector& pos) const override {
        return eta; 
    }
};

}

