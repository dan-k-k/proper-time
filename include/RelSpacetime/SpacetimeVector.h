// include/RelSpacetime/SpacetimeVector.h
#pragma once
#include <Eigen/Dense>

namespace RelSpacetime {

    // 4-vector (t, x, y, z). c, M = 1.
    class SpacetimeVector {
    public:
        SpacetimeVector();
        SpacetimeVector(double t, double x, double y, double z);

        double t() const { return data[0]; } // Assessors
        double x() const { return data[1]; }
        double y() const { return data[2]; }
        double z() const { return data[3]; }

        const Eigen::Vector4d& getEigenVector() const { return data; }

    private:
        Eigen::Vector4d data; 
    };

} // namespace RelSpacetime

