// src/SpacetimeVector.cpp
#include "RelSpacetime/SpacetimeVector.h"

namespace RelSpacetime {

    SpacetimeVector::SpacetimeVector() : data(Eigen::Vector4d::Zero()) {}

    SpacetimeVector::SpacetimeVector(double t, double x, double y, double z) {
        data << t, x, y, z;
    }

} // namespace RelSpacetime

