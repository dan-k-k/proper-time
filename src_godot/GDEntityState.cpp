// src_godot/GDEntityState.cpp
#include "GDEntityState.h"
#include <godot_cpp/core/class_db.hpp>

using namespace godot;

void GDEntityState::_bind_methods() {
    ClassDB::bind_method(D_METHOD("get_position"), &GDEntityState::get_position);
    ClassDB::bind_method(D_METHOD("set_position", "pos"), &GDEntityState::set_position);

    ClassDB::bind_method(D_METHOD("get_four_velocity"), &GDEntityState::get_four_velocity);
    ClassDB::bind_method(D_METHOD("set_four_velocity", "vel"), &GDEntityState::set_four_velocity);

    ClassDB::bind_method(D_METHOD("get_proper_time"), &GDEntityState::get_proper_time);
    ClassDB::bind_method(D_METHOD("set_proper_time", "t"), &GDEntityState::set_proper_time);

    ClassDB::bind_method(D_METHOD("get_coordinate_time"), &GDEntityState::get_coordinate_time);
    ClassDB::bind_method(D_METHOD("set_coordinate_time", "t"), &GDEntityState::set_coordinate_time);

    ClassDB::bind_method(D_METHOD("get_time_scale"), &GDEntityState::get_time_scale);
    ClassDB::bind_method(D_METHOD("set_time_scale", "s"), &GDEntityState::set_time_scale);

    ADD_PROPERTY(PropertyInfo(Variant::VECTOR4, "position"), "set_position", "get_position");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR4, "four_velocity"), "set_four_velocity", "get_four_velocity");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "proper_time"), "set_proper_time", "get_proper_time");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "coordinate_time"), "set_coordinate_time", "get_coordinate_time");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "time_scale"), "set_time_scale", "get_time_scale");
}

Vector4 GDEntityState::get_position() const {
    const auto& x = core_state.position;
    return Vector4(x.t(), x.x(), x.y(), x.z()); // (t, x, y, z) 
}

void GDEntityState::set_position(const Vector4& pos) {
    core_state.position = RelSpacetime::SpacetimeVector(pos.x, pos.y, pos.z, pos.w); // x=t, y=x, z=y, w=z mapped to SpacetimeVector
}

Vector4 GDEntityState::get_four_velocity() const {
    const auto& v = core_state.fourVelocity;
    return Vector4(v[0], v[1], v[2], v[3]);
}

void GDEntityState::set_four_velocity(const Vector4& vel) {
    core_state.fourVelocity = Eigen::Vector4d(vel.x, vel.y, vel.z, vel.w);
}

double GDEntityState::get_proper_time() const { return core_state.properTime; }
void GDEntityState::set_proper_time(double t) { core_state.properTime = t; }

double GDEntityState::get_coordinate_time() const { return core_state.coordinateTime; }
void GDEntityState::set_coordinate_time(double t) { core_state.coordinateTime = t; }

double GDEntityState::get_time_scale() const { return core_state.timeScale; }
void GDEntityState::set_time_scale(double s) { core_state.timeScale = s; }

