// src_godot/GDEntityState.h
#pragma once
#include <godot_cpp/classes/ref_counted.hpp>
#include "RelSpacetime/RelativityTypes.h" 

using namespace godot;

class GDEntityState : public RefCounted {
    GDCLASS(GDEntityState, RefCounted)

private:
    RelSpacetime::EntityState core_state; 

protected:
    static void _bind_methods(); 

public:
    GDEntityState() {}
    
    // Position
    Vector4 get_position() const; 
    void set_position(const Vector4& pos);
    
    // Four-Velocity
    Vector4 get_four_velocity() const;
    void set_four_velocity(const Vector4& vel);

    // Times
    double get_proper_time() const;
    void set_proper_time(double t);

    double get_coordinate_time() const;
    void set_coordinate_time(double t);

    double get_time_scale() const;
    void set_time_scale(double s);

    // C++ Interop
    const RelSpacetime::EntityState& get_core_state() const { return core_state; }
    void set_core_state(const RelSpacetime::EntityState& state) { core_state = state; }
};

