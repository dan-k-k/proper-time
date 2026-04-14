// src_godot/GDRelativityAPI.h
#pragma once
#include <godot_cpp/classes/node.hpp>
#include <memory>
#include "RelSpacetime/RelativityVisualsAPI.h"
#include "RelSpacetime/WorldlineIntegrator.h"
#include "RelSpacetime/MetricTensor.h" 
#include "GDEntityState.h" 

using namespace godot;

class GDRelativityAPI : public Node {
    GDCLASS(GDRelativityAPI, Node)

private:
    RelSpacetime::RelativityVisualsAPI* core_api;
    std::unique_ptr<RelSpacetime::MetricTensor> active_metric;
    Vector3 black_hole_position;

protected:
    static void _bind_methods();

public:
    GDRelativityAPI();
    ~GDRelativityAPI();

    void initialize_universe(Dictionary universe_params);

    Ref<GDEntityState> compute_next_state(
        Ref<GDEntityState> current_state,
        Dictionary engine_inputs,
        double server_dt
    );

    Dictionary get_visual_state_for_opponent(
        Ref<GDEntityState> observer_state,
        TypedArray<GDEntityState> target_history,
        double previous_travel_time = -1.0 
    );

    Dictionary get_visual_state_for_black_hole(
        Ref<GDEntityState> observer_state
    );
};

