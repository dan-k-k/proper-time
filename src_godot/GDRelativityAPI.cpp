#include "GDRelativityAPI.h"
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/projection.hpp>

using namespace godot;
using namespace RelSpacetime;

void GDRelativityAPI::_bind_methods() {
    ClassDB::bind_method(D_METHOD("initialize_universe", "universe_params"), &GDRelativityAPI::initialize_universe);

    ClassDB::bind_method(
        D_METHOD("compute_next_state", "current_state", "engine_inputs", "server_dt"),
        &GDRelativityAPI::compute_next_state
    );

    ClassDB::bind_method(
        D_METHOD("get_visual_state_for_opponent", "observer_state", "target_history", "previous_travel_time"),
        &GDRelativityAPI::get_visual_state_for_opponent,
        DEFVAL(-1.0)
    );

    ClassDB::bind_method(
        D_METHOD("get_visual_state_for_black_hole", "observer_state"),
        &GDRelativityAPI::get_visual_state_for_black_hole
    );
}

GDRelativityAPI::GDRelativityAPI() {
    core_api = memnew(RelativityVisualsAPI(0.0));
    black_hole_position = Vector3(0, 0, 0); // Default
}

GDRelativityAPI::~GDRelativityAPI() {
    memdelete(core_api);
}

void GDRelativityAPI::initialize_universe(Dictionary universe_params) {
    double mass = universe_params.get("mass", 0.0);
    double spin = universe_params.get("spin", 0.0);
    black_hole_position = universe_params.get("position", Vector3(0, 0, 0));

    active_metric = std::make_unique<KerrMetric>(mass, spin);
}

Ref<GDEntityState> GDRelativityAPI::compute_next_state(
    Ref<GDEntityState> current_state,
    Dictionary engine_inputs,
    double server_dt
) {
    if (current_state.is_null() || !active_metric) {
        return current_state; 
    }

    EngineInputs inputs;
    Vector3 thrust = engine_inputs.get("local_proper_thrust", Vector3(0,0,0));
    Vector3 fwd = engine_inputs.get("forward", Vector3(1,0,0));
    Vector3 rgt = engine_inputs.get("right", Vector3(0,1,0));
    Vector3 up = engine_inputs.get("up", Vector3(0,0,1));

    inputs.localProperThrust = Eigen::Vector3d(thrust.x, thrust.y, thrust.z);
    inputs.forward = Eigen::Vector3d(fwd.x, fwd.y, fwd.z);
    inputs.right = Eigen::Vector3d(rgt.x, rgt.y, rgt.z);
    inputs.up = Eigen::Vector3d(up.x, up.y, up.z);

    EntityState core_state = current_state->get_core_state();

    // Godot into metric space
    double t = core_state.position.t();
    double x_metric = core_state.position.x() - black_hole_position.x;
    double y_metric = core_state.position.y() - black_hole_position.y;
    double z_metric = core_state.position.z() - black_hole_position.z;
    core_state.position = RelSpacetime::SpacetimeVector(t, x_metric, y_metric, z_metric);

    double r_squared = (x_metric * x_metric) + (y_metric * y_metric) + (z_metric * z_metric);
    if (r_squared < 0.01) { // divide by zero issues
        return current_state; // crushed
    }

    EntityState next_core = WorldlineIntegrator::computeNextState(
        core_state, 
        inputs, 
        server_dt, 
        *active_metric
    );

    // back into Godot Space
    double t_new = next_core.position.t();
    double x_godot = next_core.position.x() + black_hole_position.x;
    double y_godot = next_core.position.y() + black_hole_position.y;
    double z_godot = next_core.position.z() + black_hole_position.z;
    next_core.position = RelSpacetime::SpacetimeVector(t_new, x_godot, y_godot, z_godot);

    Ref<GDEntityState> next_state;
    next_state.instantiate();
    next_state->set_core_state(next_core);
    
    return next_state;
}

Dictionary GDRelativityAPI::get_visual_state_for_opponent(
    Ref<GDEntityState> observer_state,
    TypedArray<GDEntityState> target_history, 
    double previous_travel_time
) {
    Dictionary result;
    Dictionary primary_dict;
    Dictionary secondary_dict;

    if (observer_state.is_null() || target_history.size() == 0 || !active_metric) {
        primary_dict["is_visible"] = false;
        secondary_dict["is_visible"] = false;
        result["primary"] = primary_dict;
        result["secondary"] = secondary_dict;
        return result;
    }

    // Godot array into C++ std::vector
    std::vector<EntityState> history_buffer;
    history_buffer.reserve(target_history.size());
    for (int i = 0; i < target_history.size(); ++i) {
        Ref<GDEntityState> gd_state = target_history[i];
        if (gd_state.is_valid()) {
            EntityState state = gd_state->get_core_state();
            state.position = RelSpacetime::SpacetimeVector(
                state.position.t(),
                state.position.x() - black_hole_position.x,
                state.position.y() - black_hole_position.y,
                state.position.z() - black_hole_position.z
            );
            history_buffer.push_back(state);
        }
    }

    // Translate Observer
    EntityState observer = observer_state->get_core_state();
    observer.position = RelSpacetime::SpacetimeVector(
        observer.position.t(),
        observer.position.x() - black_hole_position.x,
        observer.position.y() - black_hole_position.y,
        observer.position.z() - black_hole_position.z
    );

    // Camera-independent WORLD FRAME 
    LocalFrame world_frame(
        active_metric->getMetricAt(observer.position), 
        observer.fourVelocity
    );

    VisualRenderState vis = core_api->getVisualStateForOpponent(
        observer, world_frame, history_buffer, *active_metric, previous_travel_time
    );

    // Package for .gd
    Projection godot_tetrad;
    godot_tetrad.columns[0] = Vector4(vis.localTetrad(0,0), vis.localTetrad(1,0), vis.localTetrad(2,0), vis.localTetrad(3,0));
    godot_tetrad.columns[1] = Vector4(vis.localTetrad(0,1), vis.localTetrad(1,1), vis.localTetrad(2,1), vis.localTetrad(3,1));
    godot_tetrad.columns[2] = Vector4(vis.localTetrad(0,2), vis.localTetrad(1,2), vis.localTetrad(2,2), vis.localTetrad(3,2));
    godot_tetrad.columns[3] = Vector4(vis.localTetrad(0,3), vis.localTetrad(1,3), vis.localTetrad(2,3), vis.localTetrad(3,3));

    // Package Primary Image
    primary_dict["is_visible"] = vis.primary.isVisible;
    primary_dict["direction"] = Vector3(vis.primary.apparentDirection.x(), vis.primary.apparentDirection.y(), vis.primary.apparentDirection.z());
    primary_dict["travel_time"] = vis.primary.travelTime; 
    primary_dict["doppler"] = vis.primary.dopplerShift; 
    primary_dict["intensity"] = vis.primary.beamingFactor; 
    primary_dict["arc_axis"] = Vector3(vis.primary.arcTangentAxis.x(), vis.primary.arcTangentAxis.y(), vis.primary.arcTangentAxis.z());
    primary_dict["stretch"] = vis.primary.tangentialStretch;
    primary_dict["squish"] = vis.primary.radialSquish;
    primary_dict["ring_factor"] = vis.primary.ringFactor;

    // Package Secondary Image
    secondary_dict["is_visible"] = vis.secondary.isVisible;
    secondary_dict["direction"] = Vector3(vis.secondary.apparentDirection.x(), vis.secondary.apparentDirection.y(), vis.secondary.apparentDirection.z());
    secondary_dict["travel_time"] = vis.secondary.travelTime; 
    secondary_dict["doppler"] = vis.secondary.dopplerShift; 
    secondary_dict["intensity"] = vis.secondary.beamingFactor; 
    secondary_dict["arc_axis"] = Vector3(vis.secondary.arcTangentAxis.x(), vis.secondary.arcTangentAxis.y(), vis.secondary.arcTangentAxis.z());
    secondary_dict["stretch"] = vis.secondary.tangentialStretch;
    secondary_dict["squish"] = vis.secondary.radialSquish;
    secondary_dict["ring_factor"] = vis.secondary.ringFactor;

    // Final
    result["primary"] = primary_dict;
    result["secondary"] = secondary_dict;
    result["local_tetrad"] = godot_tetrad;
    
    return result;
}

Dictionary GDRelativityAPI::get_visual_state_for_black_hole(
    Ref<GDEntityState> observer_state
) {
    Dictionary result;

    if (observer_state.is_null() || !active_metric) {
        result["is_visible"] = false;
        return result;
    }

    // Observer into Metric space
    EntityState observer = observer_state->get_core_state();
    observer.position = RelSpacetime::SpacetimeVector(
        observer.position.t(),
        observer.position.x() - black_hole_position.x,
        observer.position.y() - black_hole_position.y,
        observer.position.z() - black_hole_position.z
    );

    // Observer frame using metric and 4-velocity
    LocalFrame observer_frame(
        active_metric->getMetricAt(observer.position), 
        observer.fourVelocity
    );

    BlackHoleRenderState bh_state = core_api->getBlackHoleVisualState(
        observer, observer_frame, *active_metric
    );

    // Package
    result["is_visible"] = true; 
    result["direction"] = Vector3(bh_state.apparentDirection.x(), bh_state.apparentDirection.y(), bh_state.apparentDirection.z());
    result["angular_radius"] = bh_state.shadowAngularRadius;
    result["doppler_scale"] = bh_state.dopplerScale;
    result["squish_axis"] = Vector3(bh_state.squishAxis.x(), bh_state.squishAxis.y(), bh_state.squishAxis.z());

    return result;
}

