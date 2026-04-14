// src_godot/register_types.cpp
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/godot.hpp>

#include "GDEntityState.h"
#include "GDRelativityAPI.h"

using namespace godot;

void initialize_relativity_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
    }

    ClassDB::register_class<GDEntityState>();
    ClassDB::register_class<GDRelativityAPI>();
}

void uninitialize_relativity_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
    }
}

extern "C" {

GDExtensionBool GDE_EXPORT rel_spacetime_library_init(
    GDExtensionInterfaceGetProcAddress p_get_proc_address,
    const GDExtensionClassLibraryPtr p_library,
    GDExtensionInitialization *r_initialization
) {
    godot::GDExtensionBinding::InitObject init_obj(p_get_proc_address, p_library, r_initialization);

    init_obj.register_initializer(initialize_relativity_module);
    init_obj.register_terminator(uninitialize_relativity_module);
    init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SCENE);

    return init_obj.init();
}

}

