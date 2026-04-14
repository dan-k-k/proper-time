# tests/CMakeLists.txt
# Black Hole Tests
add_executable(RunBlackHole test_blackhole.cpp)
target_link_libraries(RunBlackHole PRIVATE RelSpacetimeCore)
add_test(NAME BlackHoleTest COMMAND RunBlackHole)

# Integrator Tests
add_executable(RunIntegratorTests test_integrator.cpp)
target_link_libraries(RunIntegratorTests PRIVATE RelSpacetimeCore)
add_test(NAME IntegratorTest COMMAND RunIntegratorTests)

# Kerr Tests
add_executable(RunKerrTests test_kerr.cpp)
target_link_libraries(RunKerrTests PRIVATE RelSpacetimeCore)
add_test(NAME KerrTest COMMAND RunKerrTests)

# Metrics Tests
add_executable(RunMetricsTests test_metrics.cpp)
target_link_libraries(RunMetricsTests PRIVATE RelSpacetimeCore)
add_test(NAME MetricsTest COMMAND RunMetricsTests)

# Optics Tests
add_executable(RunOpticsTests test_optics.cpp)
target_link_libraries(RunOpticsTests PRIVATE RelSpacetimeCore)
add_test(NAME OpticsTest COMMAND RunOpticsTests)

# Shooting Tests
add_executable(RunShootingTests test_shooting.cpp)
target_link_libraries(RunShootingTests PRIVATE RelSpacetimeCore)
add_test(NAME ShootingTest COMMAND RunShootingTests)

# Solver Tests
add_executable(RunSolverTests test_solver.cpp)
target_link_libraries(RunSolverTests PRIVATE RelSpacetimeCore)
add_test(NAME SolverTest COMMAND RunSolverTests)

# Tetrad Tests
add_executable(RunTetradTests test_tetrad.cpp)
target_link_libraries(RunTetradTests PRIVATE RelSpacetimeCore)
add_test(NAME TetradTest COMMAND RunTetradTests)

# Transform Tests
add_executable(RunTransformTests test_transform.cpp)
target_link_libraries(RunTransformTests PRIVATE RelSpacetimeCore)
add_test(NAME TransformTest COMMAND RunTransformTests)

# Integration Tests
add_executable(RunIntegrationTests test_integration.cpp)
target_link_libraries(RunIntegrationTests PRIVATE RelSpacetimeCore)
add_test(NAME IntegrationTest COMMAND RunIntegrationTests)

