# Each test runs its cases and their analysis together.
find_package(Python COMPONENTS Interpreter REQUIRED)

add_warpx_test(
    test_2d_load_balance_split_merge  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_load_balance_split_merge  # inputs
    OFF  # analysis: included in the run wrapper
    OFF  # checksum: check field, particle, and layout invariants
    OFF  # dependency
    RUN_WRAPPER ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_test.py split_merge
)

add_warpx_test(
    test_2d_load_balance_costs  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_load_balance_costs  # inputs
    OFF  # analysis: included in the run wrapper
    OFF  # checksum: check field, particle, and layout invariants
    OFF  # dependency
    RUN_WRAPPER ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_test.py costs
)
