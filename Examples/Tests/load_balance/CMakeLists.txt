# Heuristic costs and stationary particles keep these MPI tests deterministic.
# Use literal names so the source CI checker can find every test input.
add_warpx_test(
    test_2d_load_balance_costs  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_load_balance_costs  # inputs
    "analysis.py costs"  # analysis
    OFF  # checksum: check exact field, particle, and layout invariants
    OFF  # dependency
)

add_warpx_test(
    test_2d_load_balance_default  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_load_balance_default  # inputs
    "analysis.py default"  # analysis
    OFF  # checksum: check exact field, particle, and layout invariants
    OFF  # dependency
)

add_warpx_test(
    test_2d_load_balance_merge  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_load_balance_merge  # inputs
    "analysis.py merge"  # analysis
    OFF  # checksum: check exact field, particle, and layout invariants
    OFF  # dependency
)

add_warpx_test(
    test_2d_load_balance_refined  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_load_balance_refined  # inputs
    "analysis.py refined"  # analysis
    OFF  # checksum: check exact field, particle, and layout invariants
    OFF  # dependency
)

add_warpx_test(
    test_2d_load_balance_split  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_load_balance_split  # inputs
    "analysis.py split"  # analysis
    OFF  # checksum: check exact field, particle, and layout invariants
    OFF  # dependency
)
