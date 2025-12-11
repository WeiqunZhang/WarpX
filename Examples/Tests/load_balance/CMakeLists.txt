# Heuristic costs and stationary particles keep these MPI tests deterministic.
foreach(case IN ITEMS costs default merge refined split)
    add_warpx_test(
        test_2d_load_balance_${case}  # name
        2  # dims
        2  # nprocs
        inputs_test_2d_load_balance_${case}  # inputs
        "analysis.py ${case}"  # analysis
        OFF  # checksum: check exact field, particle, and layout invariants
        OFF  # dependency
    )
endforeach()
