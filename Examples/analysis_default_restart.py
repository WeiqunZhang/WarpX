#!/usr/bin/env python3

import argparse
import os

import numpy as np
import yt


def check_restart(filename, tolerance=1e-12):
    """
    Compare output data generated from initial run with output data generated after restart.

    Parameters
    ----------
    filename : str
        Name of the plotfile containing the output data generated after restart.
    tolerance : float, optional (default = 1e-12)
        Relative error between restart and original data must be smaller than tolerance.
    """
    # Load output data generated after restart
    ds_restart = yt.load(filename)

    # Load output data generated from initial run
    benchmark = os.path.join(os.getcwd().replace("_restart", ""), filename)
    ds_benchmark = yt.load(benchmark)

    # yt 4.0+ has rounding issues with our domain data:
    # RuntimeError: yt attempted to read outside the boundaries
    # of a non-periodic domain along dimension 0.
    for ds in (ds_restart, ds_benchmark):
        if "force_periodicity" in dir(ds):
            ds.force_periodicity()

    # Compare the data on every mesh refinement level
    for lev in range(ds_benchmark.max_level + 1):
        print(f"\nlevel {lev}")
        dims = ds_benchmark.domain_dimensions * ds_benchmark.refine_by**lev
        ad_restart = ds_restart.covering_grid(
            level=lev, left_edge=ds_restart.domain_left_edge, dims=dims
        )
        ad_benchmark = ds_benchmark.covering_grid(
            level=lev, left_edge=ds_benchmark.domain_left_edge, dims=dims
        )
        compare_fields(ad_restart, ad_benchmark, ds_benchmark.field_list, tolerance)


def compare_fields(ad_restart, ad_benchmark, field_list, tolerance):
    """
    Compare all fields of two yt data containers.

    Parameters
    ----------
    ad_restart : yt data container (e.g. covering grid) of the restarted run
    ad_benchmark : yt data container of the initial run
    field_list : list of fields to compare
    tolerance : float
        Relative error between restart and original data must be smaller than tolerance.
    """
    # Loop over all fields (all particle species, all particle attributes, all grid fields)
    # and compare output data generated from initial run with output data generated after restart
    print(f"\ntolerance = {tolerance}")
    print()
    for field in field_list:
        dr = ad_restart[field].squeeze().v
        db = ad_benchmark[field].squeeze().v
        assert dr.shape == db.shape
        if db.size == 0:
            print(f"field: {field}; no data")
            continue
        error = np.amax(np.abs(dr - db))
        if np.amax(np.abs(db)) != 0.0:
            error /= np.amax(np.abs(db))
        print(f"field: {field}; error = {error}")
        assert error < tolerance
    print()


if __name__ == "__main__":
    # define parser
    parser = argparse.ArgumentParser()

    # add arguments: output file path
    parser.add_argument(
        "--path",
        help="path to output file",
        type=str,
        required=True,
    )

    # add arguments: relative tolerance
    default_tolerance = 1e-12
    parser.add_argument(
        "--rtol",
        help="relative tolerance between restart and original",
        type=float,
        required=False,
        default=default_tolerance,
    )

    # parse arguments
    args = parser.parse_args()

    # compare restart results against original results
    check_restart(filename=args.path, tolerance=args.rtol)
