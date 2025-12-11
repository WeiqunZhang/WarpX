#!/usr/bin/env python3

# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check fields, particles, layouts, and efficiency after runtime load balancing."""

import math
import sys
from pathlib import Path

case = sys.argv[1]
path = Path("diags/reducedfiles")


def read_diagnostic(name):
    lines = (path / f"{name}.txt").read_text().splitlines()
    return lines[0].split(), [line.split() for line in lines[1:] if line.strip()]


_, energy = read_diagnostic("FE")
assert [int(row[0]) for row in energy] == list(range(5))
for row in energy:
    # No initial fields or particle velocities: these must remain exactly zero.
    assert all(float(value) == 0.0 for value in row[2:]), row

_, efficiency = read_diagnostic("LBE")
assert [int(row[0]) for row in efficiency] == list(range(5))
for row in efficiency:
    if int(row[0]) >= 2:
        assert all(0.0 < float(value) <= 1.0 for value in row[2:]), row
        if case == "default":
            assert float(row[2]) == 1.0, row

_, particles = read_diagnostic("PN")
assert [int(row[0]) for row in particles] == list(range(5))
particle_counts = [int(float(row[2])) for row in particles]
assert particle_counts[0] > 0
assert all(count == particle_counts[0] for count in particle_counts), particle_counts

header, costs = read_diagnostic("LBC")
# Locate each box using the header, allowing the extra GPU ID column on GPUs
# and NaN padding when the number of boxes changes during the run.
box_columns = [i for i, name in enumerate(header) if "cost_box_" in name]
layouts = {}
for row in costs:
    boxes = []
    for start in box_columns:
        if math.isfinite(float(row[start])):
            boxes.append(
                {
                    "level": int(float(row[start + 2])),
                    "lo": (int(float(row[start + 3])), int(float(row[start + 4]))),
                    "particles": int(float(row[start + 7])),
                }
            )
    layouts[int(row[0])] = boxes

assert sorted(layouts) == list(range(5)), layouts.keys()

before, after = layouts[1], layouts[2]
assert len(layouts[4]) == len(after), (
    "The second load balance changed the settled layout"
)

if case == "default":
    assert len(before) == len(after) == 4
elif case == "costs":
    # Initial constant-density injection populates only the left of two boxes.
    # Runtime balancing must split that expensive box and retain the empty one.
    assert len(before) == 2
    assert {box["lo"] for box in after} == {(0, 0), (0, 8), (16, 0)}, after
    assert all(box["particles"] > 0 for box in after if box["lo"][0] == 0)
    assert all(box["particles"] == 0 for box in after if box["lo"][0] == 16)
elif case == "merge":
    assert len(after) < len(before), (before, after)
elif case == "split":
    # A low threshold forces splitting until the default minimum of 8 cells.
    # Splitting at equality would produce 4-cell sides and additional origins.
    assert len(before) == 4
    assert {box["lo"] for box in after} == {
        (i, j) for i in range(0, 32, 8) for j in range(0, 32, 8)
    }, after
elif case == "refined":
    for level in (0, 1):
        old_boxes = [box for box in before if box["level"] == level]
        new_boxes = [box for box in after if box["level"] == level]
        assert old_boxes
        assert len(new_boxes) == 4 * len(old_boxes), (old_boxes, new_boxes)
else:
    raise ValueError(f"Unknown load-balance case: {case}")

print(f"{case}: fields, particles, layout, and efficiency checks passed")
