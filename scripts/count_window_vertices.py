#!/usr/bin/env python3
"""Vertex counts of the displacement-histogram polytopes of the centered five-point window (offsets -2..2).

For m = 1..10 (default; --mmax to change), enumerate the displacement histograms (x_{-2}, ..., x_2) of all
permutations tau of {0..m-1} with |i - tau(i)| <= 2 and count the vertices of their convex hull exactly (integer
arithmetic in the 3-dimensional coordinates (x_{-2}, x_{-1}, x_1)). Expected 1, 2, 5, 5, 9, 6, 9, 9, 9, 10.
"""
import argparse

from common import Timer, finish, out_dir, report

from dsxray.windows import centered_window_histograms, five_point_window_vertex_count

EXPECTED = [1, 2, 5, 5, 9, 6, 9, 9, 9, 10]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/window_vertices")
    ap.add_argument("--mmax", type=int, default=10)
    ap.add_argument("--full", action="store_true", help="same as the default")
    a = ap.parse_args()
    t = Timer()
    counts = [five_point_window_vertex_count(m) for m in range(1, a.mmax + 1)]
    hist = [len(centered_window_histograms(m)) for m in range(1, a.mmax + 1)]
    res = {"vertex_counts": counts, "distinct_histograms": hist}
    res["checks"] = {"counts": report("five-point window vertex counts", counts[:10] == EXPECTED[:len(counts)],
                                      f"{counts} ({t()} s)")}
    finish(res, out_dir(a.out), "window_vertices")


if __name__ == "__main__":
    main()
