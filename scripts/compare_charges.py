#!/usr/bin/env python3
from matplotlib import pyplot
import sys


def parse_charges(path):
    with open(path) as file:
        for line in file.readlines():
            try:
                yield float(line)
            except ValueError:
                pass


diffs = [abs(a - b) for (a, b) in zip(parse_charges(sys.argv[1]),
                                      parse_charges(sys.argv[2]))]
print(
    f"min error: {min(diffs)}, max error: {max(diffs)}, mean error: {sum(diffs) / len(diffs)}")
if len(sys.argv) >= 4:
    pyplot.scatter(range(len(diffs)), diffs, 0.1)
    pyplot.savefig(sys.argv[3])
if max(diffs) > 1e-9:
    exit(1)
else:
    exit(0)
