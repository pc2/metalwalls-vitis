#!/usr/bin/env python3
import json
import sys

out_path = sys.argv[1]
n_ranks = int(sys.argv[2])

with open("metrics.json", mode="r") as metrics_file:
    metrics_json = json.load(metrics_file)

with open("charges.out", mode="r") as charges_file:
    charges = [float(line) for line in charges_file.readlines()
               if len(line) > 0 and not line.startswith("#")]

with open(out_path, mode="r") as out_file:
    out_json = json.load(out_file)

out_json.append({
    "n_ranks": n_ranks,
    "metrics": metrics_json,
    "charges": charges
})

with open(out_path, mode="w") as out_file:
    json.dump(out_json, out_file)
