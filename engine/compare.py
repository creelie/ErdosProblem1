"""compare.py -- compares a fresh engine run with its record (engine/anaxia_regions.json or
engine/anaxia_centre.json).

Usage: python3 compare.py <fresh.json> <record.json>
The verdict, the reason and the measurement of every composite, and the summary fields other than
the memory statistics, must agree; floating-point numbers to a relative 1e-9.
"""
import json
import math
import sys


def same(x, y):
    if isinstance(x, float) or isinstance(y, float):
        return (isinstance(x, (int, float)) and isinstance(y, (int, float))
                and math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-12))
    if isinstance(x, dict) and isinstance(y, dict):
        return x.keys() == y.keys() and all(same(x[k], y[k]) for k in x)
    if isinstance(x, list) and isinstance(y, list):
        return len(x) == len(y) and all(same(u, v) for u, v in zip(x, y))
    return x == y


fresh, record = (json.load(open(f)) for f in sys.argv[1:3])
keys = [k for k in record["summary"] if k != "memory"]
bad = [k for k in keys if not same(fresh["summary"].get(k), record["summary"].get(k))]
a = {c["composite"]: c for c in fresh["composites"]}
b = {c["composite"]: c for c in record["composites"]}
if set(a) != set(b):
    bad.append("composite identifiers")
else:
    bad += [i for i in sorted(a) if not all(same(a[i][f], b[i][f]) for f in ("parts", "outcome", "reason", "measurement"))]
if bad:
    raise SystemExit(f"engine run differs from the record: {bad[:10]}")
print(f"engine run agrees with the record: {len(a)} composites, best c = {fresh['summary']['best']['c']}")
