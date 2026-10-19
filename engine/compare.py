"""compare.py -- compares a fresh engine run with the record engine/anaxia_regions.json.

Usage: python3 compare.py <fresh.json> <record.json>
The verdict, the reason and the measurement of every composite, and the summary fields other than
the memory statistics, must agree.
"""
import json
import sys

fresh, record = (json.load(open(f)) for f in sys.argv[1:3])
keys = ("engine", "library", "closed", "accepted", "rejected", "conditional", "best", "pmax",
        "accepted_constants")
bad = [k for k in keys if fresh["summary"].get(k) != record["summary"].get(k)]
a = {c["composite"]: c for c in fresh["composites"]}
b = {c["composite"]: c for c in record["composites"]}
if set(a) != set(b):
    bad.append("composite identifiers")
else:
    bad += [i for i in sorted(a) if any(a[i][f] != b[i][f] for f in ("parts", "outcome", "reason", "measurement"))]
if bad:
    raise SystemExit(f"engine run differs from the record: {bad[:10]}")
print(f"engine run agrees with the record: {len(a)} composites, best c = {fresh['summary']['best']['c']}")
