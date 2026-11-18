#!/usr/bin/env python3
"""Join original field metadata to per-library read recovery and observed ASV diversity."""
import argparse
import csv
import math
from collections import Counter
from pathlib import Path


def read_rows(path):
    with path.open() as file:
        yield from csv.DictReader(file, delimiter="\t")


def summarize(group, out):
    tracking = {r["library_id"]: r for r in read_rows(group / "read_tracking.tsv")}
    richness, depth, n_log_n = Counter(), Counter(), Counter()
    for row in read_rows(group / "asv_counts.tsv"):
        key, count = row["library_id"], int(row["count"])
        if key not in tracking or count < 1:
            raise ValueError(f"Invalid sparse count for {key}")
        depth[key] += count
        richness[key] += 1
        n_log_n[key] += count * math.log(count)
    records = []
    for row in read_rows(group / "library_metadata.tsv"):
        key = row["library_id"]
        retained = depth[key]
        tr = tracking.pop(key)
        if retained != int(tr["nonchim"]):
            raise ValueError(f"{key}: sum of ASV counts disagrees with tracking")
        read_in = int(tr["input"])
        extra = {"reads_after_primer_handling":read_in, "nonchimeric_reads":retained,
                 "nonchim_fraction_postprimer": format(retained/read_in, ".6g") if read_in else "",
                 "observed_asvs":richness[key],
                 "shannon_ln":format(math.log(retained)-n_log_n[key]/retained, ".8g") if retained else ""}
        collisions = set(row) & set(extra)
        if collisions:
            raise ValueError(f"Reserved output column in manifest: {collisions}")
        records.append({**row, **extra})
    if not records or tracking:
        raise ValueError("Metadata/tracking library mismatch")
    with out.open("x", newline="") as file:
        writer = csv.DictWriter(file, records[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("group", type=Path)
    p.add_argument("output_tsv", type=Path)
    args = p.parse_args()
    summarize(args.group, args.output_tsv)
