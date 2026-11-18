#!/usr/bin/env python3
"""Check read accounting, stable sequence identity and recovery of known demo templates."""
import argparse
import csv
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def require(test, message):
    if not test:
        raise ValueError(message)


def verify(demo, prepared, output):
    expected = defaultdict(set)
    for item in rows(demo / "truth.tsv"):
        expected[item["analysis_group"]].add(item["sequence"])
    raw = {x["library_id"]: x for x in rows(prepared / "primer_tracking.tsv")}
    for group, truth in expected.items():
        path = output/"groups"/group
        lines = (path/"representative_sequences.fasta").read_text().splitlines()
        sequences = dict(zip([line[1:] for line in lines[::2]], lines[1::2], strict=True))
        require(truth <= set(sequences.values()), f"{group}: expected templates missing")
        for key, seq in sequences.items():
            require(key == "ASV_" + sha256(seq.encode()).hexdigest(), f"{key}: unstable identifier")
        counts = Counter()
        for row in rows(path/"asv_counts.tsv"):
            require(row["asv_id"] in sequences, "unrecognized ASV")
            counts[row["library_id"]] += int(row["count"])
        tracked = rows(path/"read_tracking.tsv")
        expected_libs = {r["library_id"] for r in rows(demo/"manifest.tsv") if r["analysis_group"] == group}
        require({r["library_id"] for r in tracked} == expected_libs, f"{group}: library coverage")
        require(set(counts) <= expected_libs, f"{group}: unexpected library")
        for row in tracked:
            lid = row["library_id"]
            before = int(raw[lid]["primer_free_reads_or_pairs"])
            require(before == int(row["input"]), f"{lid}: preprocessing count mismatch")
            require(0 < counts[lid] == int(row["nonchim"]) <= int(row["filtered"]) <= before,
                    f"{lid}: invalid count balance")
        print(f"{group}: recovered {len(truth)}/{len(truth)} known sequences; "
              f"{len(sequences)} total ASVs; {len(tracked)} libraries; {sum(counts.values())} final reads/pairs")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("demo", "prepared", "denoised"):
        p.add_argument(key, type=Path)
    a = p.parse_args()
    verify(a.demo, a.prepared, a.denoised)
