#!/usr/bin/env python3
"""Build the documented BAQ2420 example from the pinned public QIIME 2 artifact."""
import argparse
import csv
import gzip
import hashlib
import io
import re
import subprocess
import zipfile
from pathlib import Path

ARTIFACT_SHA256 = "f38e1a86c593f64c1391937e67d3f6e3d7b361b82f80fdd5b79b556e4be9fe10"
DATA = "210fed7d-8f02-4cee-b888-166ca00d6581/data/"
SAMPLES = ("BAQ2420.1.1", "BAQ2420.1.2", "BAQ2420.1.3", "BAQ2420.2", "BAQ2420.3")
GROUP = "Atacama_V4_13bp"
RUN = "M00176_65_A41FR"


def write_table(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare(artifact, metadata, output):
    if hashlib.sha256(artifact.read_bytes()).hexdigest() != ARTIFACT_SHA256:
        raise ValueError("Artifact hash differs from the tested QIIME 2 release")
    with metadata.open(newline="", encoding="utf-8-sig") as f:
        site = {r["sample-id"]: r for r in csv.DictReader(f, delimiter="\t")}
    output.mkdir(parents=True, exist_ok=False)
    manifests = []
    with zipfile.ZipFile(artifact) as archive:
        records = list(csv.DictReader(io.StringIO(archive.read(DATA+"MANIFEST").decode())))
        index = {(r["sample-id"], r["direction"]): r["filename"] for r in records}
        for sid in SAMPLES:
            info = site[sid]
            if info["site-name"] != "BAQ2420":
                raise ValueError("Unexpected sample metadata")
            raw, trimmed = [], []
            for direction, mate, oligo in (("forward", "R1", "GTG[CT]CAGC[AC]GCCGCGGTAA"),
                                           ("reverse", "R2", "GGACTAC[ACGT][ACG]GGGT[AT]TCTAAT")):
                p = output / index[sid, direction]
                p.write_bytes(archive.read(DATA+p.name))
                with gzip.open(p, "rt") as fastq:
                    for i, line in enumerate(fastq):
                        if i % 4 == 0 and not line.startswith("@M00176:65:000000000-A41FR:"):
                            raise ValueError("Unexpected instrument/run/flow cell")
                        if i % 4 == 1 and re.match(oligo, line):
                            raise ValueError("Amplification-primer sequence starts a read")
                raw.append(p)
                trimmed.append(output / f"{sid}_{mate}.fastq.gz")
            # Reproduce the tutorial's QIIME 2 left trims and 150-base endpoints.
            subprocess.run(["cutadapt", "-u", "13", "-U", "13", "-l", "137", "-L", "137",
                            "-m", "137:137", "-o", trimmed[0], "-p", trimmed[1], *raw], check=True)
            manifests.append(dict(library_id=sid, sample_id=sid, study_id="ERP019482", run_id=RUN,
                                  analysis_group=GROUP, marker="16S", layout="paired", locus="V4",
                                  primer_f="GTGCCAGCMGCCGCGGTAA", primer_r="GGACTACHVGGGTWTCTAAT",
                                  orientation="forward", endpoint="merged_V4_minus_13bp_each_end",
                                  platform="illumina", quality_source="full_phred33", primer_state="removed",
                                  read1=trimmed[0].name, read2=trimmed[1].name, sample_type="soil",
                                  transect=info["transect-name"], site_name=info["site-name"],
                                  source_depth_code=info["depth"]))
    write_table(output/"manifest.tsv", manifests)
    write_table(output/"run_settings.tsv",
                [dict(study_id="ERP019482", run_id=RUN, analysis_group=GROUP,
                      trunc_len_f=0, trunc_len_r=0, max_ee_f=2, max_ee_r=2,
                      min_len=137, min_overlap=12, max_mismatch=0)])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for v in ("artifact", "metadata", "output"):
        p.add_argument(v, type=Path)
    args = p.parse_args()
    prepare(args.artifact, args.metadata, args.output)
