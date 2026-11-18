#!/usr/bin/env python3
"""Validate demultiplexed Illumina FASTQs and prepare oriented, primer-free libraries."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from itertools import zip_longest
from pathlib import Path

IDS = ("library_id", "sample_id", "study_id", "run_id", "analysis_group")
GROUP = ("marker", "locus", "primer_f", "primer_r", "layout", "orientation", "endpoint")
REQUIRED = {*IDS, *GROUP, "platform", "quality_source", "primer_state", "read1", "read2"}
IUPAC = "ACGTRYSWKMBDHVN"
COMPLEMENT = str.maketrans(IUPAC, "TGCAYRSWMKVHDBN")


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def read_manifest(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or REQUIRED - set(fields):
            raise ValueError("Duplicate columns or missing: " + ", ".join(sorted(REQUIRED - set(fields))))
        rows = list(reader)
    if not rows:
        raise ValueError("Manifest is empty")
    seen_ids, seen_paths, groups = set(), set(), {}
    for row in rows:
        if None in row or None in row.values():
            raise ValueError("Manifest row does not match its header")
        lid = row["library_id"]
        if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", row[k]) for k in IDS):
            raise ValueError(f"{lid}: invalid ID")
        if lid in seen_ids:
            raise ValueError(f"Duplicate library_id: {lid}")
        seen_ids.add(lid)
        if row["platform"] != "illumina" or row["quality_source"] != "full_phred33":
            raise ValueError(f"{lid}: use original Illumina Phred+33 FASTQ; see docs/decisions.md")
        if row["marker"] not in ("16S", "18S", "ITS") or row["layout"] not in ("single", "paired"):
            raise ValueError(f"{lid}: unsupported marker or layout")
        if row["primer_state"] not in ("present", "removed") or row["orientation"] != "forward":
            raise ValueError(f"{lid}: declare oriented R1 and known primer state")
        if not row["locus"].strip() or not row["endpoint"].strip():
            raise ValueError(f"{lid}: provide locus and retained endpoint definition")
        for primer in ("primer_f", "primer_r"):
            row[primer] = row[primer].upper()
            if not re.fullmatch(f"[{IUPAC}]+", row[primer]):
                raise ValueError(f"{lid}: specify the assay oligo in {primer} (5' to 3')")
        signature = tuple(row[k] for k in GROUP)
        if groups.setdefault(row["analysis_group"], signature) != signature:
            raise ValueError(f"{lid}: incompatible assays in analysis_group {row['analysis_group']}")
        if bool(row["read2"]) != (row["layout"] == "paired"):
            raise ValueError(f"{lid}: read2 disagrees with layout")
        for field in ("read1", "read2") if row["read2"] else ("read1",):
            read = (path.parent / row[field]).resolve()
            if not read.is_file() or read in seen_paths:
                raise ValueError(f"{lid}: FASTQ missing or repeated: {read}")
            seen_paths.add(read)
            row[field] = str(read)
    return fields, rows


def fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as stream:
        while header := stream.readline():
            sequence = stream.readline().rstrip(b"\r\n")
            plus = stream.readline()
            quality = stream.readline().rstrip(b"\r\n")
            if (not header.startswith(b"@") or not plus.startswith(b"+") or not sequence or
                    len(sequence) != len(quality) or min(quality) < 33 or max(quality) > 93):
                raise ValueError(f"Invalid 4-line Phred+33 FASTQ record in {path}: {header[:70]!r}")
            tokens = header[1:].split()
            if not tokens:
                raise ValueError(f"Empty read identifier in {path}")
            read_id = tokens[0]
            end = None
            if read_id.endswith((b"/1", b"/2")):
                read_id, end = read_id[:-2], int(read_id[-1:])
            elif len(tokens) > 1 and re.match(rb"^[12]:", tokens[1]):
                end = int(tokens[1][:1])
            yield read_id, end, quality


def audit(row, qc_table, sample_limit):
    paths = [row["read1"]] + ([row["read2"]] if row["read2"] else [])
    cycles = [[] for _ in paths]
    varies = [False] * len(paths)
    distinct = [set() for _ in paths]
    count = 0
    for count, pair in enumerate(zip_longest(*(fastq(p) for p in paths)), 1):
        if len(pair) == 2 and (None in pair or pair[0][0] != pair[1][0] or
                              pair[0][1] == 2 or pair[1][1] == 1):
            raise ValueError(f"{row['library_id']}: unsynchronized mates at record {count}")
        for index, (_, _, q) in enumerate(pair):
            if not varies[index]:
                varies[index] = q.count(q[:1]) != len(q)
                distinct[index].add(q[0] - 33)
            if count <= sample_limit:
                while len(cycles[index]) < len(q):
                    cycles[index].append(Counter())
                for hist, score in zip(cycles[index], q):
                    hist[score - 33] += 1
    if not count and row.get("sample_type") in ("blank", "negative_control"):
        return 0
    if not count:
        raise ValueError(f"{row['library_id']}: empty FASTQ")
    for i, hists in enumerate(cycles):
        if not varies[i] and (len(distinct[i]) == 1 or distinct[i] <= {3, 30}):
            raise ValueError(f"{row['library_id']}: uniform qualities on R{i+1}; verify original base qualities")
        for pos, hist in enumerate(hists, 1):
            total = sum(hist.values())
            cumulative, quantiles = 0, []
            for q, num in sorted(hist.items()):
                cumulative += num
                while len(quantiles) < 3 and cumulative >= total * (0.1, 0.5, 0.9)[len(quantiles)]:
                    quantiles.append(q)
            qc_table.writerow([row["library_id"], i + 1, pos, total, *quantiles])
    return count


def trim_command(row, stage, threads):
    f, r = row["primer_f"], row["primer_r"]
    lid = row["library_id"]
    cmd = ["cutadapt", "-j", str(threads), "-e", "0.1", "-O", "8", "-m", "1",
           "--discard-untrimmed", "--json", str(stage / f"{lid}.cutadapt.json"),
           "-a", f"^{f}...{revcomp(r)}", "-o", str(stage / f"{lid}_R1.fastq.gz")]
    if row["layout"] == "paired":
        cmd += ["--pair-filter=any", "-A", f"^{r}...{revcomp(f)}", "-p", str(stage / f"{lid}_R2.fastq.gz")]
    return cmd + [row["read1"]] + ([row["read2"]] if row["read2"] else [])


def prepare(args):
    fields, rows = read_manifest(args.manifest.resolve())
    if args.threads < 1 or args.qc_reads < 1:
        raise ValueError("threads and qc-reads must be positive")
    if any(r["primer_state"] == "present" for r in rows) and not shutil.which("cutadapt"):
        raise ValueError("cutadapt executable required for primer_state=present")
    out = args.out.resolve()
    if out.exists():
        raise ValueError(f"Output exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{out.name}-", dir=out.parent))
    try:
        if args.hash_inputs:
            with (stage/"input_sha256.tsv").open("w") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["library_id", "mate", "sha256"])
                for row in rows:
                    for mate in ("read1", "read2"):
                        if row[mate]:
                            with open(row[mate], "rb") as raw:
                                checksum = hashlib.file_digest(raw, "sha256").hexdigest()
                            writer.writerow([row["library_id"], mate, checksum])
        with (stage/"quality_by_cycle.tsv").open("w") as qc, (stage/"primer_tracking.tsv").open("w") as tracking:
            quality, track = [csv.writer(h, delimiter="\t", lineterminator="\n") for h in (qc, tracking)]
            quality.writerow(["library_id", "mate", "position", "sampled_bases", "q10", "q50", "q90"])
            track.writerow(["library_id", "input_reads_or_pairs", "primer_free_reads_or_pairs", "primer_state_input"])
            for row in rows:
                n = audit(row, quality, args.qc_reads)
                kept, state = n, row["primer_state"]
                if state == "present" and n:
                    cmd = trim_command(row, stage, args.threads)
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
                    report = json.loads((stage/f"{row['library_id']}.cutadapt.json").read_text())
                    kept = report["read_counts"]["output"]
                    if not kept and row.get("sample_type") not in ("blank", "negative_control"):
                        raise ValueError(f"{row['library_id']}: zero primer matches; verify orientation, primers and spacers")
                    for i in (1, 2) if row["layout"] == "paired" else (1,):
                        row[f"read{i}"] = f"{row['library_id']}_R{i}.fastq.gz"
                row["primer_state"] = "removed"
                track.writerow([row["library_id"], n, kept, state])
        with (stage/"trimmed_manifest.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        shutil.copyfile(args.manifest, stage/"input_manifest.tsv")
        stage.rename(out)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--hash-inputs", action="store_true", help="write SHA-256 of original file bytes (extra I/O)")
    parser.add_argument("--qc-reads", type=int, default=10000,
                        help="first reads per library and mate for per-cycle quantiles; all records checked")
    try:
        prepare(parser.parse_args())
    except (ValueError, OSError, EOFError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"{parser.prog}: {exc}\n")
