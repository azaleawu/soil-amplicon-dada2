#!/usr/bin/env python3
"""Generate reproducible synthetic FASTQs and their known insert sequences. No biological data."""
import argparse
import csv
import gzip
import random
from contextlib import ExitStack
from pathlib import Path

FWD, REV = "ACGTTGACTGGACGTA", "TTACGAGCCATGATGC"  # artificial fixture oligos
BASES = "ACGT"


def rc(seq):
    return seq.translate(str.maketrans(BASES, "TGCA"))[::-1]


def tsv(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def with_errors(seq, rng):
    qual, called = [], []
    for base in seq:
        q = rng.choice((22, 25, 28, 30, 32, 34, 37))
        qual.append(chr(33 + q))
        called.append(rng.choice(BASES.replace(base, "")) if rng.random() < 10 ** (-q / 10) else base)
    return "".join(called), "".join(qual)


def generate(out, reads):
    out.mkdir(parents=True, exist_ok=False)
    rng = random.Random(5701)
    configs = [("16S", "V4", "paired", [165, 166, 165]),
               ("18S", "V4", "single", [185, 190, 189]),
               ("ITS", "ITS2", "paired", [82, 133, 174])]
    libraries, settings, truths = [], [], []
    for marker, locus, layout, lengths in configs:
        group = "demo_" + marker + "_" + layout
        templates = ["".join(rng.choices(BASES, k=size)) for size in lengths]
        truth = [t if layout == "paired" else t[:100] for t in templates]
        truths += [dict(analysis_group=group, sequence=seq) for seq in truth]
        for batch in (1, 2) if marker == "16S" else (1,):
            run = "Flowcell" + str(batch)
            settings.append(dict(study_id="Demo", run_id=run, analysis_group=group,
                                 trunc_len_f=100 if layout == "single" else 0, trunc_len_r=0,
                                 max_ee_f=2, max_ee_r=2, min_len=60, min_overlap=12))
            for rep in (1, 2):
                lid = f"{group}_B{batch}_L{rep}"
                sample = "CoreA_Top" if rep == 1 else "CoreA_Deep"
                stem = f"{lid}_R"
                state = "removed" if marker == "18S" else "present"
                libraries.append(dict(library_id=lid, sample_id=sample, study_id="Demo", run_id=run,
                                      analysis_group=group, marker=marker, layout=layout, locus=locus,
                                      primer_f=FWD, primer_r=REV, orientation="forward",
                                      endpoint="merged_insert" if layout=="paired" else "R1_100",
                                      platform="illumina", quality_source="full_phred33", primer_state=state,
                                      read1=stem+"1.fastq.gz", read2=stem+"2.fastq.gz" if layout=="paired" else "",
                                      profile_id="CoreA", depth_top_cm=(0 if rep==1 else 30),
                                      depth_bottom_cm=(10 if rep==1 else 40), sample_type="soil",
                                      source_accession="SIMULATED"))
                with ExitStack() as stack:
                    mates = [stack.enter_context(gzip.open(out/(stem+str(i)+".fastq.gz"), "wt"))
                             for i in ((1,2) if layout == "paired" else (1,))]
                    for n in range(reads):
                        seq = rng.choices(templates, weights=[5, 3, 2])[0]
                        for direction, file in enumerate(mates, 1):
                            insert = seq if direction == 1 else rc(seq)
                            tag, tail = (FWD, rc(REV)) if direction==1 else (REV, rc(FWD))
                            raw = (tag+insert+tail)[:124] if state=="present" else insert[:116]
                            bases, qual = with_errors(raw, rng)
                            file.write(f"@{lid}.{n}/{direction}\n{bases}\n+\n{qual}\n")
    tsv(out/"manifest.tsv", libraries)
    tsv(out/"run_settings.tsv", settings)
    tsv(out/"truth.tsv", truths)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--reads", type=int, default=4000, help="reads/pairs per library")
    args = p.parse_args()
    if args.reads < 1:
        p.error("positive count required")
    generate(args.out, args.reads)
