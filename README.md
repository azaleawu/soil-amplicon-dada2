# Soil amplicon analysis with DADA2

**Shuai Wu · 武帅**

A working pipeline for soil 16S, fungal ITS and eukaryotic 18S surveys. It links sequence variants to libraries, physical samples, soil profiles and sampling depths. Across studies it learns sequencing errors within each batch and combines exact sequences only within a declared, compatible assay group.

[中文方法指南](docs/decisions.md) · [Research profile / 个人简介](docs/profile.md) · [Public code catalogue / 代码索引](references/README.md) · [Validation](docs/validation.md) · [Atacama soil example](docs/soil-example.md)

## Run a complete example

From an environment with micromamba and git, run the following commands. The demo generates eight libraries locally with variable per-base qualities, known sequence variants, an ITS read-through case, and two independent 16S sequencing runs.

```bash
git clone https://github.com/azaleawu/soil-amplicon-dada2.git
cd soil-amplicon-dada2
micromamba create -y -f environment.yml
micromamba run -n soil-dada2 python examples/make_demo.py demo
micromamba run -n soil-dada2 python scripts/prepare.py demo/manifest.tsv prepared --threads 2 --hash-inputs
micromamba run -n soil-dada2 Rscript scripts/denoise.R prepared/trimmed_manifest.tsv demo/run_settings.tsv denoised 2
micromamba run -n soil-dada2 python examples/check_demo.py demo prepared denoised
micromamba run -n soil-dada2 python scripts/group_summary.py denoised/groups/demo_16S_paired demo/16S_summary.tsv
```

The [Atacama Desert example](docs/soil-example.md) runs the same workflow on five published soil libraries and traces 10,442 input pairs to 4,261 nonchimeric pairs and 128 ASVs. The example includes the public data download, its checksum, editable inputs and per-library results.

## Apply the pipeline to a study

1. Fill one manifest row per library. Record the source accession, soil profile, sampling depth and the relationship between libraries and physical samples. Check the locus, actual primer sequences, orientation, sequence boundaries and quality provenance. [Field definitions and a quality-plot example](docs/inputs.md) · [Manifest example](examples/manifest.example.tsv)
2. Run `scripts/prepare.py` on the manifest. It checks the full files and mate order, summarizes original per-cycle quality, and uses Cutadapt when primers are present. Inspect primer recovery and post-primer quality before setting truncation. [Marker decisions](docs/decisions.md)
3. Set your study/run/group keys and filtering parameters, then run `scripts/denoise.R`. It filters, learns errors, infers ASVs, merges paired reads, removes chimeras within a batch and writes one sparse count table per analysis group. Single-end input skips merging. [Run-setting profiles](profiles/README.md)
4. Assign taxonomy to each group's representative FASTA with a versioned SILVA, PR2 or UNITE reference. The classifier checks a supplied SHA-256 and saves bootstrap support. [Reference downloads, hashes and commands](docs/taxonomy.md)
5. Join the count audit to physical-sample metadata with `scripts/group_summary.py`. It reports retained reads, observed ASVs and Shannon diversity (natural log) for every library, including controls. These unrarefied library summaries support QC. Ecological comparisons need an explicit treatment of sequencing depth, profile pairing and study effects.

## Files to inspect after a run

| File | What to check |
| --- | --- |
| `prepared/primer_tracking.tsv`, `quality_by_cycle.tsv` | Raw input versus primer-free retention; per-cycle Q10/median/Q90 from up to the first 10,000 input reads per mate by default |
| `denoised/runs/<study>/<run>/<group>/` | Learned error curves and models per direction, read tracking and the nonchimeric sequence table |
| `denoised/groups/<group>/read_tracking.tsv` | Input to DADA2, filtered, denoised, merged (paired only) and nonchimeric counts per library |
| `asv_counts.tsv`, `representative_sequences.fasta`, `library_metadata.tsv` inside a group | Sparse integer counts (`library_id, asv_id, count`); stable sequence-based ASV IDs; study and soil sampling metadata |
| `denoised/provenance/` and optional taxonomy output | Run parameters, input table fingerprints, R environment, database version and checksum |

## Coverage and attribution

The curated catalogue contains 22 soil projects or data-based tutorials, 4 general implementations and 12 methodological or database resources. Code entries link to inspected files at fixed revisions. A separate GitHub search ledger records all retrieved candidates. Coverage is dated 18 November 2026; indexing, language and code availability limit recall. [Annotated catalogue and search strategy](references/README.md)

The executable path targets oriented, demultiplexed Illumina reads with full Phred+33 qualities. The catalogue also points to specialist workflows for long reads and functional markers. The repository contains original integration code and links to the upstream authors. Cite the DADA2 algorithm, the reference database release and relevant study-specific protocols when reporting an analysis. [DADA2](https://doi.org/10.1038/nmeth.3869) · [Cutadapt](https://doi.org/10.14806/ej.17.1.200) · [Author's GitHub](https://github.com/azaleawu)
