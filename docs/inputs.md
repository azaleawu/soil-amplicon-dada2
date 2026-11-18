# FASTQ manifest and run settings

Start with demultiplexed, oriented, four-line Illumina FASTQ or FASTQ.GZ. Each record contains full per-base Phred+33 qualities. Keep one row per pair of mate files or single-end file. Paths are relative to the input TSV or absolute. For a small working example generate the demo, which writes its own matching FASTQs, manifest and run settings.

| Column | Meaning |
| --- | --- |
| `study_id`, `run_id` | Source study and physical sequencing batch. Many SRA run accessions can belong to one sequencing batch. If the batch cannot be reconstructed, record the evidence and use a conservative independent unit; check that it supplies enough reads to learn errors. |
| `sample_id`, `library_id` | Physical sample identifier within study, globally unique file/library identifier. Technical replicates can share sample_id. For split lanes of the same library, document lane handling and preserve their mapping. |
| `analysis_group` | Declared comparison group. The script verifies marker, locus, oligo sequences, orientation, read layout and endpoint label match within a group. You supply the evidence that the labels reflect the reads. |
| `marker`, `locus` | Marker is `16S`, `18S` or `ITS`. Name the actual region e.g. `V4`, `V3-V4`, `ITS1` or `ITS2`. Keep unresolved assays in their own groups until verified. |
| `primer_f`, `primer_r` | Both assay oligos, written 5′→3′, excluding barcodes and spacers. Supply them even for primer-free inputs as evidence of targeting. IUPAC degeneracies are supported. |
| `platform`, `quality_source` | This implementation accepts `illumina` and `full_phred33`. Check the archive's description and quality distribution. |
| `layout`, `orientation`, `endpoint` | `paired` or `single`, R1 in `forward` orientation, and a definition such as `merged_insert` or `R1_180`. Document upstream direction normalization. |
| `primer_state` | `present` requires the expected anchored 5′ primer at the start of each retained mate. The opposite primer at the 3′ end is optional. `removed` declares both ends clean, including possible short-insert read-through. Independently check the source treatment first. |
| `read1`, `read2` | Keep read2 empty for single-end data. Do not place two studies or two libraries in the same FASTQ. |
| Extra columns | Preserved. Suggested fields are `profile_id`, `site_id`, `depth_top_cm`, `depth_bottom_cm`, `horizon`, `sample_type`, `source_accession`, `protocol_source`, coordinates and their source, pooling and replicate IDs. For blanks use `blank` or `negative_control`. Missing metadata can remain empty. |

Portable identifiers use letters, digits, dots, dashes and underscores and start with a letter or digit. The demo's oligos are artificial. For real data look up the two sequences used in that study, even if primer names look familiar. For example, the original and modified EMP 515F/806R oligos differ. [EMP sequences and revisions](https://earthmicrobiome.ucsd.edu/protocols-and-standards/16s/)

## Configuring denoising

Use one TSV row per study_id + run_id + analysis_group combination. The [profile sheet](../profiles/run_settings.example.tsv) can be edited in a spreadsheet or a text editor. Copy the three key fields from your manifest. Use the quality plot on **primer-free** reads to choose truncation. The Python QC table summarizes the *input* reads; primer positions can shift the cycle numbers. Inspect the appropriate quality distribution in R, for example:

```r
library(dada2)
m <- read.delim("prepared/trimmed_manifest.tsv", colClasses = "character")
f <- m$read1[1]
if (!startsWith(f, "/")) f <- file.path("prepared", f)
plotQualityProfile(f)
```

| Setting | Interpretation |
| --- | --- |
| `trunc_len_f`, `trunc_len_r` | Fixed length after primer removal; 0 disables fixed length cutting. Single-end groups with fixed truncation must use the same length across runs. |
| `max_ee_f`, `max_ee_r`, `min_len` | Maximum expected errors per mate and minimum retained length. R2 fields are unused on single-end data. |
| `min_overlap`, `max_mismatch` | Required merged overlap and tolerated differences for paired reads. The provided profile starts with 12 and 0. |
| `pool`, `chimera_method`, `learn_nbases` | Optional. Defaults are false, consensus and 100,000,000 bases for error learning. Both read directions have their own error model. |

Cutadapt uses a required anchored 5′ primer, an optional 3′ read-through primer, 10% error tolerance and an 8-base minimum overlap for an unanchored match. Both mates must meet the required-primer condition in a paired library. A protocol with a leading spacer, partial 5′ primer, removed 5′ primer plus retained opposite 3′ primer, or mixed orientations needs its own upstream trimming rule; pass the validated, normalized files as `removed`. Cutadapt's JSON report records linked-adapter matches. With `--hash-inputs`, the preparation step stores SHA-256 hashes for the original FASTQ bytes. Empty blank controls are retained in tracking. Inference uses maxN 0, truncQ 2, PhiX removal, and the worker count supplied on the command line.

If preprocessing fails, the output path stays available for a corrected attempt. A failed R run preserves diagnostics in OUTDIR.partial. Inspect that directory and give the next run a new output name. A completed output directory is never reused by these entry points.
