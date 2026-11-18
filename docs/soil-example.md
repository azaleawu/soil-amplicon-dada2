# Atacama Desert soil 16S example

Five soil libraries collected at BAQ2420 on the Baquedano transect in northern Chile yielded 128 nonchimeric ASVs from 10,442 paired reads in a test of this pipeline. The sequence input is the public, demultiplexed **QIIME 2 2024.10 teaching subset** of Neilson et al. (2017). The original study is archived as ENA **ERP019482** and Qiita **10360**. [Study](https://doi.org/10.1128/mSystems.00195-16) · [ENA](https://www.ebi.ac.uk/ena/browser/view/ERP019482) · [QIIME 2 release and sample selection history](https://docs.qiime2.org/2024.10/tutorials/atacama-soils/)

## Selection and endpoints

The example takes every BAQ2420 library present in the released artifact. All five have instrument M00176, run 65 and flow cell 000000000-A41FR in their FASTQ headers. Sample names and depth codes come from the QIIME 2 metadata. The original depth categories are retained as codes. [Sample metadata](https://data.qiime2.org/2024.10/tutorials/atacama-soils/sample_metadata.tsv)

The paper reports Illumina MiSeq sequencing of V4 with the EMP 515F/806R pair and identifies the assay as the same amplification protocol used by Fierer et al. (2012). The example records the original Caporaso EMP oligos as a protocol-based annotation. The cited Methods do not print the individual oligo sequences. Neither the original nor the updated full-length primer pattern occurs at its expected 5′ read start in this subset (exact IUPAC matching). The preprocessing script follows the QIIME 2 tutorial’s trimming decision and retains bases 14–150 of each mate. These five libraries share a separate analysis group because the retained endpoints are part of the ASV definition. [Paper Methods](https://doi.org/10.1128/mSystems.00195-16) · [Original and updated EMP oligos](https://earthmicrobiome.org/protocols-and-standards/16s/) · [Tutorial trim settings](https://docs.qiime2.org/2024.10/tutorials/atacama-soils/)

## Reproduce from the repository root

Use a fresh output directory and an activated environment with Python, Cutadapt and DADA2. The example was executed with R 4.5.3, DADA2 1.38.0 and Cutadapt 5.2 under Linux. The preparation script verifies the release artifact with SHA-256 before extracting any reads. QIIME 2 is unnecessary for reading this artifact.

```bash
mkdir -p work/atacama
curl --fail --location --output work/atacama/demux.qza \
  https://docs.qiime2.org/2024.10/data/tutorials/atacama-soils/demux.qza
curl --fail --location --output work/atacama/metadata.tsv \
  https://data.qiime2.org/2024.10/tutorials/atacama-soils/sample_metadata.tsv
python examples/make_atacama.py work/atacama/demux.qza \
  work/atacama/metadata.tsv work/atacama/input
python scripts/prepare.py work/atacama/input/manifest.tsv \
  work/atacama/prepared --hash-inputs
Rscript scripts/denoise.R work/atacama/prepared/trimmed_manifest.tsv \
  work/atacama/input/run_settings.tsv work/atacama/asvs 1
```

The artifact UUID is `210fed7d-8f02-4cee-b888-166ca00d6581`. The tested archive has SHA-256 `f38e1a86c593f64c1391937e67d3f6e3d7b361b82f80fdd5b79b556e4be9fe10`. The metadata file has SHA-256 `7cff810ad86a621ebc78a16b690255d839ea58e2622ad45a11a838d0f8c5b3bd`.

## Observed counts

Each library retained every pair through the 13-base endpoint trim. DADA2 used maxEE 2 per mate, independent denoising, a minimum exact overlap of 12 bp and consensus chimera removal. The denoised columns count individual mates. Input, filtered, merged and nonchimeric columns count pairs.

| Library | Source depth code | Input pairs | Filtered | Denoised F | Denoised R | Merged | Nonchimeric | ASVs present |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BAQ2420.1.1 | 1 | 2,343 | 2,115 | 1,700 | 1,701 | 575 | 575 | 19 |
| BAQ2420.1.2 | 2 | 2,103 | 1,905 | 1,639 | 1,679 | 770 | 770 | 30 |
| BAQ2420.1.3 | 3 | 2,059 | 1,887 | 1,670 | 1,731 | 983 | 983 | 43 |
| BAQ2420.2 | 2 | 1,972 | 1,784 | 1,583 | 1,619 | 1,037 | 1,037 | 38 |
| BAQ2420.3 | 2 | 1,965 | 1,741 | 1,607 | 1,611 | 896 | 896 | 39 |
| All five libraries | | **10,442** | **9,432** | **8,199** | **8,341** | **4,261** | **4,261** | **128 distinct** |

The final table retained 40.8% of input pairs. The number of merged pairs was substantially lower than either denoised mate count. No inferred variant was flagged as chimeric in this subset. ASV lengths ranged from 226 to 250 bases. The subset supports a check of real FASTQ handling, read accounting and marker endpoints. Ecological comparisons require the broader sampling design and full sequencing depth.
