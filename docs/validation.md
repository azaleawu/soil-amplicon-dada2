# Validation record · 18 November 2026

## FASTQ to nonchimeric counts

The commands in the repository README were executed using R 4.5.3, DADA2 1.38.0 and Cutadapt 5.2. The generator uses random nucleotide templates, per-base quality variation and substitutions sampled according to the simulated quality. Each library starts with 4,000 reads or read pairs. The paired 16S libraries span two error-learning runs. The paired ITS fixture contains a short insert that reads through the opposite primer. The 18S input is primer-free and single-end.

| Synthetic fixture | Libraries | Input | Nonchimeric reads or pairs | Known templates recovered | Additional ASVs |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16S paired | 4 | 16,000 | 15,988 | 3/3 | 0 |
| 18S single | 2 | 8,000 | 8,000 | 3/3 | 0 |
| ITS paired, variable length | 2 | 8,000 | 7,996 | 3/3 | 0 |

The verifier checks the actual representative sequences, their SHA-256 identities, the membership of each library, the join from raw-read handling to DADA2, and the sum of nonchimeric counts. The same fixture also passed with R 4.5.0 and DADA2 1.34.0. Small simulated libraries are useful for detecting implementation errors. Sensitivity to rare organisms, contaminants, primer bias and realistic sequencing errors is specific to the assay and sample matrix.

## Published soil data

Five libraries from the public QIIME 2 Atacama soil dataset were processed end to end with DADA2 1.38.0. From 10,442 paired reads, 9,432 passed filtering and 4,261 merged successfully. The final nonchimeric table contained 128 distinct ASVs. [Dataset, exact commands and library-level counts](soil-example.md)

## Classification back ends

For SILVA 138.2, PR2 5.1.1 and UNITE dynamic all-eukaryotes 2025-02-19, three 250-base fragments from a deterministic 256-record subset of each reference were run through the classifier with R 4.5.3 and DADA2 1.38.0. All nine returned the genus of the reference record and preserved their original sequence identifiers and row order. The returned rank counts were 6, 9 and 7 respectively. The full 156,820-record UNITE FASTA was also converted successfully. This self-reference check tests reference parsing and data plumbing. An independent estimate of classification accuracy would require different query sequences or a mock community. [Reference files, verified hashes and output schema](taxonomy.md)
