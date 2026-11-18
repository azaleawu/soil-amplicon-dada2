# Taxonomic assignment of one ASV group

After the core run, classify a group using its `representative_sequences.fasta` and `library_metadata.tsv` (the latter includes a `marker` column). Supply the same marker on the command line. The expected SHA-256 is the digest of the exact file passed to the classifier, including its compression. Download the reference once, verify the publisher's checksum when provided, and retain the accepted SHA-256 for later runs. The R environment needs DADA2, Biostrings and digest.

| Marker argument | Reference to use | Rank schema |
| --- | --- | --- |
| `16S` | DADA2-formatted SILVA SSU 138.2 **toGenus** training FASTA for bacteria and archaea. [Versioned deposit](https://zenodo.org/records/14169026) | Domain, Phylum, Class, Order, Family, Genus |
| `18S` | PR2 version 5.x, file ending in `_SSU_dada2.fasta.gz`. [Releases](https://github.com/pr2database/pr2database/releases) · [PR2 rank definitions](https://pr2-database.org/documentation/pr2-taxonomy-9-levels/) | Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species |
| `ITS` | UNITE General FASTA dynamic release, all eukaryotes, or a seven-rank DADA2-formatted fungal ITS training FASTA from UNITE. Use the release date, taxon coverage and clustering variant to identify the reference. [Versioned example](https://doi.org/10.15156/BIO/3301231) · [Reference formats](https://unite.ut.ee/repository.php) | Kingdom, Phylum, Class, Order, Family, Genus, Species |

`16S` supports the six-rank training set. SILVA `toSpecies` and `assignSpecies` serve other workflows. Older PR2 v4 files contain eight ranks. For UNITE General FASTA the script extracts the seven taxonomic fields, retains named species as given and treats unidentified ranks and `sp` labels as missing. A trailing `sh__` tag is excluded from the rank hierarchy. UNITE UTAX and SINTAX files carry a different header syntax and the script requests the General FASTA if one of those files is supplied. Extract any `.zip` or `.tgz` archive before use. [DADA2 training reference guide](https://benjjneb.github.io/dada2/training.html)

## Verified reference files

These SHA-256 values were calculated from the versioned downloads. For UNITE the digest identifies the extracted, uncompressed FASTA, for PR2 and SILVA it identifies the distributed `.gz` file.

| File | SHA-256 |
| --- | --- |
| `silva_nr99_v138.2_toGenus_trainset.fa.gz` | `e0af8246b900a922ea3ea68e7cb1a5e70ceed4db84246b6ea9865eb132452093` |
| `pr2_version_5.1.1_SSU_dada2.fasta.gz` | `0c8728abcbb2126eed2c7e587f820cbce39c138cdfdb51239bbf18621498462d` |
| `sh_general_release_dynamic_all_19.02.2025.fasta` | `29ec3e4232c3ae7b0e95bf30d376a0c92e74d0834ff6f2212df71c9fcf7f9a91` |

## Linux example

In this example `results/groups/soil_16s_v4` is the output group to classify. Change that path to your actual group. Zenodo reports MD5 `1764e2a36b4500ccb1c7d5261948a414` for the indicated six-rank SILVA download. The following commands verify the download and save its SHA-256 for reuse. [SILVA 138.2 deposit](https://zenodo.org/records/14169026)

```bash
reference=refs/silva_nr99_v138.2_toGenus_trainset.fa.gz
mkdir -p refs
curl -fL -o "$reference" \
  'https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1'
printf '%s  %s\n' 1764e2a36b4500ccb1c7d5261948a414 "$reference" | md5sum --check - && \
  sha256sum "$reference" > refs/silva_138.2_toGenus.sha256
```

```bash
Rscript scripts/taxonomy.R \
  results/groups/soil_16s_v4 16S "$reference" '138.2_toGenus' \
  "$(cut -d' ' -f1 refs/silva_138.2_toGenus.sha256)" \
  results/taxonomy/soil_16s_v4_silva138.2 4 50
```

The positional arguments are `GROUP MARKER REFERENCE_FASTA RELEASE SHA256 OUTPUT_DIR [THREADS] [MIN_BOOT]`. Threads default to 1 and minimum bootstrap to 50. On macOS `shasum -a 256 FILE` reports SHA-256. Each run writes a new output directory, which allows locally curated databases or different thresholds to be compared side by side. The classifier tests both strand orientations.

| Output | Contents |
| --- | --- |
| `taxonomy.tsv` | One row per `asv_id`, in representative FASTA order. Unsupported taxonomic calls are `NA`. |
| `bootstrap.tsv` | Per-rank DADA2 bootstrap support on the 0–100 scale for the same ASV IDs. |
| `provenance.tsv` | Marker, database family, supplied release, reference filename and checksum, input FASTA checksum, DADA2 version, threshold, threads and seed. |
| `sessionInfo.txt` | Versions of R and the loaded packages. |

Inspect the unfiltered assignments before deriving clade-specific tables. An 18S library may include fungi, animals and plants. Report how the protist subset was defined and retain the source table. In a short-read ITS or 18S assay, a species label and its bootstrap value express the classifier's support within the chosen reference and amplicon. For repeated analyses, hold the database release and target sequence interval constant and rerun classification when the reference changes.

## Validation of the three reference formats

The verification run used R 4.5.3 and DADA2 1.38.0. For each of the three files above, a deterministic sample of 256 reference records provided a small training set. Three distinct 250-base fragments with known genus labels from that set were classified with `minBoot=50`. All nine ASV identifiers matched the input order in both output tables. Each marker recovered all three expected genus labels with genus-level bootstrap support of 100. The SILVA, PR2 and UNITE outputs had 6, 9 and 7 rank columns, respectively. The UNITE header conversion was also run over all 156,820 records in the complete reference. Because query fragments were drawn from the training subsets, these checks establish file compatibility and correct output mapping. Independent mock communities would be needed to measure taxonomic accuracy.
