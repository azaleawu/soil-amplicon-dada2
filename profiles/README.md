# Marker profiles

The master sheet below is a starting point for editing the run-specific settings. Its group and run names are **examples**. Copy the rows relevant to your manifest and change their identifiers. A zero truncation length disables a fixed end cut. The example read lengths for single-end data must be replaced after inspecting each assay.

| Profile in the sheet | What to set before running |
|---|---|
| 16S V4 paired | Pick both truncation positions from the post-primer quality profile and insert-length range. Compare the predicted overlap with the chosen `min_overlap`. |
| Fungal ITS2 paired | Begin with no fixed truncation. Review read-through, merge loss, and lengths. Change locus and group for ITS1. |
| 18S V4 paired | Confirm orientation, real insert range and coverage of the taxa of interest. Choose truncation lengths using measured quality and merge retention. |
| R1 only | Use a common primer, direction and empirically chosen retained endpoint across the entire single-end group. Make a separate group for each marker and locus. |

The default `pool=false` infers variants independently per library using shared run-specific error estimates. Set `pool=pseudo` only as a documented change, e.g. when recovery of low-abundance variants is material. Record a fresh run-settings snapshot for a sensitivity analysis. For paired reads with no fixed truncation, both entries are zero. For single-end data, the reverse parameters are unused.

[Run-settings example](run_settings.example.tsv) · [Profile decisions 中文](../docs/decisions.md) · [Reference classifier selection](../docs/taxonomy.md)
