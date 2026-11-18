#!/usr/bin/env Rscript
# Assign taxonomy to the unmodified ASVs of a single output group.
usage <- paste("Rscript scripts/taxonomy.R GROUP MARKER REFERENCE_FASTA RELEASE SHA256",
               "OUTPUT_DIR [THREADS=1] [MIN_BOOT=50]")

integer_arg <- function(text, label, minimum, maximum) {
  value <- suppressWarnings(as.integer(text))
  if (!grepl("^[0-9]+$", text) || is.na(value) || value < minimum || value > maximum) {
    stop("Invalid ", label, ": ", text, call. = FALSE)
  }
  value
}

reference_header <- function(path) {
  input <- gzfile(path, "rt")  # gzfile reads both ordinary and compressed FASTA.
  on.exit(close(input))
  header <- readLines(input, n = 1L, warn = FALSE)
  if (length(header) != 1L || !startsWith(header, ">")) {
    stop("Expected a reference FASTA. Extract any ZIP or TAR archive first.")
  }
  if (grepl(";tax=", header, ignore.case = TRUE)) {
    stop("UTAX/SINTAX reference detected. Use the DADA2 training or UNITE General FASTA release.")
  }
  header
}

# UNITE General headers carry the taxonomy in the fifth pipe-delimited field.
# Explicit parsing also covers all-eukaryote editions that start with nonfungal SHs.
prepare_unite <- function(input_path, output_path) {
  input <- gzfile(input_path, "rt")
  on.exit(close(input))
  output <- gzfile(output_path, "wt")
  on.exit(close(output), add = TRUE)
  ranks <- "^k__[^;]*;p__[^;]*;c__[^;]*;o__[^;]*;f__[^;]*;g__[^;]*;s__[^;]*;?$"
  repeat {
    lines <- readLines(input, n = 30000L, warn = FALSE)
    if (!length(lines)) break
    at <- which(startsWith(lines, ">"))
    if (length(at)) {
      fields <- strsplit(lines[at], "|", fixed = TRUE)
      if (any(lengths(fields) < 5L)) stop("Malformed UNITE General FASTA header.")
      taxa <- trimws(vapply(fields, `[[`, character(1), 5L))
      taxa <- sub(";sh__[^;]+;?$", "", taxa)
      if (any(!grepl(ranks, taxa))) stop("Unexpected UNITE taxonomic hierarchy. Check the release format.")
      taxa <- sub(";$", "", taxa)
      taxa <- gsub("(^|;)[kpcofgs]__(unidentified)?(?=;|$)", "\\1_DADA2_UNSPECIFIED", taxa, perl = TRUE)
      taxa <- sub(";s__([^;]*_)?sp[.]?$", ";_DADA2_UNSPECIFIED", taxa)
      lines[at] <- paste0(">", taxa, ";")
    }
    writeLines(lines, output)
  }
}

main <- function(args) {
  if (identical(args, "--help")) { cat(usage, "\n"); return(invisible(NULL)) }
  if (length(args) < 6L || length(args) > 8L) stop(usage, call. = FALSE)
  group <- normalizePath(args[1], mustWork = TRUE)
  marker <- toupper(args[2])
  reference <- normalizePath(args[3], mustWork = TRUE)
  release <- args[4]
  expected <- tolower(args[5])
  out <- args[6]
  threads <- if (length(args) >= 7L) integer_arg(args[7], "threads", 1L, .Machine$integer.max) else 1L
  min_boot <- if (length(args) == 8L) integer_arg(args[8], "minimum bootstrap", 0L, 100L) else 50L
  domains <- list(`16S` = "SILVA", `18S` = "PR2", ITS = "UNITE")
  if (!marker %in% names(domains)) stop("Marker must be 16S, 18S or ITS.")
  if (!nzchar(trimws(release)) || grepl("[\r\n\t]", release)) stop("Specify a valid database release.")
  if (!grepl("^[0-9a-f]{64}$", expected)) stop("Expected a 64-character SHA-256 of the reference file.")
  if (file.exists(out)) stop("Output exists. Choose a new output directory for this assignment.")

  counts_meta <- utils::read.delim(file.path(group, "library_metadata.tsv"), colClasses = "character",
                                   quote = "", check.names = FALSE, na.strings = "")
  if (!identical(unique(toupper(counts_meta[["marker"]])), marker)) {
    stop("Requested marker disagrees with library_metadata.tsv or the group contains mixed markers.")
  }
  sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
  actual <- sha256(reference)
  if (!identical(expected, actual)) stop("Reference SHA-256 mismatch. Observed ", actual)
  reads_path <- file.path(group, "representative_sequences.fasta")
  fasta <- Biostrings::readBStringSet(reads_path)
  ids <- names(fasta)
  seqs <- unname(as.character(fasta))
  if (!length(seqs) || length(ids) != length(seqs) || any(!nzchar(ids)) || any(grepl("[[:space:]]", ids)) ||
      anyDuplicated(ids) || anyDuplicated(seqs) || any(!grepl("^[ACGT]+$", seqs))) {
    stop("Representative FASTA must contain unique ASV identifiers and unique uppercase A/C/G/T sequences.")
  }

  ranks <- switch(marker,
    `16S` = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
    `18S` = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species"),
    ITS = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  header <- reference_header(reference)
  ref_format <- "DADA2 training FASTA"
  classifier_input <- reference
  if (grepl("|", header, fixed = TRUE)) {
    if (marker != "ITS" || !grepl("\\|re[pf]s\\|k__", header)) stop("Unsupported reference header format.")
    ref_format <- "UNITE General FASTA"
    classifier_input <- tempfile(fileext = ".fasta.gz")
    on.exit(unlink(classifier_input), add = TRUE)
    prepare_unite(reference, classifier_input)
  } else if (!grepl(";", header, fixed = TRUE)) {
    stop("Expected semicolon-separated taxonomy in DADA2 training FASTA headers.")
  }

  set.seed(1L)
  calls <- dada2::assignTaxonomy(seqs, classifier_input, minBoot = min_boot, tryRC = TRUE,
                                outputBootstraps = TRUE, taxLevels = ranks, multithread = threads)
  # Some DADA2 documentation describes the taxonomy component as 'taxa', while
  # released versions return 'tax'. Use exact indexing for either spelling.
  taxa <- calls[["tax"]]
  if (is.null(taxa)) taxa <- calls[["taxa"]]
  boots <- calls[["boot"]]
  if (!identical(colnames(taxa), ranks) || !identical(colnames(boots), ranks)) {
    stop("Reference rank count disagrees with the selected marker. See docs/taxonomy.md.")
  }
  if (!setequal(rownames(taxa), seqs) || !setequal(rownames(boots), seqs)) stop("Unmatched output sequences.")
  tax_table <- data.frame(asv_id = ids, taxa[seqs, , drop = FALSE], check.names = FALSE, row.names = NULL)
  boot_table <- data.frame(asv_id = ids, boots[seqs, , drop = FALSE], check.names = FALSE, row.names = NULL)

  parent_dir <- dirname(out)
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(pattern = ".taxonomy-", tmpdir = parent_dir)
  if (!dir.create(stage)) stop("Could not create output staging directory.")
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  write_tsv <- function(data, filename) {
    utils::write.table(data, file.path(stage, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  write_tsv(tax_table, "taxonomy.tsv")
  write_tsv(boot_table, "bootstrap.tsv")
  write_tsv(data.frame(group = basename(group), marker = marker, reference_family = domains[[marker]],
                       release = release, reference_file = basename(reference), reference_sha256 = actual,
                       reference_format = ref_format, representative_fasta_sha256 = sha256(reads_path),
                       asvs = length(ids), min_boot = min_boot, try_rc = TRUE, seed = 1L, threads = threads,
                       dada2_version = as.character(utils::packageVersion("dada2"))), "provenance.tsv")
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(stage, "sessionInfo.txt"))
  if (file.exists(out) || !file.rename(stage, out)) stop("Could not publish completed output directory.")
  message("Wrote taxonomy for ", length(ids), " ASVs to ", out)
}

main(commandArgs(trailingOnly = TRUE))
