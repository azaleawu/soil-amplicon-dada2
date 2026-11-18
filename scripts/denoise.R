#!/usr/bin/env Rscript
# Infer ASVs from demultiplexed, primer-free, synchronized Phred+33 Illumina FASTQs.
# Usage: Rscript scripts/denoise.R manifest.tsv settings.tsv results/denoised [threads]
# A group is a curator-defined unit with matching locus, orientation and endpoints.

fail <- function(...) stop(paste0(...), call. = FALSE)
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
read_tsv <- function(path, required) {
  x <- read.delim(path, colClasses = "character", na.strings = character(),
                  check.names = FALSE, comment.char = "", quote = "", fill = FALSE)
  if (!nrow(x) || anyDuplicated(names(x))) fail("Empty or duplicate-header TSV: ", path)
  missing <- setdiff(required, names(x))
  if (length(missing)) fail("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  x
}
run_key <- function(x) paste(x$study_id, x$run_id, x$analysis_group, sep = "/")
id_columns <- c("study_id", "run_id", "analysis_group", "library_id", "sample_id")

read_manifest <- function(path) {
  m <- read_tsv(path, c(id_columns, "marker", "layout", "read1", "read2"))
  for (nm in id_columns) {
    if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", m[[nm]]))) fail("Invalid ID in ", nm)
  }
  if (anyDuplicated(m$library_id)) fail("library_id must be unique across the manifest")
  if (any(!m$marker %in% c("16S", "18S", "ITS"))) fail("marker must be 16S, 18S or ITS")
  if (any(!m$layout %in% c("paired", "single"))) fail("layout must be paired or single")
  if (any(!nzchar(m$read1) | (m$layout == "paired") != nzchar(m$read2))) {
    fail("Every library needs read1, and only paired libraries need read2")
  }
  for (g in split(m, m$analysis_group)) {
    if (length(unique(g$marker)) != 1L || length(unique(g$layout)) != 1L) {
      fail("Mixed markers or layouts in group ", g$analysis_group[1])
    }
  }
  for (nm in c("read1", "read2")) {
    i <- nzchar(m[[nm]])
    raw <- path.expand(m[[nm]][i])
    relative <- !grepl("^(/|[[:alpha:]]:[/\\\\])", raw)
    raw[relative] <- file.path(dirname(normalizePath(path)), raw[relative])
    m[[nm]][i] <- normalizePath(raw, mustWork = TRUE)
  }
  paths <- c(m$read1, m$read2[nzchar(m$read2)])
  if (anyDuplicated(paths)) fail("FASTQ path assigned more than once")
  m
}

read_settings <- function(path, m) {
  required <- c("study_id", "run_id", "analysis_group", "trunc_len_f", "trunc_len_r",
                "max_ee_f", "max_ee_r", "min_len", "min_overlap")
  s <- read_tsv(path, required)
  defaults <- list(max_mismatch = "0", pool = "false", chimera_method = "consensus",
                   learn_nbases = "100000000")
  for (nm in names(defaults)) {
    if (!nm %in% names(s)) s[[nm]] <- defaults[[nm]]
    s[[nm]][!nzchar(s[[nm]])] <- defaults[[nm]]
  }
  keys <- run_key(s)
  if (anyDuplicated(keys) || !setequal(keys, run_key(m))) {
    fail("Settings must contain exactly one row per study_id/run_id/analysis_group in manifest")
  }
  integers <- c("trunc_len_f", "trunc_len_r", "min_len", "min_overlap", "max_mismatch", "learn_nbases")
  for (nm in c(integers, "max_ee_f", "max_ee_r")) {
    values <- suppressWarnings(as.numeric(s[[nm]]))
    if (any(!is.finite(values) | values < 0)) fail("Expected finite nonnegative values for ", nm)
    if (nm %in% integers && any(values != floor(values))) fail("Expected integers for ", nm)
    s[[nm]] <- values
  }
  is_paired <- m$layout[match(keys, run_key(m))] == "paired"
  if (any(s$min_len < 1 | s$learn_nbases < 1 | (is_paired & s$min_overlap < 1))) {
    fail("min_len, learn_nbases and paired min_overlap must be positive")
  }
  if (any((s$trunc_len_f > 0 & s$trunc_len_f < s$min_len) |
          (is_paired & s$trunc_len_r > 0 & s$trunc_len_r < s$min_len))) {
    fail("Fixed truncation is shorter than min_len")
  }
  for (group in unique(s$analysis_group[!is_paired])) {
    if (length(unique(s$trunc_len_f[s$analysis_group == group])) > 1L) {
      fail("Inconsistent single-end truncation across runs in ", group)
    }
  }
  s$pool <- tolower(s$pool)
  if (any(!s$pool %in% c("false", "pseudo", "true"))) fail("pool must be false, pseudo or true")
  if (any(!s$chimera_method %in% c("consensus", "pooled", "per-sample"))) {
    fail("Unsupported chimera_method")
  }
  s
}

has_reads <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  length(readLines(con, n = 1, warn = FALSE)) == 1L
}

fit_errors <- function(paths, nbases, threads, folder, direction) {
  model <- dada2::learnErrors(paths, nbases = nbases, randomize = TRUE,
                              multithread = threads, verbose = TRUE)
  saveRDS(model, file.path(folder, paste0("error_model_", direction, ".rds")))
  grDevices::pdf(file.path(folder, paste0("error_model_", direction, ".pdf")),
                 width = 9, height = 7, useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  print(dada2::plotErrors(model, nominalQ = TRUE))
  model
}

# DADA2 can return an individual S4 object when only one library survives.
get_calls <- function(files, model, pool, threads) {
  calls <- dada2::dada(unname(files), err = model, pool = pool, multithread = threads)
  if (methods::is(calls, "dada")) calls <- list(calls)
  if (length(calls) != length(files)) fail("Unexpected number of denoised libraries")
  setNames(calls, names(files))
}

infer_block <- function(f, r, err_f, err_r, s, threads, pool) {
  fw <- get_calls(f, err_f, pool, threads)
  f_counts <- vapply(fw, function(x) sum(dada2::getUniques(x)), numeric(1))
  if (is.null(r)) {
    return(list(variants = lapply(fw, dada2::getUniques), f = f_counts, r = NULL))
  }
  rv <- get_calls(r, err_r, pool, threads)
  r_counts <- vapply(rv, function(x) sum(dada2::getUniques(x)), numeric(1))
  merged <- lapply(names(f), function(id) {
    x <- dada2::mergePairs(fw[[id]], f[[id]], rv[[id]], r[[id]],
                           minOverlap = s$min_overlap, maxMismatch = s$max_mismatch)
    setNames(x$abundance, x$sequence)
  })
  list(variants = setNames(merged, names(f)), f = f_counts, r = r_counts)
}

process_run <- function(m, s, folder, threads) {
  if (!dir.create(folder, recursive = TRUE)) fail("Cannot create run output ", folder)
  filt <- file.path(folder, "filtered")
  if (!dir.create(filt)) fail("Cannot create filter directory")
  ids <- m$library_id
  paired <- m$layout[1] == "paired"
  f <- setNames(file.path(filt, paste0(ids, "_R1.fastq.gz")), ids)
  r <- setNames(file.path(filt, paste0(ids, "_R2.fastq.gz")), ids)
  counts <- data.frame(library_id = ids, input = 0, filtered = 0, denoised_f = 0,
                       denoised_r = if (paired) 0 else NA_real_,
                       merged = if (paired) 0 else NA_real_, nonchim = 0,
                       row.names = ids)
  present <- vapply(m$read1, has_reads, logical(1))
  if (paired && any(present != vapply(m$read2, has_reads, logical(1)))) {
    fail("Empty forward/reverse mismatch in ", run_key(s))
  }
  if (any(present)) {
    a <- list(fwd = m$read1[present], filt = unname(f[present]), truncLen = s$trunc_len_f,
              maxEE = s$max_ee_f, minLen = s$min_len, maxN = 0, truncQ = 2, rm.phix = TRUE,
              qualityType = "FastqQuality", compress = TRUE, multithread = threads)
    if (paired) {
      a$rev <- m$read2[present]
      a$filt.rev <- unname(r[present])
      a$truncLen <- c(s$trunc_len_f, s$trunc_len_r)
      a$maxEE <- c(s$max_ee_f, s$max_ee_r)
    }
    filtered <- do.call(dada2::filterAndTrim, a)
    counts$input[present] <- filtered[, "reads.in"]
    counts$filtered[present] <- filtered[, "reads.out"]
  }

  keep <- ids[counts$filtered > 0]
  seqtab <- matrix(0L, nrow = 0, ncol = 0)
  if (length(keep)) {
    err_f <- fit_errors(unname(f[keep]), s$learn_nbases, threads, folder, "F")
    err_r <- if (paired) fit_errors(unname(r[keep]), s$learn_nbases, threads, folder, "R")
    # Independent inference processes one library in memory at a time.
    pool <- switch(s$pool, "false" = FALSE, "pseudo" = "pseudo", "true" = TRUE)
    blocks <- if (identical(pool, FALSE)) as.list(keep) else list(keep)
    variants <- list()
    for (block in blocks) {
      x <- infer_block(f[block], if (paired) r[block], err_f, err_r, s, threads, pool)
      counts[block, "denoised_f"] <- x$f
      if (paired) {
        counts[block, "denoised_r"] <- x$r
        counts[block, "merged"] <- vapply(x$variants, sum, numeric(1))
      }
      variants[block] <- x$variants
    }
    positive <- lengths(variants) > 0L
    if (any(positive)) {
      seqtab <- dada2::makeSequenceTable(variants[positive])
      seqtab <- dada2::removeBimeraDenovo(seqtab, method = s$chimera_method,
                                          multithread = threads, verbose = TRUE)
      if (ncol(seqtab)) counts[rownames(seqtab), "nonchim"] <- rowSums(seqtab)
    }
  }
  if (any(counts$nonchim > counts$denoised_f | counts$denoised_f > counts$filtered |
          counts$filtered > counts$input)) fail("Invalid read balance in ", run_key(s))
  if (paired && any(counts$nonchim > counts$merged | counts$merged > counts$denoised_r |
                    counts$merged > counts$denoised_f | counts$denoised_r > counts$filtered)) {
    fail("Invalid paired read balance in ", run_key(s))
  }
  if (!sum(counts$nonchim)) warning("No nonchimeric reads retained in ", run_key(s))
  write_tsv(counts, file.path(folder, "read_tracking.tsv"))
  saveRDS(seqtab, file.path(folder, "sequence_table.rds"))
  unlink(filt, recursive = TRUE)
  counts
}

export_group <- function(m, all_tracking, run_dirs, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  track <- all_tracking[match(m$library_id, all_tracking$library_id), , drop = FALSE]
  write_tsv(m, file.path(out, "library_metadata.tsv"))
  write_tsv(track, file.path(out, "read_tracking.tsv"))
  # Keep group exports sparse. A global study can have millions of all-zero cells.
  target <- file.path(out, "asv_counts.tsv")
  write_tsv(data.frame(library_id = character(), asv_id = character(), count = integer()), target)
  dictionary <- new.env(hash = TRUE, parent = emptyenv())
  for (folder in run_dirs) {
    tab <- readRDS(file.path(folder, "sequence_table.rds"))
    seqs <- colnames(tab)
    if (!length(seqs)) next
    asvs <- vapply(seqs, function(seq) {
      feature <- dictionary[[seq]]
      if (is.null(feature)) {
        feature <- paste0("ASV_", digest::digest(seq, algo = "sha256", serialize = FALSE))
        dictionary[[seq]] <- feature
      }
      feature
    }, character(1), USE.NAMES = FALSE)
    for (i in seq_len(nrow(tab))) {
      hits <- which(tab[i, ] > 0L)
      if (!length(hits)) next
      rows <- data.frame(library_id = rownames(tab)[i], asv_id = asvs[hits],
                         count = unname(tab[i, hits]))
      write.table(rows, target, append = TRUE, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
    }
  }
  sequences <- ls(dictionary, all.names = TRUE, sorted = TRUE)
  features <- vapply(sequences, function(seq) dictionary[[seq]], character(1), USE.NAMES = FALSE)
  if (anyDuplicated(features)) fail("Sequence hash collision")
  writeLines(as.vector(rbind(paste0(">", features, recycle0 = TRUE), sequences)),
             file.path(out, "representative_sequences.fasta"))
}

main <- function(args) {
  if (!length(args) || args[1] %in% c("-h", "--help")) {
    cat("Usage: Rscript scripts/denoise.R trimmed_manifest.tsv run_settings.tsv OUTDIR [THREADS]\n")
    return(invisible(NULL))
  }
  if (!length(args) %in% 3:4) fail("Expected manifest, settings, output directory and optional threads")
  threads <- if (length(args) == 4L) suppressWarnings(as.numeric(args[4])) else 1L
  if (!is.finite(threads) || threads < 1 || threads != floor(threads)) fail("Invalid thread count")
  m <- read_manifest(args[1])
  s <- read_settings(args[2], m)
  for (pkg in c("dada2", "digest")) {
    if (!requireNamespace(pkg, quietly = TRUE)) fail("Missing R package: ", pkg)
  }
  target <- sub("/+$", "", args[3])
  if (!nzchar(target) || file.exists(target)) fail("Choose a new output directory: ", target)
  partial <- paste0(target, ".partial")
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(partial)) fail("Output occupied, or an earlier run left diagnostics at ", partial)
  message("Staged outputs: ", partial)
  provenance <- file.path(partial, "provenance")
  dir.create(provenance)
  write_tsv(m, file.path(provenance, "trimmed_manifest.tsv"))
  write_tsv(s, file.path(provenance, "run_settings.tsv"))
  checks <- tools::md5sum(args[1:2])
  write_tsv(data.frame(file = names(checks), md5 = unname(checks)),
            file.path(provenance, "input_tables_md5.tsv"))
  writeLines(c(paste("UTC start", format(Sys.time(), tz = "UTC", usetz = TRUE)),
               paste("threads", threads), "random seed reset to 1 for each run",
               capture.output(sessionInfo())), file.path(provenance, "session.txt"))

  run_dirs <- file.path(partial, "runs", s$study_id, s$run_id, s$analysis_group)
  tracking <- vector("list", nrow(s))
  for (i in seq_len(nrow(s))) {
    required_key <- run_key(s[i, ])
    message("Processing ", required_key)
    rows <- m[run_key(m) == required_key, , drop = FALSE]
    rows <- rows[order(rows$library_id), , drop = FALSE]
    set.seed(1)
    tracking[[i]] <- process_run(rows, s[i, ], run_dirs[i], threads)
  }
  tracking <- do.call(rbind, tracking)
  for (group in unique(m$analysis_group)) {
    export_group(m[m$analysis_group == group, , drop = FALSE], tracking,
                 run_dirs[s$analysis_group == group], file.path(partial, "groups", group))
  }
  if (file.exists(target) || !file.rename(partial, target)) fail("Cannot publish output ", target)
  message("Completed: ", target)
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
