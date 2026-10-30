#!/usr/bin/env Rscript
# Bidirectional MR at p < 5e-8, unpruned and with greedy 1 Mb distance thinning.
# Instruments are selected using the exposure GWAS in each direction.
# All alleles must be normalized to the same reference assembly and forward
# strand. Harmonization aligns effects using this fixed strand convention.
# Physical distance thinning does not establish a particular LD r2 threshold.

usage <- function() {
  cat(paste0(
    "Usage: Rscript /main/scripts/mr/run_mr.R [options]\n",
    "  --input-dir PATH         Prepared inputs [/main/results/mr/inputs]\n",
    "  --output-dir PATH        Results [/main/results/mr]\n",
    "  --include-relaxed        Also analyze p < 5e-6 sensitivity inputs\n",
    "  --plots                  Save scatterplots from the fitted harmonized data\n",
    "  --seed INTEGER           Reproducible bootstrap seed [20261030]\n"
  ))
}
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) { usage(); quit(status = 0) }
config <- list(input_dir = "/main/results/mr/inputs", output_dir = "/main/results/mr",
               include_relaxed = FALSE, plots = FALSE,
               seed = 20261030L)
i <- 1L
while (i <= length(args)) {
  option <- args[[i]]
  if (option %in% c("--include-relaxed", "--plots")) {
    config[[gsub("-", "_", substring(option, 3))]] <- TRUE
    i <- i + 1L
  } else if (option %in% c("--input-dir", "--output-dir", "--seed")) {
    if (i == length(args)) stop("Missing value for ", option)
    name <- gsub("-", "_", substring(option, 3))
    config[[name]] <- args[[i + 1L]]
    i <- i + 2L
  } else stop("Unknown option: ", option)
}
config$seed <- suppressWarnings(as.integer(config$seed))
if (is.na(config$seed) || config$seed < 0) stop("--seed must be a nonnegative R integer")
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop("Install TwoSampleMR before running this script; see /main/README.md.")
}
if (config$plots && !requireNamespace("ggplot2", quietly = TRUE)) stop("--plots requires ggplot2")

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
}
write_tsv <- function(data, path) {
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}
as_mr_data <- function(data, role) {
  TwoSampleMR::format_data(
    data, type = role, phenotype_col = role, snp_col = "SNP",
    beta_col = paste0("beta.", role), se_col = paste0("se.", role),
    eaf_col = paste0("eaf.", role), effect_allele_col = paste0("effect_allele.", role),
    other_allele_col = paste0("other_allele.", role), pval_col = paste0("pval.", role),
    id_col = paste0("id.", role), chr_col = paste0("chr.", role),
    pos_col = paste0("pos.", role)
  )
}
thin_1mb <- function(data) {
  # Sort by exposure P; reject every new variant within <= 1,000,000 bp
  # of any already-kept variant on the same chromosome.
  data <- data[order(data$pval.exposure, data$chr.exposure, data$pos.exposure, data$SNP), ]
  kept <- integer()
  for (j in seq_len(nrow(data))) {
    nearby <- data$chr.exposure[kept] == data$chr.exposure[j] &
      abs(data$pos.exposure[kept] - data$pos.exposure[j]) <= 1000000
    if (!any(nearby)) kept <- c(kept, j)
  }
  data$SNP[kept]
}

manifest <- read_tsv(file.path(config$input_dir, "analysis_inputs.tsv"))
required <- c("direction", "exposure_file", "outcome_file", "candidate_p_cutoff")
if (!all(required %in% names(manifest)) || !nrow(manifest)) stop("Invalid or empty input manifest")
if (anyDuplicated(manifest$direction) || any(!manifest$direction %in% c("chip_to_ad", "ad_to_chip"))) {
  stop("Input manifest must contain each requested direction exactly once")
}
thresholds <- if (config$include_relaxed) c(5e-8, 5e-6) else 5e-8
if (any(!is.finite(manifest$candidate_p_cutoff)) ||
    any(manifest$candidate_p_cutoff < max(thresholds))) {
  stop("Inputs do not cover the requested threshold. Prepare again with --include-relaxed if needed.")
}
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
harmonised_dir <- file.path(config$output_dir, "harmonised")
dir.create(harmonised_dir, showWarnings = FALSE)
if (config$plots) dir.create(file.path(config$output_dir, "plots"), showWarnings = FALSE)
all_results <- list(); all_retention <- list(); all_counts <- list()
analysis_number <- 0L
# Clear the two estimate tables with explicit headers before computing. A run
# with no estimable instruments must never leave estimates from an older run.
empty_results <- data.frame(analysis_id = character(), direction = character(),
                            method = character(), nsnp = integer(), b = numeric(),
                            se = numeric(), pval = numeric(), OR = numeric(),
                            OR_low_95 = numeric(), OR_high_95 = numeric())
write_tsv(empty_results, file.path(config$output_dir, "mr_results.tsv"))
write_tsv(empty_results, file.path(config$output_dir, "ivw_summary.tsv"))

for (m in seq_len(nrow(manifest))) {
  direction <- manifest$direction[m]
  exp_raw <- read_tsv(file.path(config$input_dir, manifest$exposure_file[m]))
  out_raw <- read_tsv(file.path(config$input_dir, manifest$outcome_file[m]))
  if (anyDuplicated(exp_raw$SNP) || anyDuplicated(out_raw$SNP)) stop("Duplicate SNP IDs in ", direction)
  for (threshold in thresholds) {
    significant <- exp_raw[exp_raw$pval.exposure < threshold, , drop = FALSE]
    # Restrict candidates to available outcome alleles before distance thinning.
    # Only EXPOSURE P-values select instruments; outcome P-values never do.
    common <- significant[significant$SNP %in% out_raw$SNP, , drop = FALSE]
    retained_1mb <- thin_1mb(common)
    for (thinning in c("unpruned", "1mb")) {
      analysis_number <- analysis_number + 1L
      tag <- paste(direction, format(threshold, scientific = TRUE), thinning, sep = "_")
      tag <- gsub("e-0", "e-", tag, fixed = TRUE)
      chosen <- if (thinning == "unpruned") common$SNP else retained_1mb
      audit <- data.frame(
        analysis_id = rep(tag, nrow(significant)), direction = rep(direction, nrow(significant)),
        p_threshold = rep(threshold, nrow(significant)),
        thinning = rep(thinning, nrow(significant)), SNP = significant$SNP,
        chr = significant$chr.exposure, pos = significant$pos.exposure,
        pval_exposure = significant$pval.exposure,
        matched_outcome = significant$SNP %in% out_raw$SNP,
        selected_for_harmonisation = significant$SNP %in% chosen,
        mr_keep = rep(FALSE, nrow(significant)), stringsAsFactors = FALSE
      )
      audit$exclusion_reason <- ifelse(!audit$matched_outcome, "outcome_unavailable_or_alleles_differ",
                                     ifelse(!audit$selected_for_harmonisation, "within_1mb_of_lower_p_variant", ""))
      counts <- data.frame(analysis_id = tag, direction = direction, p_threshold = threshold,
                           thinning = thinning, n_significant_exposure = nrow(significant),
                           n_matched_outcome = nrow(common), n_selected = length(chosen),
                           n_harmonised_kept = 0L, mean_F = NA_real_, min_F = NA_real_,
                           status = "no_matching_instruments", error_message = "",
                           stringsAsFactors = FALSE)
      if (length(chosen)) {
        exposure <- as_mr_data(common[match(chosen, common$SNP), , drop = FALSE], "exposure")
        outcome <- as_mr_data(out_raw[out_raw$SNP %in% chosen, , drop = FALSE], "outcome")
        dat <- TwoSampleMR::harmonise_data(exposure, outcome, action = 1L)
        # Persist the authoritative data for the estimates AND optional plots.
        write_tsv(dat, file.path(harmonised_dir, paste0(tag, ".tsv")))
        keep <- !is.na(dat$mr_keep) & dat$mr_keep
        kept_snps <- dat$SNP[keep]
        audit$mr_keep <- audit$SNP %in% kept_snps
        rejected <- audit$selected_for_harmonisation & !audit$mr_keep
        audit$exclusion_reason[rejected] <- "not_retained_by_harmonisation"
        counts$n_harmonised_kept <- length(unique(kept_snps))
        counts$status <- "no_instruments_after_harmonisation"
        if (length(kept_snps)) {
          strength <- (dat$beta.exposure[keep] / dat$se.exposure[keep])^2
          counts$mean_F <- mean(strength)
          counts$min_F <- min(strength)
          methods <- if (length(kept_snps) == 1L) "mr_wald_ratio" else "mr_ivw"
          if (length(kept_snps) >= 3L) methods <- c(methods, "mr_weighted_median", "mr_egger_regression")
          # Each scenario has its own fixed seed for bootstrap-based estimates.
          set.seed(as.integer((as.double(config$seed) + analysis_number) %% .Machine$integer.max))
          fit <- tryCatch(TwoSampleMR::mr(dat, method_list = methods), error = identity)
          if (inherits(fit, "error")) {
            counts$status <- "estimation_failed"
            counts$error_message <- gsub("[\t\r\n]", " ", conditionMessage(fit))
            fit <- NULL
          }
          if (!is.null(fit) && nrow(fit)) {
            fit$analysis_id <- tag
            fit$direction <- direction
            fit$p_threshold <- threshold
            fit$thinning <- thinning
            fit$n_harmonised_kept <- counts$n_harmonised_kept
            fit$OR <- exp(fit$b)
            fit$OR_low_95 <- exp(fit$b - 1.96 * fit$se)
            fit$OR_high_95 <- exp(fit$b + 1.96 * fit$se)
            all_results[[length(all_results) + 1L]] <- fit
            counts$status <- if (any(is.finite(fit$pval))) "estimated" else "estimation_returned_no_finite_p"
            if (config$plots && counts$status == "estimated") {
              plots <- TwoSampleMR::mr_scatter_plot(fit, dat[keep, , drop = FALSE])
              if (length(plots)) ggplot2::ggsave(
                file.path(config$output_dir, "plots", paste0(tag, ".pdf")),
                plot = plots[[1]], width = 7, height = 5
              )
            }
          } else if (counts$status != "estimation_failed") counts$status <- "no_estimate_returned"
        }
      }
      all_retention[[length(all_retention) + 1L]] <- audit
      all_counts[[length(all_counts) + 1L]] <- counts
      message(tag, ": ", counts$n_harmonised_kept, " retained; ", counts$status)
    }
  }
}
write_tsv(do.call(rbind, all_retention), file.path(config$output_dir, "instrument_retention.tsv"))
write_tsv(do.call(rbind, all_counts), file.path(config$output_dir, "analysis_summary.tsv"))
settings <- data.frame(harmonise_action = 1L, seed = config$seed,
                       distance_bp = 1000000, selection_order = "outcome_match_then_distance_thinning",
                       TwoSampleMR_version = as.character(utils::packageVersion("TwoSampleMR")))
write_tsv(settings, file.path(config$output_dir, "run_settings.tsv"))
writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "session_info.txt"))
if (!length(all_results)) stop("No MR estimates were produced; inspect analysis_summary.tsv and input_qc.tsv")
results <- do.call(rbind, all_results)
write_tsv(results, file.path(config$output_dir, "mr_results.tsv"))
write_tsv(results[results$method == "Inverse variance weighted", , drop = FALSE],
          file.path(config$output_dir, "ivw_summary.tsv"))
message("Results written to ", config$output_dir)
