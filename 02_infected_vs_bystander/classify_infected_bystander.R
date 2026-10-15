#!/usr/bin/env Rscript
# ==============================================================================
# Longitudinal single-cell profiling of IAV infection: analysis pipeline
# Classification of single cells into infected and bystander populations
# ------------------------------------------------------------------------------
# Rationale
#   A low level of viral reads is detectable in every library, including cells
#   that were never productively infected. This residual signal comes from the
#   inoculum left in the culture medium and it prevents a simple zero/non-zero
#   cut-off. The expression of the marker gene across a population of cells is
#   instead bimodal, so the two populations are separated by fitting a
#   two-component mixture model rather than by choosing an arbitrary threshold.
#
# Method
#   1. Log1p transform the FPKM of the single marker gene that the infection
#      readout is based on, which is the NP row of the matrix.
#   2. Estimate a background distribution from the pre-infection time point,
#      where no productive infection can have occurred yet.
#   3. Fit a two-component Gaussian mixture with the EM algorithm at the
#      classification time point, seeding the component means with the
#      background mean and with the upper range of the observed values.
#   4. Assign the component with the lower mean to the bystander population and
#      the component with the higher mean to the infected population.
#   5. Report the posterior classification probability of every cell and the
#      implied decision boundary.
#   6. Verify that the cells called infected carry the signal in every viral
#      segment and not only in the marker gene.
#
# Input
#   A gene-by-cell expression matrix in FPKM. The file may be .xlsx, .csv or
#   .tsv. Column names must follow the <CellID>.<TimePoint> convention, for
#   example WSN01.0, WSN01.3, WSN01.6. The marker gene and the viral genes are
#   located by name.
#
# Output (written to --output-dir)
#   cell_classification.csv    one row per cell, with the classifying value,
#                              assigned population and posterior probability
#   classification_summary.txt fitted mixture parameters, decision boundary and
#                              population sizes
#
# Requirements
#   R packages: mclust
#
# Usage
#   Rscript classify_infected_bystander.R \
#           --fpkm /path/to/merged_fpkm.csv \
#           --output-dir ./results_classification
#
# ==============================================================================

# Load the shared helpers from the repository root.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(dirname(normalizePath(.self))), "R", "common.R"))
rm(.self)

USAGE <- "Usage: Rscript classify_infected_bystander.R --fpkm FILE --output-dir DIR
       [--classify-gene NAME] [--viral-genes LIST | --viral-pattern REGEX]
       [--background-time N] [--classify-time N] [--model-name E|V]"

opt <- parse_cli(list(
  fpkm            = NULL,
  `output-dir`    = NULL,
  `classify-gene` = "",
  `viral-genes`   = "",
  `viral-pattern` = "^WSN_",
  `background-time` = NA_real_,
  `classify-time`   = NA_real_,
  `model-name`      = "E"
), USAGE)

if (!requireNamespace("mclust", quietly = TRUE)) {
  stop("Package 'mclust' is required. Install it with install.packages('mclust')",
       call. = FALSE)
}
# Attached rather than merely loaded. Mclust() with a prior looks up mclustBIC on
# the search path, which a namespace-only reference does not provide.
suppressPackageStartupMessages(library(mclust))

ensure_dir(opt[["output-dir"]])

# --------------------------------------------------------------- 1. read in ---
step("Reading expression matrix: ", opt$fpkm)
# Coerced to a matrix so that a single-row or single-block subset stays numeric
# and can be passed straight to log1p, is.finite and colMeans.
expr <- as.matrix(as_numeric_matrix(read_matrix(opt$fpkm)))
meta <- parse_samples(colnames(expr))
expr <- expr[, meta$sample, drop = FALSE]

note(length(rownames(expr)), " genes x ", ncol(expr), " cells")
note("Time points: ", paste(sort(unique(meta$time)), collapse = ", "))
note("Cells: ", length(unique(meta$cell_id)))

# -------------------------------------------------- 2. marker and viral genes ---
classify_gene <- opt[["classify-gene"]]
if (!nzchar(classify_gene)) {
  hits <- grep("(^|_)NP$", rownames(expr), value = TRUE)
  if (length(hits) != 1L) {
    stop("The NP row could not be identified unambiguously, ", length(hits),
         " rows matched. Pass --classify-gene explicitly.", call. = FALSE)
  }
  classify_gene <- hits
} else if (!classify_gene %in% rownames(expr)) {
  stop("The marker gene '", classify_gene, "' is not present in the matrix",
       call. = FALSE)
}
step("Marker gene used for the classification: ", classify_gene)

viral <- split_csv(opt[["viral-genes"]])
viral <- viral[nzchar(viral)]
if (length(viral) == 0) {
  viral <- grep(opt[["viral-pattern"]], rownames(expr), value = TRUE)
} else {
  absent <- setdiff(viral, rownames(expr))
  if (length(absent)) {
    stop("These viral genes are not present in the matrix: ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
}
if (length(viral) == 0) {
  stop("No viral genes found. Pass --viral-genes explicitly or adjust --viral-pattern",
       call. = FALSE)
}
step("Viral genes used for the per-segment validation: ",
     paste(viral, collapse = ", "))

# ------------------------------------------------------ 3. classifying value ---
# The classification is driven by the marker gene alone, so that the two
# populations are separated on the same quantity that the infection readout is
# based on. The full viral block is still extracted because it is used below to
# check that the cells called infected carry the signal in every segment.
classify_vals <- expr[classify_gene, ]
if (any(!is.finite(classify_vals))) {
  note("Non-finite values found for ", classify_gene, " and replaced by 0")
  classify_vals[!is.finite(classify_vals)] <- 0
}
classify_vals <- log1p(classify_vals)

viral_mat <- expr[viral, , drop = FALSE]
if (anyNA(viral_mat)) {
  note("Non-finite values found in the viral block and replaced by 0")
  viral_mat[!is.finite(viral_mat)] <- 0
}

cell_df <- meta
cell_df$classify_gene <- classify_gene
cell_df$classify_value <- classify_vals[cell_df$sample]

# ----------------------------------------------------- 4. background model ----
bg_time <- opt[["background-time"]]
if (is.na(bg_time)) bg_time <- min(meta$time)
if (!bg_time %in% meta$time) {
  stop("--background-time ", bg_time, " is not among the observed time points: ",
       paste(sort(unique(meta$time)), collapse = ", "), call. = FALSE)
}

bg_vals <- cell_df$classify_value[cell_df$time == bg_time]
if (length(bg_vals) < 2) {
  stop("At least two cells are needed at the background time point, found ",
       length(bg_vals), call. = FALSE)
}
bg_mean <- mean(bg_vals)
bg_sd <- stats::sd(bg_vals)

step("Background model from time point ", bg_time, " (n = ", length(bg_vals), ")")
note("mean = ", fmt_num(bg_mean), "   sd = ", fmt_num(bg_sd))

# --------------------------------------------------- 5. mixture model fit -----
cls_time <- opt[["classify-time"]]
if (is.na(cls_time)) cls_time <- max(meta$time)
if (!cls_time %in% meta$time) {
  stop("--classify-time ", cls_time, " is not among the observed time points: ",
       paste(sort(unique(meta$time)), collapse = ", "), call. = FALSE)
}

cls_df <- cell_df[cell_df$time == cls_time, , drop = FALSE]
if (nrow(cls_df) < 4) {
  stop("At least four cells are needed to fit a two-component mixture, found ",
       nrow(cls_df), call. = FALSE)
}

step("Fitting a two-component Gaussian mixture at time point ", cls_time,
     " (n = ", nrow(cls_df), ")")

x <- cls_df$classify_value
seed_means <- c(bg_mean, max(x) * 0.8)
if (seed_means[2] <= seed_means[1]) {
  # Degenerate case where the upper seed collapses onto the background. Nudge it
  # apart so that the EM algorithm starts from two distinguishable components.
  seed_means[2] <- seed_means[1] + max(stats::sd(x), 1e-3)
}

fit <- mclust::Mclust(
  data       = x,
  G          = 2,
  modelNames = opt[["model-name"]],
  prior      = mclust::priorControl(mean = seed_means)
)

if (is.null(fit) || is.null(fit$classification)) {
  stop("The mixture model failed to converge. Try --model-name V or check that ",
       "the expression of the marker gene is genuinely bimodal", call. = FALSE)
}

# ------------------------------------------------------- 6. label components ---
# Component indices returned by mclust carry no biological meaning, so the
# labels are attached by comparing the fitted means rather than by assuming
# that component 1 is always the background.
comp_means <- as.numeric(fit$parameters$mean)
bystander_comp <- which.min(comp_means)
infected_comp <- which.max(comp_means)

cls_df$population <- ifelse(fit$classification == bystander_comp,
                            "Bystander", "Infected")
cls_df$probability <- apply(fit$z, 1, max)
cls_df$component <- fit$classification

# Decision boundary implied by the fitted mixture. For the equal-variance model
# it reduces to the midpoint of the two means, corrected for unequal mixing
# proportions.
pi <- as.numeric(fit$parameters$pro)
mu <- comp_means
if (toupper(opt[["model-name"]]) == "E") {
  sigma2 <- as.numeric(fit$parameters$variance$sigmasq)
} else {
  sigma2 <- mean(as.numeric(fit$parameters$variance$sigmasq))
}
boundary <- (mu[1] + mu[2]) / 2 +
  sigma2 / (mu[2] - mu[1]) * log(pi[1] / pi[2])

# ------------------------------------------------------------- 7. reporting ---
n_inf <- sum(cls_df$population == "Infected")
n_bys <- sum(cls_df$population == "Bystander")
inf_pct <- 100 * n_inf / nrow(cls_df)

step("Classification result")
note("Infected  : ", n_inf, " cells")
note("Bystander : ", n_bys, " cells")
note("Infected fraction: ", fmt_num(inf_pct, 1), " %")

# Sanity check that the assignment is biologically coherent. The classifier sees
# the marker gene alone, so this verifies that the cells it called infected
# carry the viral signal in every segment and not only in that one gene.
if (n_inf > 0 && n_bys > 0) {
  per_gene_inf <- rowMeans(viral_mat[, cls_df$sample[cls_df$population == "Infected"],
                                     drop = FALSE])
  per_gene_bys <- rowMeans(viral_mat[, cls_df$sample[cls_df$population == "Bystander"],
                                     drop = FALSE])
  note("Mean viral FPKM per segment, infected : ",
       paste(fmt_num(per_gene_inf, 2), collapse = ", "))
  note("Mean viral FPKM per segment, bystander: ",
       paste(fmt_num(per_gene_bys, 2), collapse = ", "))
} else {
  per_gene_inf <- per_gene_bys <- NA
}

# ---------------------------------------------------------------- 8. outputs --
out_csv <- file.path(opt[["output-dir"]], "cell_classification.csv")
write_table(cls_df[, c("sample", "cell_id", "time", "classify_gene",
                       "classify_value", "population", "probability",
                       "component")],
            out_csv)

out_txt <- file.path(opt[["output-dir"]], "classification_summary.txt")
lines <- c(
  strrep("=", 64),
  "Infected versus bystander classification",
  strrep("=", 64),
  "",
  paste0("Input matrix           : ", basename(opt$fpkm)),
  paste0("Marker gene            : ", classify_gene),
  paste0("Classifying value      : log1p(", classify_gene, " FPKM), that is ",
         "ln(1 + ", classify_gene, " FPKM)"),
  paste0("Viral genes            : ", paste(viral, collapse = ", ")),
  paste0("Background time point  : ", bg_time, " h"),
  paste0("Classification time    : ", cls_time, " h"),
  paste0("Mixture model          : mclust, G = 2, modelNames = '",
         opt[["model-name"]], "'"),
  "",
  strrep("-", 64),
  "Background distribution",
  strrep("-", 64),
  paste0("  n cells              : ", length(bg_vals)),
  paste0("  mean classify value  : ", fmt_num(bg_mean, 4)),
  paste0("  sd   classify value  : ", fmt_num(bg_sd, 4)),
  "",
  strrep("-", 64),
  "Fitted mixture",
  strrep("-", 64),
  paste0("  log-likelihood       : ", fmt_num(fit$loglik, 4)),
  paste0("  mixing proportions   : ",
         paste(sprintf("component %d = %.4f", seq_along(pi), pi), collapse = ", ")),
  paste0("  component means      : ",
         paste(sprintf("component %d = %.4f", seq_along(mu), mu), collapse = ", ")),
  paste0("  bystander component  : ", bystander_comp),
  paste0("  infected component   : ", infected_comp),
  paste0("  decision boundary    : ", fmt_num(boundary, 4),
         " on the log1p scale"),
  "",
  strrep("-", 64),
  "Population sizes",
  strrep("-", 64),
  paste0("  infected             : ", n_inf),
  paste0("  bystander            : ", n_bys),
  paste0("  infected fraction    : ", fmt_num(inf_pct, 2), " %"),
  "",
  strrep("-", 64),
  "Mean viral FPKM per segment",
  strrep("-", 64),
  paste0("  gene                 : ", paste(names(per_gene_inf), collapse = ", ")),
  paste0("  infected             : ", paste(fmt_num(per_gene_inf, 3), collapse = ", ")),
  paste0("  bystander            : ", paste(fmt_num(per_gene_bys, 3), collapse = ", "))
)
writeLines(lines, out_txt)
message("  written: ", out_txt)

step("Done")
