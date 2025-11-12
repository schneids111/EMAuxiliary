# -------------------------------------------------------------------
# EMAuxiliary — lme4-style to BLIMP with auxiliaries (.mean route)
# -------------------------------------------------------------------

# Helpers: safe sd, ICC, level detection ---------------------------------------

.safe_sd <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

.icc_est <- function(x, id) {
  # One-way random effects ICC estimate from raw data (NA-tolerant)
  ok <- !is.na(x) & !is.na(id)
  x  <- x[ok]; id <- id[ok]
  if (!length(x)) return(NA_real_)
  # within-cluster variances and means
  v_within <- tapply(x, id, function(z) stats::var(z, na.rm = TRUE))
  m_by     <- tapply(x, id, function(z) mean(z, na.rm = TRUE))
  vw <- mean(v_within, na.rm = TRUE)
  vb <- stats::var(m_by, na.rm = TRUE)
  if (!is.finite(vw) || !is.finite(vb)) return(NA_real_)
  tot <- vb + vw
  if (tot <= 0) return(NA_real_)
  vb / tot
}

.is_level2_var <- function(vname, data, id) {
  if (!vname %in% names(data)) return(FALSE)
  sds <- tapply(data[[vname]], data[[id]], function(z) stats::sd(z, na.rm = TRUE))
  all(is.na(sds) | sds < 1e-8)
}

# Token utilities ---------------------------------------------------------------

# Convert raw term labels (possibly with .mean already present) to BLIMP tokens
# We keep ":" while validating, then later convert ":" -> "*".
.tokens_from_terms <- function(terms_vec) {
  # terms_vec is like c("x", "x:x.mean", "m", "x.mean")
  terms_vec[is.na(terms_vec)] <- ""
  terms_vec[nzchar(terms_vec)]
}

# Extract "base" predictor names that appear as main effects (no interactions)
.base_terms_from_formula <- function(fml) {
  tl_raw <- attr(stats::terms(lme4::nobars(fml)), "term.labels")
  main   <- unique(unlist(strsplit(tl_raw, ":", fixed = TRUE)))
  main   <- main[nzchar(main)]
  # strip ".mean" suffix from bases
  sub("\\.mean$", "", main)
}

# Standardize selected columns in-place (z-score)
.standardize_inplace <- function(df, cols) {
  zcols <- intersect(cols, names(df))
  if (!length(zcols)) return(list(data = df, changed = character(0)))
  changed <- character(0)
  for (v in zcols) {
    x <- df[[v]]
    if (is.numeric(x)) {
      sdv <- stats::sd(x, na.rm = TRUE)
      if (is.finite(sdv) && sdv > 0) {
        df[[v]] <- as.numeric(scale(x)[,1])
        changed <- c(changed, v)
      }
    }
  }
  list(data = df, changed = unique(changed))
}

# Mean token chooser for auxiliaries:
# - if predictor is pure L2 => use raw name (no ".mean")
# - else => use "name.mean"
.mean_token <- function(xname, data, id) {
  if (.is_level2_var(xname, data, id)) xname else paste0(xname, ".mean")
}

# BLIMP last PSR reader (final row's Highest PSR)
.blimp_last_psr <- function(blimp_fit) {
  out <- utils::capture.output(rblimp::output(blimp_fit))
  start <- grep("^\\s*BURN-IN POTENTIAL SCALE REDUCTION \\(PSR\\) OUTPUT:", out)
  if (!length(start)) return(NA_real_)
  # Find the last "Comparing iterations" line in the block
  block <- out[(start[1]):length(out)]
  cmp   <- grep("^\\s*Comparing iterations", block)
  if (!length(cmp)) return(NA_real_)
  last  <- block[tail(cmp, 1)]
  # Extract the number just before "Parameter"
  m <- regmatches(last, regexpr("\\s([0-9]+\\.[0-9]+)\\s+Parameter", last))
  as.numeric(trimws(sub("Parameter.*$", "", sub("^\\s*", "", m))))
}

# ------------------------------------------------------------------------------
# EMAuxiliary(): main wrapper
# ------------------------------------------------------------------------------

EMAuxiliary <- function(
  formula,
  data,
  id,
  aux = NULL,
  ordinal = NULL,
  nominal = NULL,
  center_group = NULL,   # variables to group-mean center (enables x.mean)
  center_grand = NULL,   # variables to grand-mean center (not related to .mean)
  burn = 5000,
  iter = 10000,
  chains = 2,
  seed = 12345
) {
  stopifnot(is.data.frame(data))
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Install lme4.")
  if (!requireNamespace("rblimp", quietly = TRUE)) stop("Install rblimp.")

  # Fixed part and random part
  fixed_form <- lme4::nobars(formula)
  fb         <- lme4::findbars(formula)      # random terms
  rand_slopes <- character(0)
  if (length(fb)) {
    inside_terms <- unique(unlist(lapply(fb, function(b) {
      left <- deparse(b[[2]])
      attr(stats::terms(as.formula(paste("~", left))), "term.labels")
    })))
    rand_slopes <- setdiff(inside_terms, "1")
  }

  # Outcome name
  y <- deparse(formula[[2L]])

  # ---------------- Build RHS with .mean validation ---------------------------

  # Raw term labels with ":" to split later; keep as-is
  tl_raw <- attr(stats::terms(fixed_form), "term.labels")
  tl_raw <- tl_raw[nzchar(tl_raw)]
  rhs_tokens_colon <- .tokens_from_terms(tl_raw)

  # Split interactions into atomic bits to find any *.mean in use
  rhs_atoms <- unique(unlist(strsplit(rhs_tokens_colon, ":", fixed = TRUE)))
  rhs_atoms <- rhs_atoms[nzchar(rhs_atoms)]
  means_in_rhs <- unique(sub("\\.mean$", "", rhs_atoms[grepl("\\.mean$", rhs_atoms)]))

  # Validation: every base that appears as *.mean must be in center_group (or be y)
  cg_allow <- unique(c(center_group, y)) # y.mean is allowed even if user didn't list it
  cg_allow <- cg_allow[!is.na(cg_allow) & nzchar(cg_allow)]
  missing_bases <- setdiff(means_in_rhs, cg_allow)
  if (length(missing_bases)) {
    stop(
      "You used ", paste0(paste0(missing_bases, ".mean"), collapse = ", "),
      " in the formula but did not include the base name(s) in `center_group`.\n",
      "Add: center_group = c(", paste(sprintf('\"%s\"', missing_bases), collapse = ", "), ", ...).",
      call. = FALSE
    )
  }

  # Now produce the '*' BLIMP RHS tokens
  blimp_terms <- if (length(rhs_tokens_colon)) gsub(":", "*", rhs_tokens_colon, fixed = TRUE) else character(0)
  rhs_tokens  <- unique(blimp_terms)
  fixed_rhs   <- if (length(rhs_tokens)) paste(rhs_tokens, collapse = " ") else "1"

  # ---------------- Random slopes: forbid ordinal/nominal and *.mean ----------

  slope_vars <- rand_slopes
  if (length(slope_vars)) {
    # strip possible ".mean" if someone wrote it (random slopes must be L1)
    if (any(grepl("\\.mean\\b", slope_vars))) {
      stop("Random slopes may only be on Level-1 predictors; do not use '*.mean' in random-slope terms.")
    }
    # Disallow ordinal/nominal as random slopes
    bad_on_slopes <- intersect(slope_vars, unique(c(ordinal %||% character(0), nominal %||% character(0))))
    if (length(bad_on_slopes)) {
      stop("Ordinal/nominal variables cannot appear as random-slope terms: ",
           paste(bad_on_slopes, collapse = ", "))
    }
  }

  # ---------------- Build CENTER blocks --------------------------------------

  # We pass exactly what the user requested for predictors:
  # - center_group → groupmean
  # - center_grand → grandmean
  # Additionally, if auxiliaries are used, we add y to groupmean so
  # we can legally use y.mean in auxiliary RHS (no effect on focal RHS scaling).
  center_group_to_pass <- unique(c(center_group, if (length(aux)) y))
  center_group_to_pass <- center_group_to_pass[!is.na(center_group_to_pass) & nzchar(center_group_to_pass)]
  center_grand_to_pass <- unique(center_grand)
  center_grand_to_pass <- center_grand_to_pass[!is.na(center_grand_to_pass) & nzchar(center_grand_to_pass)]

  center_lines <- character(0)
  if (length(center_group_to_pass))
    center_lines <- c(center_lines, paste0("groupmean = ", paste(center_group_to_pass, collapse = " ")))
  if (length(center_grand_to_pass))
    center_lines <- c(center_lines, paste0("grandmean = ", paste(center_grand_to_pass, collapse = " ")))
  center_block <- if (length(center_lines)) paste0("CENTER: ", paste(center_lines, collapse = "; "), ";") else NULL

  # ---------------- Auxiliary model construction ------------------------------

  # Base predictors (no interactions): their raw names (strip .mean suffix)
  base_preds <- .base_terms_from_formula(formula)

  # Determine L1 vs L2 for auxiliaries
  aux <- aux %||% character(0)
  aux <- unique(aux[nzchar(aux)])

  # ICC and "near pure L1" test: ICC < .05
  aux_icc <- setNames(rep(NA_real_, length(aux)), aux)
  for (a in aux) {
    if (a %in% names(data)) aux_icc[a] <- .icc_est(data[[a]], data[[id]])
  }
  near_L1 <- names(aux_icc)[which(is.finite(aux_icc) & aux_icc < 0.05)]

  # Predictor pools for aux RHS:
  # Mixed-level aux → use full RHS (including interactions) + y + y.mean
  full_rhs_for_aux <- rhs_tokens_colon
  # L1-only aux (near_L1) → y + base x’s only (no .mean)
  # Pure L2 aux → y.mean + L2 predictors / means of mixed predictors

  # Helper to detect if an aux is L2-only
  aux_is_L2 <- function(a) .is_level2_var(a, data, id)

  # Build RHS per aux variable
  aux_lines <- character(0)
  for (a in aux) {
    if (aux_is_L2(a)) {
      # Pure L2 aux: only y.mean and L2 predictors (raw L2 names or x.mean forms)
      L2_preds <- character(0)
      for (p in base_preds) {
        if (.is_level2_var(p, data, id)) {
          L2_preds <- c(L2_preds, p)  # pure L2 stays raw
        } else {
          # mixed predictor: use p.mean
          L2_preds <- c(L2_preds, paste0(p, ".mean"))
        }
      }
      rhs <- unique(c(paste0(y, ".mean"), L2_preds))
      rhs <- rhs[nzchar(rhs)]
      if (length(rhs))
        aux_lines <- c(aux_lines, paste0(a, " ~ ", paste(rhs, collapse = " "), ";"))

    } else if (a %in% near_L1) {
      # Near pure L1 aux: only y and base x’s (no .mean, no interactions)
      rhs <- unique(c(y, base_preds))
      rhs <- rhs[rhs != a & nzchar(rhs)]
      if (length(rhs))
        aux_lines <- c(aux_lines, paste0(a, " ~ ", paste(rhs, collapse = " "), ";"))

    } else {
      # Mixed-level aux: y + y.mean + full focal RHS (incl. interactions)
      rhs <- unique(c(y, paste0(y, ".mean"), full_rhs_for_aux))
      rhs <- rhs[rhs != a & nzchar(rhs)]
      if (length(rhs))
        aux_lines <- c(aux_lines, paste0(a, " ~ ", paste(rhs, collapse = " "), ";"))
    }
  }

  aux_block <- if (length(aux_lines)) paste0("auxiliary.model:\n  ", paste(aux_lines, collapse = "\n  ")) else NULL
  aux_model_arg <- if (length(aux_lines)) paste(aux_lines, collapse = "\n  ") else NULL

  # ---------------- VARIABLES list (observed only) ----------------------------

  # All variables referenced as observed in BLIMP (exclude .mean tokens)
  rhs_observed_atoms <- unique(unlist(strsplit(rhs_tokens_colon, ":", fixed = TRUE)))
  rhs_observed_bases <- unique(sub("\\.mean$", "", rhs_observed_atoms))
  all_vars <- unique(c(id, y, rhs_observed_bases, aux, ordinal, nominal))
  all_vars <- all_vars[nzchar(all_vars)]

  # ---------------- Auto-standardize continuous auxiliaries -------------------

  data_run <- data
  cont_aux <- aux[aux %in% names(data_run)]
  # avoid scaling ordinal/nominal auxiliaries
  cat_vars <- unique(c(ordinal %||% character(0), nominal %||% character(0)))
  cont_aux <- setdiff(cont_aux, cat_vars)
  # keep only numeric
  cont_aux <- cont_aux[vapply(cont_aux, function(v) is.numeric(data_run[[v]]), logical(1))]
  if (length(cont_aux)) {
    std_res <- .standardize_inplace(data_run, cont_aux)
    if (length(std_res$changed)) {
      message("Standardized auxiliaries (z): ", paste(std_res$changed, collapse = ", "))
    }
    data_run <- std_res$data
  }

  # ---------------- Compose BLIMP code (for printing only) --------------------

  model_line <- if (!length(slope_vars)) {
    paste0("MODEL:\n  ", y, " ~ ", fixed_rhs, ";")
  } else {
    paste0("MODEL:\n  ", y, " ~ ", fixed_rhs, " | ", paste(slope_vars, collapse = " "), ";")
  }

  blimp_code <- paste(
    "DATA: <in-memory by rblimp>;",
    paste0("VARIABLES: ", paste(all_vars, collapse = " "), ";"),
    paste0("CLUSTERID: ", id, ";"),
    if (length(ordinal %||% character(0))) paste0("ORDINAL: ", paste(ordinal, collapse = " "), ";") else NULL,
    if (length(nominal %||% character(0))) paste0("NOMINAL: ", paste(nominal, collapse = " "), ";") else NULL,
    if (!is.null(center_block)) center_block else NULL,
    model_line,
    if (!is.null(aux_block)) aux_block else NULL,
    paste0("SEED: ", seed, ";"),
    paste0("BURN: ", burn, ";"),
    paste0("ITER: ", iter, ";"),
    paste0("CHAINS: ", chains, ";"),
    sep = "\n"
  )

  # ---------------- Run BLIMP -------------------------------------------------

  args <- list(
    data      = data_run,
    clusterid = id,
    model     = paste0(y, " ~ ", fixed_rhs, if (length(slope_vars)) paste0(" | ", paste(slope_vars, collapse = " ")) else ""),
    seed      = seed,
    burn      = burn,
    iter      = iter,
    chains    = chains
  )
  if (length(center_group_to_pass) || length(center_grand_to_pass)) {
    center_pieces <- character(0)
    if (length(center_group_to_pass))
      center_pieces <- c(center_pieces, paste0("groupmean = ", paste(center_group_to_pass, collapse = " ")))
    if (length(center_grand_to_pass))
      center_pieces <- c(center_pieces, paste0("grandmean = ", paste(center_grand_to_pass, collapse = " ")))
    args$center <- paste(center_pieces, collapse = "; ")
  }
  if (length(ordinal %||% character(0))) args$ordinal <- paste(ordinal, collapse = " ")
  if (length(nominal %||% character(0))) args$nominal <- paste(nominal, collapse = " ")

  # pass auxiliary model through a recognized arg name if available
  if (!is.null(aux_model_arg)) {
    fm <- names(formals(rblimp::rblimp))
    aux_names <- c("auxiliary.model", "auxmodel", "auxiliarymodel", "auxiliary")
    matched <- intersect(aux_names, fm)
    if (length(matched)) {
      args[[matched[1]]] <- aux_model_arg
    } else {
      # inline fallback
      args$model <- paste0(args$model, ";\nauxiliary.model:\n  ", aux_model_arg, ";")
    }
  }

  fit <- rblimp::rblimp(
    data      = args$data,
    clusterid = args$clusterid,
    model     = args$model,
    seed      = args$seed,
    burn      = args$burn,
    iter      = args$iter,
    chains    = args$chains,
    center    = args$center %||% NULL,
    ordinal   = args$ordinal %||% NULL,
    nominal   = args$nominal %||% NULL,
    `auxiliary.model` = args[["auxiliary.model"]] %||% NULL
  )

  # ---------------- PSR + scale audit (only if PSR > 1.05) -------------------

  psr <- tryCatch(.blimp_last_psr(fit), error = function(e) NA_real_)
  if (is.finite(psr) && psr > 1.05) {
    header <- sprintf(
      "Final PSR is %.3f. The final PSR should be < 1.05.\nResults may be unreliable. Consider longer burn/iter, more chains, or simplifying auxiliaries.\n",
      psr
    )
    # scale audit on focal variables (y + base predictors)
    base_preds_for_audit <- base_preds[base_preds %in% names(data)]
    focal_cols <- unique(c(y, base_preds_for_audit))
    sds <- vapply(focal_cols, function(v) .safe_sd(data[[v]]), numeric(1))
    sd_ok <- sds[is.finite(sds) & sds > 0]
    ratio_thresh <- 6
    abs_thresh   <- 3
    addendum <- NULL
    if (length(sd_ok) >= 1) {
      max_sd <- max(sd_ok); min_sd <- min(sd_ok)
      ratio  <- if (length(sd_ok) >= 2) max_sd/min_sd else 1
      trigger <- (length(sd_ok) >= 2 && ratio > ratio_thresh) || any(sd_ok > abs_thresh)
      if (trigger) {
        sd_tbl <- paste(sprintf("  %s: %.3f", names(sds), sds), collapse = "\n")
        addendum <- paste0(
          "Potential cause for poor PSR: large scale differences among focal variables.\n",
          sprintf("  max SD = %.3f, min SD = %.3f, max/min = %.2f (thresholds: ratio > %g OR any SD > %g)\n",
                  max_sd, min_sd, ratio, ratio_thresh, abs_thresh),
          "Focal SDs:\n", sd_tbl, "\n\n",
          "Tip: Consider rescaling (e.g., z-scaling) outcome and largest-SD predictors BEFORE calling EMAuxiliary()."
        )
      }
    }
    warning(paste0(header, if (!is.null(addendum)) paste0("\n", addendum) else ""), call. = FALSE)
  }

  # Return
  structure(
    list(
      blimp_code = blimp_code,
      fit        = fit,
      meta       = list(y = y)
    ),
    class = "EMAuxiliary_fit"
  )
}

# ------------------------------------------------------------------------------
# Printers
# ------------------------------------------------------------------------------

# Accept either wrapper object or raw blimp_obj
.get_blimp_fit <- function(obj) if (!is.null(obj$fit)) obj$fit else obj

blimp_print_psr <- function(obj) {
  fit <- .get_blimp_fit(obj)
  out <- utils::capture.output(rblimp::output(fit))

  start <- grep("^\\s*BURN-IN POTENTIAL SCALE REDUCTION \\(PSR\\) OUTPUT:", out)
  if (!length(start)) {
    cat("PSR section not found.\n")
    return(invisible())
  }
  stop1 <- grep("^\\s*METROPOLIS-HASTINGS ACCEPTANCE RATES", out)
  stop2 <- grep("^\\s*$", out)
  stops <- c(stop1[1], stop2[stop2 > start[1]][1], length(out))
  end <- min(stops[is.finite(stops)], na.rm = TRUE)
  cat(paste(out[start[1]:end], collapse = "\n"), "\n")
}

blimp_print_focal <- function(obj) {
  fit <- .get_blimp_fit(obj)
  out <- utils::capture.output(rblimp::output(fit))
  y <- if (!is.null(obj$meta) && !is.null(obj$meta$y)) obj$meta$y else {
    m <- regmatches(out, regexpr("^\\s*Outcome Variable:\\s+\\S+", out))
    if (length(m) && nchar(m[1])) sub("^\\s*Outcome Variable:\\s+", "", m[1]) else ""
  }
  if (!nchar(y)) {
    cat("Could not infer focal outcome name.\n")
    return(invisible())
  }
  y_esc <- gsub("([][{}()+*.^$|?\\\\])", "\\\\\\1", y, perl = TRUE)

  start <- grep(paste0("^\\s*Outcome Variable:\\s+", y_esc, "\\b"), out)
  if (!length(start)) {
    cat("Focal section not found for outcome '", y, "'.\n", sep = "")
    return(invisible())
  }
  next_outcome <- grep("^\\s*Outcome Variable:\\s+\\S+", out)
  next_outcome <- next_outcome[next_outcome > start[1]]
  pred_est     <- grep("^\\s*PREDICTOR MODEL ESTIMATES:", out)
  pred_est     <- pred_est[pred_est > start[1]]
  aux_block    <- grep("^\\s+auxiliary\\.model block:", out)

  candidates <- c(next_outcome[1], pred_est[1], aux_block[1], length(out))
  end <- min(candidates[is.finite(candidates)], na.rm = TRUE) - 1L
  if (!is.finite(end) || end < start[1]) end <- length(out)

  cat(paste(out[start[1]:end], collapse = "\n"), "\n")
}

# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b
