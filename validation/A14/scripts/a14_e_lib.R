
A14_SEED <- 20261118L

a14_set_libpaths <- function(root) {
  rl <- file.path(root, "tools", "Rlib")
  .libPaths(c(rl, .libPaths()))
  invisible(.libPaths())
}

a14_load_packages <- function() {
  suppressPackageStartupMessages({
    library(ape); library(phylolm); library(data.table)
  })
}

a14_bh <- function(p) {            # Benjamini-Hochberg adjusted p-values, input order kept
  m <- length(p); o <- order(p); q <- numeric(m); prev <- 1
  for (r in m:1) { i <- o[r]; prev <- min(prev, p[i] * m / r); q[i] <- prev }
  q
}
a14_z <- function(v) (v - mean(v)) / sd(v)

a14_mde_or <- function(n, prevalence, neff_frac = 1, alpha = 0.05, power = 0.8, vif = 1) {
  se <- sqrt(vif / (n * neff_frac * prevalence * (1 - prevalence)))
  exp((qnorm(1 - alpha / 2) + qnorm(power)) * se)
}

a14_read_s2 <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    x <- as.data.frame(readxl::read_excel(path, sheet = "Table S2", .name_repair = "minimal"), stringsAsFactors = FALSE)
    x <- x[, 1:12]
  } else {
    x <- as.data.frame(data.table::fread(path, sep = "\t", data.table = FALSE))
  }
  names(x) <- sub("^X_FC_resistance$", "FC5_resistance", names(x))
  keep <- !is.na(x$assembly) & !grepl("^=", x$assembly) & !grepl("^=", x$species)
  x <- x[keep, , drop = FALSE]
  drugs <- grep("_resistance$", names(x), value = TRUE)
  for (d in drugs) x[[d]] <- as.integer(as.character(x[[d]]))
  x$species <- as.character(x$species); x$assembly <- as.character(x$assembly); x$Order <- as.character(x$Order)
  x
}

a14_read_names <- function(path) {
  x <- if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE) else
    as.data.frame(data.table::fread(path, sep = "\t", data.table = FALSE))
  x[!duplicated(x$assembly_fullID_updated), ]
}

a14_count_pumps <- function(retrieval_path, og_member_path = NULL, frag_pairs_path = NULL) {
  rt <- data.table::fread(retrieval_path, sep = "\t", data.table = FALSE, select = c("species_code", "genome_kind", "passes_rule_final_seq", "cdhit99_relation"))
  rt <- rt[rt$genome_kind == "y1000_genome", ]
  pass <- rt[rt$passes_rule_final_seq == "Y", ]
  cnt <- as.data.frame(table(species_code = pass$species_code), stringsAsFactors = FALSE)
  names(cnt)[2] <- "n_pass"
  allg <- data.frame(species_code = unique(rt$species_code), stringsAsFactors = FALSE)
  out <- merge(allg, cnt, by = "species_code", all.x = TRUE)
  out$n_pass[is.na(out$n_pass)] <- 0L
  if (!is.null(og_member_path)) {
    og <- data.table::fread(og_member_path, sep = "\t", data.table = FALSE, select = c("species_code", "family", "pass_flag"))
    og <- og[og$pass_flag == "PASS", ]
    fam <- as.data.frame.matrix(table(og$species_code, og$family))
    colnames(fam) <- c("n_pdr11aus1like", "n_pdr12like", "n_pdr5like", "n_snq2like")[match(colnames(fam), c("PDR11/AUS1-like", "PDR12-like", "PDR5-like", "SNQ2-like"))]
    fam$species_code <- rownames(fam)
    out <- merge(out, fam, by = "species_code", all.x = TRUE)
    for (k in c("n_pdr11aus1like", "n_pdr12like", "n_pdr5like", "n_snq2like")) if (k %in% names(out)) out[[k]][is.na(out[[k]])] <- 0L
    out$n_families <- rowSums(out[, c("n_pdr11aus1like", "n_pdr12like", "n_pdr5like", "n_snq2like")] > 0)
  }
  if (!is.null(frag_pairs_path)) {
    fp <- data.table::fread(frag_pairs_path, sep = "\t", data.table = FALSE)
    fp <- fp[fp$has_PASS == "N", ]
    fo <- as.data.frame(table(species_code = fp$species_code), stringsAsFactors = FALSE)
    names(fo)[2] <- "n_frag_only_families"
    out <- merge(out, fo, by = "species_code", all.x = TRUE)
    out$n_frag_only_families[is.na(out$n_frag_only_families)] <- 0L
  }
  out
}

a14_read_index <- function(path) {
  ix <- data.table::fread(path, sep = "\t", data.table = FALSE)
  comp <- ix[["cds_len_equals_3L+3(complete,stop)"]]
  data.frame(species_code = ix$speciesCode, n_proteins = ix$n_proteins,
             frac_partial = 1 - comp / ix$n_cds, stringsAsFactors = FALSE)
}

a14_build_frame <- function(s2, names_tab, counts, index, tree, drug = "Fluconazole_resistance") {
  m <- merge(s2, names_tab[, c("assembly_fullID_updated", "FINAL_Tip_ID")], by.x = "assembly", by.y = "assembly_fullID_updated", all.x = TRUE)
  if (any(is.na(m$FINAL_Tip_ID))) stop("S2 assemblies without a name-table row: ", sum(is.na(m$FINAL_Tip_ID)))
  if (any(duplicated(m$FINAL_Tip_ID))) stop("duplicated tree tips among S2 rows")
  m <- merge(m, counts, by.x = "FINAL_Tip_ID", by.y = "species_code", all.x = TRUE)
  m <- merge(m, index, by.x = "FINAL_Tip_ID", by.y = "species_code", all.x = TRUE)
  if (any(is.na(m$n_proteins))) stop("S2 rows without a proteome-index row: ", sum(is.na(m$n_proteins)))
  for (k in intersect(c("n_pass", "n_pdr5like", "n_snq2like", "n_pdr12like", "n_pdr11aus1like", "n_families", "n_frag_only_families"), names(m))) m[[k]][is.na(m[[k]])] <- 0L
  need <- c("n_pass", "n_proteins", "frac_partial")
  if (any(!complete.cases(m[, need]))) stop("rows without pump count or genome covariates: ", sum(!complete.cases(m[, need])))
  if (any(!(m$FINAL_Tip_ID %in% tree$tip.label))) stop("S2 species not on the tree: ", sum(!(m$FINAL_Tip_ID %in% tree$tip.label)))
  tr <- drop.tip(tree, setdiff(tree$tip.label, m$FINAL_Tip_ID))
  rownames(m) <- m$FINAL_Tip_ID
  m <- m[tr$tip.label, ]
  list(df = m, tree = tr)
}

a14_build_frame_chx <- function(s3, counts, index, tree, chx_cols) {
  d <- s3[, c("FINAL_Tip_ID", chx_cols)]
  d <- d[rowSums(!is.na(d[, chx_cols, drop = FALSE]) & d[, chx_cols, drop = FALSE] != "") > 0, ]
  for (cc in chx_cols) d[[cc]] <- suppressWarnings(as.integer(d[[cc]]))
  m <- merge(d, counts, by.x = "FINAL_Tip_ID", by.y = "species_code", all.x = TRUE)
  m <- merge(m, index, by.x = "FINAL_Tip_ID", by.y = "species_code", all.x = TRUE)
  m$n_pass[is.na(m$n_pass)] <- 0L
  m <- m[complete.cases(m[, c("n_pass", "n_proteins", "frac_partial")]) & m$FINAL_Tip_ID %in% tree$tip.label, ]
  tr <- drop.tip(tree, setdiff(tree$tip.label, m$FINAL_Tip_ID))
  rownames(m) <- m$FINAL_Tip_ID
  m <- m[tr$tip.label, ]
  m$species <- m$FINAL_Tip_ID
  list(df = m, tree = tr)
}

a14_add_predictors <- function(df) {
  df$x_raw <- df$n_pass
  df$x_z <- a14_z(df$n_pass)
  df$x_sqrt_z <- a14_z(sqrt(df$n_pass))
  df$x_log_z <- a14_z(log2(df$n_pass + 1))
  df$prot_z <- a14_z(log2(df$n_proteins))
  df$frag_z <- a14_z(df$frac_partial)
  if ("n_frag_only_families" %in% names(df)) df$x_fragfam_z <- a14_z(df$n_pass + df$n_frag_only_families)
  if ("n_families" %in% names(df)) df$x_rich_z <- a14_z(df$n_families)
  for (k in c("n_pdr5like", "n_snq2like", "n_pdr12like", "n_pdr11aus1like")) {
    if (k %in% names(df) && sd(df[[k]]) > 0) df[[paste0(sub("^n_", "", k), "_z")]] <- a14_z(df[[k]])
  }
  df
}

a14_tree_eigvecs <- function(tree, k = 2) {   # David 2025 template: leading eigenvectors of vcv(tree)
  ev <- eigen(vcv(tree), symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  colnames(ev) <- paste0("PC", seq_len(k)); rownames(ev) <- tree$tip.label
  ev
}

a14_vif <- function(df, vars) {
  if (length(vars) < 2) return(setNames(1, vars))
  sapply(vars, function(v) { r2 <- summary(lm(reformulate(setdiff(vars, v), v), data = df))$r.squared; 1 / (1 - r2) })
}

a14_pglm <- function(df, tree, response, predictors, btol = 35, log.alpha.bound = 4, method = "logistic_MPLE", boot = 0) {
  f <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  warns <- character(0)
  fit <- withCallingHandlers(
    tryCatch(phyloglm(f, data = df, phy = tree, method = method, btol = btol, log.alpha.bound = log.alpha.bound, boot = boot),
             error = function(e) e),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (inherits(fit, "error")) return(list(ok = FALSE, error = conditionMessage(fit), warnings = warns, n = nrow(df)))
  co <- summary(fit)$coefficients
  tab <- data.frame(term = rownames(co), est = co[, "Estimate"], se = co[, "StdErr"], z = co[, "z.value"], p = co[, "p.value"],
                    stringsAsFactors = FALSE, row.names = NULL)
  list(ok = TRUE, fit = fit, tab = tab, alpha = unname(fit$alpha), convergence = fit$convergence, alphaWarn = fit$alphaWarn,
       warnings = warns, n = fit$n, events = sum(df[[response]] == 1), tree_height = fit$mean.tip.height,
       bound_hit = any(grepl("bound", warns, ignore.case = TRUE)) || isTRUE(fit$alphaWarn != 0 && !is.null(fit$alphaWarn)))
}

A14_FIT_CHAIN <- list(list(btol = 35, method = "logistic_MPLE"), list(btol = 100, method = "logistic_MPLE"), list(btol = 100, method = "logistic_IG10"))
a14_pglm_chain <- function(df, tree, response, predictors, chain = A14_FIT_CHAIN, log.alpha.bound = 4) {
  last <- NULL
  for (k in seq_along(chain)) {
    r <- a14_pglm(df, tree, response, predictors, btol = chain[[k]]$btol, log.alpha.bound = log.alpha.bound, method = chain[[k]]$method)
    r$attempt <- k; r$chain_method <- chain[[k]]$method; r$chain_btol <- chain[[k]]$btol
    r$clean <- r$ok && !any(grepl("converge|boundary", r$warnings, ignore.case = TRUE)) && isTRUE(r$convergence == 0)
    if (r$clean) return(r)
    last <- r
  }
  last
}

a14_or <- function(res, term, sd_raw = NA_real_) {
  r <- res$tab[res$tab$term == term, ]
  out <- data.frame(term = term, beta = r$est, se = r$se, z = r$z, p = r$p, OR = exp(r$est),
                    OR_lo = exp(r$est - qnorm(0.975) * r$se), OR_hi = exp(r$est + qnorm(0.975) * r$se))
  if (is.finite(sd_raw)) { out$OR_per_unit <- exp(r$est / sd_raw); out$OR_per_unit_lo <- exp((r$est - qnorm(0.975) * r$se) / sd_raw); out$OR_per_unit_hi <- exp((r$est + qnorm(0.975) * r$se) / sd_raw) }
  out
}

A14_E1_PRIMARY <- c("x_z", "prot_z", "frag_z")

a14_run_e1 <- function(df, tree, response = "fluc", loo = TRUE, boot = 0, order_col = "Order") {
  sdx <- sd(df$x_raw)
  prim <- a14_pglm_chain(df, tree, response, A14_E1_PRIMARY)
  if (!prim$ok) stop("primary model failed: ", prim$error)
  out <- list(primary = prim, primary_clean = isTRUE(prim$clean), primary_attempt = prim$attempt, primary_or = a14_or(prim, "x_z", sdx), sd_x = sdx, n = prim$n, events = prim$events,
              vif = a14_vif(df, A14_E1_PRIMARY))
  prev <- mean(df[[response]] == 1)
  out$mde <- data.frame(neff_frac = c(1, 0.3, 0.1), MDE_OR = sapply(c(1, 0.3, 0.1), function(f) a14_mde_or(nrow(df), prev, f, vif = max(out$vif["x_z"]))))
  out$mde_realised_OR <- exp((qnorm(0.975) + qnorm(0.8)) * prim$tab$se[prim$tab$term == "x_z"])
  sens <- list()
  add <- function(name, res, term = "x_z") sens[[name]] <<- list(res = res, term = term)
  add("S_unadjusted", a14_pglm_chain(df, tree, response, "x_z"))
  add("S_noprot", a14_pglm_chain(df, tree, response, c("x_z", "frag_z")))
  add("S_nofrag", a14_pglm_chain(df, tree, response, c("x_z", "prot_z")))
  add("S_sqrt", a14_pglm_chain(df, tree, response, c("x_sqrt_z", "prot_z", "frag_z")), "x_sqrt_z")
  add("S_log2p1", a14_pglm_chain(df, tree, response, c("x_log_z", "prot_z", "frag_z")), "x_log_z")
  if ("x_fragfam_z" %in% names(df)) add("S_fragfam", a14_pglm_chain(df, tree, response, c("x_fragfam_z", "prot_z", "frag_z")), "x_fragfam_z")
  if ("x_rich_z" %in% names(df)) add("S_richness", a14_pglm_chain(df, tree, response, c("x_rich_z", "prot_z", "frag_z")), "x_rich_z")
  d0 <- df[df$n_pass > 0, ]       # genomes without any PASS sequence (likely failed retrieval for a universal family) removed
  if (nrow(d0) < nrow(df)) {
    t0 <- drop.tip(tree, setdiff(tree$tip.label, rownames(d0)))
    d0$x_z <- a14_z(d0$n_pass); d0$prot_z <- a14_z(log2(d0$n_proteins)); d0$frag_z <- a14_z(d0$frac_partial)
    add("S_nozero", a14_pglm_chain(d0, t0, response, A14_E1_PRIMARY))
  }
  add("S_btol10", a14_pglm(df, tree, response, A14_E1_PRIMARY, btol = 10))
  add("S_alphabound8", a14_pglm_chain(df, tree, response, A14_E1_PRIMARY, log.alpha.bound = 8))
  ev <- a14_tree_eigvecs(tree, 2)
  d2 <- cbind(df, ev[rownames(df), , drop = FALSE])
  add("S_eigvec2", a14_pglm_chain(d2, tree, response, c(A14_E1_PRIMARY, "PC1", "PC2")))
  out$sens <- sens
  g <- glm(as.formula(paste(response, "~", paste(A14_E1_PRIMARY, collapse = " + "))), data = df, family = binomial())
  gs <- summary(g)$coefficients
  out$glm_reference <- data.frame(beta = gs["x_z", 1], se = gs["x_z", 2], p = gs["x_z", 4], OR = exp(gs["x_z", 1]))
  if (loo && order_col %in% names(df)) {
    lo <- list()
    for (o in sort(unique(df[[order_col]]))) {
      d <- df[df[[order_col]] != o, ]
      if (sum(d[[response]] == 1) < 25 || sum(d[[response]] == 0) < 25) { lo[[o]] <- list(status = "NOT_EVALUABLE"); next }
      t2 <- drop.tip(tree, setdiff(tree$tip.label, rownames(d)))
      d$x_z <- a14_z(d$n_pass); d$prot_z <- a14_z(log2(d$n_proteins)); d$frag_z <- a14_z(d$frac_partial)
      r <- a14_pglm_chain(d, t2, response, A14_E1_PRIMARY)
      lo[[o]] <- if (r$ok) list(status = "OK", beta = r$tab$est[r$tab$term == "x_z"], p = r$tab$p[r$tab$term == "x_z"], n = r$n) else list(status = "FAILED")
    }
    out$loo <- lo
  }
  if (boot > 0) {
    set.seed(A14_SEED)
    out$boot <- a14_pglm(df, tree, response, A14_E1_PRIMARY, boot = boot)
  }
  out$label <- a14_label_e1(out)
  out
}

a14_label_e1 <- function(out) {
  if (!is.null(out$primary_clean) && !out$primary_clean)
    return(list(label = "not evaluable (the primary model did not converge cleanly in the whole fallback chain)", failed = "0_convergence"))
  b <- out$primary_or$beta; p <- out$primary_or$p
  if (!(p < 0.05 && b > 0)) {
    return(list(label = if (p < 0.05 && b < 0) "not supported (significant in the opposite direction)" else "not supported", failed = "1_primary"))
  }
  fails <- character(0)
  for (nm in names(out$sens)) {
    s <- out$sens[[nm]]
    if (!s$res$ok) { fails <- c(fails, paste0("2_", nm, "_failed")); next }
    be <- s$res$tab$est[s$res$tab$term == s$term]
    if (!(be > 0)) fails <- c(fails, paste0("2_", nm))
  }
  if (!is.null(out$loo)) {
    ok <- Filter(function(z) identical(z$status, "OK"), out$loo)
    n_pos <- sum(sapply(ok, function(z) z$beta > 0))
    if (length(ok) < length(out$loo) - 1 || n_pos < length(ok)) fails <- c(fails, "3_leave_one_order_out")
  }
  list(label = if (length(fails) == 0) "robust" else "supported, not robust", failed = fails)
}

a14_run_e2_drugs <- function(df, tree, drugs, min_events = 25) {
  rows <- lapply(drugs, function(d) {
    y <- df[[d]]
    if (sum(y == 1, na.rm = TRUE) < min_events) return(data.frame(drug = d, status = "NOT_TESTABLE", events = sum(y == 1)))
    dd <- df; dd$y <- y
    r <- a14_pglm(dd, tree, "y", A14_E1_PRIMARY)
    if (!r$ok) return(data.frame(drug = d, status = "FAILED", events = sum(y == 1)))
    o <- a14_or(r, "x_z", sd(df$x_raw))
    data.frame(drug = d, status = "OK", events = sum(y == 1), n = r$n, o[, c("beta", "se", "p", "OR", "OR_lo", "OR_hi")])
  })
  res <- do.call(rbind, lapply(rows, function(r) { for (k in c("beta", "se", "p", "OR", "OR_lo", "OR_hi", "n")) if (!k %in% names(r)) r[[k]] <- NA; r }))
  res$q <- NA_real_
  ok <- res$status == "OK"
  res$q[ok] <- a14_bh(res$p[ok])
  res
}

a14_run_e2_family <- function(df, tree, response = "fluc") {
  fam <- intersect(c("pdr5like_z", "snq2like_z", "pdr12like_z", "pdr11aus1like_z"), names(df))
  r <- a14_pglm(df, tree, response, c(fam, "prot_z", "frag_z"))
  if (!r$ok) return(list(ok = FALSE, error = r$error))
  tab <- r$tab[r$tab$term %in% fam, ]
  tab$q <- a14_bh(tab$p)
  tab$OR <- exp(tab$est)
  list(ok = TRUE, tab = tab, vif = a14_vif(df, c(fam, "prot_z", "frag_z")), alpha = r$alpha, n = r$n)
}

a14_code_targetsite <- function(ts) {
  split_aa <- function(x) strsplit(as.character(x), ";", fixed = TRUE)
  e <- split_aa(ts$eL42_56); u <- split_aa(ts$uL15_38)
  ts$eL42_mut_any <- sapply(e, function(v) if (all(v %in% c("", "NA", NA))) NA_integer_ else as.integer(any(v[!v %in% c("", "NA")] != "P")))
  ts$eL42_mut_all <- sapply(e, function(v) if (all(v %in% c("", "NA", NA))) NA_integer_ else as.integer(all(v[!v %in% c("", "NA")] != "P")))
  ts$uL15_mut_any <- sapply(u, function(v) if (all(v %in% c("", "NA", NA))) NA_integer_ else as.integer(any(v[!v %in% c("", "NA")] != "Q")))
  ts$ts_mut_any <- ifelse(is.na(ts$eL42_mut_any) | is.na(ts$uL15_mut_any), NA_integer_, as.integer(ts$eL42_mut_any == 1 | ts$uL15_mut_any == 1))
  ts$eL42_mixed <- ifelse(is.na(ts$eL42_mut_any), NA_integer_, as.integer(ts$eL42_mut_any != ts$eL42_mut_all))   # copies disagree
  ts
}

a14_run_e2_chx <- function(df, tree, ts, chx_cols = c("O1_Cycloheximide_0.01_", "O2_Cycloheximide_0.1_"), min_events = 25) {
  ts <- a14_code_targetsite(ts)
  d <- merge(df, ts[, c("species", "eL42_mut_any", "uL15_mut_any", "ts_mut_any", "eL42_mixed")], by.x = "FINAL_Tip_ID", by.y = "species", all.x = TRUE)
  rownames(d) <- d$FINAL_Tip_ID
  rows <- list()
  for (cc in chx_cols) {
    for (variant in c("adjusted", "restricted", "adjusted_nomixed")) {
      dd <- d[!is.na(d[[cc]]) & !is.na(d$ts_mut_any), ]
      if (variant == "restricted") dd <- dd[dd$ts_mut_any == 0, ]
      if (variant == "adjusted_nomixed") dd <- dd[!is.na(dd$eL42_mixed) & dd$eL42_mixed == 0, ]
      dd$y <- as.integer(dd[[cc]])
      if (nrow(dd) < 50 || sum(dd$y == 1) < min_events || sum(dd$y == 0) < min_events) { rows[[length(rows) + 1]] <- data.frame(cycloheximide = cc, variant = variant, status = "NOT_TESTABLE", n = nrow(dd), events = sum(dd$y == 1)); next }
      t2 <- drop.tip(tree, setdiff(tree$tip.label, rownames(dd)))
      dd <- dd[t2$tip.label, ]
      dd$x_z <- a14_z(dd$n_pass); dd$prot_z <- a14_z(log2(dd$n_proteins)); dd$frag_z <- a14_z(dd$frac_partial)
      preds <- c("x_z", "prot_z", "frag_z", if (variant != "restricted") "ts_mut_any")
      r <- a14_pglm(dd, t2, "y", preds)
      if (!r$ok) { rows[[length(rows) + 1]] <- data.frame(cycloheximide = cc, variant = variant, status = "FAILED", n = nrow(dd), events = sum(dd$y == 1)); next }
      o <- a14_or(r, "x_z", sd(dd$n_pass))
      rows[[length(rows) + 1]] <- data.frame(cycloheximide = cc, variant = variant, status = "OK", n = r$n, events = sum(dd$y == 1), o[, c("beta", "se", "p", "OR", "OR_lo", "OR_hi")])
    }
  }
  res <- do.call(rbind, lapply(rows, function(r) { for (k in c("beta", "se", "p", "OR", "OR_lo", "OR_hi")) if (!k %in% names(r)) r[[k]] <- NA; r }))
  res$q <- NA_real_; ok <- res$status == "OK" & res$variant %in% c("adjusted", "restricted")   # family FE2c; adjusted_nomixed is a sensitivity
  res$q[ok] <- a14_bh(res$p[ok])
  res
}
