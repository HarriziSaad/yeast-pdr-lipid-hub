
args <- commandArgs(trailingOnly = TRUE)
kv <- function(key, default = NULL) { h <- grep(paste0("^", key, "="), args, value = TRUE); if (length(h)) sub(paste0("^", key, "="), "", h[1]) else default }
root <- kv("root"); out <- kv("out"); toy <- identical(kv("toy", "0"), "1"); boot <- as.integer(kv("boot", "0"))
if (is.null(root) || is.null(out)) stop("root= and out= are required")
this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this)), "a14_e_lib.R"))
a14_set_libpaths(root); a14_load_packages()
dir.create(out, showWarnings = FALSE, recursive = TRUE)

p <- function(...) file.path(root, ...)
s2_path   <- kv("s2",         p("phase0_verify", "data_access", "downloads", "harrison_phenotypes", "harrison2026_S2_resistance_532species.xlsx"))
names_p   <- kv("names",      p("phase0_verify", "data_access", "downloads", "y1000", "ALL_NAMES_FILE.xlsx"))
tree_p    <- kv("tree",       p("phase0_verify", "data_access", "downloads", "y1000", "1154yeasts_1403OGs_ml_timetree.tree"))
retr_p    <- kv("retrieval",  p("phase1_orthology", "phase1a", "outputs", "retrieval_table.tsv"))
og_p      <- kv("og_members", p("phase1_orthology", "phase1a", "outputs", "fragment_audit", "og_member_paralog_assignment.tsv"))
frag_p    <- kv("frag_pairs", p("phase1_orthology", "phase1a", "outputs", "fragment_audit", "audit_pairs_genome_family.tsv"))
index_p   <- kv("index",      p("phase1_orthology", "data_staging", "y1000_merged_tables", "y1000_proteome_index.tsv"))
s3_p      <- kv("s3",         p("phase0_verify", "data_access", "downloads", "harrison_phenotypes", "harrison2024_S3_metabolic_traits.tsv"))
ts_p      <- kv("targetsite", NULL)

set.seed(A14_SEED)
tree <- read.tree(tree_p)
stopifnot(Ntip(tree) == 1154 || toy)
s2 <- a14_read_s2(s2_path)
nm <- a14_read_names(names_p)
cnt <- a14_count_pumps(retr_p, og_p, if (file.exists(frag_p)) frag_p else NULL)
idx <- a14_read_index(index_p)
fr <- a14_build_frame(s2, nm, cnt, idx, tree)
df <- a14_add_predictors(fr$df); tr <- fr$tree
df$fluc <- df$Fluconazole_resistance
if (!toy && nrow(df) != 532) stop("expected n = 532, got ", nrow(df))
cat("analysed species:", nrow(df), " tree tips:", Ntip(tr), " fluconazole resistant:", sum(df$fluc == 1), "\n")
cat("marginal pump count (n_pass): min", min(df$n_pass), "median", median(df$n_pass), "max", max(df$n_pass), " SD", round(sd(df$n_pass), 3), "\n")
write.table(df[, c("FINAL_Tip_ID", "species", "assembly", "Order", "n_pass", "n_proteins", "frac_partial")], file.path(out, "A14_E_analysis_frame.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

e1 <- a14_run_e1(df, tr, "fluc", loo = TRUE, boot = boot)
sens_tab <- do.call(rbind, lapply(names(e1$sens), function(n) {
  s <- e1$sens[[n]]
  if (!s$res$ok) return(data.frame(model = n, status = "FAILED"))
  r <- s$res$tab[s$res$tab$term == s$term, ]
  data.frame(model = n, status = "OK", beta = r$est, se = r$se, p = r$p, OR = exp(r$est), alpha = s$res$alpha, bound_hit = s$res$bound_hit)
}))
loo_tab <- if (!is.null(e1$loo)) do.call(rbind, lapply(names(e1$loo), function(n) { z <- e1$loo[[n]]; data.frame(order_left_out = n, status = z$status, beta = if (is.null(z$beta)) NA else z$beta, p = if (is.null(z$p)) NA else z$p) })) else NULL
write.table(e1$primary_or, file.path(out, "A14_E1_primary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(e1$primary$tab, file.path(out, "A14_E1_primary_full_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sens_tab, file.path(out, "A14_E1_sensitivities.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(loo_tab)) write.table(loo_tab, file.path(out, "A14_E1_leave_one_order_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(sprintf("label: %s", e1$label$label), sprintf("failed checks: %s", paste(e1$label$failed, collapse = ", ")),
             sprintf("n = %d, events = %d, SD of the pump count = %.4f", e1$n, e1$events, e1$sd_x),
             sprintf("alpha (OU, per tree-height units) = %.5g, convergence = %s, bound hit = %s", e1$primary$alpha, e1$primary$convergence, e1$primary$bound_hit),
             sprintf("realised minimum detectable OR per SD (80%% power, two-sided 0.05, from the realised SE) = %.3f", e1$mde_realised_OR),
             "a priori minimum detectable OR per SD at n_eff = 100 / 30 / 10 % of n:", paste(sprintf("  %.0f%%: %.3f", 100 * e1$mde$neff_frac, e1$mde$MDE_OR), collapse = "\n"),
             "VIF:", paste(sprintf("  %s: %.2f", names(e1$vif), e1$vif), collapse = "\n"),
             sprintf("non-phylogenetic logistic regression (reference only): beta = %.4f, p = %.4g", e1$glm_reference$beta, e1$glm_reference$p),
             "primary-model warnings:", paste0("  ", e1$primary$warnings)), file.path(out, "A14_E1_summary.txt"))
saveRDS(e1[c("primary_or", "label", "mde", "mde_realised_OR", "vif", "glm_reference")], file.path(out, "A14_E1_result.rds"))
cat("E1 label:", e1$label$label, " | OR per SD:", sprintf("%.3f (%.3f-%.3f)", e1$primary_or$OR, e1$primary_or$OR_lo, e1$primary_or$OR_hi), " p =", signif(e1$primary_or$p, 3), "\n")

drugs <- setdiff(grep("_resistance$", names(df), value = TRUE), "Fluconazole_resistance")
e2a <- a14_run_e2_drugs(df, tr, drugs)
write.table(e2a, file.path(out, "A14_E2a_other_drugs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
e2b <- a14_run_e2_family(df, tr, "fluc")
if (e2b$ok) write.table(e2b$tab, file.path(out, "A14_E2b_family_counts_fluconazole.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
e2b_other <- do.call(rbind, lapply(drugs, function(d) {
  if (sum(df[[d]] == 1, na.rm = TRUE) < 25) return(NULL)
  r <- a14_run_e2_family(df, tr, d)
  if (!r$ok) return(NULL)
  data.frame(drug = d, term = r$tab$term, beta = r$tab$est, se = r$tab$se, OR = exp(r$tab$est),
             OR_lo = exp(r$tab$est - qnorm(0.975) * r$tab$se), OR_hi = exp(r$tab$est + qnorm(0.975) * r$tab$se))
}))
if (!is.null(e2b_other)) write.table(e2b_other, file.path(out, "A14_E2b_family_counts_other_drugs_descriptive.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(ts_p) && file.exists(ts_p)) {
  s3 <- data.table::fread(s3_p, sep = "\t", data.table = FALSE)
  names(s3)[names(s3) == "species name"] <- "FINAL_Tip_ID"
  chx <- grep("Cycloheximide", names(s3), value = TRUE)
  allg <- a14_build_frame_chx(s3, cnt, idx, tree, chx)
  ts <- data.table::fread(ts_p, sep = "\t", data.table = FALSE)
  e2c <- a14_run_e2_chx(allg$df, allg$tree, ts, chx)
  write.table(e2c, file.path(out, "A14_E2c_cycloheximide.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out, "A14_E_sessionInfo.txt"))
cat("done\n")
