# Generates cert/data/*.csv inputs and the R reference values in sgmm_reference.dta
# (via reference.csv) using staggeredGMM 0.1.0 from CRAN.
# Usage: Rscript make_reference.R <output_dir>

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(staggeredGMM))

make_test_panel <- function(cohorts = c(3L, 4L), n_per = 2L, n_never = 2L,
                            T_total = 4L, effects = NULL, sd = 0, seed = 1L) {
  set.seed(seed)
  g_vec <- c(rep(as.integer(cohorts), each = n_per), rep(0L, n_never))
  N <- length(g_vec)
  unit <- rep(seq_len(N), each = T_total)
  time <- rep(seq_len(T_total), times = N)
  g <- g_vec[unit]
  alpha <- stats::rnorm(N)
  lambda <- stats::rnorm(T_total)
  tau <- numeric(length(unit))
  if (!is.null(effects)) {
    m <- match(paste(g, time, sep = "\r"), paste(effects$g, effects$t, sep = "\r"))
    tau[!is.na(m)] <- effects$value[m[!is.na(m)]]
  }
  eps <- if (sd > 0) stats::rnorm(length(unit), sd = sd) else 0
  data.frame(unit_id = unit, year = time, cohort = g,
             y = alpha[unit] + lambda[time] + tau + eps)
}
tiny_effects <- function() data.frame(g = c(3L, 3L, 4L), t = c(3L, 4L, 4L), value = c(1, 2, 3))
inject_mcar <- function(dat, frac = 0.10, seed = 7L) {
  set.seed(seed)
  n <- as.integer(nrow(dat) * frac)
  dat$y[sample(nrow(dat), n)] <- NA_real_
  dat
}

f17 <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
write_data <- function(dat, name) {
  d <- dat
  for (nm in names(d)) if (is.numeric(d[[nm]])) d[[nm]] <- f17(d[[nm]])
  utils::write.csv(d, file.path(out_dir, "data", paste0(name, ".csv")),
                   row.names = FALSE, quote = TRUE)
}

rows <- list()
add <- function(case, key, value, idx = seq_along(value)) {
  if (length(value) == 0L) return(invisible())
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case, key = key, idx = idx,
    value = if (is.character(value)) value else f17(as.numeric(value)),
    stringsAsFactors = FALSE)
}

run_case <- function(case, dat, data_name, ..., jtests = list()) {
  w <- character(0)
  fit <- withCallingHandlers(
    gmm_staggered(dat, ...),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  add(case, "data", data_name, 1L)
  add(case, "coef", unname(fit$coefficients))
  add(case, "se", fit$catt$std_error)
  add(case, "g", fit$catt$g); add(case, "t", fit$catt$t)
  add(case, "identified", as.integer(fit$catt$identified))
  add(case, "n_comparisons", fit$catt$n_comparisons)
  add(case, "vcov", as.numeric(fit$vcov))
  for (a in c("CW", "EW")) {
    A <- fit$aggregate[[a]]
    add(case, paste0(a, "_estimate"), A$estimate, 1L)
    add(case, paste0(a, "_std_error"), A$std_error, 1L)
    add(case, paste0(a, "_estimate_identified"), A$estimate_identified, 1L)
    add(case, paste0(a, "_std_error_identified"), A$std_error_identified, 1L)
    add(case, paste0(a, "_share"), A$identified_weight_share, 1L)
    add(case, paste0("w_", a), fit$weights[[a]])
  }
  cv <- fit$convergence
  add(case, "converged", as.integer(cv$converged), 1L)
  add(case, "solve_ok", as.integer(cv$solve_ok), 1L)
  add(case, "termination", cv$termination, 1L)
  add(case, "n_iter", cv$n_iter, 1L)
  add(case, "max_change", cv$max_change, 1L)
  d <- fit$design
  for (k in c("N_units", "T_periods", "n_cohorts", "N_beta", "N_identified",
              "N_2x2", "n_dropped_2x2", "moment_rank", "n_kept", "n_zero_lag"))
    add(case, k, d[[k]], 1L)
  add(case, "N_g", unname(d$N_g))
  add(case, "has_never", as.integer(fit$controls$has_never), 1L)
  add(case, "n_never_zero", fit$controls$n_never_zero, 1L)
  add(case, "n_never_late", fit$controls$n_never_late, 1L)
  add(case, "cov_n_fallback", fit$controls$cov_n_fallback, 1L)
  add(case, "balanced", as.integer(fit$validation$balanced), 1L)
  add(case, "n_missing_cells", fit$validation$n_missing_cells, 1L)
  add(case, "n_obs", fit$validation$n_obs, 1L)
  add(case, "n_warnings", length(w), 1L)
  if (length(w)) add(case, "warning", w)
  add(case, "W_trace", sum(diag(fit$.internal$W)), 1L)
  add(case, "W_sum", sum(fit$.internal$W), 1L)
  for (jt in jtests) {
    jw <- character(0)
    j <- withCallingHandlers(
      gmm_j_test(fit, type = jt$type, window = jt$window),
      warning = function(cnd) { jw <<- c(jw, conditionMessage(cnd)); invokeRestart("muffleWarning") })
    p <- paste0("J_", jt$type, jt$window, "_")
    add(case, paste0(p, "stat"), j$statistic, 1L)
    add(case, paste0(p, "df"), j$df, 1L)
    add(case, paste0(p, "p"), j$p_value, 1L)
    add(case, paste0(p, "base"), j$base_period, 1L)
    add(case, paste0(p, "n_moments"), j$n_moments, 1L)
    add(case, paste0(p, "n_dropped"), j$n_dropped, 1L)
    add(case, paste0(p, "moments"), j$moments$moment)
    add(case, paste0(p, "n_warnings"), length(jw), 1L)
  }
  invisible(fit)
}

sp <- sim_panel
write_data(sp, "sim_panel")
JT <- list(list(type = "local", window = 3L), list(type = "full", window = 3L),
           list(type = "local", window = 5L))

run_case("S1", sp, "sim_panel", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", weighting = "pooled_toeplitz", jtests = JT)
run_case("S2", sp, "sim_panel", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", weighting = "cohort_toeplitz", jtests = JT)
run_case("S3", sp, "sim_panel", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", weighting = "unrestricted", jtests = JT)
run_case("S4", sp, "sim_panel", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", covar = c("x1", "x2"), jtests = JT)

s5 <- inject_mcar(sp, 0.10, 7L)
s5$y[s5$cohort == 13L & s5$year == 15L] <- NA_real_
write_data(s5, "s5")
run_case("S5", s5, "s5", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", jtests = JT)
run_case("S5b", s5, "s5", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", weighting = "cohort_toeplitz", jtests = JT)
run_case("S5c", s5, "s5", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", covar = c("x1", "x2"))

bb <- beck_banks
write_data(bb, "beck_banks")
bb6 <- bb[bb$branch_reform > 1976, ]
write_data(bb6, "beck_banks_1977")
run_case("S6", bb6, "beck_banks_1977", yname = "ln_gini", tname = "wrkyr",
         idname = "state", gname = "branch_reform")
run_case("S6b", bb6, "beck_banks_1977", yname = "ln_gini", tname = "wrkyr",
         idname = "state", gname = "branch_reform", weighting = "cohort_toeplitz")

t1 <- make_test_panel(effects = tiny_effects(), sd = 0)
write_data(t1, "t1")
for (w in c("pooled_toeplitz", "cohort_toeplitz", "unrestricted"))
  run_case(paste0("T1_", w), t1, "t1", yname = "y", tname = "year",
           idname = "unit_id", gname = "cohort", weighting = w)

t2 <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 0L, T_total = 4L,
                      effects = tiny_effects(), sd = 0.05, seed = 51L)
write_data(t2, "t2")
run_case("T2", t2, "t2", yname = "y", tname = "year", idname = "unit_id", gname = "cohort")

t3 <- make_test_panel(cohorts = c(3L, 4L), n_per = 2L, n_never = 2L, T_total = 5L,
                      sd = 0.2, seed = 77L)
set.seed(78); nu <- length(unique(t3$unit_id))
t3$xa <- stats::rnorm(nu)[t3$unit_id]; t3$xb <- stats::rnorm(nu)[t3$unit_id]
write_data(t3, "t3")
run_case("T3", t3, "t3", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", covar = c("xa", "xb"))

t4 <- make_test_panel(cohorts = c(3L, 5L), n_per = 5L, n_never = 5L, T_total = 7L,
                      sd = 0.4, seed = 91L)
write_data(t4, "t4")
for (w in c("pooled_toeplitz", "cohort_toeplitz", "unrestricted"))
  run_case(paste0("T4_", w), t4, "t4", yname = "y", tname = "year",
           idname = "unit_id", gname = "cohort", weighting = w,
           jtests = list(list(type = "full", window = 3L)))

t5 <- make_test_panel(cohorts = c(3L, 4L, 5L), n_per = 3L, n_never = 3L,
                      T_total = 6L, sd = 0.05)
write_data(t5, "t5")
run_case("T5", t5, "t5", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", never_treated = FALSE)

t6 <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 3L, T_total = 5L,
                      effects = tiny_effects(), sd = 0.05, seed = 31L)
t6$cohort[t6$cohort == 0L] <- 99L
write_data(t6, "t6")
run_case("T6", t6, "t6", yname = "y", tname = "year", idname = "unit_id", gname = "cohort")

run_case("T7", sp, "sim_panel", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", max_iter = 1L)

t8 <- make_test_panel(cohorts = c(3L, 4L), n_per = 6L, n_never = 0L, T_total = 6L,
                      sd = 0.3, seed = 53L)
write_data(t8, "t8")
run_case("T8", t8, "t8", yname = "y", tname = "year", idname = "unit_id", gname = "cohort")

t9 <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 3L, T_total = 5L,
                      sd = 0.2, seed = 71L)
set.seed(721); xu <- stats::rnorm(length(unique(t9$unit_id)))
t9$x <- xu[t9$unit_id] + 0.1 * t9$year
t9 <- t9[order(t9$unit_id, -t9$year), ]
write_data(t9, "t9")
run_case("T9", t9, "t9", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", covar = "x")

t10 <- make_test_panel(cohorts = c(3L, 4L), n_per = 5L, n_never = 5L, T_total = 4L,
                       sd = 0.3, seed = 55L)
write_data(t10, "t10")
run_case("T10", t10, "t10", yname = "y", tname = "year", idname = "unit_id",
         gname = "cohort", weighting = "unrestricted", max_iter = 1L)

res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out_dir, "reference.csv"), row.names = FALSE)
cat("cases:", length(unique(res$case)), " rows:", nrow(res), "\n")
