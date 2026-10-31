## ----setup, include = FALSE, echo = FALSE, message = FALSE--------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>")
library(mize)

result_summary <- function(res) {
  objective <- if (is.null(res$f)) NA_real_ else res$f
  data.frame(
    status = res$status,
    terminate = res$terminate$what,
    iter = res$iter,
    nf = res$nf,
    ng = res$ng,
    objective = signif(objective, 6),
    row.names = NULL
  )
}

## ----example-problem----------------------------------------------------------
rb_fg <- list(
  fn = function(x) {
    100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
  },
  gr = function(x) {
    c(
      -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
      200 * (x[2] - x[1] * x[1])
    )
  }
)
rb0 <- c(-1.2, 1)

## ----termination-by-max-iter--------------------------------------------------
res <- mize(
  rb0,
  rb_fg,
  max_iter = 10,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)

## ----assert-max-iter, include = FALSE-----------------------------------------
stopifnot(
  identical(res$terminate$what, "max_iter"),
  identical(res$status, "budget_exhausted"),
  identical(res$converged, FALSE)
)

## ----absolute-tolerance-------------------------------------------------------
res <- mize(
  rb0,
  rb_fg,
  abs_tol = 1e-8,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)

## ----assert-absolute-tolerance, include = FALSE-------------------------------
stopifnot(
  identical(res$terminate$what, "abs_tol"),
  identical(res$status, "converged"),
  identical(res$converged, TRUE)
)

## ----relative-tolerance-------------------------------------------------------
slow_minimum <- 10
slow_fg <- list(
  fn = function(x) (x - slow_minimum)^2,
  gr = function(x) 2 * (x - slow_minimum)
)

relative_result <- mize(
  0,
  slow_fg,
  method = "SD",
  line_search = "constant",
  step0 = 1e-4,
  max_iter = 5,
  abs_tol = NULL,
  rel_tol = 1e-3,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)

## ----relative-tolerance-table, echo=FALSE-------------------------------------
relative_summary <- data.frame(
  Reason = relative_result$terminate$what,
  `Relative change` = relative_result$terminate$val,
  Parameter = relative_result$par,
  `Distance from 10` = abs(relative_result$par - slow_minimum),
  check.names = FALSE
)
knitr::kable(relative_summary, digits = 4, align = c("l", "r", "r", "r"))

## ----assert-relative-tolerance, include = FALSE-------------------------------
stopifnot(
  identical(relative_result$terminate$what, "rel_tol"),
  identical(relative_result$status, "converged"),
  identical(relative_result$converged, TRUE),
  relative_result$terminate$val < 1e-3,
  abs(relative_result$par - slow_minimum) > 9
)

## ----gradient-2-norm-tolerance------------------------------------------------
res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = 1e-3,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)

## ----assert-gradient-2-norm-tolerance, include = FALSE------------------------
stopifnot(
  identical(res$terminate$what, "grad_tol"),
  identical(res$status, "converged"),
  identical(res$converged, TRUE)
)

## ----gradient-infinity-norm-tolerance-----------------------------------------
res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = 1e-3,
  step_tol = NULL
)
result_summary(res)

## ----assert-gradient-infinity-norm-tolerance, include = FALSE-----------------
stopifnot(
  identical(res$terminate$what, "ginf_tol"),
  identical(res$status, "converged"),
  identical(res$converged, TRUE)
)

## ----step-tolerance-----------------------------------------------------------
res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = 1e-4
)
result_summary(res)

## ----assert-step-tolerance, include = FALSE-----------------------------------
stopifnot(
  identical(res$terminate$what, "step_tol"),
  identical(res$status, "converged"),
  identical(res$converged, TRUE)
)

## ----evaluation-budget-comparison---------------------------------------------
budget_names <- c("max_fn", "max_gr", "max_fg")

run_with_budget <- function(budget) {
  arguments <- list(
    par = rb0,
    fg = rb_fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  arguments[[budget]] <- 10
  do.call(mize, arguments)
}

budget_runs <- setNames(
  lapply(budget_names, run_with_budget),
  budget_names
)

## ----evaluation-budget-table, echo=FALSE--------------------------------------
budget_summary <- do.call(
  rbind,
  lapply(names(budget_runs), function(budget) {
    run <- budget_runs[[budget]]
    data.frame(
      Control = budget,
      Reason = run$terminate$what,
      `Function calls` = run$nf,
      `Gradient calls` = run$ng,
      `Total calls` = run$nf + run$ng,
      check.names = FALSE
    )
  })
)
rownames(budget_summary) <- NULL
knitr::kable(
  budget_summary,
  align = c("l", "l", "r", "r", "r")
)

## ----assert-evaluation-budgets, include=FALSE---------------------------------
budget_reasons <- vapply(
  budget_runs,
  function(run) run$terminate$what,
  character(1)
)
budget_statuses <- vapply(budget_runs, function(run) run$status, character(1))
budget_convergence <- vapply(
  budget_runs,
  function(run) run$converged,
  logical(1)
)

stopifnot(
  identical(unname(budget_reasons), budget_names),
  all(budget_statuses == "budget_exhausted"),
  !any(budget_convergence),
  budget_runs$max_fn$nf == 10,
  budget_runs$max_gr$ng == 10,
  budget_runs$max_fg$nf + budget_runs$max_fg$ng == 10
)

## ----non-finite-gradient------------------------------------------------------
bad_fg <- list(
  fn = function(x) sum(x^2),
  gr = function(x) rep(NaN, length(x))
)

failure_result <- mize(c(1, 1), bad_fg)
failure_result[c("status", "message", "terminate")]

## ----assert-non-finite-gradient, include=FALSE--------------------------------
stopifnot(
  identical(failure_result$status, "failed"),
  identical(failure_result$converged, FALSE),
  identical(failure_result$terminate$what, "gr_inf")
)

## ----checking-convergence-less-often------------------------------------------
res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = 1e-3,
  ginf_tol = NULL,
  step_tol = NULL,
  check_conv_every = 5
)
result_summary(res)

## ----assert-checking-convergence-less-often, include = FALSE------------------
stopifnot(
  identical(res$terminate$what, "grad_tol"),
  identical(res$status, "converged"),
  res$iter %% 5 == 0
)

## ----best-and-last-results----------------------------------------------------
separating_fg <- list(
  fn = function(x) -3 * x^3 + 5 * x^2 - x,
  gr = function(x) -9 * x^2 + 10 * x - 1
)

separated_result <- mize(
  0,
  separating_fg,
  method = "SD",
  line_search = "constant",
  step0 = 1,
  max_iter = 1,
  abs_tol = 0,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = 0,
  step_tol = NULL,
  store_progress = TRUE
)

point_comparison <- data.frame(
  Point = c("Returned best", "Terminal"),
  Parameter = c(separated_result$best_par, separated_result$last_par),
  Objective = c(separated_result$best_f, separated_result$last_f),
  `Gradient infinity norm` = c(
    separated_result$ginfn,
    separated_result$terminate$val
  ),
  check.names = FALSE
)
knitr::kable(point_comparison, align = c("l", "r", "r", "r"))

## ----assert-best-and-last-results, include = FALSE----------------------------
stopifnot(
  identical(separated_result$terminate$what, "ginf_tol"),
  identical(separated_result$status, "converged"),
  identical(separated_result$best_par, separated_result$par),
  separated_result$best_f < separated_result$last_f,
  identical(separated_result$best_par, 0),
  identical(separated_result$last_par, 1),
  identical(separated_result$terminate$val, 0),
  identical(separated_result$ginfn, 1)
)

