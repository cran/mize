## ----setup, include = FALSE, echo = FALSE, message = FALSE--------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>")
library(mize)

## ----define-rosenbrock--------------------------------------------------------
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

## ----check-gradient-----------------------------------------------------------
rb0 <- c(-1.2, 1)
gradient_check <- check_mize_gradient(rb_fg, rb0)

gradient_summary <- data.frame(
  "Maximum absolute error" = formatC(
    gradient_check$max_abs_error,
    format = "e",
    digits = 2
  ),
  "Maximum relative error" = formatC(
    gradient_check$max_rel_error,
    format = "e",
    digits = 2
  ),
  check.names = FALSE
)
knitr::kable(gradient_summary, align = c("r", "r"))

## ----assert-gradient, include=FALSE-------------------------------------------
stopifnot(
  is.finite(gradient_check$max_rel_error),
  gradient_check$max_rel_error < 1e-6
)

## ----run-mize-----------------------------------------------------------------
res <- mize(rb0, rb_fg, store_progress = TRUE)

## ----result-summary-----------------------------------------------------------
res[c("par", "f", "status", "message")]

## ----callback-summary---------------------------------------------------------
res[c("nf", "ng", "nh", "nhi")]

## ----result-termination-------------------------------------------------------
res[c("converged", "terminate")]

## ----assert-result, include=FALSE---------------------------------------------
initial_f <- rb_fg$fn(rb0)
convergence_reasons <- c(
  "abs_tol", "rel_tol", "grad_tol", "ginf_tol", "step_tol"
)

stopifnot(
  is.finite(res$best_f),
  res$best_f < initial_f,
  identical(res$status, "converged"),
  isTRUE(res$converged),
  res$terminate$what %in% convergence_reasons
)

## ----progress-plot, echo=FALSE, fig.width=6.5, fig.height=4.2, fig.cap="Rosenbrock objective values over cumulative function evaluations. The vertical axis uses a logarithmic scale.", fig.alt="Line plot of the Rosenbrock objective decreasing from about 24 to near zero as cumulative function evaluations increase. The objective axis is logarithmic."----
progress_for_plot <- subset(res$progress, is.finite(f) & f > 0)

plot(
  progress_for_plot$nf,
  progress_for_plot$f,
  type = "o",
  pch = 16,
  cex = 0.6,
  log = "y",
  xlab = "Cumulative function evaluations",
  ylab = "Objective"
)

stopifnot(
  nrow(progress_for_plot) > 1,
  tail(progress_for_plot$f, 1) < progress_for_plot$f[1]
)

## ----compact-progress---------------------------------------------------------
progress_view <- res$progress[, c("f", "step", "nf", "ng"), drop = FALSE]

## ----progress-table, echo=FALSE-----------------------------------------------
progress_rows <- unique(c(
  seq_len(min(3L, nrow(progress_view))),
  seq.int(max(1L, nrow(progress_view) - 2L), nrow(progress_view))
))
progress_excerpt <- progress_view[progress_rows, , drop = FALSE]

progress_table <- data.frame(
  Iteration = as.integer(rownames(progress_excerpt)),
  Objective = trimws(formatC(progress_excerpt$f, format = "g", digits = 4)),
  Step = trimws(formatC(progress_excerpt$step, format = "g", digits = 4)),
  "Function calls" = progress_excerpt$nf,
  "Gradient calls" = progress_excerpt$ng,
  check.names = FALSE
)
knitr::kable(
  progress_table,
  align = rep("r", ncol(progress_table))
)

