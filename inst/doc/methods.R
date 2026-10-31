## ----setup, include = FALSE, echo = FALSE, message = FALSE--------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>")
library(mize)

rb_fg <- list(
  fn = function(x) {
    100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
  },
  gr = function(x) {
    c(
      -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
      200 * (x[2] - x[1] * x[1])
    )
  },
  fg = function(x) {
    a <- x[2] - x[1] * x[1]
    b <- 1 - x[1]
    list(
      fn = 100 * a^2 + b^2,
      gr = c(-400 * x[1] * a - 2 * b, 200 * a)
    )
  }
)
rb0 <- c(-1.2, 1)

## ----method-surface-check, include = FALSE------------------------------------
documented_methods <- c(
  "SD", "BFGS", "SR1", "L-BFGS", "CG", "TN", "NEWTON", "PHESS", "NAG",
  "Momentum", "DBD"
)
stopifnot(all(vapply(
  documented_methods,
  function(method) is.list(make_mize(method = method)),
  logical(1)
)))

## ----c2-diagnostics-----------------------------------------------------------
c2_runs <- lapply(c(`0.1` = 0.1, `0.5` = 0.5), function(c2_value) {
  mize(
    rb0,
    rb_fg,
    method = "CG",
    cg_update = "HZ+",
    line_search = "More-Thuente",
    c2 = c2_value,
    step_next_init = "hz",
    max_iter = 4,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL,
    store_progress = TRUE
  )
})

c2_diagnostics <- do.call(rbind, lapply(seq_along(c2_runs), function(i) {
  progress <- c2_runs[[i]]$progress
  progress <- progress[!is.na(progress$ls_reason), , drop = FALSE]
  data.frame(
    `c2` = names(c2_runs)[i],
    Iteration = seq_len(nrow(progress)),
    `Proposed alpha` = progress$alpha_init,
    `Selected alpha` = progress$alpha,
    `Reason / outcome` = paste(progress$ls_reason, progress$ls_outcome, sep = " / "),
    `Local work (fn / gr)` = paste(progress$ls_nf, progress$ls_ng, sep = " / "),
    check.names = FALSE
  )
}))
rownames(c2_diagnostics) <- NULL
c2_display <- c2_diagnostics
c2_display[["Proposed alpha"]] <- formatC(
  c2_display[["Proposed alpha"]],
  format = "g",
  digits = 4
)
c2_display[["Selected alpha"]] <- formatC(
  c2_display[["Selected alpha"]],
  format = "g",
  digits = 4
)
knitr::kable(c2_display)

## ----c2-check, include = FALSE------------------------------------------------
stopifnot(
  all(vapply(c2_runs, function(result) {
    identical(result$status, "budget_exhausted") &&
      identical(result$terminate$what, "max_iter")
  }, logical(1))),
  all(is.finite(c2_diagnostics[["Proposed alpha"]])),
  all(c2_diagnostics[["Proposed alpha"]] > 0),
  all(is.finite(c2_diagnostics[["Selected alpha"]])),
  all(c2_diagnostics[["Selected alpha"]] > 0),
  !isTRUE(all.equal(
    c2_runs[[1]]$progress$alpha,
    c2_runs[[2]]$progress$alpha
  ))
)

## ----armijo-and-fixed-steps---------------------------------------------------
backtracking_res <- mize(
  rb0,
  rb_fg,
  method = "SD",
  line_search = "Backtracking",
  step0 = 1,
  step_down = 0.5,
  max_iter = 10
)

constant_res <- mize(
  rb0,
  rb_fg,
  method = "SD",
  line_search = "Constant",
  norm_direction = TRUE,
  step0 = 0.01,
  max_iter = 10
)

bold_res <- mize(
  rb0,
  rb_fg,
  method = "SD",
  line_search = "Bold Driver",
  step_down = 0.5,
  step_up = 1.1,
  max_iter = 10
)

## ----bounded-comparison-helpers, include = FALSE------------------------------
bounded_controls <- list(
  max_iter = 100,
  abs_tol = 0,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL,
  store_progress = TRUE
)

run_bounded <- function(method, ...) {
  do.call(
    mize,
    c(list(par = rb0, fg = rb_fg, method = method), list(...), bounded_controls)
  )
}

summarize_runs <- function(results) {
  rows <- lapply(seq_along(results), function(i) {
    result <- results[[i]]
    data.frame(
      Run = names(results)[[i]],
      `Best objective` = result$best_f,
      `Last objective` = result$last_f,
      `Function callbacks` = result$nf,
      `Gradient callbacks` = result$ng,
      Status = result$status,
      Termination = result$terminate$what,
      check.names = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary
}

assert_bounded_runs <- function(results) {
  summary <- summarize_runs(results)
  stopifnot(
    all(summary$Status == "budget_exhausted"),
    all(summary$Termination == "max_iter"),
    all(is.finite(summary[["Best objective"]])),
    all(is.finite(summary[["Last objective"]])),
    all(summary[["Best objective"]] <= summary[["Last objective"]])
  )
}

plot_current_objective <- function(results) {
  traces <- lapply(results, function(result) {
    keep <- is.finite(result$progress$nf) &
      is.finite(result$progress$f) & result$progress$f > 0
    result$progress[keep, c("nf", "f"), drop = FALSE]
  })
  stopifnot(all(vapply(traces, nrow, integer(1)) > 0L))

  colors <- rep(
    c("#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7"),
    length.out = length(traces)
  )
  line_types <- rep(c(1, 2, 3, 4, 5), length.out = length(traces))
  plot(
    NA_real_,
    NA_real_,
    xlim = range(unlist(lapply(traces, `[[`, "nf"))),
    ylim = range(unlist(lapply(traces, `[[`, "f"))),
    log = "y",
    xlab = "Cumulative function callbacks",
    ylab = "Current objective (log scale)"
  )
  grid()
  for (i in seq_along(traces)) {
    lines(traces[[i]]$nf, traces[[i]]$f, col = colors[[i]], lty = line_types[[i]])
  }
  legend(
    "topright",
    legend = names(results),
    col = colors,
    lty = line_types,
    bty = "n",
    cex = 0.75,
    ncol = if (length(results) > 3L) 2 else 1
  )
}

## ----nag-q-comparison, fig.width = 7, fig.height = 4.5, fig.cap = "Current objective against cumulative function callbacks for two NAG schedules. Both runs use the same start and fixed 100-iteration budget.", fig.alt = "Log-scale line plot comparing current objective over function callbacks for NAG with q equal to zero and q equal to 0.001."----
nag_runs <- list(
  `q = 0` = run_bounded("NAG", nest_q = 0),
  `q = 0.001` = run_bounded("NAG", nest_q = 0.001)
)
assert_bounded_runs(nag_runs)
knitr::kable(summarize_runs(nag_runs), digits = 6)
plot_current_objective(nag_runs)

## ----deterministic-momentum-schedule------------------------------------------
cosine_schedule <- function(iter, max_iter) {
  phase <- min(iter / max_iter, 1)
  0.4 + 0.2 * (1 - cos(pi * phase))
}

schedule_values <- data.frame(
  Iteration = c(1, 5, 10),
  Momentum = vapply(
    c(1, 5, 10),
    cosine_schedule,
    numeric(1),
    max_iter = 10
  )
)
knitr::kable(schedule_values, digits = 3)

## ----momentum-type-comparison, fig.width = 7, fig.height = 4.5, fig.cap = "Current objective for classical and Nesterov-style momentum under one fixed schedule and 100-iteration budget.", fig.alt = "Log-scale line plot comparing current objective over function callbacks for classical and Nesterov-style momentum."----
momentum_runs <- list(
  Classical = run_bounded("Momentum", mom_schedule = 0.9),
  Nesterov = run_bounded(
    "Momentum",
    mom_schedule = 0.9,
    mom_type = "nesterov"
  )
)
assert_bounded_runs(momentum_runs)
knitr::kable(summarize_runs(momentum_runs), digits = 6)
plot_current_objective(momentum_runs)

## ----additive-dbd-------------------------------------------------------------
dbd_res <- mize(
  rb0,
  rb_fg,
  method = "DBD",
  step0 = "rasmussen",
  step_down = 0.8,
  step_up = 0.2,
  step_up_fun = "+",
  dbd_weight = 0.5,
  max_iter = 10
)

## ----adaptive-restart-comparison, fig.width = 7, fig.height = 4.5, fig.cap = "Current objective for classical momentum with no restart and four adaptive-restart configurations under the same 100-iteration budget.", fig.alt = "Log-scale line plot comparing current objective over function callbacks for momentum without restart and with function- or gradient-based restart at waiting periods of ten or one iterations."----
restart_runs <- list(
  None = momentum_runs[["Classical"]],
  `fn, wait = 10` = run_bounded(
    "Momentum",
    mom_schedule = 0.9,
    restart = "fn",
    restart_wait = 10
  ),
  `gr, wait = 10` = run_bounded(
    "Momentum",
    mom_schedule = 0.9,
    restart = "gr",
    restart_wait = 10
  ),
  `fn, wait = 1` = run_bounded(
    "Momentum",
    mom_schedule = 0.9,
    restart = "fn",
    restart_wait = 1
  ),
  `gr, wait = 1` = run_bounded(
    "Momentum",
    mom_schedule = 0.9,
    restart = "gr",
    restart_wait = 1
  )
)
assert_bounded_runs(restart_runs)
knitr::kable(summarize_runs(restart_runs), digits = 6)
plot_current_objective(restart_runs)

