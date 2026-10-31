## ----setup, include = FALSE, echo = FALSE, message = FALSE--------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>")
library(mize)

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

## ----reinitializing-optim-----------------------------------------------------
first_batch <- stats::optim(
  par = rb0,
  fn = rb_fg$fn,
  gr = rb_fg$gr,
  method = "BFGS",
  control = list(maxit = 10)
)
second_batch <- stats::optim(
  par = first_batch$par,
  fn = rb_fg$fn,
  gr = rb_fg$gr,
  method = "BFGS",
  control = list(maxit = 10)
)
c(first = first_batch$value, second = second_batch$value)

## ----create-and-initialize----------------------------------------------------
opt <- make_mize(method = "BFGS")
opt <- mize_init(
  opt,
  par = rb0,
  fg = rb_fg,
  max_iter = 3,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)

## ----optional-observation-witness---------------------------------------------
observation_opt <- make_mize(
  method = "SD",
  line_search = "constant",
  step0 = 0.0001,
  par = rb0,
  fg = rb_fg,
  max_iter = 3,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
observation_step <- mize_step(observation_opt, rb0, rb_fg)
names(observation_step)

## ----assert-optional-observation-witness, include = FALSE---------------------
stopifnot(
  is.null(observation_step$f),
  is.null(observation_step$g),
  identical(observation_step$nh, 0),
  identical(observation_step$nhi, 0)
)

## ----minimal-loop-------------------------------------------------------------
par <- rb0
opt <- make_mize(
  method = "BFGS",
  par = par,
  fg = rb_fg,
  max_iter = 3,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)

while (!opt$is_terminated) {
  par_old <- par
  step <- mize_step(opt, par, rb_fg)
  opt <- step$opt
  par <- step$par
  if (opt$is_terminated) {
    break
  }

  step_info <- mize_step_summary(
    opt,
    par,
    rb_fg,
    par_old = par_old,
    calc_fn = FALSE,
    calc_gr = FALSE
  )
  opt <- step_info$opt
  if (opt$is_terminated) {
    break
  }

  opt <- check_mize_convergence(step_info)
}

opt[c("status", "converged", "message", "terminate")]

## ----assert-minimal-loop, include = FALSE-------------------------------------
stopifnot(
  isTRUE(opt$is_terminated),
  identical(opt$status, "budget_exhausted"),
  identical(opt$converged, FALSE),
  identical(opt$terminate$what, "max_iter")
)

## ----summary-schema, echo=FALSE-----------------------------------------------
summary_schema <- data.frame(
  `Field group` = c(
    "Core",
    "Requested observations",
    "Line search",
    "Direction diagnostics",
    "Momentum"
  ),
  Examples = c(
    "`iter`, `nf`, `ng`, `nh`, `nhi`, `step`",
    "`f`, `g2n`, `ginfn`",
    "`alpha`, `ls_reason`, `ls_outcome`, `ls_nf`, `ls_ng`",
    "`direction_reason`",
    "`mu`"
  ),
  Availability = c(
    "Generally available",
    "When calculated, cached, or requested",
    "Searches that report these diagnostics",
    "Methods that report direction provenance",
    "Methods with a momentum stage"
  ),
  check.names = FALSE
)
knitr::kable(summary_schema)

## ----production-loop-definition-----------------------------------------------
numeric_or_na <- function(x) {
  if (is.null(x)) NA_real_ else as.numeric(x)
}

character_or_na <- function(x) {
  if (is.null(x)) NA_character_ else as.character(x)
}

progress_record <- function(step_info) {
  data.frame(
    iter = as.integer(step_info$iter),
    f = as.numeric(step_info$f),
    nf = as.integer(step_info$nf),
    ng = as.integer(step_info$ng),
    nh = as.integer(step_info$nh),
    nhi = as.integer(step_info$nhi),
    step = as.numeric(step_info$step),
    alpha = numeric_or_na(step_info$alpha),
    ls_outcome = character_or_na(step_info$ls_outcome)
  )
}

bind_progress <- function(records) {
  if (length(records) == 0L) {
    return(data.frame(
      iter = integer(),
      f = double(),
      nf = integer(),
      ng = integer(),
      nh = integer(),
      nhi = integer(),
      step = double(),
      alpha = double(),
      ls_outcome = character()
    ))
  }
  do.call(rbind, records)
}

initialize_stateful <- function(opt, par, fg) {
  initial_info <- mize_step_summary(
    opt,
    par,
    fg,
    calc_fn = TRUE,
    calc_gr = FALSE
  )
  state <- list(
    opt = initial_info$opt,
    par = par,
    best_par = NULL,
    best_f = Inf,
    progress = list()
  )

  if (state$opt$is_terminated || is.null(initial_info$f)) {
    return(state)
  }
  if (!is.finite(initial_info$f)) {
    state$opt <- check_mize_convergence(initial_info)
    return(state)
  }

  state$best_par <- par
  state$best_f <- initial_info$f
  state$progress[[1L]] <- progress_record(initial_info)
  state
}

advance_stateful <- function(state, fg, max_steps = Inf) {
  steps <- 0L

  while (!state$opt$is_terminated && steps < max_steps) {
    par_old <- state$par
    step <- mize_step(state$opt, state$par, fg)
    state$opt <- step$opt
    state$par <- step$par
    if (state$opt$is_terminated) {
      break
    }

    step_info <- mize_step_summary(
      state$opt,
      state$par,
      fg,
      par_old = par_old,
      calc_fn = TRUE
    )
    state$opt <- step_info$opt
    if (state$opt$is_terminated) {
      break
    }

    if (!is.null(step_info$f) && is.finite(step_info$f)) {
      if (step_info$f < state$best_f) {
        state$best_f <- step_info$f
        state$best_par <- state$par
      }

      next_record <- length(state$progress) + 1L
      state$progress[[next_record]] <- progress_record(step_info)
    }

    state$opt <- check_mize_convergence(step_info)
    steps <- steps + 1L
  }

  state
}

## ----production-run-----------------------------------------------------------
controls <- list(
  max_iter = 30,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)

new_bfgs_state <- function() {
  opt <- do.call(
    make_mize,
    c(list(method = "BFGS", par = rb0, fg = rb_fg), controls)
  )
  initialize_stateful(opt, rb0, rb_fg)
}

initial_state <- new_bfgs_state()
initial_observation <- data.frame(
  `Starting objective` = initial_state$best_f,
  `Function calls` = initial_state$opt$counts$fn,
  `Progress records` = length(initial_state$progress),
  check.names = FALSE
)
knitr::kable(initial_observation)

uninterrupted <- advance_stateful(initial_state, rb_fg)
uninterrupted_progress <- bind_progress(uninterrupted$progress)
knitr::kable(
  tail(uninterrupted_progress, 3),
  digits = 6,
  row.names = FALSE
)
uninterrupted$opt[c("status", "converged", "message", "terminate")]

## ----assert-production-run, include = FALSE-----------------------------------
initial_recheck <- mize_step_summary(
  initial_state$opt,
  initial_state$par,
  rb_fg,
  calc_fn = TRUE,
  calc_gr = FALSE
)
stopifnot(
  identical(initial_state$best_par, rb0),
  identical(initial_state$best_f, rb_fg$fn(rb0)),
  identical(initial_state$opt$counts$fn, 1),
  identical(length(initial_state$progress), 1L),
  identical(initial_recheck$nf, initial_state$opt$counts$fn),
  identical(initial_recheck$f, initial_state$best_f),
  isTRUE(uninterrupted$opt$is_terminated),
  identical(uninterrupted$opt$status, "budget_exhausted"),
  identical(uninterrupted$opt$converged, FALSE),
  identical(uninterrupted$opt$terminate$what, "max_iter"),
  is.finite(uninterrupted$best_f),
  !is.null(uninterrupted$best_par)
)

## ----exact-function-budget----------------------------------------------------
budget_opt <- make_mize(
  method = "SD",
  line_search = "constant",
  step0 = 0.0001,
  par = rb0,
  fg = rb_fg,
  max_fn = 0,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
budget_state <- initialize_stateful(
  budget_opt,
  rb0,
  rb_fg
)
budget_result <- data.frame(
  `Best recorded` = !is.null(budget_state$best_par),
  `Function calls` = budget_state$opt$counts$fn,
  Status = budget_state$opt$status,
  Reason = budget_state$opt$terminate$what,
  check.names = FALSE
)
knitr::kable(budget_result)

## ----assert-exact-function-budget, include = FALSE----------------------------
stopifnot(
  is.null(budget_state$best_par),
  identical(budget_state$best_f, Inf),
  identical(length(budget_state$progress), 0L),
  identical(budget_state$opt$counts$fn, 0),
  identical(budget_state$opt$status, "budget_exhausted"),
  identical(budget_state$opt$terminate$what, "max_fn")
)

## ----checkpoint-and-resume----------------------------------------------------
partial <- advance_stateful(
  new_bfgs_state(),
  rb_fg,
  max_steps = 10
)

checkpoint_file <- tempfile(fileext = ".rds")
saveRDS(partial, checkpoint_file)
rm(partial)

resumed <- readRDS(checkpoint_file)
unlink(checkpoint_file)
checkpoint_progress_is_list <- is.list(resumed$progress)
checkpoint_records_are_typed <- all(vapply(
  resumed$progress,
  is.data.frame,
  logical(1)
))
resumed <- advance_stateful(resumed, rb_fg)

## ----compare-lifecycles-------------------------------------------------------
one_shot <- do.call(
  mize,
  c(list(par = rb0, fg = rb_fg, method = "BFGS"), controls)
)

comparison <- data.frame(
  run = c("one-shot", "stateful", "checkpoint/resume"),
  objective = c(one_shot$f, uninterrupted$best_f, resumed$best_f),
  nf = c(one_shot$nf, uninterrupted$opt$counts$fn, resumed$opt$counts$fn),
  ng = c(one_shot$ng, uninterrupted$opt$counts$gr, resumed$opt$counts$gr),
  status = c(
    one_shot$status,
    uninterrupted$opt$status,
    resumed$opt$status
  ),
  terminate = c(
    one_shot$terminate$what,
    uninterrupted$opt$terminate$what,
    resumed$opt$terminate$what
  )
)
knitr::kable(
  comparison,
  digits = 8,
  col.names = c(
    "Run",
    "Objective",
    "Function calls",
    "Gradient calls",
    "Status",
    "Reason"
  )
)

comparison_tolerance <- 1e-10
equivalent <- c(
  one_shot_best = isTRUE(all.equal(
    one_shot$par,
    uninterrupted$best_par,
    tolerance = comparison_tolerance
  )),
  one_shot_objective = isTRUE(all.equal(
    one_shot$f,
    uninterrupted$best_f,
    tolerance = comparison_tolerance
  )),
  one_shot_last = isTRUE(all.equal(
    one_shot$last_par,
    uninterrupted$par,
    tolerance = comparison_tolerance
  )),
  resumed_current = isTRUE(all.equal(
    resumed$par,
    uninterrupted$par,
    tolerance = comparison_tolerance
  )),
  resumed_best_par = isTRUE(all.equal(
    resumed$best_par,
    uninterrupted$best_par,
    tolerance = comparison_tolerance
  )),
  resumed_best = isTRUE(all.equal(
    resumed$best_f,
    uninterrupted$best_f,
    tolerance = comparison_tolerance
  ))
)

## ----assert-compare-lifecycles, include = FALSE-------------------------------
stopifnot(
  all(equivalent),
  checkpoint_progress_is_list,
  checkpoint_records_are_typed,
  all(comparison$status == "budget_exhausted"),
  all(comparison$terminate == "max_iter"),
  !file.exists(checkpoint_file)
)

## ----reinitialization-lifecycle-----------------------------------------------
counts_before <- uninterrupted$opt$counts
reinitialized <- mize_init(
  uninterrupted$opt,
  rb0,
  rb_fg
)

data.frame(
  iter_before = uninterrupted$opt$iter,
  iter_after = reinitialized$iter,
  counts_preserved = identical(counts_before, reinitialized$counts),
  terminated_after = reinitialized$is_terminated
)

## ----assert-reinitialization-lifecycle, include = FALSE-----------------------
stopifnot(
  identical(reinitialized$iter, 0),
  identical(reinitialized$counts, counts_before),
  identical(reinitialized$is_terminated, FALSE),
  is.null(reinitialized$status),
  is.null(reinitialized$terminate)
)

