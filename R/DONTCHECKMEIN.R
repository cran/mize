# Rough code to run mize with various settings over all the test functions in
# the funconstrain package https://github.com/jlmelville/mize and
# https://github.com/jlmelville/funconstrain

# Example
# devtools::install_github("jlmelville/mize")
# devtools::install_github("jlmelville/funconstrain")
# library("mize")
# library("funconstrain")

# Test L-BFGS on just Rosenbrock logging info to screen
# res <- benchmark_funco_df(funconstrain::rosen(), verbose = TRUE, method = "L-BFGS")

# If dataset can be of variable dimensionality, set parameter n.
# Extended Rosenbrock with 1000 variables
# res <- benchmark_funco_df(funconstrain::ex_rosen(), n = 1000, verbose = TRUE, method = "L-BFGS")

# Test BFGS with loose line search on all datasets
# res <- benchmark_funco_all(method = "BFGS", scale_hess = TRUE, abs_tol = 0, rel_tol = 0, step_tol = .Machine$double.eps, ginf_tol = 0, grad_tol = 0, line_search = "mt", step_next_init = "quad", step0 = "scipy", c2 = 0.9,  max_iter = 10000)

# names: vector of funconstrain dataset names to test
# hess_fd: if TRUE create a Hessian function via finite difference approximation
# restart: if TRUE repeat optimization using optimized par as new starting point
# verbose: if TRUE log progress to console
# ...: options to pass to mize
# returns a dataframe with final function value, number of function evaluations, number
#   of gradient evaluations, number of iterations, infinity norm of g and 2-norm of g,
#   for each dataset in names.
benchmark_funco_all <-
  function(names = c(
    "rosen",
    "freud_roth",
    "powell_bs",
    "brown_bs",
    "beale",
    "jenn_samp",
    "helical",
    "bard",
    "gauss",
    "meyer",
    "gulf",
    "box_3d",
    "powell_s",
    "wood",
    "kow_osb",
    "brown_den",
    "osborne_1",
    "biggs_exp6",
    "osborne_2",
    "watson",
    "ex_rosen",
    "ex_powell",
    "penalty_1",
    "penalty_2",
    "var_dim",
    "trigon",
    "brown_al",
    "disc_bv",
    "disc_ie",
    "broyden_tri",
    "broyden_band",
    "linfun_fr",
    "linfun_r1",
    "linfun_r1z",
    "chebyquad"
  ),
  hess_fd = FALSE,
  restart = FALSE,
  verbose = FALSE,
  ...) {
    res <- NULL

    for (name in names) {
      dataset <- get(name)()
      bm <-
        benchmark_funco_df(
          dataset,
          name = name,
          hess_fd = hess_fd,
          restart = restart,
          log_name = verbose,
          ...
        )
      if (!is.null(res)) {
        res <- rbind(res, bm)
      }
      else {
        res <- bm
      }
    }
    res
  }

# dataset: one of the funconstrain datasets (or something with its format)
# n: dimensionality if needed for the dataset
# name: name to associate with the results
# hess_fd: if TRUE create a Hessian function via finite difference approximation
# restart: if TRUE repeat optimization using optimized par as new starting point
# log_name: if TRUE log progress to console
# ...: options to pass to mize
benchmark_funco_df <-
  function(dataset,
           n = NULL,
           name = NULL,
           hess_fd = FALSE,
           restart = FALSE,
           log_name = FALSE,
           ...) {
    if (class(dataset$x0) == "numeric") {
      x0 <- dataset$x0
    }
    else {
      if (!is.null(n)) {
        x0 <- dataset$x0(n = n)
      }
      else {
        x0 <- dataset$x0()
      }
    }

    if (hess_fd) {
      dataset$hs <- make_hfd(dataset$fn)
    }

    res <- mize(fg = dataset, par = x0, ...)
    if (restart) {
      res <- mize(fg = dataset, par = res$par, ...)
    }
    if (log_name) {
      message(name)
    }

    df <-
      data.frame(
        f = res$f,
        nf = res$nf,
        ng = res$ng,
        iter = res$iter
      )

    if (!is.null(res$ginfn)) {
      df$ginf <- res$ginfn
    }
    else {
      df$ginf <- norm_inf(dataset$gr(res$par))
    }
    if (!is.null(res$g2n)) {
      df$gr2 <- res$g2n
    }
    else {
      df$gr2 <- norm2(dataset$gr(res$par))
    }
    row.names(df) <- name
    df
  }

# A rough and ready finite difference Hessian approach
hfd <- function(par, fn, rel_eps = sqrt(.Machine$double.eps)) {
  hs <- matrix(0, nrow = length(par), ncol = length(par))
  for (i in 1:length(par)) {
    for (j in i:length(par)) {
      oldxi <- par[i]
      oldxj <- par[j]

      if (oldxi != 0 && oldxj != 0) {
        eps <- min(oldxi, oldxj) * rel_eps
      }
      else {
        eps <- 1e-3
      }
      if (i != j) {
        par[i] <- par[i] + eps
        par[j] <- par[j] + eps
        fpp <- fn(par)

        par[j] <- oldxj - eps
        fpm <- fn(par)

        par[i] <- oldxi - eps
        par[j] <- oldxj + eps
        fmp <- fn(par)

        par[j] <- oldxj - eps
        fmm <- fn(par)

        par[i] <- oldxi
        par[j] <- oldxj

        val <- (fpp - fpm - fmp + fmm) / (4 * eps * eps)

        hs[i, j] <- val
        hs[j, i] <- val
      }
      else {
        f <- fn(par)
        oldxi <- par[i]

        par[i] <- oldxi + 2 * eps
        fpp <- fn(par)

        par[i] <- oldxi + eps
        fp <- fn(par)

        par[i] <- oldxi - 2 * eps
        fmm <- fn(par)

        par[i] <- oldxi - eps
        fm <- fn(par)

        par[i] <- oldxi

        hs[i, i] <-
          (-fpp + 16 * fp - 30 * f + 16 * fm - fmm) / (12 * eps * eps)
      }
    }
  }
  hs
}

# Given an objective function, return the fd Hessian function
make_hfd <- function(fn, eps = 1.e-3) {
  function(par) {
    hfd(par, fn, eps)
  }
}

comp_res <- function(res1, res2) {
  sum(sign(res1$ng - res2$ng))
}

