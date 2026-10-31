## ----setup, include = FALSE, echo = FALSE, message = FALSE--------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>")
library(mize)

## ----callback-factory---------------------------------------------------------
distance_epsilon <- 1e-10

make_mmds_callbacks <- function(
  distances,
  combined = FALSE,
  epsilon = distance_epsilon
) {
  target <- as.matrix(distances)
  preparation_count <- 0L

  prepare <- function(par) {
    preparation_count <<- preparation_count + 1L
    coordinates <- matrix(par, ncol = 2, byrow = TRUE)
    list(
      coordinates = coordinates,
      distances = as.matrix(stats::dist(coordinates))
    )
  }

  stress <- function(prepared) {
    0.5 * sum((target - prepared$distances)^2)
  }

  gradient <- function(prepared) {
    weights <- (target - prepared$distances) /
      (prepared$distances + epsilon)
    coordinates <- prepared$coordinates
    result <- matrix(
      nrow = nrow(coordinates),
      ncol = ncol(coordinates)
    )

    for (i in seq_len(nrow(coordinates))) {
      differences <- sweep(-coordinates, 2, -coordinates[i, ])
      result[i, ] <- colSums(differences * weights[, i])
    }

    as.vector(t(result)) * -2
  }

  callbacks <- list(
    fn = function(par) stress(prepare(par)),
    gr = function(par) gradient(prepare(par))
  )
  if (combined) {
    callbacks$fg <- function(par) {
      prepared <- prepare(par)
      list(
        fn = stress(prepared),
        gr = gradient(prepared)
      )
    }
  }

  list(
    fg = callbacks,
    preparations = function() preparation_count,
    reset_preparations = function() {
      preparation_count <<- 0L
      invisible(NULL)
    }
  )
}

## ----seeded-start-------------------------------------------------------------
target_distances <- datasets::eurodist
set.seed(42)
initial_par <- rnorm(attr(target_distances, "Size") * 2)

separate_callbacks <- make_mmds_callbacks(target_distances)
combined_callbacks <- make_mmds_callbacks(
  target_distances,
  combined = TRUE
)

## ----gradient-check-----------------------------------------------------------
gradient_check <- check_mize_gradient(
  combined_callbacks$fg,
  initial_par
)

gradient_summary <- data.frame(
  max_abs_error = gradient_check$max_abs_error,
  max_rel_error = gradient_check$max_rel_error,
  fg_fn_abs_error = gradient_check$fg_consistency$fn$abs_error,
  fg_gr_max_abs_error =
    gradient_check$fg_consistency$gr$max_abs_error
)
knitr::kable(
  signif(gradient_summary, digits = 4),
  row.names = FALSE,
  col.names = c(
    "Maximum absolute error",
    "Maximum relative error",
    "fg/fn absolute error",
    "fg/gr maximum absolute error"
  )
)

## ----assert-gradient-check, include = FALSE-----------------------------------
stopifnot(
  gradient_check$max_rel_error < 1e-3,
  isTRUE(gradient_check$fg_consistency$checked),
  gradient_check$fg_consistency$fn$abs_error < 1e-10,
  gradient_check$fg_consistency$gr$max_abs_error < 1e-10
)

## ----optimization-------------------------------------------------------------
mmds_controls <- list(
  method = "L-BFGS",
  max_iter = 200,
  abs_tol = NULL,
  rel_tol = 1e-8,
  grad_tol = NULL,
  step_tol = NULL
)

separate_callbacks$reset_preparations()
combined_callbacks$reset_preparations()

separate_result <- do.call(
  mize,
  c(
    list(par = initial_par, fg = separate_callbacks$fg),
    mmds_controls
  )
)
combined_result <- do.call(
  mize,
  c(
    list(par = initial_par, fg = combined_callbacks$fg),
    mmds_controls
  )
)

results <- list(
  separate = separate_result,
  combined = combined_result
)
initial_stress <- gradient_check$fn
pair_count <- choose(attr(target_distances, "Size"), 2)
preparation_counts <- c(
  separate = separate_callbacks$preparations(),
  combined = combined_callbacks$preparations()
)

initial_rmse_km <- sqrt(initial_stress / pair_count)
final_rmse_km <- sqrt(combined_result$f / pair_count)

callback_comparison <- data.frame(
  callbacks = c("separate fn/gr", "combined fg"),
  iterations = vapply(results, `[[`, numeric(1), "iter"),
  objective_calls = vapply(results, `[[`, numeric(1), "nf"),
  gradient_calls = vapply(results, `[[`, numeric(1), "ng"),
  distance_preparations = unname(preparation_counts),
  check.names = FALSE
)
knitr::kable(
  callback_comparison,
  row.names = FALSE,
  col.names = c(
    "Callbacks",
    "Iterations",
    "Objective calls",
    "Gradient calls",
    "Distance preparations"
  )
)

## ----assert-optimization, include = FALSE-------------------------------------
stopifnot(
  all(vapply(results, `[[`, character(1), "status") == "converged"),
  all(vapply(results, `[[`, logical(1), "converged")),
  all(vapply(
    results,
    function(result) result$terminate$what,
    character(1)
  ) == "rel_tol"),
  all(vapply(results, `[[`, numeric(1), "f") < initial_stress),
  isTRUE(all.equal(
    separate_result$par,
    combined_result$par,
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    separate_result$f,
    combined_result$f,
    tolerance = 1e-10
  )),
  identical(separate_result$nf, combined_result$nf),
  identical(separate_result$ng, combined_result$ng),
  unname(preparation_counts["combined"]) <
    unname(preparation_counts["separate"])
)

## ----geographic-reference-----------------------------------------------------
geography_path <- system.file(
  "extdata",
  "eurodist-geography.csv",
  package = "mize"
)
map_path <- system.file(
  "extdata",
  "europe-map.csv",
  package = "mize"
)
city_geography <- utils::read.csv(geography_path)
europe_map <- utils::read.csv(map_path)

project_geography <- function(longitude, latitude) {
  earth_radius_km <- 6371.009
  reference_longitude <- 9
  reference_latitude <- 48
  degrees_to_radians <- pi / 180

  cbind(
    x = earth_radius_km * cos(reference_latitude * degrees_to_radians) *
      (longitude - reference_longitude) * degrees_to_radians,
    y = earth_radius_km *
      (latitude - reference_latitude) * degrees_to_radians
  )
}

## ----align-to-geography-------------------------------------------------------
center_coordinates <- function(par) {
  coordinates <- matrix(par, ncol = 2, byrow = TRUE)
  sweep(coordinates, 2, colMeans(coordinates))
}

align_to_reference <- function(configuration, reference) {
  configuration <- scale(
    configuration,
    center = TRUE,
    scale = FALSE
  )
  reference_center <- colMeans(reference)
  centered_reference <- sweep(reference, 2, reference_center)
  decomposition <- svd(crossprod(configuration, centered_reference))
  orthogonal_transform <- decomposition$u %*% t(decomposition$v)
  scale_factor <- sum(decomposition$d) / sum(configuration^2)

  list(
    coordinates = sweep(
      scale_factor * configuration %*% orthogonal_transform,
      2,
      reference_center,
      "+"
    ),
    scale = scale_factor,
    determinant = det(orthogonal_transform)
  )
}

centered_coordinates <- center_coordinates(combined_result$par)
stopifnot(identical(city_geography$city, labels(target_distances)))
geographic_coordinates <- project_geography(
  city_geography$longitude,
  city_geography$latitude
)
alignment <- align_to_reference(
  centered_coordinates,
  geographic_coordinates
)
aligned_coordinates <- alignment$coordinates
geographic_residual_km <- sqrt(rowSums(
  (aligned_coordinates - geographic_coordinates)^2
))
geographic_rmse_km <- sqrt(mean(geographic_residual_km^2))

## ----plot-geographic-alignment, echo = FALSE, fig.width = 7, fig.height = 7, fig.alt = sprintf("Map of Europe comparing the geographic and similarity-aligned metric MDS locations of 21 cities. Black open circles mark geographic locations, blue triangles mark aligned MDS locations, and orange segments show pointwise discrepancies; the residual RMSE is about %.0f kilometres.", geographic_rmse_km)----
map_coordinates <- project_geography(
  europe_map$longitude,
  europe_map$latitude
)
old_par <- graphics::par(
  mar = c(3.6, 3.6, 3.4, 1),
  mgp = c(2.1, 0.65, 0),
  tcl = -0.25
)

graphics::plot(
  range(map_coordinates[, 1]),
  range(map_coordinates[, 2]),
  type = "n",
  asp = 1,
  axes = FALSE,
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Recovered road-distance geometry versus geographic Europe",
  cex.main = 0.95
)
longitude_breaks <- seq(-10, 30, by = 10)
latitude_breaks <- seq(35, 60, by = 5)
grid_x <- project_geography(
  longitude_breaks,
  rep(48, length(longitude_breaks))
)[, 1]
grid_y <- project_geography(
  rep(9, length(latitude_breaks)),
  latitude_breaks
)[, 2]
graphics::abline(
  v = grid_x,
  h = grid_y,
  col = "#E5E7EB",
  lty = 3
)
for (indices in split(seq_len(nrow(europe_map)), europe_map$group)) {
  graphics::polygon(
    map_coordinates[indices, 1],
    map_coordinates[indices, 2],
    col = "#F3F4F6",
    border = "#B8C0C8",
    lwd = 0.55
  )
}
graphics::axis(
  1,
  at = grid_x,
  labels = c("10°W", "0°", "10°E", "20°E", "30°E"),
  cex.axis = 0.75
)
graphics::axis(
  2,
  at = grid_y,
  labels = paste0(latitude_breaks, "°N"),
  las = 1,
  cex.axis = 0.75
)
graphics::box(col = "#8B949E")
graphics::segments(
  geographic_coordinates[, 1],
  geographic_coordinates[, 2],
  aligned_coordinates[, 1],
  aligned_coordinates[, 2],
  col = grDevices::adjustcolor("#D55E00", alpha.f = 0.78),
  lwd = 1.35
)
graphics::points(
  aligned_coordinates,
  pch = 17,
  col = "#0072B2",
  cex = 0.8
)
graphics::points(
  geographic_coordinates,
  pch = 21,
  bg = "white",
  col = "#111827",
  cex = 0.75
)
label_positions <- c(
  4, 2, 1, 2, 2, 4, 4, 4, 2, 4, 2,
  2, 2, 2, 4, 4, 4, 2, 4, 4, 4
)
graphics::text(
  geographic_coordinates[, 1],
  geographic_coordinates[, 2],
  labels = city_geography$city,
  pos = label_positions,
  offset = 0.28,
  cex = 0.55,
  col = "#111827"
)
graphics::legend(
  "topleft",
  legend = c(
    "Geographic location",
    "Aligned MDS location",
    "Residual"
  ),
  pch = c(21, 17, NA),
  pt.bg = c("white", NA, NA),
  col = c("#111827", "#0072B2", "#D55E00"),
  lty = c(NA, NA, 1),
  lwd = c(NA, NA, 1.35),
  pt.cex = c(0.75, 0.8, NA),
  cex = 0.68,
  bg = grDevices::adjustcolor("white", alpha.f = 0.9),
  box.col = "#D1D5DB",
  inset = 0.015
)
graphics::mtext(
  sprintf(
    "Similarity alignment; pointwise residual RMSE = %.0f km",
    geographic_rmse_km
  ),
  side = 3,
  line = 0.25,
  cex = 0.72,
  col = "#4B5563"
)
graphics::par(old_par)

## ----assert-geographic-alignment, include = FALSE-----------------------------
stopifnot(
  nrow(city_geography) == attr(target_distances, "Size"),
  all(is.finite(aligned_coordinates)),
  alignment$scale > 0,
  abs(abs(alignment$determinant) - 1) < 1e-12,
  geographic_rmse_km < 350
)

## ----assert-stress-invariance, include = FALSE--------------------------------
angle <- pi / 3
rotation <- matrix(
  c(cos(angle), -sin(angle), sin(angle), cos(angle)),
  nrow = 2,
  byrow = TRUE
)
translated <- sweep(
  centered_coordinates,
  2,
  c(1000, -500),
  "+"
)
rotated <- centered_coordinates %*% rotation
reflected <- centered_coordinates
reflected[, 1] <- -reflected[, 1]

stress_at <- function(coordinates) {
  combined_callbacks$fg$fn(as.vector(t(coordinates)))
}
invariant_stress <- vapply(
  list(translated, rotated, reflected),
  stress_at,
  numeric(1)
)
relative_invariance_error <- max(
  abs(invariant_stress - combined_result$f)
) / max(1, abs(combined_result$f))
stopifnot(relative_invariance_error < 1e-12)

