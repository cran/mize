cm_weights <- function(mu) { c(1, mu, mu * mu, mu * mu * mu) }

nag_weights <- function(mu) { c(1 + mu, mu * mu, mu * mu * mu, 0) }

mumat <- function(wfun, mus = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99)) {
  rel_weights <- matrix(nrow = length(mus), ncol = 5)
  for (i in 1:length(mus)) {
    mu <- mus[i]
    weights <- wfun(mu)
    rel_weights[i, ] <- c(mu, weights / sum(weights))
  }
  colnames(rel_weights) <- c("mu", "s4", "s3", "s2", "s1")
  
  rel_weights
}

nagd_weights <- function(beta, mu) { 
  # ubu <- 1 - mu + beta * mu
  # ubub <- beta
  # c(ubu, ubu * mu, ubu * mu * mu, (1 - beta) * mu * mu * mu)
  
}

cmd_weights <- function(mu) {
  mu1 <- 1 - mu
  c(mu1, mu * mu1, mu * mu * mu1, mu * mu * mu)
}