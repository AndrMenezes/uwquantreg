#' @importFrom stats runif
#' @name uweibull
#' @aliases uweibull duweibull puweibull quweibull ruweibull
#'
#' @title unit-Weibull distribution
#'
#' @description Density function, distribution function, quantile function, random number generation function
#' for the unit-Weibull distribution re-parametrized in terms of the \eqn{\tau}th quantile.
#'
#' @author André Felipe B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references Mazucheli, J., Menezes, A. F. B and Ghitany, M. E. (2018). The unit-Weibull distribution and
#' associated inference. \emph{Journal of Applied Probability and Statistics} \bold{13} (2), 1--22.
#'
#' Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P., Ghitany, M. E., 2020.
#' The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
#' modelling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics} \bold{47} (6), 954--974.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the \eqn{\tau}-quantile \eqn{\tau \in (0, 1)}.
#' @param phi nonnegative shape parameter.
#' @param tau the parameter to specify which quantile use in the parametrization.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{duweibull} gives the density, \code{puweibull} gives the distribution function,
#' \code{quweibull} gives the quantile function and \code{ruweibull} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' set.seed(6969)
#' x <- ruweibull(n = 1000, mu = 0.5, phi = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' hist(x, prob = TRUE, main = 'unit-Weibull distribution')
#' lines(S, duweibull(x = S, mu = 0.5, phi = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, puweibull(q = S, mu = 0.5, phi = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(quweibull(p = S, mu = 0.5, phi = 1.5, tau = 0.5), col = 2)



# Density -----------------------------------------------------------------
#' @rdname uweibull
#' @export
duweibull <- function(x, mu, phi, tau, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  lxmu <- (log(x) / log(mu))
  pdf  <- phi / x * log(tau) / log(mu) * lxmu^(phi - 1) * tau^(lxmu^phi)
  if(log) {
    log(pdf)
  } else {
    pdf
  }
}

# CDF ---------------------------------------------------------------------
#' @rdname uweibull
#' @export
puweibull <- function(q, mu, phi, tau, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q >= 0, q <= 1, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  p <- tau^((log(q) / log(mu))^phi)
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

# Quantile ----------------------------------------------------------------
#' @rdname uweibull
#' @export
quweibull <- function(p, mu, phi, tau, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p >= 0, p <= 1, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- exp(p)
  }
  mu ^ exp((log(-log(p)) - log(-log(tau))) / phi)
}

# Random numbers ----------------------------------------------------------
#' @rdname uweibull
#' @export
ruweibull <- function(n, mu, phi, tau)
{
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  u <- runif(n)
  quweibull(u, mu, phi, tau)
}

#' @importFrom stats rbinom
#' @name zouweibull
#' @aliases zouweibull dzouweibull pzouweibull qzouweibull rzouweibull
#'
#' @title Zero-or-One unit-Weibull distribution
#'
#' @description Density function, distribution function, quantile function, random number generation function
#' for the zero-or-one unit-Weibull distribution re-parametrized in terms of the \eqn{\tau}th quantile.
#'
#' @author André Felipe B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references Mazucheli, J., Menezes, A. F. B and Ghitany, M. E. (2018). The unit-Weibull distribution and
#' associated inference. \emph{Journal of Applied Probability and Statistics} \bold{13} (2), 1--22.
#'
#' Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P., Ghitany, M. E., 2020.
#' The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
#' modelling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics} \bold{47} (6), 954--974.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the \eqn{tau}-quantile \eqn{\tau \in (0, 1)}.
#' @param phi nonnegative shape parameter.
#' @param tau the parameter to specify which quantile use in the parametrization.
#' @param nu probability of inflation.
#' @param inflation specify the zero (0) or the one (1) inflation.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{dzouweibull} gives the density, \code{pzouweibull} gives the distribution function,
#' \code{qzouweibull} gives the quantile function and \code{rzouweibull} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link{uweibull}}.
#'
#' @examples
#' set.seed(6969)
#' x <- rzouweibull(n = 1e4, mu = 0.5, phi = 1.5, tau = 0.5, nu = 0.3, inflation = 1)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], l = 1000)
#' plot(ecdf(x))
#' lines(S, pzouweibull(q = S, mu = 0.5, phi = 1.5, tau = 0.5, nu = 0.3, inflation = 1), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qzouweibull(p = S, mu = 0.5, phi = 1.5, tau = 0.5, nu = 0.3, inflation = 1), col = 2)
#'
#' probs <- c(0.1, 0.2, 0.5, 0.7, 0.9)
#' (p <- pzouweibull(q = probs, mu = 0.1, phi = 2, tau = 0.5, nu = 0.5, inflation = 1))
#' qzouweibull(p, mu = 0.1, phi = 2, tau = 0.5, nu = 0.5, inflation = 1)
#'
#' (p <- pzouweibull(q = probs, mu = 0.1, phi = 2, tau = 0.5, nu = 0.5, inflation = 0))
#' qzouweibull(p, mu = 0.1, phi = 2, tau = 0.5, nu = 0.5, inflation = 0)
#'



# Density -----------------------------------------------------------------
#' @rdname zouweibull
#' @export
dzouweibull <- function(x, mu, phi, tau, nu, inflation, log = FALSE) {
  id <- x == inflation
  d_c <- (1 - nu) * duweibull(x[!id], mu, phi, tau)
  d <- ifelse(id, nu, d_c)
  if (log) {
    log(d)
  } else {
    d
  }
}


# CDF ---------------------------------------------------------------------
#' @rdname zouweibull
#' @export
pzouweibull <- function(q, mu, phi, tau, nu, inflation, lower.tail = TRUE, log.p = FALSE) {
  p <- nu * (q >= inflation) + (1 - nu) * puweibull(q, mu, phi, tau)
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}


# Quantile ----------------------------------------------------------------
#' @rdname zouweibull
#' @export
qzouweibull <- function(p, mu, phi, tau, nu, inflation, lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- exp(p)
  }
  if (length(mu) != length(p)) mu <- rep(mu, length(p))
  if (length(phi) != length(p)) phi <- rep(phi, length(p))
  if (length(nu) != length(p)) nu <- rep(nu, length(p))

  if (inflation == 1) {
    t <- p[p <= 1 - nu] / (1 - nu[p <= 1 - nu])
    w <- quweibull(t, mu[p <= 1 - nu], phi[p <= 1 - nu], tau)
    q <- ifelse(p <= 1 - nu, w, 1)
  } else if (inflation == 0) {
    t <- p[nu < p]
    w <- quweibull((t - nu[nu < p])/(1 - nu[nu < p]), mu[nu < p], phi[nu < p], tau)
    q <- ifelse(nu >= p, 0, w)
  }
  q
}

# Random numbers ----------------------------------------------------------
#' @rdname zouweibull
#' @export
rzouweibull <- function(n, mu, phi, tau, nu, inflation) {
  rd  <- rbinom(n, 1, prob = nu)
  rc  <- ruweibull(length(mu), mu, phi, tau)
  ifelse(rd == 1, inflation, rc)
}
