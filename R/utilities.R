# Compute the negative of the log-like ------------------------------------

lluwquantreg <- function(par, tau, linkinv, linkinv_phi, X, Z, y) {

  # location parameter (mu)
  p           <- ncol(X)
  beta        <- par[seq.int(length.out = p)]
  mu          <- linkinv(X %*% beta)

  # shape parameter (phi)
  if (is.null(Z)) {
    phi <- par[-seq.int(length.out = p)]
  }
  else {
    q           <- ncol(Z)
    gamma       <- par[seq.int(length.out = p) + q]
    phi         <- linkinv_phi(Z %*% gamma)
  }

  # Log-likelihood

  # ll <- suppressWarnings(duweibull(y, mu, phi, tau, log = TRUE))
  # if (any(!is.finite(ll))) NaN else -sum(ll)

  ll <- -suppressWarnings(
    sum(log(phi / y)) + sum(log( log(tau) / log(mu) )) +
      sum(log(((log(y) / log(mu))) ^ (phi - 1))) +
      log(tau) * sum((log(y) / log(mu)) ^ phi)
  )
  ll
}

# Estimate the shape parameter solving a non-linear equation --------------

est_phi <- function(y)
{
  g_phi <- function(phi, y) {
    n <- length(y)
    n / phi + sum(log(-log(y))) - n / sum((-log(y))^phi) * sum((-log(y))^phi * log(-log(y)))
  }
  out <- tryCatch(stats::uniroot(g_phi, interval = c(1e-04, 100), y=y)[["root"]], error = function(e) NULL)
  ifelse(is.null(out), 1.5, out)
}


# Format output ------------------------------------------------------------

format.perc <- function(probs, digits) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}

