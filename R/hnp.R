#' @title Half-Normal Plot with Simulation Envelopes for \code{uwquantreg} and \code{zouwquantreg} Objects
#'
#' @description Produces a (half-)normal plot from a fitted model object of class \code{\link{uwquantreg}} and \code{\link{uwquantreg}}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object fitted model object of class \code{\link{uwquantreg}}.
#' @param nsim number of simulations used to compute envelope. Default is 99.
#' @param resid.type type of residuals to be used: The currently options are \code{cox-snell} or \code{quantile}.
#' @param halfnormal logical. If \code{TRUE}, a half-normal plot is produced. If \code{FALSE}, a normal plot is produced.
#' @param plot Should the (half-)normal plot be plotted? Default is \code{TRUE}.
#' @param level confidence level of the simulated envelope. Default is 0.95.
#' @param ... currently not used.
#'
#' @references Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P., Ghitany, M. E., 2020.
#' The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
#' modelling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics} \bold{47} (6), 954--974.
#'
#'
#' @importFrom stats qnorm qexp residuals quantile median
#' @importFrom graphics plot lines
#'
#' @name hnp
NULL

#' @rdname hnp
#' @export
hnp <- function(object, ...){
  UseMethod("hnp", object)
}


#' @rdname hnp
#' @export
hnp.uwquantreg <- function(object, nsim = 99,
                           halfnormal = TRUE, plot = TRUE, level = 0.95,
                           resid.type = c("cox-snell", "quantile"), ...) {

  mu   <- object$fitted.values$mu
  phi  <- object$fitted.values$phi
  tau  <- object$tau
  init <- object$coefficients

  control <- object$control

  data     <- object$df
  formula  <- object$formula
  link     <- object$link$name
  link.phi <- object$link.phi$name
  n        <- nrow(data)
  id       <- as.character(formula)[2]

  ysim    <- matrix(
    data = ruweibull(n = n * nsim, mu = mu, phi = phi, tau = tau),
    nrow = n,
    ncol = nsim
  )

  res_sim <- sapply(1:nsim, function(j) {
    df      <- data
    df[,id] <- ysim[,j]
    obj     <- tryCatch(expr = uwquantreg(
      formula = formula,
      tau = tau,
      link = link,
      link.phi = link.phi,
      data = df,
      start = init,
      control = control
      ),
      error = function(e) NULL)
    if(!is.null(obj)) {
      residuals(obj, type = resid.type)
    }
    else {
      rep(NaN, n)
    }
  })

  is.na(res_sim) <- sapply(res_sim, is.infinite)

  if (halfnormal) {
    res_obs <- sort(abs(residuals(object, type = resid.type)))
    res_sim <- apply(res_sim, 2, function(x) sort(abs(x), na.last = TRUE))
    if (resid.type == "quantile") {
      res_teo <- qnorm((1:n + n - 1/8) / (2 * n + 0.5))
    }
    if (resid.type == "cox-snell") {
      res_teo <- qexp((1:n + n - 1/8) / (2 * n + 0.5))
    }
  }

  else {
    res_obs <- sort(residuals(object, type = resid.type))
    res_sim <- apply(res_sim, 2, function(x) sort(x, na.last = TRUE))
    if (resid.type == "quantile") {
      res_teo <- qnorm((1:n - 3 / 8) / (n + 1 / 4))
    }
    if (resid.type == "cox-snell") {
      res_teo <- qexp((1:n - 3 / 8) / (n + 1 / 4))
    }
  }

  alpha   <- (1 - level)/2
  res_lwr <- apply(res_sim, 1, quantile, probs = alpha, na.rm = T)
  res_upr <- apply(res_sim, 1, quantile, probs = 1 - alpha, na.rm = T)
  res_mid <- apply(res_sim, 1, median, na.rm = T)

  if (plot) {
    Ry <- c(min(res_lwr), max(res_upr))
    Rx <- range(res_teo)

    plot(
      x = res_teo,
      y = res_obs,
      xlab = 'Theoretical quantiles',
      ylab = 'Residuals',
      xlim = Rx,
      ylim = Ry,
      bty = 'o',
      pch = 3
    )
    lines(x = res_teo, y = res_lwr)
    lines(x = res_teo, y = res_upr)
    lines(x = res_teo, y = res_mid, lty = 2)
  }

  list(
    res_obs = res_obs,
    res_teo = res_teo,
    lower   = res_lwr,
    median  = res_mid,
    upper   = res_upr
  )
}

#' @rdname hnp
#' @export

hnp.zouwquantreg <- function(object, nsim = 99, halfnormal = TRUE, plot = TRUE, level = 0.95, ...) {

  mu  <- fitted(object, type = "mu")
  phi <- fitted(object, type = "phi")
  nu  <- fitted(object, type = "nu")
  tau <- object$tau

  data       <- object$fit_infl$data
  formula    <- object$formula
  formula.nu <- object$formula.nu
  link       <- object$link
  link.phi   <- object$link.phi
  link.nu    <- object$link.nu
  control    <- object$control
  inflation  <- object$inflation
  init       <- list(object$fit_quant$coefficients, object$fit_infl$coefficients)

  n  <- nrow(data)
  id <- as.character(Formula::as.Formula(formula))[2]

  ysim    <- matrix(rzouweibull(n * nsim, mu, phi, tau, nu, inflation), nrow = n, ncol = nsim)

  res_sim <- sapply(1:nsim, function(j) {
    df      <- data
    df[,id] <- ysim[,j]
    obj     <- tryCatch(zouwquantreg(
      formula = formula,
      tau  = tau,
      inflation  = inflation,
      link = link,
      link.phi = link.phi,
      formula.nu = formula.nu,
      link.nu = link.nu,
      data   = df,
      start = init,
      control = control
    ),
    error = function(e) NULL)
    if (!is.null(obj)) {
      residuals(obj)
    }
    else {
      rep(NaN, n)
    }
  })

  is.na(res_sim) <- sapply(res_sim, is.infinite)

  if (halfnormal) {
    res_obs <- sort(abs(residuals(object)))
    res_sim <- apply(res_sim, 2, function(x) sort(abs(x), na.last = TRUE))
    res_teo <- qnorm((1:n + n - 1/8) / (2 * n + 0.5))
  }

  else {
    res_obs <- sort(residuals(object))
    res_sim <- apply(res_sim, 2, function(x) sort(x, na.last = TRUE))
    res_teo <- qnorm((1:n - 3 / 8) / (n + 1 / 4))
  }

  alpha   <- (1 - level)/2
  res_lwr <- apply(res_sim, 1, quantile, probs = alpha, na.rm = T)
  res_upr <- apply(res_sim, 1, quantile, probs = 1 - alpha, na.rm = T)
  res_mid <- apply(res_sim, 1, median, na.rm = T)

  if (plot) {
    Ry <- c(min(res_lwr), max(res_upr))
    Rx <- range(res_teo)

    plot(x = res_teo, y = res_obs,
         xlab = 'Theoretical quantiles',
         ylab = 'Residuals', xlim = Rx, ylim = Ry,
         bty = 'o', pch = 3)
    lines(x = res_teo, y = res_lwr)
    lines(x = res_teo, y = res_upr)
    lines(x = res_teo, y = res_mid, lty = 2)
  }

  list(
    res_obs = res_obs,
    res_teo = res_teo,
    lower   = res_lwr,
    median  = res_mid,
    upper   = res_upr
  )

}


