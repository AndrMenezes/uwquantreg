#' @title Prediction Method for \code{uwquantreg} Objects
#'
#' @description Extract various types of predictions from unit-Weibull quantile regression models.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object fitted model object of class \code{\link{uwquantreg}}.
#' @param newdata optionally, a data frame in which to look for variables with which to predict.
#' If omitted, the original observations are used.
#' @param type character indicating type of predictions. The options are \code{mu}, \code{phi} and \code{all}.
#' @param interval type of interval desired. The options are \code{none} and \code{confidence}.
#' @param level coverage probability for the confidence intervals. Default is \code{0.95}.
#' @param se.fit logical. If \code{TRUE} return the asymptotic standard errors.
#' @param ... currently not used.
#'
#'
#' @importFrom stats qnorm
#' @importFrom Formula as.Formula


#' @rdname predict.uwquantreg
#' @export
predict.uwquantreg <- function(object, newdata, type, interval = "confidence", level = 0.95,
                               se.fit = FALSE, ...) {
  Vcov <- object$vcov
  p    <- length(object$mu_coefficients)
  q    <- length(object$phi_coefficients)
  ff   <- Formula::as.Formula(object$formula)

  if (length(type) != 1) type <- "mu"

  if (type == "phi" && is.null(object$data$Z)) {
    stop("It is not possible predict phi, because the model has constant shape")
  }

  if (missing(newdata)) {
    mats <- switch(type,
                   "mu"  = list(X = object$data$X),
                   "phi" = list(Z = object$data$Z),
                   "all" = list(X = object$data$X,
                                Z = object$data$Z))
  }

  else {
    mats <- switch(type,
                   "mu"  = list(X = model.matrix(ff[-2], newdata, rhs = 1)),
                   "phi" = list(Z = model.matrix(ff[-2], newdata, rhs = 2)),
                   "all" = list(X = model.matrix(ff[-2], newdata, rhs = 1),
                                Z = model.matrix(ff[-2], newdata, rhs = 2)))
  }


  if (type == "all" || type == "mu") {
    beta       <- object$mu_coefficients
    linkobj.mu <- object$link
    Vcov.beta  <- Vcov[1:p, 1:p]
    g.mu.hat   <- as.numeric(mats$X %*% beta)
  }

  if (type == "all" || type == "phi") {
    gamma       <- object$phi_coefficients
    linkobj.phi <- object$link.phi
    Vcov.gamma  <- Vcov[1:p + q, 1:p + q]
    g.phi.hat   <- as.numeric(mats$Z %*% gamma)
  }

  if (type == "mu") {
    mu.hat <- linkobj.mu$linkinv(g.mu.hat)
    out    <- cbind(fit = mu.hat)
    if (interval == "confidence") {
      J        <- linkobj.mu$mu.eta(mu.hat) * mats$X
      variance <- tcrossprod(J %*% Vcov.beta, J)
      stderror <- sqrt(diag(variance))
      qa       <- qnorm(1 - level / 2)
      out      <- cbind(out,
                        lower = out[, "fit"] - qa * stderror,
                        upper = out[, "fit"] + qa * stderror)
    }
    if(se.fit) {
      out <- cbind(out, se.fit = stderror)
    }
  }

  if (type == "phi") {
    phi.hat <- linkobj.phi$linkinv(g.phi.hat)
    out     <- cbind(fit = phi.hat)
    if (interval == "confidence") {
      J        <- linkobj.phi$mu.eta(phi.hat) * Z
      variance <- tcrossprod(J %*% Vcov.gamma, J)
      stderror <- sqrt(diag(variance))
      qa       <- qnorm(1 - level / 2)
      out      <- cbind(out,
                        lower = out[, "fit"] - qa * stderror,
                        upper = out[, "fit"] + qa * stderror)
    }
    if (se.fit) {
      out <- cbind(out, se.fit = stderror)
    }
  }

  if (type == "all") {
    mu.hat  <- linkobj.mu$linkinv(g.mu.hat)
    phi.hat <- linkobj.phi$linkinv(g.phi.hat)
    out     <- cbind(fit.mu = mu.hat, fit.phi = phi.hat)
    if (interval == "confidence") {
      J.mu   <- linkobj.mu$mu.eta(mu.hat) * X
      var.mu <- tcrossprod(J.mu %*% Vcov.beta, J.mu)
      std.mu <- sqrt(diag(var.mu))

      J.phi   <- linkobj.phi$mu.eta(phi.hat) * Z
      var.phi <- tcrossprod(J.phi %*% Vcov.gamma, J.phi)
      std.phi <- sqrt(diag(var.phi))

      qa       <- qnorm(1 - level / 2)
      out      <- cbind(out,
                        lower.mu = out[, "fit.mu"] - qa * std.mu,
                        upper.mu = out[, "fit.mu"] + qa * std.mu,
                        lower.phi = out[, "fit.phi"] - qa * std.phi,
                        upper.phi = out[, "fit.phi"] + qa * std.phi)
    }
    if (se.fit) {
      out <- cbind(out, se.fit.mu = std.mu, se.fit.phi = std.phi)
    }
  }
  if (missing(newdata)) {
    out
  }  else {
    as.data.frame(cbind(newdata, out))
  }
}

#' @title Prediction Method for \code{zouwquantreg} Objects
#'
#' @description Extract various types of predictions from Zero-or-One unit-Weibull quantile regression models.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object fitted model object of class \code{\link{uwquantreg}}.
#' @param newdata optionally, a data frame in which to look for variables with which to predict.
#' If omitted, the original observations are used.
#' @param type character indicating type of predictions. The currently option is \code{quantile}.
#' @param ... currently not used.
#'
#'
#' @importFrom stats make.link model.frame model.matrix


#' @rdname predict.zouwquantreg
#' @export
predict.zouwquantreg <- function(object, newdata, type = "quantile", ...) {

  ff1 <- Formula::as.Formula(object$formula)
  ff2 <- Formula::as.Formula(object$formula.nu)
  tau <- object$tau
  inf <- object$inflation

  if (missing(newdata)) {
    mf <- model.frame(ff1, data = object$fit_infl$data)
    X  <- model.matrix(ff1, data = object$fit_infl$data, rhs = 1)
    if (length(ff1)[2] == 1) {Z <- NULL}
    else {Z  <- model.matrix(ff1, data = object$fit_infl$data, rhs = 1)}
    mats <- list(
      X = X,
      Z = Z,
      W = model.matrix(object$fit_infl)
    )
  } else {
    mats <- list(X = model.matrix(ff1[-2], newdata, rhs = 1),
                 W = model.matrix(ff2[-2], newdata, rhs = 1))
    if(length(ff1)[2] == 2) {
      mats$Z <- model.matrix(ff1[-2], newdata, rhs = 2)
    }
  }

  beta       <- object$fit_quant$mu_coefficients
  linkobj.mu <- make.link(object$link)
  g.mu.hat   <- c(mats$X %*% beta)
  mu.hat     <- linkobj.mu$linkinv(g.mu.hat)

  alpha      <- object$fit_infl$coefficients
  linkobj.nu <- make.link(object$link.nu)
  g.nu.hat   <- c(mats$W %*% alpha)
  nu.hat     <- linkobj.nu$linkinv(g.nu.hat)

  if(length(ff1)[2] == 2) {
    gamma       <- object$fit_quant$phi_coefficients
    linkobj.phi <- make.link(object$linkobj.phi)
    g.phi.hat    <- c(mats$Z %*% gamma)
    phi.hat     <- linkobj.phi$linkinv(g.phi.hat)

  } else {
    phi.hat <- object$fit_quant$phi_coefficients
  }

  out <- qzouweibull(
    p = rep(tau, length(mu.hat)),
    mu = mu.hat,
    phi = phi.hat,
    tau = tau,
    nu = nu.hat,
    inflation = inf
  )
  out <- as.data.frame(cbind(ypred=out, tau=tau))
  if (missing(newdata)) {
    out
  }  else {
    as.data.frame(cbind(newdata, out))
  }

}

