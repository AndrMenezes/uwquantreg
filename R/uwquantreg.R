#' @name uwquantreg.control
#' @aliases uwquantreg.control
#'
#' @title Control Parameters for unit-Weibull quantile regression
#'
#' @description Parameters that control fitting of unit-weibull quantile regression models using \code{uwquantreg}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param method characters string specifying the method argument passed to \code{optim}.
#' @param maxit integer specifying the \code{maxit} argument (maximal number of iterations) passed to \code{optim}.
#' @param reltol relative convergence tolerance passed to \code{optim}.
#' @param ... arguments passed to \code{optim}.
#'
#' @references Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P., Ghitany, M. E., 2020.
#' The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
#' modelling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics} \bold{47} (6), 954--974.
#'

#' @rdname uwquantreg.control
#' @export
uwquantreg.control <- function(method = "BFGS", maxit = 5000, reltol = 1e-8, ...) {
  out <- list(
    method = method,
    maxit = maxit,
    reltol = reltol
  )
  out <- c(out, list(...))
  if (!is.null(out$fnscale)) warning("fnscale must not be modified")
  out$fnscale <- 1
  out
}

#' @title unit-Weibull quantile regression
#'
#' @description Fit unit-Weibull quantile regression.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param formula symbolic description of the quantile model.
#' @param tau the quantile(s) to be estimated, number between 0 and 1.
#' @param data data.frame contain the variables in the model.
#' @param link character specification of the link function in the quantile model. Default is \code{logit}.
#' @param link.phi character specification of the link function in the shape model. Default is \code{log}.
#' @param start an optional vector with starting values for all parameters
#' @param control a list of control arguments specified via \code{uwquantreg.control}.
#' @param y numeric vector of response variable.
#' @param X numeric matrix. Regressor matrix for the quantile model.
#' @param Z numeric matrix. Regressor matrix for the shape model. Default is constant shape model.
#'
#'
#' @importFrom stats optim make.link model.frame model.matrix model.response
#' @importFrom quantreg rq.fit
#' @importFrom Formula as.Formula
#' @name uwquantreg
NULL

#' @rdname uwquantreg
#' @export
uwquantreg <- function(formula, tau, data, link, link.phi = NULL, start = NULL,
                       control = uwquantreg.control()) {


  if (missing(link)) link <- "logit"

  # building model matrices
  ff <- Formula::as.Formula(formula)
  mf <- model.frame(ff, data = data)
  y  <- model.response(mf, "numeric")
  X  <- model.matrix(ff, data = data, rhs = 1)
  p  <- ncol(X)

  if (length(ff)[2] == 1) {
    Z <- NULL
    q <- 1
  } else {
    Z  <- model.matrix(ff, data = data, rhs = 2)
    q  <- ncol(Z)
  }

  ## check response variable
  if(!(min(y) > 0 & max(y) < 1)) {
    stop("invalid dependent variable, all observations must be in (0, 1)")
  }

  fit <- uwquantreg.fit(
    y = y,
    X = X,
    Z = Z,
    tau = tau,
    link = link,
    link.phi = link.phi,
    start = start,
    control = control
  )

  mu_coefficients  <- fit$par[1:p]
  phi_coefficients <- fit$par[1:q + p]
  vcov             <- solve(fit$hessian)

  # fitted values
  linkobj.mu  <- make.link(link)
  fitted_mu   <- linkobj.mu$linkinv(X %*% mu_coefficients)
  fitted_phi  <- phi_coefficients
  linkobj.phi <- NULL
  if (!is.null(Z)) {
    linkobj.phi <- make.link(link.phi)
    fitted_phi  <- linkobj.phi$linkinv(Z %*% phi_coefficients)
  }

  # Output
  out <- list(
    call             = match.call(),
    formula          = ff,
    control          = control,
    link             = linkobj.mu,
    link.phi         = linkobj.phi,
    # link             = link,
    # link.phi         = link.phi,
    tau              = tau,
    loglik           = -fit$value,
    vcov             = vcov,
    coefficients     = fit$par,
    mu_coefficients  = mu_coefficients,
    phi_coefficients = phi_coefficients,
    fitted.values    = list(mu = as.numeric(fitted_mu),
                            phi = as.numeric(fitted_phi)),
    nobs             = length(y),
    npar             = length(fit$par),
    df.residual      = length(y) - length(fit$par),
    data             = list(X = X, Z = Z, y = y),
    df               = data
  )
  class(out) <- "uwquantreg"
  out
}

#' @rdname uwquantreg
#' @export
uwquantreg.fit <- function(y, X, Z = NULL, tau, link, link.phi, start, control = uwquantreg.control()) {
  n <- length(y)
  p <- ncol(X)

  linkobj <- make.link(link)
  method  <- control$method
  control$method <- NULL

  # Initial guess
  if (is.null(start)) {
    # For beta
    ystar        <- linkobj$linkfun(y)
    reg_ini      <- suppressWarnings(quantreg::rq.fit(X, ystar, tau = tau))
    start        <- reg_ini$coefficients
    names(start) <- colnames(X)

    # For phi/gamma
    if (is.null(Z)) {
      q     <- 1
      phi   <- est_phi(y)
      start <- c(start, "phi" = phi)
      linkobj.phi <- NULL
    }
    else {
      q     <- ncol(Z)
      gamma <- rep(0, q)
      start <- c(start, gamma)

      names(start)[1:q + p] <- colnames(Z)
      linkobj.phi   <- make.link(link.phi)
    }
  }

  # Maximization
  opt <- optim(
    par     = start,
    fn      = lluwquantreg,
    method  = method,
    hessian = TRUE,
    control = control,
    X       = X,
    Z       = Z,
    y       = y,
    tau     = tau,
    linkinv = linkobj$linkinv,
    linkinv_phi = linkobj.phi$linkinv
  )
  ## check if the optim converged
  if (opt$convergence > 0) {
    opt$converged <- FALSE
    warning("optimization failed to converge")
  } else {
    opt$converged <- TRUE
  }
  opt
}
