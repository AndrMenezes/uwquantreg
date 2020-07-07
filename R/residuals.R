#' @title Residuals Method for \code{uwquantreg} Objects
#'
#' @description Extract various types of residuals from unit-Weibull quantile regression models:
#' cox-snell and quantile residuals for now.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object fitted model object of class \code{\link{uwquantreg}}.
#' @param type character indicating type of residuals.
#' @param ... currently not used.
#'
#' @references Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P., Ghitany, M. E., 2020.
#' The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
#' modelling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics} \bold{47} (6), 954--974.
#'
#'

#' @rdname residuals.uwquantreg
#' @export
residuals.uwquantreg <- function(object, type = c("cox-snell", "quantile"), ...) {

  mu  <- object$fitted.values$mu
  phi <- object$fitted.values$phi
  tau <- object$tau
  y   <- object$data$y

  res <- switch(type,
                "cox-snell" = {
                  -log(puweibull(y, mu, phi, tau, lower.tail = FALSE))
                },
                "quantile" = {
                  qnorm(puweibull(y, mu, phi, tau))
                })

  res
}


#' @title Residuals Method for \code{zouwquantreg} Objects
#'
#' @description Extract quantile residuals from Zero-or-One unit-Weibull quantile regression models.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object fitted model object of class \code{\link{zouwquantreg}}.
#' @param ... currently not used.

#' @rdname residuals.zouwquantreg
#' @export
#'
residuals.zouwquantreg <- function(object, ...) {

  mu  <- fitted(object, type = "mu")
  phi <- fitted(object, type = "phi")
  nu  <- fitted(object, type = "nu")
  tau <- object$tau
  y   <- object$y
  inf <- object$inflation
  id  <- which(y != inf)
  n   <- length(y)

  if (inf == 0) {
    res_inf <- sapply(1:n, function(i) runif(1, min = 0, max = nu[i]))
  }
  else {
    res_inf <- sapply(1:n, function(i) runif(1, min = nu[i], max = 1))
  }

  res_cont <- nu[id] + (1 - nu[id]) * puweibull(q =  y[id], mu = mu, phi = phi, tau = tau)

  res <- ifelse(y == inf, res_inf, res_cont)
  qnorm(res)
}
