#' @title Methods for \code{zouwquantreg} Objects
#'
#' @description Methods for extracting information from fitted unit-Weibull quantile regression
#' objects of class \code{\link{zouwquantreg}}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object,x fitted model object of class \code{\link{zouwquantreg}}.
#' @param digits  minimal number of _significant_ digits
#'   the estimated parameters is returned and printed. Default is \code{FALSE}.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param type character indicating type of fitted values to return.
#' @param what character indicating what component return the confidence intervals.
#' If \code{quant} return coefficients from quantile model,
#' if \code{infl} return coefficients from inflation model. If \code{both} a list with parameters from
#' quantile and inflation models.
#' @param ... additional argument(s) for methods. Currently not used.
#'
#' @importFrom stats pnorm cov2cor coef vcov printCoefmat confint
#'
#' @name methods-zouwquantreg
NULL


#' @rdname methods-zouwquantreg
#' @export

print.zouwquantreg <- function(x, digits = 4, ...) {
  inf <- ifelse(x$inflation == 0, "Zero", "One")

  cat("\n", inf, " Inflated unit-Weibull quantile regression model \n", sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Mu coefficients (quantile model with ", x$link, " link and tau = ", x$tau,  "): \n", sep = "")

  print.default(FF(x$fit_quant$mu_coefficients, Digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")

  if(!is.null(x$link.phi)) {
    cat("Phi coefficients (shape model with ", x$link.phi, " link):", "\n", sep = "")

    print.default(FF(x$fit_quant$phi_coefficients, Digits = digits), print.gap = 2, quote = FALSE)

    cat("\n")
  }
  else {
    cat("Model with constant shape parameter:", "\n", sep = "")

    print.default(FF(x$fit_quant$phi_coefficients, Digits = digits), print.gap = 2, quote = FALSE)

    cat("\n")
  }

  cat("Nu coefficients (degenerate model at ", inf, " with ", x$link.nu, " link function)\n", sep = "")

  print.default(FF(x$fit_infl$coefficients, Digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")

  invisible(x)
}


# Summary -----------------------------------------------------------------

#' @rdname methods-zouwquantreg
#' @export
summary.zouwquantreg <- function(object, ...) {

  s1 <- summary(object$fit_quant)
  s2 <- summary(object$fit_infl)

  table <- list(quant = s1$coeftable,
                infl = s2$coefficients)

  out <- list(coeftable   = table,
              call        = object$call,
              tau         = object$tau,
              link        = object$link,
              link.phi    = object$link.phi,
              link.nu     = object$link.nu,
              inflation   = object$inflation,
              npar        = c(mu  = length(object$fit_quant$mu_coefficients),
                              phi = length(object$fit_quant$phi_coefficients),
                              nu  = length(object$fit_infl$coefficients)))

  class(out) <- "summary.zouwquantreg"
  out
}

# Print output summary ----------------------------------------------------


#' @export
print.summary.zouwquantreg <- function(x, digits = 4, ...) {
  p <- x$npar[1]
  q <- x$npar[2]

  inf <- ifelse(x$inflation == 0, "Zero", "One")

  cat("\n Wald-tests for ", inf, " Inflated unit-Weibull quantile regression model", "\n" ,sep = "")

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Mu coefficients: (quantile model with ", x$link, " link and tau = ", x$tau,  "): \n", sep = "")
  printCoefmat(x$coeftable$quant[1:p, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
  cat("\n")

  if(q > 1) {
    cat("Phi coefficients (shape model with ", x$link.phi, " link):", "\n", sep = "")
    printCoefmat(x$coeftable$quant[1:p + q, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
    cat("\n")
  }

  else {
    id <- 1:p + q
    cat("Model with constant shape:", "\n", sep = "")
    printCoefmat(x$coeftable$quant[id[length(id)], , drop = FALSE], digits = digits, has.Pvalue = TRUE)
    cat("\n")
  }

  cat("Nu coefficients (degenerate model at ", inf, " with ", x$link.nu, " link function): \n", sep = "")
  printCoefmat(x$coeftable$infl[, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
  cat("\n")

  invisible(x)
}


# coef function -----------------------------------------------------------

#' @rdname methods-zouwquantreg
#' @export
coef.zouwquantreg <- function(object, ...) {
  if(!missing(...)) {
    warning("Extra arguments discarded")
  }
  list(
    quant = coef(object$fit_quant),
    infl = coef(object$fit_infl)
  )
}


# vcov function -----------------------------------------------------------

#' @rdname methods-zouwquantreg
#' @export
vcov.zouwquantreg <- function(object, ...) {
  if(!missing(...)) {
    warning("Extra arguments discarded")
  }
  list(
    quant = vcov(object$fit_quant),
    infl = vcov(object$fit_infl)
  )
}


# logLik function ---------------------------------------------------------

#' @rdname methods-zouwquantreg
#' @export
logLik.zouwquantreg <- function(object, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  ll <- object$loglik
  attr(ll, "df")   <- object$npar
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}


# confint function --------------------------------------------------------

#' @rdname methods-zouwquantreg
#' @export
confint.zouwquantreg <- function(object, parm, level = 0.95, what = c("both", "quant", "infl"), ...)
{
  if (length(what) > 1) what <- "both"

  pnames <- list(
    names(object$fit_quant$coefficients),
    names(object$fit_infl$coefficients)
  )
  if (missing(parm)) {
    parm <- pnames
  }
  else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }

  switch (what,
    "both" = list(
      quant = confint(object$fit_quant, parm=parm[[1]], level=level),
      infl = confint(object$fit_infl, parm=parm[[2]], level=level)
    ),
    "quant" = confint(object$fit_quant, parm=parm[[1]], level=level),
    "infl" = confint(object$fit_infl, parm=parm[[2]], level=level)
  )
}

#' @rdname methods-zouwquantreg
#' @export
fitted.zouwquantreg <- function(object, type = c("mu", "phi", "nu", "all"),  ...) {

  if (!missing(...)) {
    warning("Extra arguments discarded")
  }

  if (length(type) != 1) type <- "mu"

  switch (type,
          "all" = list(mu = fitted(object$fit_quant, type = "mu"),
                       phi = fitted(object$fit_quant, type = "phi"),
                       nu = fitted(object$fit_infl)),
          "mu"  = fitted(object$fit_quant, type = "mu"),
          "phi" = fitted(object$fit_quant, type = "phi"),
          "nu"  = fitted(object$fit_infl)
  )
}

