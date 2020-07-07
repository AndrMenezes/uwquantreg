#' @title Zero-or-One unit-Weibull quantile regression
#'
#' @description Fit Zero-or-One unit-Weibull quantile regression.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param formula symbolic description of the quantile model.
#' @param tau the quantile(s) to be estimated, number between 0 and 1.
#' @param inflation specify the zero (0) or the one (1) inflation.
#' @param data data.frame contain the variables in the model.
#' @param link character specification of the link function in the quantile model. Default is \code{logit}.
#' @param link.phi character specification of the link function in the shape model. Default is \code{log}.
#' @param formula.nu a symbolic description of the model to be fitted for the inflation part.
#' @param link.nu character specification of the link function inflation model. Default is \code{logit}.
#' @param start a list with length equal specifying starting values for the parameters in the linear predictor.
#' The first entrie refers to continous part and the second for the inflation. Default is \code{NULL}
#' @param control a list of control arguments specified via \code{uwquantreg.control}.
#'
#'
#' @importFrom stats glm binomial fitted vcov logLik
#'
#' @rdname zouwquantreg
#' @export

zouwquantreg <- function(formula, tau, inflation, data, link, link.phi = NULL, formula.nu, link.nu, start = NULL,
                         control = uwquantreg.control())
{

  if (!is.null(start)) {
    start_quant <- start[[1]]
    start_infl <- start[[2]]
  } else {
    start_quant <- NULL
    start_infl <- NULL
  }

  if (missing(link)) link <- "logit"
  if (missing(link.nu)) link.nu <- "logit"

  ff <- Formula::as.Formula(formula)

  # adjusting data
  id              <- as.character(ff)[2]
  data_quant      <- subset(data, data[, id] != inflation)
  data_infl       <- data
  data_infl[, id] <- ifelse(data[, id] == inflation, 1, 0)

  # fitting models
  fit_quant <- uwquantreg(
    formula = ff,
    tau = tau,
    link = link,
    link.phi = link.phi,
    data = data_quant,
    start = start_quant,
    control = control
  )

  fit_infl <- glm(
    formula = formula.nu,
    family = binomial(link = link.nu),
    data = data_infl,
    start = start_infl
  )

  # Output
  out <- list(
    call       = match.call(),
    formula    = formula,
    formula.nu = formula.nu,
    link       = link,
    link.phi   = link.phi,
    link.nu    = link.nu,
    tau        = tau,
    inflation  = inflation,
    loglik     = c(logLik(fit_quant)) + c(logLik(fit_infl)),
    nobs       = nrow(data),
    npar       = fit_quant$npar + length(fit_infl$coefficients),
    fit_quant  = fit_quant,
    fit_infl   = fit_infl,
    y          = data[, id]

  )
  class(out) <- "zouwquantreg"
  out
}
