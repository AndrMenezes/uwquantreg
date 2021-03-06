% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hnp.R
\name{hnp}
\alias{hnp}
\alias{hnp.uwquantreg}
\alias{hnp.zouwquantreg}
\title{Half-Normal Plot with Simulation Envelopes for \code{uwquantreg} and \code{zouwquantreg} Objects}
\usage{
hnp(object, ...)

\method{hnp}{uwquantreg}(
  object,
  nsim = 99,
  halfnormal = TRUE,
  plot = TRUE,
  level = 0.95,
  resid.type = c("cox-snell", "quantile"),
  ...
)

\method{hnp}{zouwquantreg}(object, nsim = 99, halfnormal = TRUE, plot = TRUE, level = 0.95, ...)
}
\arguments{
\item{object}{fitted model object of class \code{\link{uwquantreg}}.}

\item{...}{currently not used.}

\item{nsim}{number of simulations used to compute envelope. Default is 99.}

\item{halfnormal}{logical. If \code{TRUE}, a half-normal plot is produced. If \code{FALSE}, a normal plot is produced.}

\item{plot}{Should the (half-)normal plot be plotted? Default is \code{TRUE}.}

\item{level}{confidence level of the simulated envelope. Default is 0.95.}

\item{resid.type}{type of residuals to be used: The currently options are \code{cox-snell} or \code{quantile}.}
}
\description{
Produces a (half-)normal plot from a fitted model object of class \code{\link{uwquantreg}} and \code{\link{uwquantreg}}.
}
\references{
Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P., Ghitany, M. E., 2020.
The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
modelling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics} \bold{47} (6), 954--974.
}
\author{
André F. B. Menezes \email{andrefelipemaringa@gmail.com}
}
