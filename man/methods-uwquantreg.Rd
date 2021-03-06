% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{methods-uwquantreg}
\alias{methods-uwquantreg}
\alias{print.uwquantreg}
\alias{summary.uwquantreg}
\alias{coef.uwquantreg}
\alias{vcov.uwquantreg}
\alias{logLik.uwquantreg}
\alias{confint.uwquantreg}
\alias{fitted.uwquantreg}
\title{Methods for \code{uwquantreg} Objects}
\usage{
\method{print}{uwquantreg}(x, digits = 4, ...)

\method{summary}{uwquantreg}(object, correlation = FALSE, ...)

\method{coef}{uwquantreg}(object, ...)

\method{vcov}{uwquantreg}(object, ...)

\method{logLik}{uwquantreg}(object, ...)

\method{confint}{uwquantreg}(object, parm, level = 0.95, ...)

\method{fitted}{uwquantreg}(object, type = c("mu", "phi", "all"), ...)
}
\arguments{
\item{digits}{minimal number of _significant_ digits}

\item{...}{additional argument(s) for methods. Currently not used.}

\item{object, x}{fitted model object of class \code{\link{uwquantreg}}.}

\item{correlation}{logical; if \code{TRUE}, the correlation matrix of
the estimated parameters is returned and printed. Default is \code{FALSE}.}

\item{parm}{a specification of which parameters are to be given confidence intervals,
either a vector of numbers or a vector of names. If missing, all parameters are considered.}

\item{level}{the confidence level required.}

\item{type}{character indicating type of fitted values to return.}
}
\description{
Methods for extracting information from fitted unit-Weibull quantile regression
objects of class \code{\link{uwquantreg}}.
}
\author{
André F. B. Menezes \email{andrefelipemaringa@gmail.com}
}
