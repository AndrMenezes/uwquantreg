% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.R
\name{residuals.uwquantreg}
\alias{residuals.uwquantreg}
\title{Residuals Method for \code{uwquantreg} Objects}
\usage{
\method{residuals}{uwquantreg}(object, type = c("cox-snell", "quantile"), ...)
}
\arguments{
\item{object}{fitted model object of class \code{\link{uwquantreg}}.}

\item{type}{character indicating type of residuals.}

\item{...}{currently not used.}
}
\description{
Extract various types of residuals from unit-Weibull quantile regression models:
cox-snell and quantile residuals for now.
}
\references{
Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P., Ghitany, M. E., 2020.
The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
modelling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics} \bold{47} (6), 954--974.
}
\author{
André F. B. Menezes \email{andrefelipemaringa@gmail.com}
}
