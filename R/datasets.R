#-----------------------------------------------------------------------
#' @title Access to Piped Water Supply
#'
#' @description the access of people in households with piped water supply in the cities of Brazil from
#' the Southeast and Northeast regions. Information obtained during the census of 2010.
#'
#' @format A \code{\link{data.frame}} with 3457 observations and 5 columns:
#'
#' \itemize{
#' \item \code{phpws }: the proportion of households with piped water supply.
#' \item \code{mhdi  }: municipal human development index.
#' \item \code{incpc }: income per capita.
#' \item \code{region }: 0 for Southeast, 1 for Northeast.
#' \item \code{pop }: total population.
#' }
#'
#' @usage data(water, package = "uwquantreg")
#'
#' @references Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P., Ghitany, M. E., 2020.
#' The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
#' modelling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics} \bold{47} (6), 954--974.
#'
"water"


#-----------------------------------------------------------------------
#' @title Transport in Campus
#'
#' @description This data set is regarding the mode of transportation in a campus. A stratified sample of
#' 60 respondents was drawn for the purpose of oversampling people who sometimes bike to campus.
#'
#' @format A \code{\link{data.frame}} with 60 observations and 7
#'   columns:
#'
#' \itemize{
#' \item \code{ntrips }: number of trips to campus in the past four weeks.
#' \item \code{nbiked  }: the number of times the respondent biked to campus.
#' \item \code{status }: respondents status (student/faculty/staff);
#' \item \code{gender }: person's gender.
#' \item \code{parking }: duration of parking permit (6, 9, or 12-month.
#' \item \code{distance }: distance to campus.
#' \item \code{propbiked }: proportion of time a respondent biked to campus.
#' }
#'
#' @usage data(transport, package = "uwquantreg")
#'
#' @references Korosteleva, O., 2019. \bold{Advanced Regression Models with SAS and R}. Taylor & Francis Group.
#'
"transport"

#-----------------------------------------------------------------------
#' @title Mortality of Young Trees
#'
#' @description Study conducting by Parks and Recreation Department on mortality of young trees planted in parks.
#'
#' @format A \code{\link{data.frame}} with 26 observations and 7
#'   columns:
#'
#' \itemize{
#' \item \code{planted }: number of planted trees.
#' \item \code{survived  }: number of trees that survived for two years.
#' \item \code{pest }: frequency of pest control;
#' \item \code{fertilization }: frequency of soil fertilization.
#' \item \code{precip }: average annual precipitation (in inches).
#' \item \code{wind }: average annual wind speed (in miles per hour).
#' \item \code{prop }: proportion of survived trees.
#' }
#'
#' @usage data(trees, package = "uwquantreg")
#'
#' @references Korosteleva, O., 2019. \bold{Advanced Regression Models with SAS and R}. Taylor & Francis Group.
#'
"trees"
