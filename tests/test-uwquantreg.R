library(uwquantreg)

rm(list = ls())
# devtools::load_all()

data(sdac, package = "simplexreg")

## check fit and methods
fit <- uwquantreg(
  formula = rcd ~ gender + chemo + ageadj,
  data = sdac,
  tau = 0.5,
  link = "logit"
)
fit
summary(fit)
coef(fit)
vcov(fit)
confint(fit)


## Check predict
head(predict(fit, type = "mu", interval = "none"))
head(predict(fit, type = "mu", interval = "confidence"))

## Check residuals
head(residuals(fit, type = "cox-snell"))
head(residuals(fit, type = "quantile"))

## Check hnp
out <- hnp(fit, resid.type = "quantile", halfnormal = TRUE)
# out <- hnp(fit, resid.type = "quantile", halfnormal = FALSE)
# out <- hnp(fit, resid.type = "cox-snell", halfnormal = TRUE)
# out <- hnp(fit, resid.type = "cox-snell", halfnormal = FALSE)



