library(uwquantreg)

rm(list = ls())
# devtools::load_all()

data("transport")

fit <- zouwquantreg(
  formula = propbiked ~ status + parking,
  tau = 0.5,
  inflation = 0,
  data = transport,
  formula.nu = propbiked ~ gender + distance
)
fit
summary(fit)
coef(fit)
vcov(fit)
confint(fit, what = "both")
confint(fit, what = "quant")
fitted(fit, type = "all")

## Check residuals
head(residuals(fit))

## Check predict
head(predict(fit))

