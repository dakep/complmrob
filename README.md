# `complmrob` R Package

This R package implements robust regression methods for compositional data.
The distribution of the estimates can be approximated with various bootstrap methods.
These bootstrap methods are available for the compositional as well as for standard robust regression estimates.
This allows for direct comparison between them.

## Usage

To fit a robust regression model to compositional predictors, use the function `complmrob`.
It takes two inputs: the formula for the model (`formula`) and the data set (`data`).

As an example, consider the linear relationship between violent crime rates and average life expectancy in 50 US states.
This data is readily available in R:
```{r}
crimes <- data.frame(lifeExp = state.x77[, "Life Exp"], USArrests[ , c("Murder", "Assault", "Rape")])
```

If it would be assumed that the crime rates are truly compositional in nature (which they are not!), `complmrob` can be used to robustly fit a linear regression model:
```{r}
lifeexp_fit <- complmrob(lifeExp ~ ., data = crimes)
```

The summary for the fit, shows the usual statistics and the confidence interval for each parameter:
```{r}
summary(lifeexp_fit)
```

The inference statistics of the parameters are based on approximations that may not be reliable for compositional data.
As an alternative, bootstrap can be used:
```{r}
bs_lifexp_fit <- bootcoefs(lifeexp_fit)
summary(bs_lifexp_fit)
plot(bs_lifexp_fit)
```
