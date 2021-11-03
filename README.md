# A Comparative Analysis of S&P 500 and Industrial Production

In this project, I analyzed two time series data, one of which comes from S&P 500 returns and one of which comes from Industrial Production Index (IPI).

This project is part of my "Advanced Time Series Analysis" course in MSc. Statistics at KU Leuven (December 2019).

# Table of Contents
- The Data
- Univariate Analysis (Model Selection)
- Diebold Mariano Test
- Forecasts
- Multivariate Analysis (Cointegration)
- Dynamic Models & Granger Causality
- VAR Model
- Comparing Predictions
- Conclusions & Comments

## The Data
The data consists of two monthly time series in USA, from January 1970 to December 2018.
1. *trsp500*: S&P 500 Total Return; Monthly Dividend Reinvest (EOP – Equity Office Properties) [^1].
2. *indpro*: Industrial Production Index: 2012 = 100 (SA – Seasonally Adjusted) [^2].

[^1]: http://www.economagic.com/em-cgi/data.exe/fedstl/trsp500
[^2]: http://www.economagic.com/em-cgi/data.exe/fedstl/indpro

- There are high autocorrelation and no seasonality effect for both series – checked with acf and monthplot. Series are persistent.
- Log-difference transformation is applied to obtain stationary time series. Augmented Dickey Fuller test is used to ensure that
with both «drift» and «trend» options.

## Univariate Analysis (Model Selection)
- In order to select a suitable model for the time series, correlograms of the transformed series are examined.
  - *trsp500*: No significant autocorrelation or partial autocorrelation occur at any lags.
  - *indpro*: Seasonal difference might be needed at lag 2 (24 months).
- After that, 2 parsimonious models are created for each time series; one of which is found with auto.arima function.
- All models are validated with Box-Ljung test on residuals. Also, correlograms of model residuals show no autocorrelation.
- For both series, different models are better in different measures. In order to decide which model is significantly better, DieboldMariano test (with squared forecast error) is used.

## Diebold-Mariano Test
- Diebold-Mariano test compares the h-step ahead forecast errors of two models. With this model, it is possible to test which of
the two models perform better.
- Paired t-test is applied. However, Newey-West correction is used on the standard error.
- According to the Diebold-Mariano tests, there are no significant differences between the models selected (α = 0.05).
- For parsimony, the models with lower error means are selected for further analyses. They are:
  - trsp500 → ARMA(1,1)
  - indpro → SARMA(1,1)(2,1)[12]

## Forecasts
- With the selected models, the values of the next 12 months (from January 2019 to December 2019) are predicted.
- However, since the predictions are for the transformed values, back-transformation is applied to get real forecasted values.

## Multivariate Analysis (Cointegration)Multivariate Analysis (Cointegration)
- Cointegration between two time series should be checked. If cointegration exists, long run equilibrium relationship between the two variables can be estimated by OLS.
- Engle-Granger approach is used to detect cointegration. Variables are regressed on each other.
- ECM (Error Correcting Model) is applied to describe short-run dynamics between two time series.
- We cannot confidently conclude that cointegration occurs.

## Dynamic Models & Granger Causality
- Using dynamic models, one stationary time series is regressed on the other one. Lagged variables are needed, because previous models were not adequate to explain the relationship between two time series.
- Distributed lag models (DLM) and autoregressive distributed lag models (ADLM) are tried with different orders.
- Although DL models seem to perform better than ADL models, the results are inconclusive. Also, DLM and ADLM are not very suitable for predictions.
- The condition of Granger-Causality is tested with lag one, since it seems to be the most adequate model tested. Granger-Causality does not occur (α = 0.05).

## VAR Model
- Two time series are combined, and a 2-dimensional VAR model is created.
- The most suitable model is found to be VAR(3).
- Estimations for production index are more accurate than those for dividend reinvestment, due to lower p-value and higher R2.
- Again, real predictions are obtained by back-transforming the log-differenced predictions.

## Comparing Predictions
- Two sets of predictions are obtained for two time series: One from univariate analysis, one from multivariate analysis.
- Luckily, real values of the series for 2019 exist until November – which may not always be the case.
- Predictions are compared, based on absolute errors and squared errors.
  -  Univariate analysis give better results for trsp500.
  - Multivariate analysis give better results for indpro.

## Conclusions & Comments
- Both time series have an increasing trend. Trend for trsp500 seems to be exponential, whereas the trend for indpro seems to be
linear.
- Sharp declines in both time series (for instance, in 2009) coincides with the economic crisis in USA.
- In order to achieve more accurate results without violating the assumptions, time series should be stationary. Stationartiy is
satisfied with log-difference transformation.
- In order to avoid computationally exhaustive analysis, the number of models tried are kept low.
- Although the plots of the time series indicate similar patterns, cointegration is not observed.
- Even though they are valid according to Box-Ljung test, error correcting model and dynamic models are inadequate to explain
the variability.
- There may be other external regressors that affects both time series, like housing prices, interest rates, GDP per capita etc.
However, a causal relationship between two time series cannot be proven. Granger causality also did not occur.
- Since the predictions from multivariate analysis perform better than those from univariate analysis for indpro, it can be said that
industrial production is affected by external factors more than dividend reinvestment is.
