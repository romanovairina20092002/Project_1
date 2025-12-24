# Time Series Forecasting of U.S. Personal Consumption Expenditures (PCE)

## Project Overview
This project applies time series forecasting techniques to predict U.S. Personal Consumption Expenditures (PCE) using monthly data.

## Requirements
- R (version 4.0 or later)
- The following R packages:
  - VIM
  - tseries
  - imputeTS
  - tsoutliers
  - forecast
  - ggplot2
  - zoo
  - knitr

## Dataset
The dataset covers the period from 1959 to 2024 and is provided in the file `PCE.csv`.

## Variables
- `PCE` — Personal Consumption Expenditures for the given date  
- `observation_date` — Date of observation (monthly frequency)

## Methodology

1. Data exploration and preprocessing
   - Loaded and inspected monthly PCE data
   - Visualised missing values
   - Converted the series into a monthly time series starting in January 1959

2. Missing data imputation
   - Compared multiple imputation methods on the training period:
     - Linear interpolation
     - Moving average
     - Kalman filter (ARIMA-based)
     - Spline interpolation
   - Artificial missing values were created to evaluate imputation performance
   - RMSE was used to compare methods
   - Linear interpolation was selected as the best approach (lowest RMSE)

3. Data transformation and outlier treatment
   - Applied a logarithmic transformation to stabilise variance
   - Detected and corrected outliers using `tsoutliers`
   - Outliers were replaced mainly during the COVID-19 disruption period (Mar 2020 – Feb 2021)
   - Used STL decomposition to inspect trend and seasonality

4. Train–test split
   - Training set: Jan 1959 – Dec 2023
   - Test set: Jan 2024 – Dec 2024

5. Forecasting models
   The following models were fitted and evaluated:
   - Drift model (random walk with drift)
   - Exponential smoothing (ETS(A,A,N))
   - ARIMA (automatic model selection)

6. Model evaluation
   - Residual diagnostics and Box–Ljung tests were conducted for all models
   - Forecast accuracy was evaluated on the test set using RMSE 
   - The best-performing model was ARIMA(0,2,2)(2,0,1)[12] 

## Results
- Best imputation method: Linear interpolation 
- Best forecasting model: ARIMA(0,2,2)(2,0,1)[12] 

## Forecast
- The final ARIMA model was trained on the full cleaned series
- PCE values were forecasted for the next 12 months
- Forecast results are visualised alongside the historical series

## License
This project is released under the MIT License. See `LICENSE`
