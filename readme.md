# Age specific forecasts of SARS-CoV-2 infection using contact data

## Introduction 

This is a framework for incorporating social contact data, measured weekly by the CoMix survey, into a semi-mechanistic forecasting model to allow interaction between age-groups as we forecast incidence of infection of SARS-CoV-2 in England. First we estimated daily infection incidence and antibody prevalence by fitting models to COVID-19 Infection Survey data. We used these estimates to fit key epidemiological parameters to augment weekly contact matrices creating a time-varying ‘next generation matrix’, which describes transmission rates within and between age-groups. 

The model is fit in 'stan' using the 'cmdstanr' package. 


## Pipeline

To run forecasts on historical dates, execute the 'historical_estimates.r' script and all models will be run for a series of forecast dates. To score the forecasts use the 'score_forecasts()' function which relys on the 'summary_preds' output generated from 'fit_model()' and combined for multiple runs within 'historical_estimates.r'.


## Methods




## Example results