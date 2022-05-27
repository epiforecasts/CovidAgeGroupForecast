library(data.table)
library(ggplot2)
library(tidyverse)

library(future.apply)
library(future.callr)
library(future)

library(cmdstanr)
library(cowplot)

source('R/fit_model.R')
source('R/plotting.R')
library(patchwork)


case_data = fread(file = '../nation_2022-05-10.csv')

select_age_groups = c("00_04",  "05_09", "10_14",  "15_19", "20_24", "25_29", "30_34", "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69", "70_74", "75_79", "80_84", "85_89", "90+")


case_data_ages = case_data[age %in% select_age_groups, ]



case_data_ages[, lower_age := as.integer(substr(age, 1, 2))]
case_data_ages[, upper_age := as.integer(substr(age, 4, 5))]

case_data_ages[ is.na(upper_age), upper_age := 90]

breaks = c(c(0),seq(9,69, 10), c(Inf))
case_data_ages[, age_group := cut(upper_age, breaks)]

age_groups = unique(cut(1:120, breaks))

select_dates = seq(as.Date('2020-01-31'), as.Date('2022-01-01'), 7)

case_data_ages = case_data_ages[ date %in% select_dates, ]


case_data_ages[, rolling_cases := sum(rollingSum), by = c('age_group', 'date') ]

case_data_age_groups = unique(case_data_ages[,c('areaCode', 'areaName', 'areaType', 'date', 'age_group', 'rolling_cases')])

case_estimates = data.table(
  date = as.Date(case_data_age_groups$date), 
  variable = case_data_age_groups$age_group, 
  q5  = case_data_age_groups$rolling_cases, 
  q10 = case_data_age_groups$rolling_cases, 
  q25 = case_data_age_groups$rolling_cases, 
  q50 = case_data_age_groups$rolling_cases, 
  q75 = case_data_age_groups$rolling_cases, 
  q90 = case_data_age_groups$rolling_cases, 
  q95 = case_data_age_groups$rolling_cases)

ggplot(case_data_age_groups) +
  geom_point(aes(x=date, y=rolling_cases, color=age_group))




case_matrix_mean  = dcast(case_data_age_groups, value.var = 'rolling_cases', date ~ age_group)

case_matrix_mean[, date := lubridate::ymd(date)]

cms_cases = readRDS('cms_cases.rds')

age_mod_cases = cmdstan_model('stan/multi-option-contact-model-cases.stan')
period=8*7 # days 
smax=4 # weeks

sr_dates = readRDS('sr_dates.rds')

dates = seq(as.Date("2020-11-27"), as.Date('2021-12-02'), 14)

plan(callr, workers = future::availableCores()-1)
all_est = list()
for(r in 1:4){ 
  print(paste0('running model ', r))
  est <- future_lapply(
    dates, fit_NGM_model_for_date_range_cases,
    age_mod = age_mod_cases,
    period = period+smax*7, 
    smax = smax,
    inf_matrix_mean = case_matrix_mean, 
    inf_matrix_sd = NULL, 
    anb_matrix_mean = NULL, 
    anb_matrix_sd = NULL, 
    cms = cms_cases, 
    pops = readRDS('population_cases.rds'),
    runindex = r, 
    quantiles=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
    contact_option=r,
    sigma_option=1,
    contact_delay = 5, 
    forecast_horizon = 4,
    future.seed=TRUE
  )
  
  all_est = append(all_est, est)
  
}



# combine the outputs from all dates into one container
summary_pars = data.table()
for(i in 1:length(all_est)){
  summary_pars = rbind(summary_pars, data.table(all_est[[i]]$summary_pars))
}

summary_preds = data.table()
for(i in 1:length(all_est)){
  summary_preds = rbind(summary_preds, data.table(all_est[[i]]$summary_preds))
}

# set dates and age groups based on time age index
summary_preds[, date := forecast_date - (period+smax*7) + (time_index-1)*7]
summary_preds[, age_group := age_groups[age_index]]
# merge with time series of median infection estimates 
summary_preds = merge(summary_preds, case_data_age_groups, by.x = c('date', 'age_group'), by.y = c('date', 'age_group'), all.x = T, all.y = F )
summary_preds[, value := rolling_cases]
summary_preds[, variable := age_group]


# calculate baselines using get_baselines function 
baselines = get_baselines(case_estimates, summary_preds, smax=4, period=8, forecast_step = 7)
# merge baselines with summary outputs
summary_preds = rbind(summary_preds[, names(baselines), with=FALSE], baselines)


summary_conts = data.table()
for(i in 1:length(all_est)){
  summary_conts = rbind(summary_conts, data.table(all_est[[i]]$summary_conts))
}
ggplot(summary_preds) +
  geom_line(data=summary_preds[ name=='forecast_gens' & date>forecast_date - 8*7], aes(x=date, y=value))+
  geom_point(data=summary_preds[ name=='forecast_gens' & date > forecast_date], aes(x=date, y=`50%`, color=run), alpha=0.8, size=0.5)+ 
  geom_ribbon(data=summary_preds[name=='forecast_gens'], aes(x=date, ymin=`5%`, ymax=`95%`, fill=run), alpha=0.2)+
  geom_line(data = summary_preds[time_index>smax & name=='next_gens'], 
            aes(x=date, y=`50%`, color=run), alpha=0.8)+
  
  geom_ribbon(data = summary_preds[time_index>smax & name=='next_gens'], 
              aes(x=date, ymin=`5%`, ymax=`95%`, fill=run), alpha=0.3)+
  facet_grid(age_group~forecast_date, scales = 'free')+
  theme(
    legend.position = 'bottom'
  )

saveRDS(summary_pars,  'outputs/summary_pars_cases.rds')
saveRDS(summary_preds, 'outputs/summary_preds_cases.rds')
saveRDS(summary_conts, 'outputs/summary_conts_cases.rds')

summary_scores = score_forecasts(summary_preds[date > as.Date('2020-10-01')], pandemic_periods, suffix = '_cases')

plot_parameters(summary_pars, d, '_cases')
