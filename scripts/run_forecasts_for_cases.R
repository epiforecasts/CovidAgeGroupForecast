library(data.table)
library(ggplot2)
library(tidyverse)

library(future.apply)
library(future.callr)
library(future)

library(cmdstanr)
library(cowplot)

source('R/fit_model_cases.R')
source('R/make_case_weekly.R')
source('R/get_polymod_cms.R')
source('R/generate_baselines.R')
library(patchwork)

weekly_cases = make_case_weekly()

case_data_age_groups = weekly_cases[[1]]
case_estimates = weekly_cases[[2]]

breaks = c(c(0),seq(9,69, 10), c(Inf))

age_groups = c('0-9',  '10-19',   '20-29',  '30-39',  '40-49',  '50-59',  '60-69', '70+')

cas_traj = ggplot(case_data_age_groups[date>lubridate::ymd(20200730)]) + 
  #geom_ribbon(aes(x=date, ymin=q10, ymax=q90))+
  geom_point(aes(x=date, y=rolling_cases, color=age_group), size=1, stroke=0)+
  facet_wrap(~age_group, nrow=1)+
  scale_y_continuous(name='Cases', limits = c(0,100000))+
  scale_x_date( name='Date', date_labels = '%b-%y')+
  scale_color_manual(name='Age group', values = as.vector(light(9)))+
  theme_minimal_hgrid()+
  ggtitle('C')+
  theme(axis.text.x = element_text(angle = 45, size=10))


saveRDS(cas_traj, 'case_traj.rds')


case_matrix_mean  = dcast(case_data_age_groups, value.var = 'rolling_cases', date ~ age_group)

case_matrix_mean[, date := lubridate::ymd(date)]

cms_cases = readRDS('data/cms_cases.rds')

age_mod_cases = cmdstan_model('stan/multi-option-contact-model-cases.stan')
period=8*7 # days 
smax=4 # weeks

sr_dates = readRDS('data/sr_dates.rds')

dates = seq(as.Date("2020-11-27")-(4*7), as.Date('2021-12-02'), 14)

plan(callr, workers = future::availableCores()-1)
all_est = list()
for(r in c(1,4,5)){ 
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
    pops = readRDS('data/population_cases.rds'),
    runindex = r, 
    quantiles=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
    contact_option=r,
    sigma_option=1,
    contact_delay = 5, 
    forecast_horizon = 4,

    ad=0.99,
    future.seed=TRUE
    
  )
  
  all_est = append(all_est, est)
  
}

cms_pmd = get_polymod_cms(breaks)

est <- future_lapply(
  dates, fit_NGM_model_for_date_range_cases,
  age_mod = age_mod_cases,
  period = period+smax*7, 
  smax = smax,
  inf_matrix_mean = case_matrix_mean, 
  inf_matrix_sd = NULL, 
  anb_matrix_mean = NULL, 
  anb_matrix_sd = NULL, 
  cms = cms_pmd, 
  pops = readRDS('data/population_cases.rds'),
  runindex = 6, 
  quantiles=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
  contact_option=1,
  sigma_option=1,
  contact_delay = 5, 
  forecast_horizon = 4,
  future.seed=TRUE,
  ad=0.99

)

all_est = append(all_est, est)


# combine the outputs from all dates into one container
summary_pars = data.table()
for(i in 1:length(all_est)){
  summary_pars = rbind(summary_pars, data.table(all_est[[i]]$summary_pars))
}



summary_diags = data.table()
for(i in 1:length(all_est)){
  summary_diags = rbind(summary_diags, data.table(all_est[[i]]$diagnostics))
}


summary_preds = data.table()
for(i in 1:length(all_est)){
  summary_preds = rbind(summary_preds, data.table(all_est[[i]]$summary_preds))
}

samples_preds = data.table()
for(i in 1:length(all_est)){
  samples_preds = rbind(samples_preds, data.table(all_est[[i]]$samples_preds))
}


samples_preds[, date := forecast_date - (period+smax*7) + (time_index-1)*7]
samples_preds[, age_group := age_groups[age_index]]
# merge with time series of median infection estimates 
samples_preds = merge(samples_preds, case_data_age_groups[,  true_value := rolling_cases][,-c('rolling_cases', 'areaCode', 'areaName', 'areaType')], by.x = c('date', 'age_group'), by.y = c('date', 'age_group'), all.x = T, all.y = F )

# calculate baseline

# set dates and age groups based on time age index
summary_preds[, date := forecast_date - (period+smax*7) + (time_index-1)*7]
summary_preds[, age_group := age_groups[age_index]]
# merge with time series of median infection estimates 
summary_preds = merge(summary_preds, case_data_age_groups, by.x = c('date', 'age_group'), by.y = c('date', 'age_group'), all.x = T, all.y = F )
summary_preds[, value := rolling_cases]
summary_preds[, variable := age_group]


samples_preds[, prediction:=value][, -c('value')]


# calculate baselines using get_baselines function 
baselines = get_baselines(case_estimates, summary_preds, smax=4, period=8, forecast_step = 7)
# merge baselines with summary outputs
summary_preds = rbind(summary_preds[, names(baselines), with=FALSE], baselines)

baselines_samples = baselines[, prediction := `50%`][, -c("5%"    ,    "10%"   ,     "25%"    ,    "50%"    ,    "75%"    ,    "90%"  , "95%")]

baselines_samples_mult = baselines_samples %>% mutate(sample=1)

for(samp in 2:100){
  baselines_samples_mult = rbind(baselines_samples_mult , baselines_samples %>% mutate(sample=samp))
}

baselines_samples_mult = baselines_samples_mult[, true_value := value][,-c('value')]

samples_preds = rbind(samples_preds[,-c('value')], baselines_samples_mult[,-c('variable')])

summary_conts = data.table()
for(i in 1:length(all_est)){
  summary_conts = rbind(summary_conts, data.table(all_est[[i]]$summary_conts))
}

saveRDS(summary_pars,  'outputs/summary_pars_cases.rds')
saveRDS(summary_preds, 'outputs/summary_preds_cases.rds')
saveRDS(summary_conts, 'outputs/summary_conts_cases.rds')
saveRDS(summary_diags, 'outputs/summary_diags_cases.rds')
saveRDS(samples_preds, 'outputs/samples_preds_cases.rds')



