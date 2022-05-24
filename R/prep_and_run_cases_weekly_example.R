library(data.table)
library(ggplot2)

case_data = fread(file = '../nation_2022-05-10.csv')

select_age_groups = c("00_04",  "05_09", "10_14",  "15_19", "20_24", "25_29", "30_34", "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69", "70_74", "75_79", "80_84", "85_89", "90+")


case_data_ages = case_data[age %in% select_age_groups, ]



case_data_ages[, lower_age := as.integer(substr(age, 1, 2))]
case_data_ages[, upper_age := as.integer(substr(age, 4, 5))]

case_data_ages[ is.na(upper_age), upper_age := 90]

breaks = c(c(0),seq(9,69, 10), c(Inf))
case_data_ages[, age_group := cut(upper_age, breaks)]

age_groups = unique(cut(1:120, breaks))

select_dates = seq(as.Date(dates[[1]])-14, as.Date(dates[[length(dates)]])+14, 7)

case_data_ages = case_data_ages[ date %in% select_dates, ]


case_data_ages[, rolling_cases := sum(rollingSum), by = c('age_group', 'date') ]

case_data_ages


case_data_age_groups = unique(case_data_ages[,c('areaCode', 'areaName', 'areaType', 'date', 'age_group', 'rolling_cases')])


ggplot(case_data_age_groups) +
  geom_point(aes(x=date, y=rolling_cases, color=age_group))




case_matrix_mean  = dcast(case_data_age_groups, value.var = 'rolling_cases', date ~ age_group)

case_matrix_mean[, date := lubridate::ymd(date)]

cms_cases = readRDS('cms_cases.rds')

age_mod_cases = cmdstan_model('stan/multi-option-contact-model-cases.stan')
period=6*7 # days 
smax=2 # weeks

fit = fit_NGM_model_for_date_range_cases(
  end_date = dates[[10]],
  age_model = age_mod_cases,
  period = period+smax*7, 
  smax = smax,
  inf_matrix_mean = case_matrix_mean, 
  inf_matrix_sd = NULL, 
  anb_matrix_mean = NULL, 
  anb_matrix_sd = NULL, 
  cms = cms_cases, 
  pops = readRDS('population_cases.rds'),
  
  runindex = 1, 
  quantiles=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
  contact_option = 4,
  sigma_option = 1, 
  contact_delay = 5, 
  forecast_horizon = 2
)


summary_preds = data.table(fit$summary_preds)
summary_preds[, age_group := age_groups[age_index]]
summary_preds[, date := forecast_date - period+smax*7 + time_index*7]
#summary_preds = merge(summary_preds, actualest_inc_long, by=c('date', 'age_group'), all.x = TRUE)


ggplot()+
  geom_point(data=summary_preds[time_index>=6 & name=='forecast_gens'], aes(x=date, y=`50%`, color=age_group), alpha=0.8)+ 
  geom_ribbon(data=summary_preds[name=='forecast_gens'], aes(x=date, ymin=`5%`, ymax=`95%`, fill=age_group), alpha=0.2)+
  geom_line(data = summary_preds[time_index>2 & time_index<6 & name=='next_gens'], 
            aes(x=date, y=`50%`, color=age_group), alpha=0.8)

  
summary_preds[name=='forecast_gens']
