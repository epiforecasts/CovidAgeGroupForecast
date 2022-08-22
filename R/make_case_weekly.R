library(data.table)
library(tidyverse)






make_case_weekly = function(){
  case_data = fread(file = 'data/nation_2022-05-10.csv')
  select_age_groups = c("00_04",  "05_09", "10_14",  "15_19", "20_24", "25_29", "30_34", "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69", "70_74", "75_79", "80_84", "85_89", "90+")
  case_data_ages = case_data[age %in% select_age_groups, ]
  case_data_ages[, lower_age := as.integer(substr(age, 1, 2))]
  case_data_ages[, upper_age := as.integer(substr(age, 4, 5))]
  
  case_data_ages[ is.na(upper_age), upper_age := 90]
  
  breaks = c(c(0),seq(9,69, 10), c(Inf))
  
  age_groups =  c('0-9',  '10-19',   '20-29',  '30-39',  '40-49',  '50-59',  '60-69', '70+')
  
  age_lookup = setNames(object=age_groups, nm = unique(cut(1:120, breaks)))

  case_data_ages[, age_group := cut(upper_age, breaks)]
  case_data_ages[, age_group := age_lookup[age_group]]
  
 
  
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
  
  
  list(case_data_age_groups, case_estimates)

}
