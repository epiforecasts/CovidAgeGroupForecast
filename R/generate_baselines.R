
samples = read.csv(file = '../samples_age_ab.csv')

same_as_last_gen = function(inf_estimates, summary_preds, generation_time=5){
  
  inf_estimates[, date:=lubridate::ymd(date)]
  inf_estimates[, date := date + 5]
  baseline_last_gen = unique(summary_preds[,c('date', 'age_group', 'time_index', 'age_index', 'forecast_date')])
  baseline_last_gen[, name:= 'forecast_gens']
  baseline_last_gen[, run := 'baseline_last_gen']
  
  
  inf_est_thinned = inf_estimates[,c('date', 'variable', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')] %>% rename(age_group = variable)

  baseline_last_gen = merge(baseline_last_gen, inf_est_thinned, by=c('date', 'age_group'), all.x = TRUE)
  
 
  baseline_last_gen = baseline_last_gen[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run')]

  colnames(baseline_last_gen) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run')
  
}


same_diff_as_last_gen = function(inf_estimates, summary_preds, generation_time=5){
  
  inf_estimates[, date:=lubridate::ymd(date)]
  inf_estimates[, date := date + 5]
  inf_estimates[, baseline = ]
  baseline_last_gen = unique(summary_preds[,c('date', 'age_group', 'time_index', 'age_index', 'forecast_date')])
  baseline_last_gen[, name:= 'forecast_gens']
  baseline_last_gen[, run := 'baseline_last_gen']
  
  
  inf_est_thinned = inf_estimates[,c('date', 'variable', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')] %>% rename(age_group = variable)
  
  baseline_last_gen = merge(baseline_last_gen, inf_est_thinned, by=c('date', 'age_group'), all.x = TRUE)
  
  
  baseline_last_gen = baseline_last_gen[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run')]
  
  colnames(baseline_last_gen) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run')
  
}

age_group = '2-10'

values_only  = inf_estimates[variable==age_group,c('q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')]

values_only

last_diff = function(values){
  
  (tail(values$q50,-5) - head(values$q50,-5)) + tail(values, -5)
  
}


