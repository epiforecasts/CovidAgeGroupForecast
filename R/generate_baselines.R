library(tidyverse)
library(data.table)

same_as_last_gen = function(inf_ests, summary_preds, generation_time=5){
  
  inf_ests[, date:=lubridate::ymd(date)]
  inf_ests_unshifted = copy(inf_ests)
  inf_ests[, date := date + 5]
  baseline_last_gen = unique(summary_preds[name=='forecast_gens',c('date', 'age_group', 'time_index', 'age_index', 'forecast_date')])
  baseline_last_gen[, name:= 'forecast_gens']
  baseline_last_gen[, run := 'baseline_last_gen']
  
  inf_est_thinned = inf_ests[,c('date', 'variable', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')] %>% rename(age_group = variable)

  baseline_last_gen = merge(baseline_last_gen, inf_est_thinned, by=c('date', 'age_group'), all.x = TRUE)
  baseline_last_gen = merge(baseline_last_gen, summary_preds[run=='1' & name=='forecast_gens', c('date', 'age_group', 'variable', 'value')], by=c('date', 'age_group'), all.x = TRUE)
  
  baseline_last_gen = baseline_last_gen[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')]
  colnames(baseline_last_gen) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')
 
  baseline_last_gen
  
}


same_diff_as_last_gen = function(inf_ests, summary_preds, generation_time=5){
  
  inf_ests[, date:=lubridate::ymd(date)]
  inf_ests[, date := date + 5]
  baseline_last_diff = unique(summary_preds[name=='forecast_gens',c('date', 'age_group', 'name', 'time_index', 'age_index', 'forecast_date')])
  baseline_last_diff[, run := 'baseline_last_diff']
  last_diff_all_groups = data.table()
  
  for(group in age_groups){
    values_only_group  = inf_ests[variable==group,c('q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')]
    last_diff_group = last_diff(values_only_group)
    last_diff_group[, age_group := group] 
    last_diff_group[, date := head(inf_ests[variable==group]$date,-5)]
    last_diff_all_groups = rbind(last_diff_all_groups, last_diff_group)
  }
  
  baseline_last_diff = merge(baseline_last_diff, last_diff_all_groups, by=c('date', 'age_group'), all.x = TRUE)
  baseline_last_diff = merge(baseline_last_diff, summary_preds[run=='1' & name=='forecast_gens', c('date', 'age_group', 'variable', 'value')], by=c('date', 'age_group'), all.x = TRUE)
  
  
  baseline_last_diff = baseline_last_diff[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')]
  colnames(baseline_last_diff) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')
  
  baseline_last_diff
  
}

same_as_last_value = function(inf_ests, summary_preds, smax=20, period=30, horizon=4){
  
  bl_values = inf_ests[date %in% unique(summary_preds$forecast_date), c('date', 'variable', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')] %>% rename(age_group = variable, forecast_date=date)
  bl_values[, forecast_date:=lubridate::ymd(forecast_date)]
  baseline_last_gen = unique(summary_preds[name=='forecast_gens'  & time_index > (smax+period),c('date', 'age_group', 'time_index', 'age_index', 'forecast_date')])
  baseline_last_gen[, name:= 'forecast_gens']
  baseline_last_gen[, run := 'baseline_last_val']
  
  baseline_last_gen = merge(baseline_last_gen, bl_values, by=c('forecast_date', 'age_group'), all.x = TRUE)
  baseline_last_gen = merge(baseline_last_gen, summary_preds[run=='1' & name=='forecast_gens' & time_index > (smax+period), c('forecast_date', 'date', 'age_group', 'variable', 'value')], by=c('forecast_date','date', 'age_group'), all.x = TRUE)
  
  baseline_last_gen = baseline_last_gen[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')]
  colnames(baseline_last_gen) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')
  
  unique(baseline_last_gen)
  
}

linex_from_last_values = function(inf_ests, summary_preds, smax=20, period=30, horizon=4, forecast_step=7){
  
  bl_values = inf_ests[date %in% c(unique(summary_preds$forecast_date), unique(summary_preds$forecast_date)-7), c('date', 'variable', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')] %>% rename(age_group = variable)
  
  for(d in as.character(unique(summary_preds$forecast_date))){
    d = as.Date(d)
    bl_values[date<=d & date >=(d-7), forecast_date:=d]
  }
  for(d in as.character(unique(summary_preds$forecast_date))){
    for(a in age_groups){
      d = as.Date(d)
      gradient = (bl_values[forecast_date==d & date==d & age_group==a,]$q50 - bl_values[forecast_date==d & date==(d-7) & age_group==a,]$q50)/7 
      
      base = bl_values[forecast_date==d & date==d & age_group==a,]
      
      for(ds in seq(forecast_step,horizon*forecast_step,forecast_step)){
        bl_values = 
          rbind(
          bl_values,
          data.table(
            date = d + ds,
            age_group = a,
            q5  = base$q5   + gradient * ds, 
            q10 = base$q10  + gradient * ds, 
            q25 = base$q25  + gradient * ds, 
            q50 = base$q50  + gradient * ds, 
            q75 = base$q75  + gradient * ds, 
            q90 = base$q90  + gradient * ds, 
            q95 = base$q95  + gradient * ds, 
            forecast_date = d
          )
              )
      }
      
    }
  }
  
  bl_values = bl_values[!(date < forecast_date)]

  baseline_last_gen = unique(summary_preds[name=='forecast_gens' & time_index > (smax+period),c('date', 'age_group', 'time_index', 'age_index', 'forecast_date', 'value')])
  baseline_last_gen[, name:= 'forecast_gens']
  baseline_last_gen[, run := 'baseline_linex_lv']
  baseline_last_gen[, variable := age_group]
  
  baseline_last_gen = merge(baseline_last_gen, bl_values, by=c('date', 'forecast_date', 'age_group'), all.x = TRUE)
  #baseline_last_gen = merge(baseline_last_gen, summary_preds[run=='1' & name=='forecast_gens' & time_index > (smax+period), c('date', 'age_group', 'variable', 'value')], by=c('date', 'age_group'), all.x = TRUE)
  
  baseline_last_gen = baseline_last_gen[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')]
  colnames(baseline_last_gen) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')
  
  unique(baseline_last_gen)
  
}



expex_from_last_values = function(inf_ests, summary_preds, smax=20, period=30, horizon=4, forecast_step=7){
  
  bl_values = inf_ests[date %in% c(unique(summary_preds$forecast_date), unique(summary_preds$forecast_date)-7), c('date', 'variable', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')] %>% rename(age_group = variable)
  
  for(d in as.character(unique(summary_preds$forecast_date))){
    d = as.Date(d)
    bl_values[date<=d & date >=(d-7), forecast_date:=d]
  }
  for(d in as.character(unique(summary_preds$forecast_date))){
    for(a in age_groups){
      d = as.Date(d)
      log_gradient = (log10(bl_values[forecast_date==d & date==d & age_group==a,]$q50) - log10(bl_values[forecast_date==d & date==(d-7) & age_group==a,]$q50))/7 
      
      base = bl_values[forecast_date==d & date==d & age_group==a,]
      
      for(ds in seq(forecast_step,horizon*forecast_step,forecast_step)){
        bl_values = 
          rbind(
            bl_values,
            data.table(
              date = d + ds,
              age_group = a,
              q5  = 10**(log10(base$q5)   + log_gradient * ds), 
              q10 = 10**(log10(base$q10)  + log_gradient * ds), 
              q25 = 10**(log10(base$q25)  + log_gradient * ds), 
              q50 = 10**(log10(base$q50)  + log_gradient * ds), 
              q75 = 10**(log10(base$q75)  + log_gradient * ds), 
              q90 = 10**(log10(base$q90)  + log_gradient * ds), 
              q95 = 10**(log10(base$q95)  + log_gradient * ds), 
              forecast_date = d
            )
          )
      }
      
    }
  }
  
  bl_values = bl_values[!(date < forecast_date)]
  
  baseline_last_gen = unique(summary_preds[name=='forecast_gens' & time_index > (smax+period),c('date', 'age_group', 'time_index', 'age_index', 'forecast_date', 'value')])
  baseline_last_gen[, name:= 'forecast_gens']
  baseline_last_gen[, run := 'baseline_expex_lv']
  baseline_last_gen[, variable := age_group]
  
  baseline_last_gen = merge(baseline_last_gen, bl_values, by=c('date', 'forecast_date', 'age_group'), all.x = TRUE)
  #baseline_last_gen = merge(baseline_last_gen, summary_preds[run=='1' & name=='forecast_gens' & time_index > (smax+period), c('date', 'age_group', 'variable', 'value')], by=c('date', 'age_group'), all.x = TRUE)
  
  baseline_last_gen = baseline_last_gen[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')]
  colnames(baseline_last_gen) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')
  
  unique(baseline_last_gen)
  
}





get_baselines = function(inf_ests, summary_preds, generation_time=5, smax=4, period=12, forecast_step=7){
  
  same_as_last_value = same_as_last_value(inf_ests, summary_preds, smax, period)
  linex_from_last_values = linex_from_last_values(inf_ests, summary_preds, smax, period, forecast_step)
  expex_from_last_values = expex_from_last_values(inf_ests, summary_preds, smax, period, forecast_step)
  
  baselines = rbind(same_as_last_value, linex_from_last_values, expex_from_last_values)
  
  baselines
}

last_diff = function(values){
  
  (tail(values$q50,-5) - head(values$q50,-5)) + tail(values, -5)
  
}

