


same_as_last_gen = function(inf_ests, summary_preds, generation_time=5){
  
  baseline_last_gen_2gens = data.table()
  
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
  
  baseline_last_gen[, generation := 1]
  
  baseline_last_gen_2gens = rbind(baseline_last_gen_2gens, baseline_last_gen)
  
  inf_ests[, date := date + 5]
  baseline_last_gen = unique(summary_preds[name=='forecast_gens',c('date', 'age_group', 'time_index', 'age_index', 'forecast_date')])
  baseline_last_gen[, name:= 'forecast_gens']
  baseline_last_gen[, run := 'baseline_last_gen']
  
  inf_est_thinned = inf_ests[,c('date', 'variable', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')] %>% rename(age_group = variable)
  
  baseline_last_gen = merge(baseline_last_gen, inf_est_thinned, by=c('date', 'age_group'), all.x = TRUE)
  baseline_last_gen = merge(baseline_last_gen, summary_preds[run=='1' & name=='forecast_gens', c('date', 'age_group', 'variable', 'value')], by=c('date', 'age_group'), all.x = TRUE)
  
  baseline_last_gen = baseline_last_gen[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')]
  colnames(baseline_last_gen) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')
  
  baseline_last_gen[, generation := 2]
  
  baseline_last_gen_2gens = rbind(baseline_last_gen_2gens, baseline_last_gen)
  
  baseline_last_gen_2gens
  
  
}


same_diff_as_last_gen = function(inf_ests, summary_preds, generation_time=5){
  
  baseline_last_diff_2gens = data.table()
  
  inf_ests[, date:=lubridate::ymd(date)]
  inf_ests[, date := date + 5]
  baseline_last_diff = unique(summary_preds[name=='forecast_gens',c('date', 'age_group', 'name', 'time_index', 'age_index', 'forecast_date')])
  baseline_last_diff[, run := 'baseline_last_diff']
  last_diff_all_groups = data.table()
  
  for(group in age_groups){
    values_only_group  = inf_ests[variable==group,c('q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')]
    last_diff_group = last_diff(values_only_group)
    last_diff_group[, age_group := group] 
    last_diff_group[, date := tail(inf_ests[variable==group]$date,-5)]
    last_diff_all_groups = rbind(last_diff_all_groups, last_diff_group)
  }
  
  baseline_last_diff = merge(baseline_last_diff, last_diff_all_groups, by=c('date', 'age_group'), all.x = TRUE)
  baseline_last_diff = merge(baseline_last_diff, summary_preds[run=='1' & name=='forecast_gens', c('date', 'age_group', 'variable', 'value')], by=c('date', 'age_group'), all.x = TRUE)
  
  
  baseline_last_diff = baseline_last_diff[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')]
  colnames(baseline_last_diff) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')
  baseline_last_diff[, generation := 1]
  
  inf_ests[, date := date + 5]
  baseline_last_diff = unique(summary_preds[name=='forecast_gens',c('date', 'age_group', 'name', 'time_index', 'age_index', 'forecast_date')])
  baseline_last_diff[, run := 'baseline_last_diff']
  last_diff_all_groups = data.table()
  
  for(group in age_groups){
    values_only_group  = inf_ests[variable==group,c('q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95')]
    last_diff_group = last_diff(values_only_group, gen=2)
    last_diff_group[, age_group := group] 
    last_diff_group[, date := tail(inf_ests[variable==group]$date,-5)]
    last_diff_all_groups = rbind(last_diff_all_groups, last_diff_group)
  }
  
  baseline_last_diff = merge(baseline_last_diff, last_diff_all_groups, by=c('date', 'age_group'), all.x = TRUE)
  baseline_last_diff = merge(baseline_last_diff, summary_preds[run=='1' & name=='forecast_gens', c('date', 'age_group', 'variable', 'value')], by=c('date', 'age_group'), all.x = TRUE)
  
  
  baseline_last_diff = baseline_last_diff[, c('date', 'age_group', 'name', 'q5', 'q10', 'q25', 'q50', 'q75', 'q90', 'q95', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')]
  colnames(baseline_last_diff) = c('date', 'age_group', 'name', '5%', '10%', '25%', '50%', '75%', '90%', '95%', 'time_index', 'age_index', 'forecast_date', 'run', 'variable', 'value')
  baseline_last_diff[, generation := 2]
  
  baseline_last_diff_2gens = rbind(baseline_last_diff_2gens, baseline_last_diff)
  
  baseline_last_diff_2gens
}

get_baselines = function(inf_ests, summary_preds, generation_time=5){
  
  last_gen_baseline = same_as_last_gen(inf_ests, summary_preds)
  last_diff_baseline = same_diff_as_last_gen(inf_ests, summary_preds)
  
  baselines = rbind(last_gen_baseline, last_diff_baseline)
  
  baselines
}

last_diff = function(values, gen=1){
  
  (tail(values$q50,-5) - head(values$q50,-5)) + tail(values, -5*gen)
  
}

