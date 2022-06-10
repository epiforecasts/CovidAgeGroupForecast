


make_infs_weekly = function(){
  
  samples = readRDS('../samples_age_ab.rds')
  estimates = read.csv('../estimates_age_ab.csv')
  
  samples = merge(samples, unique(estimates[, c('variable', 'population')]), by=c('variable'), all.x = T)
  
  inf_samples = samples[name == 'infections',]
  
  
  
  inf_samples = inf_samples[order(date)]
  inf_samples[, date:=as.Date(date)]
  
  
  
  inf_samples[, rollingsum:=roll::roll_sum(value, width = 7) * population, by=c('sample','variable')]
  
  inf_samples_weekly = inf_samples[date %in% select_dates,]
  
  
  
  inf_estimates_gen = 
    unique(data.table(
      name = 'infections',
      date = inf_samples_weekly$date, 
      variable = inf_samples_weekly$variable
    ))
  
  
  
  for(d in tail(unique(as.character(inf_samples_weekly$date)),-1)){
    for(v in unique(inf_samples_weekly$variable)){
      
      inf_estimates_gen[date == as.Date(d) & variable == v, mean := mean(inf_samples_weekly[date == as.Date(d) & variable == v,]$rollingsum, na.rm = T)]
      inf_estimates_gen[date == as.Date(d) & variable == v, sd := sd(inf_samples_weekly[date == as.Date(d) & variable == v,]$rollingsum)]
      
      ecdf_infs = ecdf(inf_samples_weekly[date == as.Date(d) & variable == v,]$rollingsum)
      inf_estimates_gen[date == as.Date(d) & variable == v, q5  := as.numeric(quantile(ecdf_infs ,0.05))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q10 := as.numeric(quantile(ecdf_infs ,0.10))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q15 := as.numeric(quantile(ecdf_infs ,0.15))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q20 := as.numeric(quantile(ecdf_infs ,0.20))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q25 := as.numeric(quantile(ecdf_infs ,0.25))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q30 := as.numeric(quantile(ecdf_infs ,0.30))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q35 := as.numeric(quantile(ecdf_infs ,0.35))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q40 := as.numeric(quantile(ecdf_infs ,0.40))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q45 := as.numeric(quantile(ecdf_infs ,0.45))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q50 := as.numeric(quantile(ecdf_infs ,0.50))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q55 := as.numeric(quantile(ecdf_infs ,0.55))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q60 := as.numeric(quantile(ecdf_infs ,0.60))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q65 := as.numeric(quantile(ecdf_infs ,0.65))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q70 := as.numeric(quantile(ecdf_infs ,0.70))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q75 := as.numeric(quantile(ecdf_infs ,0.75))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q80 := as.numeric(quantile(ecdf_infs ,0.80))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q85 := as.numeric(quantile(ecdf_infs ,0.85))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q90 := as.numeric(quantile(ecdf_infs ,0.90))]
      inf_estimates_gen[date == as.Date(d) & variable == v, q95 := as.numeric(quantile(ecdf_infs ,0.95))]
      
      inf_estimates_gen[date == as.Date(d) & variable == v, n_index := unique(inf_samples_weekly[date == as.Date(d) & variable == v,]$n_index)]
      inf_estimates_gen[date == as.Date(d) & variable == v, t_index := unique(inf_samples_weekly[date == as.Date(d) & variable == v,]$t_index)]
      inf_estimates_gen[date == as.Date(d) & variable == v, p_index := unique(inf_samples_weekly[date == as.Date(d) & variable == v,]$p_index)]
      
    }
  }
  
  inf_estimates_gen
}



