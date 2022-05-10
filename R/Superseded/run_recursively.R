library(data.table)
library(rstan)
library(rstanarm)
library(ggplot2)
library(tidyverse)

library(future.apply)
library(future.callr)
library(future)

library(cmdstanr)

age_mod = cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat.stan')

# sample from stan model

age_mods = list(cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat.stan'), cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat_fitmeans.stan'), cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat_fitmatmeans.stan'), cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat_fitnothing.stan'))

fit_NGM_model_for_date_range = function(end_date='20201201', 
                                        period=30, 
                                        age_mod = age_mods[[1]],
                                        inf_matrix_mean, 
                                        inf_matrix_sd, 
                                        anb_matrix_mean, 
                                        anb_matrix_sd, 
                                        cms, 
                                        variables = c('inf_rate', 'susceptibility', 'ab_protection'),
                                        quantiles = c(0.05, 0.5, 0.95), 
                                        samps = 100, 
                                        pops = population, 
                                        runindex=1
                                        ){
  
  end_date = lubridate::ymd(end_date)
  start_date = end_date - period

  # filter matrices to correct dates
  inf_matrix_mean = inf_matrix_mean[date < end_date & date > start_date]
  inf_matrix_sd = inf_matrix_sd[date < end_date & date > start_date]
  anb_matrix_mean = anb_matrix_mean[date < end_date & date > start_date]
  anb_matrix_sd = anb_matrix_sd[date < end_date & date > start_date]
  
  inf_matrix_mean[, sr := cut(date, breaks=c(sr_dates$min_date, max(sr_dates$max_date)), labels = sort(sr_dates$survey_round) )]
  
  sr_dates = sr_dates[survey_round %in% inf_matrix_mean$sr]
  # create vector to indicate correct contact matrix
  day_to_week = as.numeric(inf_matrix_mean$sr)
  day_to_week = day_to_week - min(day_to_week) + 1
  
  # load contact matrices into list
  contact_matrices = lapply(X=sort(unique(sr_dates$survey_round)), function(X){t(matrix(cms[sr==X]$cms, nrow=7))})
  
  
  contact_matrices_mu = lapply(X=sort(unique(sr_dates$survey_round)), function(X){t(matrix(cms[sr==X]$cms, nrow=7))/tail(pops,-1)$props})
  contact_matrices_sd = lapply(X=sort(unique(sr_dates$survey_round)), function(X){t(matrix(cms[sr==X]$csd, nrow=7))/tail(pops,-1)$props})
  
  mean_contacts_mu = lapply(contact_matrices_mu, rowSums)
  mean_contacts_sd = lapply(contact_matrices_sd, rowSums)
  
  mean_contacts_mu_mat = sapply(mean_contacts_mu, mean)
  mean_contacts_sd_mat = sapply(mean_contacts_sd, mean)
  # set up data input for stan 
  data = list( 
    T = dim(inf_matrix_mean[,-c('date', 'sr')])[1],
    W = dim(sr_dates)[1],
    A = dim(inf_matrix_mean[,-c('date', 'sr')])[2],
    inf_mu = inf_matrix_mean[,-c('date', 'sr')],
    inf_sd = inf_matrix_sd[,-'date'],
    anb_mu = anb_matrix_mean[,-'date'],
    anb_sd = anb_matrix_sd[,-'date'],
    day_to_week_converter = day_to_week,
    contact_matrices = contact_matrices,
    contact_matrices_mu = contact_matrices_mu,
    contact_matrices_sd = contact_matrices_sd, 
    population = tail(pops,-1)$props, 
    mean_contacts_mu = mean_contacts_mu,
    mean_contacts_sd = mean_contacts_sd, 
    mean_contacts_mu_mat = mean_contacts_mu_mat,
    mean_contacts_sd_mat = mean_contacts_sd_mat
  )
  # compile stan model from 'stan/age_specific_transmission.stan'

  fit = age_mod$sample(data, 
                       iter_warmup = 250,
                       iter_sampling=250, 
                       chains=1, 
                       max_treedepth = 12, 
                       adapt_delta=0.8, 
                       seed = sample.int(.Machine$integer.max, 1))
  
  
  summary_pars = fit$summary(
    variables = variables, ~ quantile(.x, probs = quantiles)
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(name = variable) %>%
    dplyr::mutate(
      index = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name),
      forecast_date = end_date, 
      run = runindex
    ) 
  
  summary_preds = fit$summary(
    variables = c('next_gens', 'forecast_gens', 'infections'), ~ quantile(.x, probs = quantiles, na.rm=TRUE)
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(name = variable) %>%
    dplyr::mutate(
      time_index = as.integer(sub("^.*\\[([0-9]+)[,][0-9]+]$", "\\1", name)),
      age_index = as.integer(sub("^.*\\[[0-9]+[,]([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name),
      forecast_date = end_date, 
      run = runindex
    )
  
  summary_conts = fit$summary(
    variables = c('contact_matrices', 'contact_matrices_aug'), ~ quantile(.x, probs = quantiles, na.rm=TRUE)
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(name = variable) %>%
    dplyr::mutate(
      time_index = as.integer(sub("^.*\\[([0-9]+)[,][0-9]+[,][0-9]+]$", "\\1", name)),
      age1_index = as.integer(sub("^.*\\[[0-9]+[,]([0-9]+)[,][0-9]+]$", "\\1", name)),
      age2_index = as.integer(sub("^.*\\[[0-9]+[,][0-9]+[,]([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name),
      forecast_date = end_date, 
      run = runindex
    )
  
  
  draws_pars = fit$draws(variables) %>%
    posterior::as_draws_df() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(sample = 1:dplyr::n()) %>%
    dplyr::filter(sample <= samps) %>%
    dplyr::select(-.chain, -.iteration, -.draw) %>%
    tidyr::pivot_longer(matches("[0-9]") & !matches('gam')) %>%
    dplyr::mutate(
      index = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name),
      forecast_date = end_date, 
      run = runindex
    )
  
  draws_preds = fit$draws(c('next_gens', 'forecast_gens', 'infections')) %>%
    posterior::as_draws_df() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(sample = 1:dplyr::n()) %>%
    dplyr::filter(sample <= samps) %>%
    dplyr::select(-.chain, -.iteration, -.draw) %>%
    tidyr::pivot_longer(matches("[0-9]") & !matches('gam')) %>%
    dplyr::mutate(
      time_index = as.integer(sub("^.*\\[([0-9]+)[,][0-9]+]$", "\\1", name)),
      age_index = as.integer(sub("^.*\\[[0-9]+[,]([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name),
      forecast_date = end_date, 
      run = runindex
    )
  
  outs = list(
            fit = list(fit),
            samples_pars = draws_pars, 
            samples_preds = draws_preds,
            summary_pars = summary_pars, 
            summary_preds  =summary_preds, 
            summary_conts = summary_conts
  )
  
  return(outs)
}




dates = list(20201201, 20210201, 20210501, 20210801, 20211201)


plan(callr, workers = future::availableCores())


actualest_inc_long = melt(inf_matrix_mean[,-'sr'], id.vars = 'date', measure.vars = colnames(inf_matrix_mean[,-c('date', 'sr')]))

actualest_inc_long[, age_group := variable]

actualest_anb_long = melt(anb_matrix_mean[,-'sr'], id.vars = 'date', measure.vars = colnames(anb_matrix_mean[,-c('date', 'sr')]))

actualest_anb_long[, age_group := variable]

ggplot(actualest_anb_long)+
  geom_point(aes(x=date, y = value, color=variable))

all_est = list()
for(r in 1:4){ 
  est <- future_lapply(
    dates, fit_NGM_model_for_date_range,
    age_mod = age_mods[[r]],
    period = 30, 
    inf_matrix_mean = inf_matrix_mean, 
    inf_matrix_sd = inf_matrix_sd, 
    anb_matrix_mean = anb_matrix_mean, 
    anb_matrix_sd = anb_matrix_sd, 
    cms = cms, 
    runindex = r
  )
  
  all_est = append(all_est, est)
  
  }



age_labs = age_groups
names(age_labs) = 1:7

plot_trajectories = function(est){
  summary_pars = data.table()
  
  for(i in 1:length(est)){
    
    summary_pars = rbind(summary_pars, data.table(est[[i]]$summary_pars))
    
  }
  
  ggplot(data.table(summary_pars)[!is.na(index)])+
    geom_rect(aes(xmin=index-0.1, xmax=index+0.1, ymin=`5%`, ymax=`95%`, fill=name),alpha=0.5)+
    geom_point(aes(x=index, y=`50%`, color=name))+
    facet_grid(name~forecast_date)+
    scale_x_continuous(name='age group', breaks = 1:7, labels = age_groups)+
    scale_fill_discrete(name='')+
    scale_color_discrete(name='')+
    theme_minimal()
  
  
  
  
  
  
  
  summary_preds = data.table()
  for(i in 1:length(est)){
    
    summary_preds = rbind(summary_preds, data.table(est[[i]]$summary_preds))
    
    
    
  }
  
  
  summary_preds[, datesf := as.Date(forecast_date) + time_index]
  summary_preds[, datest := as.Date(forecast_date) - 30 + time_index]
  summary_preds[, age_group := age_groups[age_index]]
  
  summary_preds = merge(summary_preds, actualest_inc_long, by.x = c('datest', 'age_group'), by.y = c('date', 'age_group'), all.x = T, all.y = F )
  summary_preds = merge(summary_preds, actualest_inc_long, by.x = c('datesf', 'age_group'), by.y = c('date', 'age_group'), all.x=T, all.y = F )
  age_labs = age_groups
  names(age_labs) = 1:7
  
  
  ggplot()+
    geom_rect(data=summary_preds[name == 'infections' & run==2,], aes(xmin=datest-0.2, xmax=datest+0.2, ymin=`5%`, ymax=`95%`), fill='blue', alpha=0.3)+
    geom_point(data=summary_preds[name == 'next_gens'& run==2,], aes(x=datest, y=value.x), size=0.01)+
    geom_point(data=summary_preds[name == 'forecast_gens'& run==2,], aes(x=datesf, y=value.y), size=0.01)+
    geom_rect(data=summary_preds[name == 'next_gens'& run==2,], aes(xmin=datest-0.2, xmax=datest+0.2, ymin=`5%`, ymax=`95%`), fill='limegreen', alpha=0.5)+
    geom_rect(data=summary_preds[name == 'forecast_gens'& run==2,], aes(xmin=datesf-0.2, xmax=datesf+0.2, ymin=`5%`, ymax=`95%`), fill='orange', alpha=0.5)+
    geom_point(data=summary_preds[name == 'next_gens'& run==2,], aes(x=datest, y=`50%`), color='green', size=0.01)+
    geom_point(data=summary_preds[name == 'forecast_gens'& run==2,], aes(x=datesf, y=`50%`), color='red', size=0.01)+
    scale_x_date(name='')+
    scale_y_continuous(name='infections')+
    
    facet_grid(age_index~forecast_date, scales='free', labeller =labeller(age_index=age_labs))+
    theme_minimal()
  
  summary_pars
  
  
  
}

d = c(lubridate::ymd(dates[[4]]), lubridate::ymd(dates[[4]]))
ggplot()+
  geom_rect(data=summary_preds[name == 'infections' & forecast_date %in% d,], aes(xmin=datest-0.2, xmax=datest+0.2, ymin=`5%`, ymax=`95%`), fill='blue', alpha=0.3)+
  geom_point(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(x=datest, y=value.x), size=0.01)+
  geom_point(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(x=datesf, y=value.y), size=0.01)+
  geom_rect(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(xmin=datest-0.2, xmax=datest+0.2, ymin=`5%`, ymax=`95%`), fill='limegreen', alpha=0.5)+
  geom_rect(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(xmin=datesf-0.2, xmax=datesf+0.2, ymin=`5%`, ymax=`95%`, fill=run), alpha=0.2)+
  geom_point(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(x=datest, y=`50%`), color='green', size=0.01)+
  geom_point(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(x=datesf, y=`50%`, color=run), size=0.01)+
  
  facet_grid(age_index~forecast_date, scales='free', labeller =labeller(age_index=age_labs))+
  scale_x_date(name='')+
  scale_y_continuous(name='infections')+
  scale_color_discrete()+
  theme_minimal()



ggplot(data.table(summary_pars)[!is.na(index)])+
  geom_rect(aes(xmin=index-0.4 + as.numeric(run)*0.15, xmax=index - 0.3 + as.numeric(run)*0.15, ymin=`5%`, ymax=`95%`, fill=as.character(run)),alpha=0.5)+
  geom_point(aes(x=index-0.35 + as.numeric(run)*0.15, y=`50%`, color=as.character(run)))+
  facet_grid(name~forecast_date)+
  scale_x_continuous(name='age group', breaks = 1:7, labels = age_groups)+
  scale_fill_discrete(name='')+
  scale_color_discrete(name='')+
  theme_minimal()

summary_pars$run = as.character(summary_pars$run)



summary_pars = plot_trajectories(all_est)


dates_ = lapply(dates, lubridate::ymd)


scores = data.table()

for (d in 1:5){
  for (a in 1:7){
  
  pv = dcast(data.table(est[[d]]$samples_preds)[name=='forecast_gens' & age_index==a,], value.var = 'value', time_index ~ sample)
  tv = summary_preds[forecast_date==as.character(dates_[[d]]) & name == 'forecast_gens' & age_index== a]$value.x
  scores = 
    rbind(
      scores, 
      data.table(
        forecast_date = d,
        age_index = a, 
        time_index = 1:5,
        cprs = scoringutils::crps(tv, pv),
        log_score = scoringutils::logs(tv, pv)
      ) 
    )
}}
  
ggplot(scores) + 
  geom_point(aes(x=time_index, y=cprs))+
  facet_grid(forecast_date~age_index)


scoringutils::crps(summary_preds[forecast_date==as.Date(dates[[1]])])

fit_2weeks_dec2020 = fit_NGM_model_for_date_range('20201201', 30, age_mod, inf_matrix_mean, inf_matrix_sd, anb_matrix_mean, anb_matrix_sd, cms)
fit_2weeks_feb2021 = fit_NGM_model_for_date_range('20210201', 30, inf_matrix_mean, inf_matrix_sd, anb_matrix_mean, anb_matrix_sd, cms)
fit_2weeks_may2021 = fit_NGM_model_for_date_range('20210501', 30, inf_matrix_mean, inf_matrix_sd, anb_matrix_mean, anb_matrix_sd, cms)
fit_2weeks_aug2021 = fit_NGM_model_for_date_range('20210801', 30, inf_matrix_mean, inf_matrix_sd, anb_matrix_mean, anb_matrix_sd, cms)
fit_2weeks_dec2021 = fit_NGM_model_for_date_range('20211201', 30, inf_matrix_mean, inf_matrix_sd, anb_matrix_mean, anb_matrix_sd, cms)




stan_plot(fit_2weeks_dec2020, pars = c('inf_rate[1]', 'inf_rate[2]', 'inf_rate[3]', 'inf_rate[4]', 'inf_rate[5]', 'inf_rate[6]', 'inf_rate[7]'))
stan_plot(fit_2weeks_dec2020, pars = c('susceptibility[1]', 'susceptibility[2]', 'susceptibility[3]', 'susceptibility[4]', 'susceptibility[5]', 'susceptibility[6]', 'susceptibility[7]'))

stan_plot(fit_2weeks_feb2021, pars = c('inf_rate[1]', 'inf_rate[2]', 'inf_rate[3]', 'inf_rate[4]', 'inf_rate[5]', 'inf_rate[6]', 'inf_rate[7]'))
stan_plot(fit_2weeks_feb2021, pars = c('susceptibility[1]', 'susceptibility[2]', 'susceptibility[3]', 'susceptibility[4]', 'susceptibility[5]', 'susceptibility[6]', 'susceptibility[7]'))

stan_plot(fit_2weeks_may2021, pars = c('inf_rate[1]', 'inf_rate[2]', 'inf_rate[3]', 'inf_rate[4]', 'inf_rate[5]', 'inf_rate[6]', 'inf_rate[7]'))
stan_plot(fit_2weeks_may2021, pars = c('susceptibility[1]', 'susceptibility[2]', 'susceptibility[3]', 'susceptibility[4]', 'susceptibility[5]', 'susceptibility[6]', 'susceptibility[7]'))


stan_plot(fit_2weeks_aug2021, pars = c('inf_rate[1]', 'inf_rate[2]', 'inf_rate[3]', 'inf_rate[4]', 'inf_rate[5]', 'inf_rate[6]', 'inf_rate[7]'))
stan_plot(fit_2weeks_aug2021, pars = c('susceptibility[1]', 'susceptibility[2]', 'susceptibility[3]', 'susceptibility[4]', 'susceptibility[5]', 'susceptibility[6]', 'susceptibility[7]'))



stan_plot(fit_2weeks_dec2021, pars = c('inf_rate[1]', 'inf_rate[2]', 'inf_rate[3]', 'inf_rate[4]', 'inf_rate[5]', 'inf_rate[6]', 'inf_rate[7]'))
stan_plot(fit_2weeks_dec2021, pars = c('susceptibility[1]', 'susceptibility[2]', 'susceptibility[3]', 'susceptibility[4]', 'susceptibility[5]', 'susceptibility[6]', 'susceptibility[7]'))

