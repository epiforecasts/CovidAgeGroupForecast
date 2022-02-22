library(data.table)
library(rstan)
library(rstanarm)
library(ggplot2)

library(future.apply)
library(future.callr)
library(future)

library(cmdstanr)

age_mod = cmdstan_model('stan/age_specific_transmission-fit_contacts.stan')
# sample from stan model

fit_NGM_model_for_date_range = function(end_date='20201201', 
                                        period=14, 
                                        age_mod = age_mod,
                                        inf_matrix_mean, 
                                        inf_matrix_sd, 
                                        anb_matrix_mean, 
                                        anb_matrix_sd, 
                                        cms, 
                                        variables = c('inf_rate', 'susceptibility', 'ab_protection'),
                                        quantiles = c(0.05, 0.5, 0.95), 
                                        samps = 100
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
  contact_matrices_sd = lapply(X=sort(unique(sr_dates$survey_round)), function(X){t(matrix(cms[sr==X]$csd, nrow=7))})
  
  
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
    contact_matrices_mu = contact_matrices,
    contact_matrices_sd = contact_matrices_sd
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
      forecast_date = end_date
    ) 
  
  summary_preds = fit$summary(
    variables = c('next_gens', 'forecast_gens'), ~ quantile(.x, probs = quantiles, na.rm=TRUE)
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(name = variable) %>%
    dplyr::mutate(
      time_index = as.integer(sub("^.*\\[([0-9]+)[,][0-9]+]$", "\\1", name)),
      age_index = as.integer(sub("^.*\\[[0-9]+[,]([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name),
      forecast_date = end_date
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
      forecast_date = end_date
    )
  
  draws_preds = fit$draws(c('next_gens', 'forecast_gens')) %>%
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
      forecast_date = end_date
    )
  
  outs = list(
            fit = list(fit),
            samples_pars = draws_pars, 
            samples_preds = draws_preds,
            summary_pars = summary_pars, 
            summary_preds  =summary_preds
  )
  
  return(outs)
}


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


dates = list(20201201, 20210201, 20210501, 20210801, 20211201)


plan(callr, workers = future::availableCores())


actualest_inc_long = melt(inf_matrix_mean[,-'sr'], id.vars = 'date', measure.vars = colnames(inf_matrix_mean[,-c('date', 'sr')]))

actualest_inc_long[, age_group := variable]


est <- future_lapply(
  dates, fit_NGM_model_for_date_range,
  age_mod = age_mod,
  period = 30, 
  inf_matrix_mean = inf_matrix_mean, 
  inf_matrix_sd = inf_matrix_sd, 
  anb_matrix_mean = anb_matrix_mean, 
  anb_matrix_sd = anb_matrix_sd, 
  cms = cms
)

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
  geom_point(data=summary_preds[name == 'next_gens',], aes(x=datest, y=value.x), size=0.01)+
  geom_point(data=summary_preds[name == 'forecast_gens',], aes(x=datesf, y=value.y), size=0.01)+
  geom_ribbon(data=summary_preds[name == 'next_gens',], aes(x=datest, ymin=`5%`, ymax=`95%`), fill='limegreen')+
  geom_ribbon(data=summary_preds[name == 'forecast_gens',], aes(x=datesf, ymin=`5%`, ymax=`95%`), fill='orange')+
  geom_line(data=summary_preds[name == 'next_gens',], aes(x=datest, y=`50%`), color='green')+
  geom_line(data=summary_preds[name == 'forecast_gens',], aes(x=datesf, y=`50%`), color='red')+

  facet_grid(age_index~forecast_date, scales='free', labeller =labeller(age_index=age_labs))+
  theme_minimal()

variables = c('inf_rate', 'susceptibility', 'ab_protection')
quantiles = c(0.05, 0.5, 0.95)



ggplot()+
  geom_point(data=fit_2weeks_dec2020$summary_preds[name == 'next_gens',], aes(x=datest, y=value.x), size=0.01)+
  geom_point(data=fit_2weeks_dec2020$summary_preds[name == 'forecast_gens',], aes(x=datesf, y=value.y), size=0.01)+
  geom_ribbon(data=fit_2weeks_dec2020$summary_preds[name == 'next_gens',], aes(x=datest, ymin=`5%`, ymax=`95%`), fill='limegreen')+
  geom_ribbon(data=fit_2weeks_dec2020$summary_preds[name == 'forecast_gens',], aes(x=datesf, ymin=`5%`, ymax=`95%`), fill='orange')+
  geom_line(data=fit_2weeks_dec2020$summary_preds[name == 'next_gens',], aes(x=datest, y=`50%`), color='green')+
  geom_line(data=fit_2weeks_dec2020$summary_preds[name == 'forecast_gens',], aes(x=datesf, y=`50%`), color='red')+
  
  facet_grid(age_index~forecast_date, scales='free', labeller =labeller(age_index=age_labs))+
  theme_minimal()

variables = c('inf_rate', 'susceptibility', 'ab_protection')
quantiles = c(0.05, 0.5, 0.95)



est <- rbindlist(est, use.names = TRUE, fill = TRUE)
# Add summary information to posterior summary and samples
fit_2weeks_dec2020[, summary_pars := map2(summary_pars, name, ~ as.data.table(.x)[, name := .y])]
est[, summary := map2(summary, level, ~ as.data.table(.x)[, level := .y])]
est[, samples := map2(samples, variable, ~ as.data.table(.x)[, variable := .y])]
est[, samples := map2(samples, level, ~ as.data.table(.x)[, level := .y])]

