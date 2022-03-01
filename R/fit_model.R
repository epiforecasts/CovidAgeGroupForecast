



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
