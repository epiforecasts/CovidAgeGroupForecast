



fit_NGM_model_for_date_range_cases = function(end_date='20211001', 
                                        period=30, 
                                        smax = 20,
                                        age_model = age_mod_cases,
                                        inf_matrix_mean, 
                                        inf_matrix_sd, 
                                        anb_matrix_mean, 
                                        anb_matrix_sd, 
                                        cms, 
                                        variables = c('inf_rate', 'susceptibility'),
                                        quantiles = c(0.05, 0.5, 0.95), 
                                        samps = 100, 
                                        pops = tail(population, -1), 
                                        runindex=1, 
                                        contact_option=1, 
                                        sigma_option=1, 
                                        contact_delay=0, 
                                        forecast_horizon=10,
                                        ad=0.8

){
  
  # set the end and start dates as date objects
  end_date = lubridate::ymd(end_date)
  start_date = end_date - period
  
  # filter matrices to correct dates
  inf_matrix_mean = inf_matrix_mean[date <= end_date & date >= start_date]
  if (!is.null(inf_matrix_sd)){
      inf_matrix_sd = inf_matrix_sd[date <= end_date & date <= start_date]
      anb_matrix_mean = anb_matrix_mean[date <= end_date & date <= start_date]
      anb_matrix_sd = anb_matrix_sd[date <= end_date & date <= start_date]
  
  }
  inf_matrix_mean[, sr := cut(date-contact_delay, breaks=c(sr_dates$min_date, max(sr_dates$max_date)), labels = sort(sr_dates$survey_round) )]
  
  sr_dates = sr_dates[survey_round %in% inf_matrix_mean$sr]
  # create vector to indicate correct contact matrix
  day_to_week = as.numeric(inf_matrix_mean$sr)
  day_to_week = day_to_week - min(day_to_week) + 1
  
  day_to_week = as.numeric(factor(day_to_week,levels=unique(day_to_week)))
  
  A = dim(inf_matrix_mean[,-c('date', 'sr')])[2]
  
  # load contact matrices into list
  contact_matrices = lapply(X=sort(unique(sr_dates$survey_round)), function(X){t(matrix(cms[sr==X]$cms, nrow=A))})
  
  
  # calculate symmetric contact matrix means and SD from samples
  contact_matrices_mu = lapply(X=sort(unique(sr_dates$survey_round)), function(X){t(matrix(cms[sr==X]$cms, nrow=A))/pops$props})
  contact_matrices_sd = lapply(X=sort(unique(sr_dates$survey_round)), function(X){t(matrix(cms[sr==X]$csd, nrow=A))/pops$props})
  
  # calculate  mean and sd of contacts by age group
  mean_contacts_mu = lapply(contact_matrices_mu, rowSums)
  mean_contacts_sd = lapply(contact_matrices_sd, rowSums)
  
  # calculate  mean and sd of contacts overall
  mean_contacts_mu_mat = sapply(mean_contacts_mu, mean)
  mean_contacts_sd_mat = sapply(mean_contacts_sd, mean)
  
  
  # construct data object for stan unput
  inf_matrix_mean_inp = inf_matrix_mean[,-c('date', 'sr')]
  if(!is.null(inf_matrix_sd)){ 
    
    inf_matrix_sd_inp = inf_matrix_sd[,-'date']
    anb_matrix_sd_inp = anb_matrix_sd[,-'date']
    anb_matrix_sd_inp = anb_matrix_sd[,-'date']
  }
  else{
    inf_matrix_sd_inp = NULL
    anb_matrix_sd_inp = NULL
    anb_matrix_sd_inp = NULL
  }
  
  data = list( 
    T = dim(inf_matrix_mean_inp)[1],
    W = dim(sr_dates)[1],
    A = dim(inf_matrix_mean_inp)[2],
    cases = inf_matrix_mean_inp,
    day_to_week_converter = day_to_week,
    contact_matrices = contact_matrices,
    contact_matrices_mu = contact_matrices_mu,
    contact_matrices_sd = contact_matrices_sd, 
    population = pops$props, 
    mean_contacts_mu = mean_contacts_mu,
    mean_contacts_sd = mean_contacts_sd, 
    mean_contacts_mu_mat = mean_contacts_mu_mat,
    mean_contacts_sd_mat = mean_contacts_sd_mat, 
    smax=smax, 
    horizon=forecast_horizon, 
    contact_option=contact_option, 
    sigma_option=sigma_option
  )
  
  
  init_fun <- function(...) list(w_mu=5.0, w_sigma=1.5)
  # fit the model  
  
  w_sig = 5.0/7.0; 
  w_mu = 5.0/7.0;
  
  w_prior_sig = log(((w_sig^2)/(w_mu^2))  + 1);
  w_prior_mu  = log(w_mu) - (w_prior_sig)/2;
  fit = age_model$sample(data, 
                         init = list(list(lnmu = w_prior_mu, lnsig2 = w_prior_sig)),
                       iter_warmup = 250,
                       iter_sampling=250, 
                       chains=1, 
                       max_treedepth = 12, 
                       adapt_delta=ad, 
                       seed = sample.int(.Machine$integer.max, 1))
  
  
  out <- data.table(
    fit = list(fit),
    data = list(data))
  
    diag <- fit$sampler_diagnostics(format = "df")
    diagnostics <- data.table(
      samples = nrow(diag),
      max_rhat = round(max(
        fit$summary(
          variables = NULL, posterior::rhat,
          .args = list(na.rm = TRUE)
        )$`posterior::rhat`,
        na.rm = TRUE
      ), 2),
      divergent_transitions = sum(diag$divergent__),
      per_divergent_transitions = sum(diag$divergent__) / nrow(diag),
      max_treedepth = max(diag$treedepth__)
    )
    diagnostics[, no_at_max_treedepth := sum(diag$treedepth__ == max_treedepth)]
    diagnostics[, per_at_max_treedepth := no_at_max_treedepth / nrow(diag)]
    out <- cbind(out, diagnostics)
    
    timing <- round(fit$time()$total, 1)
    out[, run_time := timing]
    out[, date := end_date]
    out[, run := runindex]

  
  # create data frames from stan outputs
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
    variables = c('next_gens', 'forecast_gens'), ~ quantile(.x, probs = quantiles, na.rm=TRUE)
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
      forecast_date = end_date, 
      run = runindex
    )
  
  
  
  # construct output container
  outs = list(
    fit = list(fit),
    samples_pars = draws_pars, 
    samples_preds = draws_preds,
    summary_pars = summary_pars, 
    summary_preds  =summary_preds, 
    summary_conts = summary_conts, 
    diagnostics = out

  )
  
  return(outs)
}
