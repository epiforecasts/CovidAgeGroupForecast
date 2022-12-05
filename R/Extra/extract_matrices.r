library(data.table)

quantiles = c(0.05, 0.5, 0.95)
summary_pars = fit$fit[[1]]$summary(
  variables = c('contact_matrices_aug'), ~ quantile(.x, probs = quantiles))%>%
  dplyr::as_tibble() %>%
  dplyr::rename(name = variable)
summary_pars = data.table(summary_pars)

summary_pars[ , first_col := substr(name, 22, 22)]

summary_pars[ , second_col := substr(name, 24, 24)]

summary_pars[ , third_col := substr(name, 26, 26)]


image(matrix(summary_pars[first_col=='5']$`50%`, ncol=7))


library(data.table)

quantiles = c(0.05, 0.5, 0.95)
summary_pars = fit$fit[[1]]$summary(
  variables = c('contact_matrices'), ~ quantile(.x, probs = quantiles))%>%
  dplyr::as_tibble() %>%
  dplyr::rename(name = variable)
summary_pars = data.table(summary_pars)

summary_pars[ , first_col := substr(name, 18, 18)]

summary_pars[ , second_col := substr(name, 24, 24)]

summary_pars[ , third_col := substr(name, 26, 26)]


image(matrix(summary_pars[first_col=='5']$`50%`, ncol=7))
