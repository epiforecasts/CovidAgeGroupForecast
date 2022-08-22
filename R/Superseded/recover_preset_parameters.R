library(data.table)
library(rstan)
library(rstanarm)
library(ggplot2)

# load samples of infections and abprev from inc2prev
samples = readRDS('../../../inc2prev/outputs/samples_combined_age_ab_A.rds')
samples = data.table(samples)

# load mean contact matrices from contactmatrixcode
cms = readRDS('cms.rds')
sr_dates = readRDS('sr_dates.rds')
population = readRDS('population.rds')

# create matrices of mean and standard deviations for infections and ab_prevalence
inf_matrix_mean = dcast(samples[name=='infections', c('sample', 'variable', 'value', 'date')], 
                        date ~ variable, 
                        value.var = 'value', 
                        fun.aggregate = mean)
inf_matrix_sd = dcast(samples[name=='infections', c('sample', 'variable', 'value', 'date')], 
                      date ~ variable, 
                      value.var = 'value', 
                      fun.aggregate = sd)

anb_matrix_mean = dcast(samples[name=='dab', c('sample', 'variable', 'value', 'date')], 
                        date ~ variable, 
                        value.var = 'value', 
                        fun.aggregate = mean)
anb_matrix_sd = dcast(samples[name=='dab', c('sample', 'variable', 'value', 'date')], 
                      date ~ variable, 
                      value.var = 'value', 
                      fun.aggregate = sd)

# set the correct order for age columns in matrices 
age_groups =  c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )
colorder= c('date', age_groups )
setcolorder(inf_matrix_mean, colorder)
setcolorder(inf_matrix_sd, colorder)
setcolorder(anb_matrix_mean, colorder)
setcolorder(anb_matrix_sd, colorder)



# calculate mean of infections (From proportions of age groups)
anb_matrix_mean[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )] = 
  anb_matrix_mean[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )] * 0.9 

inf_matrix_mean[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )] = 
  data.table(t(t(inf_matrix_mean[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )]) * 
                 tail(population$pop_total,-1)))

inf_matrix_sd[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )] = 
  data.table(t(t(inf_matrix_sd[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )]) * 
                 tail(population$pop_total,-1)))


# assign the CoMix survey round to each time step 

start_date = lubridate::ymd('20201101')
end_date = lubridate::ymd('20210501')
# filter matrices to correct dates
inf_matrix_mean = inf_matrix_mean[date < end_date & date > start_date]
inf_matrix_sd = inf_matrix_sd[date < end_date & date > start_date]
anb_matrix_mean = anb_matrix_mean[date < end_date & date > start_date]
anb_matrix_sd = anb_matrix_sd[date < end_date & date > start_date]
sr_dates = sr_dates[min_date < end_date & max_date > start_date]

inf_matrix_mean[, sr := cut(date, breaks=c(sr_dates$min_date, max(sr_dates$max_date)), labels = sort(sr_dates$survey_round) )]


# create vector to indicate correct contact matrix
day_to_week = as.numeric(inf_matrix_mean$sr)
day_to_week = day_to_week - min(day_to_week) + 1

# load contact matrices into list
contact_matrices = lapply(X=sort(unique(sr_dates$survey_round)), function(X){t(matrix(cms[sr==X]$cms, nrow=7))})

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
  contact_matrices = contact_matrices
)
# compile stan model from 'stan/age_specific_transmission.stan'
age_mod_pp = stan_model('stan/age_specific_transmission_diff_prior_pred.stan')
# sample from stan model
fit = sampling(age_mod_pp, data, warmup=50, iter=100, chains=1, control = list(max_treedepth = 12, adapt_delta=0.99),seed = sample.int(.Machine$integer.max, 1))


stan_plot(fit, pars = c('inf_rate[1]', 'inf_rate[2]', 'inf_rate[3]', 'inf_rate[4]', 'inf_rate[5]', 'inf_rate[6]', 'inf_rate[7]'))
stan_plot(fit, pars = c('susceptibility[1]', 'susceptibility[2]', 'susceptibility[3]', 'susceptibility[4]', 'susceptibility[5]', 'susceptibility[6]', 'susceptibility[7]'))

fitmat = rstan::extract(fit)

fitmat$suscept_hyper_mu

inf_matrix_mean_sim = fitmat$next_gen[1,,]
inf_matrix_mean_sim[is.na(inf_matrix_mean_sim)] = 0.0


data_ng = list( 
  T = dim(inf_matrix_mean[,-c('date', 'sr')])[1],
  W = dim(sr_dates)[1],
  A = dim(inf_matrix_mean[,-c('date', 'sr')])[2],
  inf_mu = inf_matrix_mean[,-c('date', 'sr')],
  inf_mu_ng = inf_matrix_mean_sim,
  inf_sd = inf_matrix_sd[,-'date'],
  anb_mu = anb_matrix_mean[,-'date'],
  anb_sd = anb_matrix_sd[,-'date'],
  day_to_week_converter = day_to_week,
  contact_matrices = contact_matrices
)
# compile stan model from 'stan/age_specific_transmission.stan'
age_mod_ng = stan_model('stan/age_specific_transmission_ng_inf.stan')
# sample from stan model
fit_recover = sampling(age_mod_ng, data_ng, warmup=500, iter=1000, chains=1, control = list(max_treedepth = 12, adapt_delta=0.80),seed = sample.int(.Machine$integer.max, 1))


stan_plot(fit_recover, pars = c('inf_rate[1]', 'inf_rate[2]', 'inf_rate[3]', 'inf_rate[4]', 'inf_rate[5]', 'inf_rate[6]', 'inf_rate[7]'))
stan_plot(fit_recover, pars = c('susceptibility[1]', 'susceptibility[2]', 'susceptibility[3]', 'susceptibility[4]', 'susceptibility[5]', 'susceptibility[6]', 'susceptibility[7]'))
stan_plot(fit_recover, pars = c('ab_protection'))

