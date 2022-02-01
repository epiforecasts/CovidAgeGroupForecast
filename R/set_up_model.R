library(data.table)
library(rstan)

# load samples of infections and abprev from inc2prev
samples = readRDS('../../../inc2prev/outputs/samples_combined_age_ab.rds')
samples = data.table(samples)

# load mean contact matrices from contactmatrixcode
cms = readRDS('cms.rds')
sr_dates = readRDS('sr_dates.rds')

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

# filter matrices to correct dates
inf_matrix_mean = inf_matrix_mean[date >= min(sr_dates$min_date)]
inf_matrix_sd = inf_matrix_sd[date >= min(sr_dates$min_date)]
anb_matrix_mean = anb_matrix_mean[date >= min(sr_dates$min_date)]
anb_matrix_sd = anb_matrix_sd[date >= min(sr_dates$min_date)]

# create vector to indicate correct contact matrix
inf_matrix_mean[, sr := cut(date, breaks=c(sr_dates$min_date, max(sr_dates$max_date) + 100), labels = sr_dates$survey_round )]
day_to_week = as.numeric(inf_matrix_mean$sr)
day_to_week = day_to_week - min(day_to_week) + 1

# load contact matrices into list
contact_matrices = lapply(X=unique(sr_dates$survey_round), function(X){matrix(cms[sr==X]$cms, nrow=7)})

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
age_mod = stan_model('stan/age_specific_transmission.stan')
# sample from stan model
fit = sampling(age_mod, data, warmup=100, iter=500, chains=1)


# plot susceptibility and infectiousness parameters
pairs(fit, pars = c('inf_rate[1]', 'inf_rate[2]', 'inf_rate[3]', 'inf_rate[4]', 'inf_rate[5]', 'inf_rate[6]', 'inf_rate[7]'))
pairs(fit, pars = c('susceptibility[1]', 'susceptibility[2]', 'susceptibility[3]', 'susceptibility[4]', 'susceptibility[5]', 'susceptibility[6]', 'susceptibility[7]'))

