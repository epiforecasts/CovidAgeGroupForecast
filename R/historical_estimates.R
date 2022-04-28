library(data.table)
library(ggplot2)
library(tidyverse)

library(future.apply)
library(future.callr)
library(future)

library(cmdstanr)

source('R/fit_model.R')
source('R/plotting.R')
library(patchwork)
age_groups =  c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )
age_lookup = setNames(object=1:7, nm = age_groups)

age_labs = age_groups
names(age_labs) = 1:7


# load samples of infections and abprev from inc2prev
#samples = readRDS('../../../inc2prev/outputs/samples_combined_age_ab_A.rds')
samples = readRDS('../samples_age_ab.rds')
samples = data.table(samples)

estimates_ = readRDS('../../../inc2prev/outputs/estimates_combined_age_ab_A.rds')
ab_samples = data.table(samples)[name == 'dab']
estimates = read.csv(file = '../estimates_age_ab.csv')

population = data.table(unique(estimates[,c('variable', 'population')]))
population = population[order(match(unique(variable), c("2-10", "11-15", "16-24", "25-34", "35-49", "50-69", "70+"))),]
population[, props := population/sum(population)]
estimates = data.table(estimates)
estimates[, age_group_no := age_lookup[variable]]

inf_estimates = data.table(estimates)[name == 'infections']
inf_estimates[, date := as.Date(date)]
inf_estimates_mu = inf_estimates[, c('mean', 'variable', 'date')]

inf_estimates[, mean_no := mean * population]
inf_estimates[, sd_no := sd * population]


inf_traj = ggplot(inf_estimates) + 
  geom_ribbon(aes(x=date, ymin=q10, ymax=q90))+
  geom_point(aes(x=date, y=mean_no), size=0.3, stroke=0, color='red')+
  facet_wrap(~age_group_no, nrow=1, labeller= labeller(age_group_no = age_labs))+
  scale_y_continuous(name='Infections')+
  scale_x_date(breaks = c(), name='')+
  theme_minimal_hgrid()+
  ggtitle('A')

inf_matrix_mean  = dcast(inf_estimates, value.var = 'mean_no', date ~ variable)
inf_matrix_sd = dcast(inf_estimates, value.var = 'sd_no', date ~ variable)


ab_estimates = data.table(estimates)[name == 'gen_dab']

ab_estimates[, date := min(inf_estimates$date) + (t_index-1)]

ab_estimates[, mean_no := mean * population]
ab_estimates[, sd_no := sd * population]

ab_estimates_mu = ab_estimates[, c('mean_no', 'variable', 'date')]

anb_matrix_mean  = dcast(ab_estimates, value.var = 'mean', date ~ variable)
anb_matrix_sd = dcast(ab_estimates, value.var = 'sd', date ~ variable)



anb_traj = ggplot(ab_estimates) + 
  geom_ribbon(aes(x=date, ymin=q10, ymax=q90))+
  geom_point(aes(x=date, y=mean), size=0.3, stroke=0, color='DodgerBlue')+
  facet_wrap(~age_group_no, nrow=1, labeller= labeller(age_group_no = age_labs))+
  scale_y_continuous(name='Proportion with antibodies')+
  scale_x_date(date_labels = '%b-%y')+
  theme_minimal_hgrid()+
  theme(axis.text.x = element_text(angle = 45))+
  ggtitle('B')

input_plot = inf_traj/anb_traj

ggsave('plots/inputs.pdf', input_plot, width=10, height=5)
ggsave('plots/inputs.png', input_plot, width=10, height=5, units = 'in', dpi=300)

# load mean contact matrices from contactmatrixcode
cms = readRDS('cms.rds')
sr_dates = readRDS('sr_dates.rds')
population = readRDS('population.rds')

ggplot(sr_dates) + 
  geom_rect(aes(xmin=min_date, xmax=max_date, ymax=survey_round+0.4, ymin=survey_round-0.4, fill=survey_round))+
  scale_fill_viridis()+
  scale_y_continuous(trans='reverse', name='Survey Round')+
  theme_minimal_vgrid()+
  theme(legend.position = 'none')

ggsave('plots/survey_rounds.png', width=10, height=3)



# create matrices of mean and standard deviations for infections and ab_prevalence
#inf_matrix_mean = dcast(samples[name=='infections', c('sample', 'variable', 'value', 'date')], 
#                        date ~ variable, 
#                        value.var = 'value', 
#                        fun.aggregate = mean)
#inf_matrix_sd = dcast(samples[name=='infections', c('sample', 'variable', 'value', 'date')], 
#                      date ~ variable, 
#                      value.var = 'value', 
#                      fun.aggregate = sd)
#
#anb_matrix_mean = dcast(samples[name=='dab', c('sample', 'variable', 'value', 'date')], 
#                        date ~ variable, 
#                        value.var = 'value', 
#                        fun.aggregate = mean)
#anb_matrix_sd = dcast(samples[name=='dab', c('sample', 'variable', 'value', 'date')], 
#                      date ~ variable, 
#                      value.var = 'value', 
#                      fun.aggregate = sd)
#
# set the correct order for age columns in matrices 
age_groups =  c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )
colorder= c('date', age_groups )
setcolorder(inf_matrix_mean, colorder)
setcolorder(inf_matrix_sd, colorder)
setcolorder(anb_matrix_mean, colorder)
setcolorder(anb_matrix_sd, colorder)



#inf_matrix_mean[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )] = 
#  data.table(t(t(inf_matrix_mean[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )]) * 
#                 population$population))
#
#inf_matrix_sd[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )] = 
#  data.table(t(t(inf_matrix_sd[, c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )]) * 
#                 population$population))



actualest_inc_long = melt(inf_matrix_mean[,-'sr'], id.vars = 'date', measure.vars = colnames(inf_matrix_mean[,-c('date', 'sr')]))
actualest_inc_long[, age_group := variable]

ggplot(actualest_inc_long)+
  geom_point(aes(x=date, y = value, color=variable))


age_mods = list(cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat.stan'), cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat_fitmeans.stan'), cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat_fitmatmeans.stan'), cmdstan_model('stan/age_specific_transmission-fit_contacts_symmat_fitnothing.stan'))

dates = sort(head(tail(sr_dates$min_date, -5), -10))

#dates = sapply(dates, lubridate::ymd)



period = 30 

plan(callr, workers = future::availableCores())
all_est = list()
for(r in 1:4){ 
  est <- future_lapply(
    dates, fit_NGM_model_for_date_range,
    age_mod = age_mods[[r]],
    period = period, 
    inf_matrix_mean = inf_matrix_mean, 
    inf_matrix_sd = inf_matrix_sd, 
    anb_matrix_mean = anb_matrix_mean, 
    anb_matrix_sd = anb_matrix_sd, 
    cms = cms, 
    runindex = r, 
    quantiles=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
    future.seed=TRUE
  )
  
  all_est = append(all_est, est)
  
}



summary_pars = data.table()
for(i in 1:length(all_est)){
  summary_pars = rbind(summary_pars, data.table(all_est[[i]]$summary_pars))
}

summary_preds = data.table()
for(i in 1:length(all_est)){
  summary_preds = rbind(summary_preds, data.table(all_est[[i]]$summary_preds))
}


summary_preds[name=='forecast_gens', date := as.Date(forecast_date) + time_index]
summary_preds[name=='next_gens',  date := as.Date(forecast_date) - period + time_index]
summary_preds[, age_group := age_groups[age_index]]

summary_preds = merge(summary_preds, actualest_inc_long, by.x = c('date', 'age_group'), by.y = c('date', 'age_group'), all.x = T, all.y = F )


baselines = get_baselines(inf_estimates, summary_preds)
#baselines[, variable := age_group]
#baselines [, value := NA]
summary_preds = rbind(summary_preds, baselines)


summary_conts = data.table()
for(i in 1:length(all_est)){
  summary_preds = rbind(summary_preds, data.table(all_est[[i]]$summary_conts))
}

saveRDS(summary_pars, 'outputs/summary_pars.rds')
saveRDS(summary_preds, 'outputs/summary_preds.rds')
saveRDS(summary_conts, 'outputs/summary_conts.rds')

d = lapply(dates, lubridate::ymd)
unlist(d)

d = d[10:15]


sort(dates)

plot_parameters(summary_pars, d)





plot_trajectories_one_ax(summary_preds, d)

ggsave("plots/trajectories_zoomed.pdf", width = 10, height=7)
ggsave("plots/trajectories_zoomed.png", width = 10, height=7,  units = 'in', dpi = 300 )

