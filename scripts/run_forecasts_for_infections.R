library(data.table)
library(ggplot2)
library(tidyverse)

library(future.apply)
library(future.callr)
library(future)

library(cmdstanr)
library(cowplot)
library(khroma)

source('R/fit_model.R')
source('R/make_infs_weekly.R')
source('R/get_polymod_cms.R')
source('R/generate_baselines.R')
library(patchwork)

breaks = c(2,11,16,25,35,50,70,Inf)


# set age group objects 
age_groups =  c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )
age_lookup = setNames(object=1:7, nm = age_groups)
age_labs = age_groups
names(age_labs) = 1:7

# load infection and antibody estimates from CSV file (output from inc2prev)
estimates = read.csv(file = 'data/estimates_age_ab.csv')

select_dates = seq(as.Date('2020-01-31'), as.Date('2022-01-01'), 7)

# extract population data 
population = data.table(unique(estimates[,c('variable', 'population')]))
population = population[order(match(unique(variable), c("2-10", "11-15", "16-24", "25-34", "35-49", "50-69", "70+"))),]
population[, props := population/sum(population)]

# give each age group an integer id
estimates = data.table(estimates)
estimates[, age_group_no := age_lookup[variable]]

inf_estimates = make_infs_weekly()
inf_estimates[, mean_no := mean]
inf_estimates[, sd_no := sd]

inf_estimates[, age_group_no := age_lookup[variable]]

# plot of infection time series
light = colour('light')

inf_traj = ggplot(inf_estimates) + 
  geom_pointrange(aes(x=date, y=mean, ymax=q90, ymin=q10, color=as.character(age_group_no)), size=0.3, stroke=0)+
  facet_wrap(~age_group_no, nrow=1, labeller= labeller(age_group_no = age_labs))+
  scale_y_continuous(name='Antibody prevalence')+
  scale_x_date(date_labels = '%b-%y', name='Date')+
  scale_color_manual(name='Age group', values = as.vector(light(9)), labels =age_labs)+
  theme_minimal_hgrid()+
  theme(axis.text.x = element_text(angle = 45, size=10))+
  ggtitle('B')



inf_estimates_weekly = inf_estimates[,weekly_mean := mean, by=c('variable')]
inf_estimates_weekly = inf_estimates[,weekly_sd := sd, by=c('variable')]


inf_estimates_weekly = inf_estimates_weekly[date %in% select_dates,]



inf_matrix_mean  = dcast(inf_estimates_weekly, value.var = 'weekly_mean', date ~ variable)
inf_matrix_sd = dcast(inf_estimates_weekly, value.var = 'weekly_sd', date ~ variable)

# extract daily antibody data 
ab_estimates = data.table(estimates)[name == 'gen_dab']
ab_estimates[, date := min(inf_estimates$date) + (t_index-1)]
ab_estimates[, mean_no := mean * population]
ab_estimates[, sd_no := sd * population]
ab_estimates_mu = ab_estimates[, c('mean_no', 'variable', 'date')]


# construct mean and sd of antibody prevalence matrices
ab_estimates_weekly = ab_estimates[date %in% select_dates,]
anb_matrix_mean  = dcast(ab_estimates_weekly, value.var = 'mean', date ~ variable)
anb_matrix_sd = dcast(ab_estimates_weekly, value.var = 'sd', date ~ variable)

ab_estimates_weekly[, age_group_no := age_lookup[variable]]


# plot antibody time series
anb_traj = ggplot(ab_estimates_weekly) + 
  #geom_ribbon(aes(x=date, ymin=q10, ymax=q90))+
  geom_pointrange(aes(x=date, y=mean, ymax=q90, ymin=q10, color=as.character(age_group_no)), size=0.3, stroke=0)+
  facet_wrap(~age_group_no, nrow=1, labeller= labeller(age_group_no = age_labs))+
  scale_y_continuous(name='Antibody prevalence')+
  scale_x_date(date_labels = '%b-%y', name='Date')+
  scale_color_manual(name='Age group', values = as.vector(light(9)), labels =age_labs)+
  theme_minimal_hgrid()+
  theme(axis.text.x = element_text(angle = 45, size=10))+
  ggtitle('B')


saveRDS(inf_traj, 'data/inf_traj.rds')
saveRDS(anb_traj, 'data/anb_traj.rds')


# combine and save time series plots
input_plot = inf_traj/anb_traj
ggsave('plots/inputs.pdf', input_plot, width=10, height=5)
ggsave('plots/inputs.png', input_plot, width=10, height=5, units = 'in', dpi=300)

# load mean contact matrices from contactmatrixcode
cms = readRDS('data/cms.rds')

# extract survey round dates
sr_dates = readRDS('data/sr_dates.rds')
population = readRDS('data/population.rds')


# plot survey rounds
ggplot(sr_dates) + 
  geom_rect(aes(xmin=min_date, xmax=max_date, ymax=survey_round+0.4, ymin=survey_round-0.4, fill=survey_round))+
  scale_fill_viridis()+
  scale_y_continuous(trans='reverse', name='Survey Round')+
  theme_minimal_vgrid()+
  theme(legend.position = 'none')

ggsave('plots/survey_rounds.png', width=10, height=3)


# set the correct order for age columns in matrices 
age_groups =  c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )
colorder= c('date', age_groups )
setcolorder(inf_matrix_mean, colorder)
setcolorder(inf_matrix_sd, colorder)
setcolorder(anb_matrix_mean, colorder)
setcolorder(anb_matrix_sd, colorder)

# extract time series of medians from outputs
actualest_inc_long = melt(inf_matrix_mean[,-'sr'], id.vars = 'date', measure.vars = colnames(inf_matrix_mean[,-c('date', 'sr')]))
actualest_inc_long[, age_group := variable]

ggplot(actualest_inc_long)+
  geom_point(aes(x=date, y = value, color=variable))

# compile stan model
age_mod = cmdstan_model('stan/multi-option-contact-model.stan')

# select forecast dates
dates = seq(as.Date("2020-11-27")-(4*7), as.Date('2021-12-02'), 14)
# set period over which to fit for each forecast and maximum generation interval
period = 8*7 
smax = 4 

# call the fit_NGM_model_for_date_range for each date in dates

plan(callr, workers = future::availableCores()-1)
all_est = list()


for(r in c(1,4,5)){ 

  est <- future_lapply(
    dates, fit_NGM_model_for_date_range,
    age_model = age_mod,
    period = period+smax*7, 
    smax = smax,
    inf_matrix_mean = inf_matrix_mean, 
    inf_matrix_sd = inf_matrix_sd, 
    anb_matrix_mean = anb_matrix_mean, 
    anb_matrix_sd = anb_matrix_sd, 
    cms = cms, 
    runindex = r, 
    quantiles=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
    contact_option=r,
    sigma_option=2,
    forecast_horizon = 4,
    future.seed=TRUE
  )
  
  all_est = append(all_est, est)
  
}



cms_pmd = get_polymod_cms(breaks)

age_lookup = 1:7; names(age_lookup) = levels(cms_pmd$Var1)

cms_pmd = data.table(cms_pmd)

cms_pmd[, level1:=age_lookup[Var1]]
cms_pmd[, level2:=age_lookup[Var2]]

cms_pmd = cms_pmd[order(sr, cms_pmd$level1, cms_pmd$level2),]

est <- future_lapply(
  dates, fit_NGM_model_for_date_range,
  age_model = age_mod,
  period = period+smax*7, 
  smax = smax,
  inf_matrix_mean = inf_matrix_mean, 
  inf_matrix_sd = inf_matrix_sd, 
  anb_matrix_mean = anb_matrix_mean, 
  anb_matrix_sd = anb_matrix_sd, 
  cms = cms_pmd, 
  runindex = 6, 
  quantiles=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
  contact_option=1,
  sigma_option=2,
  forecast_horizon = 4,
  future.seed=TRUE, 
  ad=0.9
)

all_est = append(all_est, est)
  



# combine the outputs from all dates into one container
summary_pars = data.table()
for(i in 1:length(all_est)){
  summary_pars = rbind(summary_pars, data.table(all_est[[i]]$summary_pars))
}


summary_diags = data.table()
for(i in 1:length(all_est)){
  summary_diags = rbind(summary_diags, data.table(all_est[[i]]$diagnostics))
}

summary_diags[run==1,]
summary_diags[run==4,]
summary_diags[run==5,]
summary_diags[run==6,]

summary_preds = data.table()
for(i in 1:length(all_est)){
  summary_preds = rbind(summary_preds, data.table(all_est[[i]]$summary_preds))
}

samples_preds = data.table()
for(i in 1:length(all_est)){
  samples_preds = rbind(samples_preds, data.table(all_est[[i]]$samples_preds))
}

samples_preds[, date := forecast_date - (period+smax*7) + (time_index-1)*7]
samples_preds[, age_group := age_groups[age_index]]
# merge with time series of median infection estimates 
samples_preds = merge(samples_preds[, prediction:=value][, -c('value')], actualest_inc_long[,  true_value := value][,-c('value')], by.x = c('date', 'age_group'), by.y = c('date', 'age_group'), all.x = T, all.y = F )

# calculate baseline


#summary_preds = data.table(fit$summary_preds)
# set dates and age groups based on time age index
summary_preds[, date := forecast_date - (period+smax*7) + (time_index-1)*7]
summary_preds[, age_group := age_groups[age_index]]
# merge with time series of median infection estimates 
summary_preds = merge(summary_preds, actualest_inc_long, by.x = c('date', 'age_group'), by.y = c('date', 'age_group'), all.x = T, all.y = F )

# calculate baselines using get_baselines function 
baselines = get_baselines(inf_estimates_weekly, summary_preds, smax = 4, period=8, forecast_step = 7)
# merge baselines with summary outputs
summary_preds = rbind(summary_preds, baselines)

baselines_samples = baselines[, prediction := `50%`][, -c("5%"    ,    "10%"   ,     "25%"    ,    "50%"    ,    "75%"    ,    "90%"  , "95%")]

baselines_samples_mult = baselines_samples %>% mutate(sample=1)

for(samp in 2:100){
  baselines_samples_mult = rbind(baselines_samples_mult , baselines_samples %>% mutate(sample=samp))
}

samples_preds

baselines_samples_mult = baselines_samples_mult[, true_value := value][,-c('value')]

samples_preds = rbind(samples_preds, baselines_samples_mult)

summary_conts = data.table()
for(i in 1:length(all_est)){
  summary_conts = rbind(summary_conts, data.table(all_est[[i]]$summary_conts))
}


# save outputs 
saveRDS(summary_pars, 'outputs/summary_pars_infweek.rds')
saveRDS(summary_preds, 'outputs/summary_preds_infweek.rds')
saveRDS(summary_conts, 'outputs/summary_conts_infweek.rds')

saveRDS(summary_diags, 'outputs/summary_diags_infweek.rds')


saveRDS(samples_preds, 'outputs/samples_preds_infweek.rds')



