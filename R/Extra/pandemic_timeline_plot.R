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
  facet_wrap(~age_group_no, ncol=1, labeller= labeller(age_group_no = age_labs))+
  scale_y_continuous(name='Infections')+
  scale_x_date()+
  theme_minimal_hgrid()+
  ggtitle('A')



ab_estimates = data.table(estimates)[name == 'gen_dab']

ab_estimates[, date := min(inf_estimates$date) + (t_index-1)]

ab_estimates[, mean_no := mean * population]
ab_estimates[, sd_no := sd * population]

ab_estimates_mu = ab_estimates[, c('mean_no', 'variable', 'date')]

anb_matrix_mean  = dcast(ab_estimates, value.var = 'mean', date ~ variable)
anb_matrix_sd = dcast(ab_estimates, value.var = 'sd', date ~ variable)
ab_estimates[, ano := as.character(age_group_no)]

ggplot(ab_estimates) + 
  geom_ribbon(aes(x=date, ymin=q10, ymax=q90, fill=ano), alpha=0.1)+
  geom_point(aes(x=date, y=mean, color=ano), size=0.3, stroke=0)+
  geom_rect(aes(x=))
  scale_y_continuous(name='Proportion with antibodies')+
  theme_minimal_hgrid()+
  scale_color_discrete(labels = age_groups)+
  scale_fill_discrete(labels = age_groups)+
  scale_x_date(date_labels = '%b %Y')

