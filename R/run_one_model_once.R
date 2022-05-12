fit = fit_NGM_model_for_date_range(
  end_date = '20201001',
  age_mod = age_mod,
  period = period+smax, 
  inf_matrix_mean = inf_matrix_mean, 
  inf_matrix_sd = inf_matrix_sd, 
  anb_matrix_mean = anb_matrix_mean, 
  anb_matrix_sd = anb_matrix_sd, 
  cms = cms, 
  runindex = 1, 
  quantiles=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
  contact_option = 1
)
summary_preds = data.table(fit$summary_preds)
summary_preds[, age_group := age_groups[age_index]]
summary_preds[, date := forecast_date - 49 + time_index]
summary_preds = merge(summary_preds, actualest_inc_long, by=c('date', 'age_group'), all.x = TRUE)

ggplot() +
  geom_point(data=summary_preds[time_index>=50 & name=='forecast_gens'], aes(x=date, y=`50%`, color=age_group), alpha=0.8)+ 
  geom_ribbon(data=summary_preds[name=='forecast_gens'], aes(x=date, ymin=`5%`, ymax=`95%`, fill=age_group), alpha=0.2)+
  geom_line(data = summary_preds[time_index>20 & time_index<50 & name=='next_gens'], 
            aes(x=date, y=`50%`, color=age_group), alpha=0.8)+
  geom_point(data = summary_preds[time_index>50 & name=='forecast_gens'], 
             aes(x=date, y=value, color=age_group), shape=18, size=3, alpha=0.8)+
  geom_vline(xintercept = 49)

fit$summary_pars


summary_pars = data.table(fit$summary_pars)
