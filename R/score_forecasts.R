library(cowplot)
pandemic_periods = 
  list(
    reduced_restrictions, lubridate::ymd('2020-07-30'), 
    schools_opened = lubridate::ymd('2020-09-04'),
    Lockdown_2, lubridate::ymd('2020-11-05'),
    L2_Easing, lubridate::ymd('2020-12-03'),
    Christmas, lubridate::ymd('2020-12-20'),
    Lockdown_3, lubridate::ymd('2021-05-01'), 
    L3_Schools_open, lubridate::ymd('2021-03-09'), 
    L3_easing, lubridate::ymd('2021-03-29'), 
    Open, lubridate::ymd('2021-10-01')
  )

pandemic_periods = data.table(
  
  periods = 
    c('reduced_restrictions', 
      'schools_opened', 
      'Lockdown_2', 
      'L2_Easing', 
      'Christmas', 
      'Lockdown_3', 
      'L3_Schools_open',
      'L3_easing', 
      'Open', 
      'END'), 
  
  date_breaks = 
    c(lubridate::ymd('2020-07-30'), 
      lubridate::ymd('2020-09-04'),
      lubridate::ymd('2020-11-05'),
      lubridate::ymd('2020-12-03'),
      lubridate::ymd('2020-12-20'),
      lubridate::ymd('2021-05-01'),
      lubridate::ymd('2021-03-09'), 
      lubridate::ymd('2021-03-29'), 
      lubridate::ymd('2021-10-01'), 
      lubridate::ymd('2022-01-20'))
  
)



score_forecasts = function(summary_preds, pandemic_periods) {
  
  preds_long = melt(summary_preds[name=='forecast_gens'], 
                    id.vars = c('date', 'age_group', 'name', 'age_index', 'time_index', 'age_index', 'forecast_date', 'run', 'value'), 
                    measure.vars = c('5%',  '10%', '25%',     '50%',   '75%', '90%',     '95%'), value.name = 'prediction', variable.name = 'quantile_chr')
  preds_long[, quantile := as.numeric(sub("%", "", quantile_chr, fixed=TRUE))/100]
  preds_long[, true_value := value]
  preds_long[, value_date := forecast_date]
  preds_long[, value_type := name]
  preds_long[, horizon := as.numeric(date - forecast_date)]
  preds_long[, model := run]
  
  preds_to_score = preds_long[,c('value_date', 'value_type', 'age_index', 'age_group', 'true_value', 'model', 'quantile', 'prediction', 'horizon')]
  
  preds_to_score[, periods:=cut(preds_to_score$value_date,
                                breaks = pandemic_periods$date_breaks, 
                                labels = head(pandemic_periods$periods,-1))]
  
  score_by_age = scoringutils::eval_forecasts(preds_to_score, 
                                              summarise_by = c("model", "age_index", "periods"))
  
  score_by_fd = scoringutils::eval_forecasts(preds_to_score, 
                                             summarise_by = c("model", "value_date"))
  
  score_by_fd_long = melt(score_by_fd, id.vars = c('model', 'value_date'), measure.vars = c('interval_score', 'bias', 'aem'), value.name = 'score', variable.name = 'score_type')
  score_by_age_long = melt(score_by_age, id.vars = c('model', 'age_index', 'periods'), measure.vars = c('interval_score', 'bias', 'aem'), value.name = 'score', variable.name = 'score_type')
  
  
  p1 =
    ggplot(score_by_fd_long) + 
    geom_point(aes(x=value_date, y=score, color=model))+
    geom_line(aes(x=value_date, y=score, color=model), linetype='dashed')+
    facet_wrap(~score_type, ncol=1, scale='free_y', labeller = labeller(score_type=setNames(c('Interval Score',  'Bias', "Absolute error of the mean"),c('interval_score', 'bias', 'aem'))))+
    scale_color_discrete(name='model', labels=c('Full contact data', 'Age-group means', 'Overall means', 'No contact data', 'Baseline - linear extrapolation', 'Baseline - last generation value'))+
    scale_x_date(name='')+
    theme_minimal_hgrid()
  
  ggsave('plots/scores_forecast_date.pdf', plot=p1, height=10, width=15)
  
  
  
  
  p2 = 
    ggplot(score_by_age_long) + 
    geom_point(aes(x=age_index, y=score, color=model))+
    geom_line(aes(x=age_index, y=score, color=model), linetype='dashed')+
    facet_grid(score_type~periods, scale='free_y', labeller = labeller(score_type=setNames(c('Interval Score',  'Bias', "Absolute error of the mean"),c('interval_score', 'bias', 'aem'))))+
    scale_color_discrete(name='model', labels=c('Full contact data', 'Age-group means', 'Overall means', 'No contact data', 'Baseline - linear extrapolation', 'Baseline - last generation value'))+
    scale_x_continuous(name='age group', breaks = 1:7, labels = age_groups, )+
    theme_minimal_hgrid()+
    theme(
      axis.text.x=element_text(angle=90)
    )
  
  ggsave('plots/scores_age_period.pdf', plot=p2, height=7, width=22)
  
  
  
  score_by_age_rel = merge(score_by_age, score_by_age[model=='baseline_last_diff', ], by = c('age_index','periods'), suffixes = c('', '_baseline'))
  score_by_age_rel[,
                   c(
                     'interval_score_rel',
                     'bias_rel', 
                     'aem_rel'
                   ) 
                   :=
                     list(
                       interval_score/interval_score_baseline,
                       bias/bias_baseline,
                       aem/aem_baseline
                     )
  ]
  
  score_by_fd_rel = merge(score_by_fd, score_by_fd[model=='baseline_last_diff', ], by = c('value_date'), suffixes = c('', '_baseline'))
  score_by_fd_rel[,
                   c(
                     'interval_score_rel',
                     'bias_rel', 
                     'aem_rel'
                   ) 
                   :=
                     list(
                       interval_score/interval_score_baseline,
                       bias/bias_baseline,
                       aem/aem_baseline
                     )
  ]
  
  
  score_by_fd_long_rel = melt(score_by_fd_rel, id.vars = c('model', 'value_date'), measure.vars = c('interval_score_rel', 'bias_rel', 'aem_rel'), value.name = 'score', variable.name = 'score_type')
  score_by_age_long_rel = melt(score_by_age_rel, id.vars = c('model', 'age_index', 'periods'), measure.vars = c('interval_score_rel', 'bias_rel', 'aem_rel'), value.name = 'score', variable.name = 'score_type')
  
  
  p1 =
    ggplot(score_by_fd_long_rel) + 
    geom_point(aes(x=value_date, y=score, color=model))+
    geom_line(aes(x=value_date, y=score, color=model), linetype='dashed')+
    facet_wrap(~score_type, ncol=1, scale='free_y', labeller = labeller(score_type=setNames(c('Interval Score',  'Bias', "Absolute error of the mean"),c('interval_score_rel', 'bias_rel', 'aem_rel'))))+
    scale_color_discrete(name='model', labels=c('Full contact data', 'Age-group means', 'Overall means', 'No contact data', 'Baseline - linear extrapolation', 'Baseline - last generation value'))+
    scale_x_date(name='')+
    theme_minimal_hgrid()
  
  ggsave('plots/scores_forecast_date_relative.pdf', plot=p1, height=10, width=15)
  
  
  
  
  p2 = 
    ggplot(score_by_age_long_rel) + 
    geom_point(aes(x=age_index, y=score, color=model))+
    geom_line(aes(x=age_index, y=score, color=model), linetype='dashed')+
    facet_grid(score_type~periods, scale='free_y', labeller = labeller(score_type=setNames(c('Interval Score',  'Bias', "Absolute error of the mean"),c('interval_score_rel', 'bias_rel', 'aem_rel'))))+
    scale_color_discrete(name='model', labels=c('Full contact data', 'Age-group means', 'Overall means', 'No contact data', 'Baseline - linear extrapolation', 'Baseline - last generation value'))+
    scale_x_continuous(name='age group', breaks = 1:7, labels = age_groups, )+
    theme_minimal_hgrid()+
    theme(
      axis.text.x=element_text(angle=90)
    )
  
  ggsave('plots/scores_age_period_relative.pdf', plot=p2, height=7, width=22)
  
  

}


score_forecasts(summary_preds[date<lubridate::ymd(20211201)], pandemic_periods)
