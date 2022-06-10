library(cowplot)
library(magrittr)
library(scoringutils)
library(ggplot2)

# set dates and names of key pandemic periods
pandemic_periods = data.table(
  
  periods = 
    c('Light \nRestrictions', 
      'Schools \nre-opened', 
      'Lockdown 2', 
      'Lockdown 2 \neasing', 
      'Christmas', 
      'Lockdown 3', 
      'Lockdown 3 \nSchools Open',
      'Lockdown 3 \neasing', 
      'Opening Up', 
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


# function to score forecasts based on single 'true' time series of infections
score_forecasts = function(samples_preds, pandemic_periods, suffix='', age_groups, 
                           labels=c('Full contact data', 
                                    'Age-group means', 
                                    'Overall means', 
                                    'No contact data',  
                                    'Baseline - last generation value', 
                                    'Baseline - linear extrapolation')) {
  
  preds_long = samples_preds[name=='forecast_gens'& date > forecast_date,]
  # manipulate outputs into correct format
  preds_long[, horizon := as.numeric(date - forecast_date)]
  preds_long[, model := run]
  preds_to_score = preds_long[,c('sample', 'age_group', 'age_index', 'date', 'prediction', 'true_value', 'model', 'forecast_date', 'horizon')]
  
  # create column to indicate pandemic periods 
  preds_to_score[, periods:=cut(preds_to_score$date,
                                breaks = pandemic_periods$date_breaks, 
                                labels = head(pandemic_periods$periods,-1))]
  
  scores = score(preds_to_score)
  
  scores[, model := as.character(model)]
  
  # calculate scores aggregated by age and pandemic period
  score_by_age =  scores %>%
    summarise_scores(by = c("model", "age_index", "periods")) 
  
  # calculate scores aggregated by forecast date
  score_by_fd = scores %>%
    summarise_scores(by = c("model", "forecast_date")) 
  
  score_by_fd_long = melt(score_by_fd, id.vars = c('model', 'forecast_date'), measure.vars = c('crps', 'ae_median'), value.name = 'score', variable.name = 'score_type')
  score_by_age_long = melt(score_by_age, id.vars = c('model', 'age_index', 'periods'), measure.vars = c('crps', 'ae_median'), value.name = 'score', variable.name = 'score_type')
  
  
  
  
  p1 =
    ggplot(score_by_fd_long) + 
    geom_point(aes(x=forecast_date, y=score, color=model))+
    geom_line(aes(x=forecast_date, y=score, color=model), linetype='dashed')+
    facet_wrap(~score_type, ncol=1, scale='free_y') + #, labeller = labeller(score_type=setNames(c('Interval Score',  'Bias', "Absolute error of the mean"),c('interval_score', 'bias', 'aem'))))+
    scale_color_discrete(name='model', labels=labels)+
    scale_x_date(name='')+
    theme_minimal_hgrid()
  
  ggsave(paste0('plots/scores_forecast_date', suffix, '.pdf'), plot=p1, height=10, width=15)
  ggsave(paste0('plots/scores_forecast_date', suffix, '.png'), plot=p1, height=10, width=15, units = 'in', dpi=300)
  
  
  
  
  
  p2 = 
    ggplot(score_by_age_long) + 
    geom_point(aes(x=age_index, y=score, color=model))+
    geom_line(aes(x=age_index, y=score, color=model), linetype='dashed')+
    facet_grid(score_type~periods, scale='free_y') + #, labeller = labeller(score_type=setNames(c('Interval Score',  'Bias', "Absolute error of the mean"),c('interval_score', 'bias', 'aem'))))+
    scale_color_discrete(name='model', labels=labels)+
    scale_x_continuous(name='age group', breaks = 1:length(age_groups), labels = age_groups, )+
    theme_minimal_hgrid()+
    theme(
      axis.text.x=element_text(angle=90)
    )
  
  ggsave(paste0('plots/scores_age_period', suffix, '.pdf'), plot=p2, height=7, width=22)
  ggsave(paste0('plots/scores_age_period', suffix, '.png'), plot=p2, height=7, width=22, units = 'in', dpi=300)
  
  
  # calculate scores aggregated by age and pandemic period relative to baseline model
  score_by_age_rel = merge(score_by_age, score_by_age[model=='baseline_expex_lv', ], by = c('age_index','periods'), suffixes = c('', '_baseline'))
  score_by_age_rel[,
                   c(
                     'crps_rel', 
                     'ae_median_rel'
                   ) 
                   :=
                     list(
                       crps/crps_baseline,
                       ae_median/ae_median_baseline
                       
                     )
  ]
  
  
  # calculate scores aggregated by forecast date relative to baseline model
  score_by_fd_rel = merge(score_by_fd, score_by_fd[model=='baseline_expex_lv', ], by = c('forecast_date'), suffixes = c('', '_baseline'))
  score_by_fd_rel[,
                   c(
                     'crps_rel', 
                     'ae_median_rel'
                   ) 
                   :=
                     list(
                       crps/crps_baseline,
                       ae_median/ae_median_baseline
                     )
  ]
  
  
  score_by_fd_long_rel = melt(score_by_fd_rel, id.vars = c('model', 'forecast_date'), measure.vars = c('crps_rel', 'ae_median_rel'), value.name = 'score', variable.name = 'score_type')
  score_by_age_long_rel = melt(score_by_age_rel, id.vars = c('model', 'age_index', 'periods'), measure.vars = c( 'ae_median_rel', 'crps_rel' ), value.name = 'score', variable.name = 'score_type')
  
  
  p1 =
    ggplot(score_by_fd_long_rel) + 
    geom_point(aes(x=forecast_date, y=score, color=model))+
    geom_line(aes(x=forecast_date, y=score, color=model), linetype='dashed')+
    facet_wrap(~score_type, ncol=1, scale='free_y') + #, labeller = labeller(score_type=setNames(c('Interval Score', "Absolute error of the mean",  'Bias'),c('interval_score_rel', 'aem_rel', 'bias'))))+
    scale_color_discrete(name='', labels=labels)+
    scale_y_continuous(trans='log2')+
    theme_minimal_hgrid()+
    theme(legend.position = 'bottom')
  
  
  ggsave(paste0('plots/scores_forecast_date_relative', suffix, '.pdf'), plot=p1, height=7, width=10)
  ggsave(paste0('plots/scores_forecast_date_relative', suffix, '.png'), plot=p1, height=7, width=10, units = 'in', dpi=300) 
  
  
  
  
  p2 = 
    ggplot(score_by_age_long_rel) + 
    geom_point(aes(y=age_index, x=score, color=model))+
    geom_path(aes(y=age_index, x=score, color=model), linetype='dashed')+
    facet_grid(periods~score_type, scale='free_x')+#, labeller = labeller(score_type=setNames(c('Interval Score', "Absolute error \nof the mean",  'Bias'),c('interval_score_rel',  'aem_rel', 'bias'))))+
    scale_color_discrete(name='', labels=labels)+
    scale_y_continuous(name='age group', breaks = 1:length(age_groups), labels = age_groups, trans='reverse')+
    scale_x_continuous(trans='log2')+
    theme_minimal_vgrid()+
    theme(
      legend.position = 'bottom'
    )
  
  ggsave(paste0('plots/scores_age_period_relative', suffix, '.pdf'), plot=p2, height=15, width=10)
  ggsave(paste0('plots/scores_age_period_relative', suffix, '.png'), plot=p2, height=15, width=10, units = 'in', dpi=300) 
  
  
  # calculate scores aggregated by pandemic period
  score_by_period = scores %>% summarise_scores(by = c("model", "periods"))
  # calculate scores aggregated over entire output
  score_overall = scores %>% summarise_scores(by = c("model"))
  # calculate scores aggregated by age 
  score_age = scores %>% summarise_scores(by = c("model", 'age_group'))
  
  
  score_age_rel = merge(score_age, score_age[model=='baseline_expex_lv', ], by = c('age_group'), suffixes = c('', '_baseline'))
  score_age_rel[,
                   c(
                     'crps_rel', 
                     'ae_median_rel'
                   ) 
                   :=
                     list(
                       crps/crps_baseline,
                       ae_median/ae_median_baseline
                     )
  ]
  
  
  score_period_rel = merge(score_by_period, score_by_period[model=='baseline_expex_lv', ], by = c('periods'), suffixes = c('', '_baseline'))
  score_period_rel[,
                c(
                  'crps_rel', 
                  'ae_median_rel'
                ) 
                :=
                  list(
                    crps/crps_baseline,
                    ae_median/ae_median_baseline
                  )
  ]
  
  
  write.csv(score_by_period, paste0('outputs/score_by_period', suffix, '.csv'))
  write.csv(score_overall,   paste0('outputs/score_overall', suffix, '.csv'))
  
  summary_scores = list(score_overall=score_overall, score_by_period=score_period_rel, score_by_age=score_age_rel)
  
  st1 = melt(summary_scores[[2]], id.vars = c('model', 'periods'), measure.vars = c('crps', 'ae_median'), variable.name = 'score_type', value.name = 'score')
  overall_table = dcast(st1, formula = periods + score_type ~ model, value.var = 'score')
  
  st2 = melt(summary_scores[[3]], id.vars = c('model', 'age_group'), measure.vars = c('crps', 'ae_median'), variable.name = 'score_type', value.name = 'score')
  age_table = dcast(st2, formula = age_group + score_type ~ model, value.var = 'score')
  
  
  
  
  write.csv(overall_table, file = paste0('outputs/overall_results', suffix, '.csv'))
  write.csv(age_table, file =     paste0('outputs/age_results',     suffix, '.csv'))
  
  summary_scores

}

# score based on summary preds (output needed before run)


# create tables of scores
