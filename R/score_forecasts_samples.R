library(cowplot)
library(magrittr)
library(scoringutils)
library(ggplot2)
library(data.table)
library(lubridate)
library(khroma)
library(stringr)
library(ggbump)


# set dates and names of key pandemic periods
pandemic_periods = data.table(
  
  periods = 
    c(
      'Lockdown 2', 
      'Lockdown 2 \neasing', 
      'Christmas + \nLockdown 3',
      'Lockdown 3 \neasing', 
      'Opening Up', 
      'END'), 
  
  date_breaks = 
    c(
      lubridate::ymd('2020-11-05'),
      lubridate::ymd('2020-12-03'),
      lubridate::ymd('2020-12-20'),
      lubridate::ymd('2021-03-09'), 
      lubridate::ymd('2021-10-01'), 
      lubridate::ymd('2022-01-20')),
  
  start_date = 
    c(
      lubridate::ymd('2020-11-05'),
      lubridate::ymd('2020-12-03'),
      lubridate::ymd('2020-12-20'),
      lubridate::ymd('2021-03-09'), 
      lubridate::ymd('2021-10-01'), 
      NA),
  
  end_date = 
    c(
      
      lubridate::ymd('2020-12-03'),
      lubridate::ymd('2020-12-20'),
      lubridate::ymd('2021-03-09'), 
      lubridate::ymd('2021-10-01'), 
      lubridate::ymd('2022-01-20'),
      NA)
    
  
)


# function to score forecasts based on single 'true' time series of infections
score_forecasts = function(samples_preds, pandemic_periods, suffix='', age_groups,
                           labels=c('Full contact data', 
                                    'Age-group means', 
                                    'Overall means', 
                                    'No contact data',  
                                    'Baseline - last generation value', 
                                    'Baseline - linear extrapolation')) {
  
  age_labs = age_groups
  names(age_labs) = 1:length(age_groups)
  
  preds_long = samples_preds[name=='forecast_gens'& date > forecast_date,]
  # manipulate outputs into correct format
  preds_long[, horizon := as.numeric(date - forecast_date)]
  preds_long[, model := as.character(run)]
  preds_to_score = preds_long[,c('sample', 'age_group', 'age_index', 'date', 'prediction', 'true_value', 'model', 'forecast_date', 'horizon')]

  preds_to_score[, horizon_long := paste0(horizon/7, ' week forecast')]
  
 
  
  
  # create column to indicate pandemic periods 
  preds_to_score[, periods:=cut(preds_to_score$forecast_date,
                                breaks = pandemic_periods$date_breaks, 
                                labels = head(pandemic_periods$periods,-1))]
  
  true_values_frame = unique(preds_to_score[sample==1,c('date', 'true_value', 'model', 'age_index', 'periods')])
  preds_ranges = scoringutils:::sample_to_range_long(preds_to_score, range = c(0,50,75), keep_quantile_col = FALSE )
  
  preds_ranges = dcast(preds_ranges, formula = age_group + age_index + date + true_value + model + forecast_date +  horizon + horizon_long +  periods  + range ~ boundary, value.var = 'prediction')  
  preds_ranges[, model := as.character(model)]

  vibrant = colour('vibrant')
  
  
 # plot raw predictions from the forecasts - for supplements - PLOT EVERYTHING
  plot_predictions = function(additional_models = c(), horizons = 7 * 1:4){
    
    preds_ranges= preds_ranges[horizon %in% horizons,]
    
    models_to_drop = setdiff(c(4,5,6), additional_models)
    pred_plot = 
    #plot_predictions(preds_to_score[model!='baseline_linex_lv'], by = c('age_group', 'horizon')) + 
    ggplot()+
    
    geom_point(data=preds_ranges[range==0 & !(model %in% models_to_drop),], aes(x=date, y=true_value), alpha=0.8, size=1)+
    geom_point(data=true_values_frame, aes(date, true_value), shape=21, fill=NA, color='black', alpha=0.2, size=1)+
    facet_grid(horizon_long~age_index, scale='free_y', labeller = labeller(age_index = age_labs))+
    geom_rect(data=preds_ranges[range!=0 & !(model %in% models_to_drop),], aes(xmin=date-7, xmax=date+7, ymin=lower, ymax=upper, fill=model, alpha=range))+
    geom_point(data=preds_ranges[range==0 & !(model %in% models_to_drop),], aes(x=date, y=lower, color=model), alpha=0.4)+
    scale_alpha_binned(breaks=c(0.0, 0.5, 0.9, 1.0), range=c(0.5, 0.2))+
    scale_color_manual(name='Model', labels=labels[sort(c(c(1,5,6), (additional_models-2)))]  , values=as.vector(vibrant(6))[sort(c(c(1,5,6), (additional_models-2)))]  )+
    scale_fill_manual(name='Model', labels=labels[sort(c(c(1,5,6), (additional_models-2)))]  , values=as.vector(vibrant(6))[sort(c(c(1,5,6), (additional_models-2)))]  )+
    scale_x_date(date_breaks = "3 months", 
                 labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                              paste(month(x, label = TRUE), "\n", year(x)), 
                                              paste(month(x, label = TRUE))))+
    scale_y_continuous(trans='log10')+
    xlab('')+
    ylab(str_to_title(substr(suffix, 2, nchar(suffix))))+
    theme_minimal_hgrid()+
    theme(
      plot.background = element_rect(fill = 'white'),
      legend.position = 'bottom',
      axis.text.x = element_text(size=9),
      legend.box = 'vertical'
    )+
    ggtitle('A')
    
    return(pred_plot)
    
  }
  
  # plot raw predictions from the forecasts - for main text
  
  plot_predictions_port = function(additional_models = c(), additional_drop = c(), horizons = 7 * 1:4, period_include=pandemic_periods$periods, plot_labs = c()){
    if(length(plot_labs) == 0 ){
      plot_labs=labels[sort(c(c(1,5,6), (additional_models-2)))]
    }
    labels=labels[sort(c(c(1,5,6), (additional_models-2)))]
    preds_ranges= preds_ranges[(horizon %in% horizons) & (periods %in% period_include),]
    
    
    models_to_drop = c(setdiff(c(4,5,6), additional_models), additional_drop)
    
    
    pred_plot = 
      #plot_predictions(preds_to_score[model!='baseline_linex_lv'], by = c('age_group', 'horizon')) + 
      ggplot()+
      
      geom_rect(data=pandemic_periods[periods %in% period_include, `:=` (plot_start =  max(start_date, min(preds_ranges$forecast_date)), plot_end = min(end_date, max(preds_ranges$date))) ], 
                aes(xmin=start_date, 
                    ymin=0, xmax=end_date, 
                    ymax=200, fill=periods), 
                alpha=0.4)+
      scale_fill_discrete(name='Period')+
      
      ggnewscale::new_scale_fill()+
      
      
      geom_point(data=preds_ranges[range==0 & !(model %in% models_to_drop),], aes(x=date, y=true_value), alpha=0.8, size=1)+
      geom_point(data=true_values_frame[(periods %in% period_include),], aes(date, true_value), shape=21, fill=NA, color='black', alpha=0.2, size=1)+
      facet_grid(age_index~horizon_long, scale='free_y', labeller = labeller(age_index = age_labs))+
      geom_rect(data=preds_ranges[range!=0 & !(model %in% models_to_drop),], aes(xmin=date-7, xmax=date+7, ymin=lower, ymax=upper, fill=model, alpha=range))+
      geom_point(data=preds_ranges[range==0 & !(model %in% models_to_drop),], aes(x=date, y=lower, color=model), alpha=0.4)+
      scale_alpha_binned(breaks=c(0.0, 0.5, 0.9, 1.0), range=c(0.5, 0.2))+
      scale_color_manual(name='Model', labels=plot_labs  , values=as.vector(vibrant(6))[sort(c(c(1,5,6), (additional_models-2)))]  )+
      scale_fill_manual(name='Model', labels=plot_labs , values=as.vector(vibrant(6))[sort(c(c(1,5,6), (additional_models-2)))]  )+
      scale_x_date(date_breaks = "3 months", 
                   labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                                paste(month(x, label = TRUE), "\n", year(x)), 
                                                paste(month(x, label = TRUE))))+
      scale_y_continuous(trans='log10', labels=scales::comma)+
      
      

      xlab('')+
      ylab(str_to_title(substr(suffix, 2, nchar(suffix))))+
      theme_minimal_hgrid()+
      theme(
        plot.background = element_rect(fill = 'white'),
        legend.position = 'bottom',
        axis.text.x = element_text(size=9), 
        legend.box = 'virtical')+
      ggtitle('A')
    
    return(pred_plot)
    
  }
  
  # plots of all models for main text
  
  
  full_model_alone_plot_2 = plot_predictions_port(additional_models = c(4,5,6), 
                                                  additional_drop = c("baseline_last_val", "baseline_linex_lv", "baseline_expex_lv"), 
                                                  horizons = c(7, 28), 
                                                  plot_labs = labels[c(1,2,3,4)])
  
 
  # plot raw predictions from the forecasts for each model against baselines 
  
  full_model_alone_plot = plot_predictions()
  ggsave(paste0('plots/prediction_plot', suffix, '.png'), width = 15, height = 7, units = 'in')
  
  full_model_alone_plot_no_data = plot_predictions(c(4))
  ggsave(paste0('plots/prediction_plot', suffix, '_no_data.png'), width = 15, height = 7, units = 'in')
  
  full_model_alone_plot_no_intr = plot_predictions(c(5))
  ggsave(paste0('plots/prediction_plot', suffix, '_no_intr.png'), width = 15, height = 7, units = 'in')
  
  full_model_alone_plot_no_intr = plot_predictions(c(6))
  ggsave(paste0('plots/prediction_plot', suffix, '_polymod.png'), width = 15, height = 7, units = 'in')

  # SCORE FORECASTS USING SCORINGUTILS
  
  scores = score(preds_to_score)
  
  ##################################
  
  
  
  # Convert samples to quntiles for coverage calculations
  quants_to_score = sample_to_quantile(preds_to_score, quantiles = seq(0.05, 0.95, 0.05), type = 7)
  
  scores_quant <- score(quants_to_score)
  scores_quant <- summarise_scores(scores_quant, by = c("model", "quantile", "horizon"))
  
  scores_quant = scores_quant[,model:=as.character(model)]
  
  scale_color_manual(name='Model', labels=labels , values=as.vector(vibrant(6)))
    
  # Coverage plots 
  covplt1 = plot_quantile_coverage(scores_quant[horizon == 7  & !(model %in% c("baseline_last_val", "baseline_linex_lv", "baseline_expex_lv")), ]) + ggtitle('1 week' )  + scale_color_manual(name='Model', labels=labels , values=as.vector(vibrant(6)), guide='none') + theme_minimal_grid()
  covplt2 = plot_quantile_coverage(scores_quant[horizon == 14 & !(model %in% c("baseline_last_val", "baseline_linex_lv", "baseline_expex_lv")), ]) + ggtitle('2 weeks')  + scale_color_manual(name='Model', labels=labels , values=as.vector(vibrant(6)), guide='none') + theme_minimal_grid() 
  covplt3 = plot_quantile_coverage(scores_quant[horizon == 21 & !(model %in% c("baseline_last_val", "baseline_linex_lv", "baseline_expex_lv")), ]) + ggtitle('3 weeks')  + scale_color_manual(name='Model', labels=labels , values=as.vector(vibrant(6)), guide='none') + theme_minimal_grid()
  covplt4 = plot_quantile_coverage(scores_quant[horizon == 28 & !(model %in% c("baseline_last_val", "baseline_linex_lv", "baseline_expex_lv")), ]) + ggtitle('4 weeks')  + scale_color_manual(name='Model', labels=labels , values=as.vector(vibrant(6)), guide='none') + theme_minimal_grid()
  
  coverage_plot = (covplt1 + covplt2)/(covplt3 + covplt4) +  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  
  coverage_plot_14 = (covplt1 / covplt4) & theme(legend.position = 'none')
  
  # Calculate the proportion of true values in forecast ranges
  wide_quantrange = dcast(quants_to_score, formula = age_group + forecast_date + horizon + periods + true_value + model ~ quantile, value.var = c('prediction'))

  prop_in_range = data.table()
  
  wide_quantrange[, in_range := (true_value >= `0.05` & true_value <= `0.95`)]
  wide_quantrange[, prop_in_range := sum(in_range)/.N, by=c('model', 'horizon', 'periods')]
  wide_quantrange[, range := '90%']
  wide_quantrange[, range_num := 90]
  
  
  prop_in_range = rbind(prop_in_range, unique(wide_quantrange[,c('age_group','horizon','model','prop_in_range', 'range', 'range_num', 'periods')]))
  
  # wide_quantrange[, in_range := (true_value >= `0.15` & true_value <= `0.85`)]
  # wide_quantrange[, prop_in_range := sum(in_range)/.N, by=c('model', 'horizon', 'periods')]
  # wide_quantrange[, range := '60%']
  # wide_quantrange[, range_num := 60]
  # 
  # prop_in_range = rbind(prop_in_range, unique(wide_quantrange[,c('age_group','horizon','model','prop_in_range', 'range', 'range_num', 'periods')]))
  
  wide_quantrange[, in_range := (true_value >= `0.25` & true_value <= `0.75`)]
  wide_quantrange[, prop_in_range := sum(in_range)/.N, by=c('model', 'horizon', 'periods')]
  wide_quantrange[, range := '50%']
  wide_quantrange[, range_num := 50]
  
  prop_in_range = rbind(prop_in_range, unique(wide_quantrange[,c('age_group','horizon','model','prop_in_range', 'range', 'range_num', 'periods')]))
  
  # plot range plots - All periods seperately for SM
  range_plot = 
    ggplot(prop_in_range[!(model %in% c("baseline_last_val", "baseline_linex_lv", "baseline_expex_lv")) & !is.na(periods)]) + 
    geom_hline(aes(yintercept=range_num/100), color='gray', size=3)+
    geom_point( map=aes(x=horizon/7, y=prop_in_range, color=model))+
    geom_bump( map=aes(x=horizon/7, y=prop_in_range, color=model))+
    scale_color_manual(name='Model', labels=labels , values=as.vector(vibrant(6)))+
    facet_grid(periods~range)+
    ylim(0,1)+
    ylab('proportion of predictions in range')+
    theme_minimal_grid()
  
  
  ggsave(paste0('plots/coverage', suffix,'.png'), coverage_plot)
  ggsave(paste0('plots/range_cov', suffix,'.png'), range_plot)
  
  # Overall coverage 
  
  prop_in_range_all = data.table()
  wide_quantrange[, in_range := (true_value >= `0.25` & true_value <= `0.75`)]
  wide_quantrange[, prop_in_range := sum(in_range)/.N, by=c('model', 'horizon')]
  wide_quantrange[, range := '50%']
  wide_quantrange[, range_num := 50]
  
  prop_in_range_all = rbind(prop_in_range_all, unique(wide_quantrange[,c('age_group','horizon','model','prop_in_range', 'range', 'range_num', 'periods')]))
  
  
  wide_quantrange[, in_range := (true_value >= `0.05` & true_value <= `0.95`)]
  wide_quantrange[, prop_in_range := sum(in_range)/.N, by=c('model', 'horizon')]
  wide_quantrange[, range := '90%']
  wide_quantrange[, range_num := 90]
  
  prop_in_range_all = rbind(prop_in_range_all, unique(wide_quantrange[,c('age_group','horizon','model','prop_in_range', 'range', 'range_num', 'periods')]))
  
  # plot range plots - overall - main text
  range_plot_all = 
    ggplot(prop_in_range_all[!(model %in% c("baseline_last_val", "baseline_linex_lv", "baseline_expex_lv")) & !is.na(periods)]) + 
    geom_hline(aes(yintercept=range_num/100), color='gray', size=3)+
    geom_point( map=aes(x=horizon/7, y=prop_in_range, color=model))+
    geom_bump( map=aes(x=horizon/7, y=prop_in_range, color=model))+
    scale_color_manual(name='Model', labels=labels , values=as.vector(vibrant(6)))+
    facet_wrap(~range, ncol=2)+
    ylim(0,1)+
    ylab('proportion of predictions in range')+
    theme_minimal_grid()
  
  

  ggsave(paste0('plots/range_cov_all', suffix,'.png'), range_plot_all, height=5, width=10)
  
  
  

  # aggregate scores by age, forecast date and period + some more detailed plots - for posterity 
  scores[, model := as.character(model)]
  
  # calculate scores aggregated by age and pandemic period
  score_by_age =  scores %>%
    summarise_scores(by = c("model", "age_index", "periods")) 
  
  score_by_fd_hz = scores %>%
    summarise_scores(by = c("model", "forecast_date", "horizon")) 
  
  score_by_age_fd_hz = scores %>%
    summarise_scores(by = c("model", "age_index", "forecast_date", "horizon")) 

  
  
  
  
  # calculate scores aggregated by forecast date
  score_by_fd = scores %>%
    summarise_scores(by = c("model", "forecast_date")) 
  
  score_by_fd_long = melt(score_by_fd, id.vars = c('model', 'forecast_date'), measure.vars = c('crps', 'ae_median'), value.name = 'score', variable.name = 'score_type')
  score_by_age_long = melt(score_by_age, id.vars = c('model', 'age_index', 'periods'), measure.vars = c('crps', 'ae_median'), value.name = 'score', variable.name = 'score_type')
  score_by_fd_hz_long = melt(score_by_fd_hz, id.vars = c('model', 'forecast_date', 'horizon'), measure.vars = c('crps', 'ae_median'), value.name = 'score', variable.name = 'score_type')
  
  
  score_by_age_fd_hz_long = melt(score_by_age_fd_hz, id.vars = c('model', 'forecast_date', 'age_index', 'horizon'), measure.vars = c('crps', 'ae_median'), value.name = 'score', variable.name = 'score_type')
  
  
  
  score_by_fd_hz_rel = merge(score_by_fd_hz, score_by_fd_hz[model=='5', ], by = c('forecast_date', 'horizon'), suffixes = c('', '_baseline'))
  score_by_fd_hz_rel[,
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
  
  score_by_age_fd_hz_rel = merge(score_by_age_fd_hz, score_by_age_fd_hz[model=='5', ], by = c('forecast_date', 'horizon', 'age_index'), suffixes = c('', '_baseline'))
  score_by_age_fd_hz_rel[,
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
  
  
  score_by_fd_hz_long_rel = melt(score_by_fd_hz_rel, id.vars = c('model', 'forecast_date', 'horizon'), measure.vars = c('crps_rel', 'ae_median_rel'), value.name = 'score', variable.name = 'score_type')
  
  hists = ggplot(score_by_fd_hz_long_rel[!(model %in% c(2, 3)) & !(model %in% c('baseline_linex_lv', 'baseline_expex_lv')) & score_type == 'crps_rel',]) + 
    geom_histogram(aes(x=score, fill=model), color='white', alpha=0.8) + 
    facet_grid(model~horizon)+
    theme_minimal_vgrid()
  dists = ggplot(score_by_fd_hz_long_rel[!(model %in% c(2, 3)) & !(model %in% c('baseline_linex_lv', 'baseline_expex_lv')) & score_type == 'crps_rel',]) + 
    geom_density(aes(x=score, fill=model), color='white', alpha=0.3) + 
    facet_wrap(~horizon, ncol=1)+
    xlim(0,2)+
    theme_minimal_vgrid()+
    theme(
      legend.position = 'bottom'
    )
  
  score_by_age_fd_hz_long_rel = melt(score_by_age_fd_hz_rel, id.vars = c('model', 'age_index',  'forecast_date', 'horizon'), measure.vars = c('crps_rel', 'ae_median_rel'), value.name = 'score', variable.name = 'score_type')
  
  
  score_by_fd_hz_long_rel[, horizon_long := paste0(horizon/7, ' week forecast')]
  
  time_series_score_age = 
    ggplot(score_by_age_fd_hz_long_rel[horizon %in% c(7, 28) & !(model %in% c(2,3,4,5,6)) & score_type == 'crps_rel']) + 
    geom_point(aes(x=forecast_date, y=score, color=model, shape=as.character(horizon)))+
    geom_line(aes(x=forecast_date, y=score, color=model), linetype='dashed')+
    facet_wrap(~score_type, ncol=1, scale='free_y') + #, labeller = labeller(score_type=setNames(c('Interval Score',  'Bias', "Absolute error of the mean"),c('interval_score', 'bias', 'aem'))))+
    scale_color_discrete(name='model', labels=labels[c(1,5,6)])+
    scale_x_date(name='')+
    facet_wrap(~age_index, nrow=1, labeller = labeller(age_index = age_labs))+
    scale_y_continuous(trans = 'log2')+
    theme_minimal_hgrid()+
    theme(
      axis.text.x = element_text(size=9)
    )
  
  time_series_score = 
    ggplot(score_by_fd_hz_long_rel[!(model %in% c(2,3,4,5,6)) & score_type == 'crps_rel']) + 
    geom_point(aes(x=forecast_date, y=score, color=model, alpha=0.6))+
    geom_line(aes(x=forecast_date, y=score, color=model), linetype='dashed', alpha=0.8)+
    scale_color_manual(name='Model', labels=labels[c(1,5,6)]  , values=as.vector(vibrant(6))[c(1,5,6)])+
    scale_x_date(name='', date_labels = '%b')+
    facet_wrap(~horizon, ncol=1, labeller = labeller(age_index = age_labs),strip.position = 'right')+
    scale_y_continuous(trans = 'log2')+
    theme_minimal_hgrid()+
    theme(legend.position = 'none') + 
    ylab('CRPS')+
    ggtitle('B')+
    theme(
      axis.text.x = element_text(size=9)
    )
  
  
  preds_and_scores = full_model_alone_plot + time_series_score + plot_layout(widths = c(4, 1))
  
  
  ggsave(paste0('plots/preds_and_scores_', suffix, '.png'), preds_and_scores,  width = 15, height = 10, units = 'in')
  
  
  time_series_score_2 = 
    ggplot(score_by_fd_hz_long_rel[!(model %in% c("baseline_last_val", "baseline_linex_lv", "baseline_expex_lv")) & score_type == 'crps_rel' & horizon %in% c(7,28)]) + 
    geom_point(aes(x=forecast_date, y=score, color=model, alpha=0.6))+
    geom_bump(aes(x=forecast_date, y=score, color=model), alpha=0.8)+
    scale_color_manual(name='Model', labels=labels[c(1,2,3,4)]  , values=as.vector(vibrant(6))[c(1,2,3,4)])+
    scale_x_date(name='', date_labels = '%b')+
    facet_wrap(~horizon_long, ncol=2, labeller = labeller(age_index = age_labs),strip.position = 'right')+
    scale_y_continuous(trans = 'log2')+
    theme_minimal_hgrid()+
    theme(legend.position = 'none') + 
    ylab('CRPS')+
    ggtitle('B')+
    theme(
      axis.text.x = element_text(size=9)
    )
  
  
  preds_and_scores = full_model_alone_plot_2 + time_series_score_2 + plot_layout(widths = c(4, 1))
  
  
  ggsave(paste0('plots/preds_and_scores_2_', suffix, '.png'), preds_and_scores,  width = 15, height = 6, units = 'in')
  
  
  
  
  
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
  score_by_age_rel = merge(score_by_age, score_by_age[model=='5', ], by = c('age_index','periods'), suffixes = c('', '_baseline'))
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
  score_by_fd_rel = merge(score_by_fd, score_by_fd[model=='5', ], by = c('forecast_date'), suffixes = c('', '_baseline'))
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
  score_by_period = scores %>% summarise_scores(by = c("model", "periods", 'horizon'))
  # calculate scores aggregated over entire output
  score_overall = scores %>% summarise_scores(by = c("model", 'horizon'))
  # calculate scores aggregated by age 
  score_age = scores %>% summarise_scores(by = c("model", 'age_group', 'horizon'))
  
  
  score_age_rel = merge(score_age, score_age[model=='5', ], by = c('age_group', 'horizon'), suffixes = c('', '_baseline'))
  score_age_rel = merge(score_age_rel, score_age[model=='5', ], by = c('age_group', 'horizon'), suffixes = c('', '_baselinelv'))
  
  score_age_rel[,
                   c(
                     'crps_rel', 
                     'ae_median_rel',
                     'crps_rellv', 
                     'ae_median_rellv'
                   ) 
                   :=
                     list(
                       crps/crps_baseline,
                       ae_median/ae_median_baseline,
                       crps/crps_baselinelv,
                       ae_median/ae_median_baselinelv
                     )
  ]
  
  
  score_period_rel = merge(score_by_period, score_by_period[model=='5', ], by = c('periods', 'horizon'), suffixes = c('', '_baseline'))
  score_period_rel = merge(score_period_rel, score_by_period[model=='5', ], by = c('periods', 'horizon'), suffixes = c('', '_baselinelv'))
  score_period_rel[,
                c(
                  'crps_rel', 
                  'ae_median_rel',
                  'crps_rellv', 
                  'ae_median_rellv'
                ) 
                :=
                  list(
                    crps/crps_baseline,
                    ae_median/ae_median_baseline,
                    crps/crps_baselinelv,
                    ae_median/ae_median_baselinelv
                  )
  ]
  
  
  write.csv(score_by_period, paste0('outputs/score_by_period', suffix, '.csv'))
  write.csv(score_overall,   paste0('outputs/score_overall', suffix, '.csv'))
  
  summary_scores = list(score_overall=score_overall, score_by_period=score_period_rel, score_by_age=score_age_rel)
  
  st1 = melt(summary_scores[[2]], id.vars = c('model', 'periods', 'horizon'), measure.vars = c('crps', 'ae_median'), variable.name = 'score_type', value.name = 'score')
  overall_table = dcast(st1, formula = periods + score_type + horizon ~ model, value.var = 'score')
  
  st2 = melt(summary_scores[[3]], id.vars = c('model', 'age_group', 'horizon'), measure.vars = c('crps', 'ae_median'), variable.name = 'score_type', value.name = 'score')
  age_table = dcast(st2, formula = age_group + score_type + horizon ~ model, value.var = 'score')
  
  
  
  
  write.csv(overall_table, file = paste0('outputs/overall_results', suffix, '.csv'))
  write.csv(age_table, file =     paste0('outputs/age_results',     suffix, '.csv'))
  
  list(summary_scores, full_model_alone_plot_2, time_series_score_2, coverage_plot_14, range_plot_all)

}

