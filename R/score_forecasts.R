
score_forecasts = function(summary_preds) {
  
  preds_long = melt(summary_preds[name=='forecast_gens'], 
                    id.vars = c('date', 'age_group', 'name', 'time_index', 'age_index', 'forecast_date', 'run', 'value'), 
                    measure.vars = c('5%',  '10%', '25%',     '50%',   '75%', '90%',     '95%'), value.name = 'prediction', variable.name = 'quantile_chr')
  
  
  preds_long[, quantile := as.numeric(sub("%", "", quantile_chr, fixed=TRUE))/100]
  
  preds_long[, true_value := value]
  
  preds_long[, value_date := forecast_date]
  preds_long[, value_type := name]
  
  preds_long[, horizon := as.numeric(date - forecast_date)]
  preds_long[, model := run]
  
  
  
  preds_to_score = preds_long[,c('value_date', 'value_type', 'age_group', 'true_value', 'model', 'quantile', 'prediction', 'horizon')]
  
  
  scores_for_table <- scoringutils::eval_forecasts(preds_to_score, 
                                         summarise_by = c("model"))
  scoringutils::score_table(scores_for_table)
  
  
  scores_coverage <- scoringutils::eval_forecasts(preds_to_score, 
                                         summarise_by = c("model", "range", "quantile"))
  scoringutils::interval_coverage(scores_coverage) + 
    ggplot2::ggtitle("Interval Coverage")
  
  scoringutils::quantile_coverage(scores_coverage) + 
    ggplot2::ggtitle("Quantile Coverage")
  
  
  scores_wis <- scoringutils::eval_forecasts(preds_to_score, 
                                         summarise_by = c("model", "range", "age_group"))
  scoringutils::wis_components(scores_wis, facet_formula = ~ age_group)
  scoringutils::range_plot(scores_wis, y = "interval_score", 
                           facet_formula = ~ age_group)
  
  
  scores_by_horizon <- scoringutils::eval_forecasts(preds_to_score, 
                                         summarise_by = c("model", "horizon"))
  scores <- scores[, horizon := as.factor(horizon)]
  scoringutils::score_heatmap(scores_by_horizon, 
                              x = "horizon", metric = "interval_score")
  
  
  
}

