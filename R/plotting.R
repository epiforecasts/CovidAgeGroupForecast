library(ggplot2)
library(tidyverse)
library(data.table)

plot_trajectories = function(summary_preds, d){
  ggplot()+
    geom_rect(data=summary_preds[name == 'infections' & forecast_date %in% d,], aes(xmin=date-0.2, xmax=date+0.2, ymin=`5%`, ymax=`95%`), fill='blue', alpha=0.3)+
    geom_point(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(x=date, y=value.x), size=0.01)+
    geom_point(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(x=date, y=value.y), size=0.01)+
    geom_rect(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(xmin=date-0.2, xmax=date+0.2, ymin=`5%`, ymax=`95%`), fill='limegreen', alpha=0.5)+
    geom_rect(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(xmin=date-0.2, xmax=date+0.2, ymin=`5%`, ymax=`95%`, fill=as.character(run)), alpha=0.2)+
    geom_point(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(x=date, y=`50%`), color='green', size=0.01)+
    geom_point(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(x=date, y=`50%`, color=as.character(run)), size=0.01)+
    
    facet_grid(age_index~forecast_date, scales='free', labeller =labeller(age_index=age_labs))+
    scale_x_date(name='')+
    scale_y_continuous(name='infections')+
    scale_color_discrete()+
    theme_minimal()
}


plot_trajectories_one_ax = function(summary_preds, d){
  ggplot()+
    #geom_rect(data=summary_preds[name == 'infections' & forecast_date %in% d,], aes(xmin=date-0.2, xmax=date+0.2, ymin=`5%`, ymax=`95%`), fill='blue', alpha=0.3)+
    geom_point(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(x=date, y=value), size=0.01)+
    geom_point(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(x=date, y=value), size=0.01)+
    #geom_rect(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(xmin=date-0.2, xmax=date+0.2, ymin=`5%`, ymax=`95%`), fill='limegreen', alpha=0.5)+
    geom_rect(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(xmin=date-0.2, xmax=date+0.2, ymin=`5%`, ymax=`95%`, fill=as.character(run)), alpha=0.2)+
    #geom_point(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(x=date, y=`50%`), color='green', size=0.01)+
    geom_point(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(x=date, y=`50%`, color=as.character(run)), size=0.01)+
    
    facet_wrap(~age_index, scales='free', labeller =labeller(age_index=age_labs), ncol=1)+
    scale_x_date(name='')+
    scale_y_continuous(name='infections')+
    scale_color_discrete(name='model', labels=c('Full contact data', 'Age-group means', 'Overall means', 'No contact data'))+
    scale_fill_discrete(name='model', labels=c('Full contact data', 'Age-group means', 'Overall means', 'No contact data'))+
    theme_minimal()
}

plot_parameters = function(summary_pars, d){
  ggplot(data.table(summary_pars[forecast_date %in% d,])[!is.na(index)])+
    geom_rect(aes(xmin=index-0.4 + as.numeric(run)*0.15, xmax=index - 0.3 + as.numeric(run)*0.15, ymin=`5%`, ymax=`95%`, fill=as.character(run)),alpha=0.5)+
    geom_point(aes(x=index-0.35 + as.numeric(run)*0.15, y=`50%`, color=as.character(run)))+
    facet_grid(name~forecast_date)+
    scale_x_continuous(name='age group', breaks = 1:7, labels = age_groups)+
    scale_fill_discrete(name='')+
    scale_color_discrete(name='')+
    theme_minimal()
  
}

