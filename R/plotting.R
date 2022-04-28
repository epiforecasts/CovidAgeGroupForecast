library(ggplot2)
library(tidyverse)
library(data.table)
library(viridis)

plot_trajectories = function(summary_preds, d){
  ggplot()+
    geom_rect(data=summary_preds[name == 'infections' & forecast_date %in% d,], aes(xmin=date-0.2, xmax=date+0.2, ymin=`5%`, ymax=`95%`), fill='blue', alpha=0.3)+
    geom_point(data=summary_preds[name == 'next_gens' & forecast_date %in% d,], aes(x=date, y=value), size=0.01)+
    geom_point(data=summary_preds[name == 'forecast_gens' & forecast_date %in% d,], aes(x=date, y=value), size=0.01)+
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
    
    facet_wrap(~age_index, labeller =labeller(age_index=age_labs), ncol=2)+
    scale_x_date(name='')+
    scale_y_continuous(name='Infections')+
    scale_color_discrete(name='model', labels=c('Full contact data', 'Age-group means', 'Overall means', 'No contact data',  'Baseline - linear extrapolation', 'Baseline - last generation value'))+
    scale_fill_discrete(name='model', labels=c('Full contact data', 'Age-group means', 'Overall means', 'No contact data', 'Baseline - linear extrapolation', 'Baseline - last generation value'))+
    theme_minimal_hgrid()+
    theme(legend.position = 'bottom')
  
  #ggsave('plots/trajectories.pdf', width=10, height=7)
  #ggsave('plots/trajectories.png', width=10, height=7, units='in', dpi=300)
}

lab_dates <- pretty(summary_pars$forecast_date)

plot_parameters = function(summary_pars, d){
  ggplot(data.table(summary_pars[forecast_date %in% d,])[!is.na(index)])+
    #geom_rect(aes(xmin=index, xmax=index, ymin=`5%`, ymax=`95%`, fill=forecast_date),alpha=0.5)+
    geom_line(aes(x=index, y=`50%`, color=forecast_date, group=forecast_date), alpha=0.2)+
    geom_point(aes(x=index, y=`50%`, color=forecast_date, group=forecast_date), alpha=0.2)+
    facet_grid(name~run, scale='free_y', labeller=labeller(run= setNames(nm=1:4, object=c('Full \ncontact data', 'Age-group \nmeans', 'Overall \nmeans', 'No \ncontact data')), name = setNames(nm=c('inf_rate', 'susceptibility'), object=c('Inectiousness', 'Susceptibility'))))+
    scale_x_continuous(name='Age group', breaks = 1:7, labels = age_groups, )+
    scale_fill_continuous(name='')+
    scale_color_viridis(name='',breaks = as.numeric(lab_dates), 
                        labels = lab_dates)+
    scale_y_continuous(name='')+
    theme_minimal_hgrid()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = 'bottom', 
          legend.key.width = unit(1, "in"), 
          legend.justification = 'centre')
  
  ggsave('plots/parameters.pdf', width=10, height =5)
  ggsave('plots/parameters.png', width=10, height =5, units='in', dpi=300)
}

