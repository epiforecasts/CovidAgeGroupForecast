library(data.table)
library(ggplot2)
library(patchwork)
library(khroma)
library(ggnewscale)
source('R/score_forecasts_samples.R')
source('R/plotting.R')


summary_preds_inf = readRDS('outputs/samples_preds_infweek.rds')
summary_preds_cas = readRDS('outputs/samples_preds_cases.rds')

inf_traj = readRDS('inf_traj.rds') 
anb_traj = readRDS('anb_traj.rds')
cas_traj = readRDS('cas_traj.rds')

inf_traj/anb_traj/cas_traj

ggsave('plots/inputs.pdf',  width=15, height=7)
ggsave('plots/inputs.png',  width=15, height=7, units = 'in', dpi=300)



age_labs_inf_reverse = 1:7
names(age_labs_inf_reverse) = c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )

age_labs_inf_forward = names(age_labs_inf_reverse)
names(age_labs_inf_forward) = age_labs_inf_reverse



outs_inf = score_forecasts(summary_preds_inf[date > as.Date('2020-10-01') & !(run %in% c(2,3) | run=='baseline_linex_lv')], 
                                 pandemic_periods,
                                 suffix = '_infections', 
                                 age_groups =  c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' ),
                                 labels = c('CoMix contact data', 
                                            'No contact data', 
                                            'No interaction',
                                            'Polymod contact data', 
                                            'Exponential baseline', 
                                            'Fixed value baseline'))


summary_scores_inf = outs_inf[[1]]

preds_plot_inf = outs_inf[[2]]
scores_plot_inf = outs_inf[[3]]
cov_plot_inf = outs_inf[[4]]
range_plot_inf = outs_inf[[5]]


overall_scores_inf = summary_scores_inf$score_overall[,c('model', 'horizon', 'crps', 'ae_median', 'bias')]

#overall_scores_inf[, model := c('CoMix contact data', 
#                                'No contact data', 
#                                'No interaction',
#                                'Polymod contact data', 
#                                'Exponential baseline', 
#                                'Fixed value baseline')]
#
ggplot(overall_scores_inf) + 
  geom_point(aes(x=model, y=crps, color=horizon))


write.csv(dcast(overall_scores_inf, formula = horizon ~ model, value.var = 'crps'), file = 'outputs/overall_infections_table.csv')
  

age_scores_inf = summary_scores_inf$score_by_age[,c('model', 'age_group', 'horizon', 'crps_rel', 'crps_rellv', 'ae_median_rel', 'ae_median_rellv', 'bias')]



age_scores_inf[, age_index := age_labs_inf_reverse[age_group]]

vibrant = colour('vibrant')



age_inf = ggplot(age_scores_inf[model != 'baseline_linex_lv',]) +
  #geom_segment(aes(x=horizon, xend=horizon, y=1, yend=crps_rel), alpha=0.3)+
  geom_point(aes(x=horizon, y=crps_rel, color=model), alpha=0.8)+
  geom_bump(aes(x=horizon, y=crps_rel, color=model), alpha=0.8)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(trans='log2', name='CRPS')+
  scale_x_discrete(labels=NULL, name='')+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)))+
  facet_wrap(~age_index, nrow=1, labeller = labeller(age_index = age_labs_inf_forward))+

  ggtitle("B")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'none'
  )

age_inf_lv = ggplot(age_scores_inf[model != 'baseline_linex_lv',]) +
  #geom_segment(aes(x=horizon, xend=horizon, y=1, yend=crps_rel), alpha=0.3)+
  geom_point(aes(x=horizon, y=crps_rellv, color=model), alpha=0.8)+
  geom_bump(aes(x=horizon, y=crps_rellv, color=model), alpha=0.8)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(trans='log2', name='CRPS')+
  scale_x_discrete(labels=NULL, name='')+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)))+
  facet_wrap(~age_index, nrow=1, labeller = labeller(age_index = age_labs_inf_forward))+
  
  ggtitle("B")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'none'
  )


period_scores_inf = summary_scores_inf$score_by_period[,c('model', 'periods', 'horizon', 'crps_rel', 'ae_median_rel', 'bias')]




period_inf = ggplot(period_scores_inf[model != 'baseline_linex_lv' & periods != 'Schools \nre-opened',]) +
  geom_bump(aes(x=horizon, y=crps_rel, color=model), alpha=0.6)+
  geom_point(aes(x=horizon, y=crps_rel, color=model), alpha=0.6)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(trans='log2', name='CRPS')+
  scale_x_discrete(labels=NULL, name='')+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)))+

  facet_wrap(~periods, nrow=1)+
  ggtitle("A")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'none'
  )


period_inf_bias = ggplot(period_scores_inf[model != 'baseline_linex_lv' & periods != 'Schools \nre-opened',]) +
  geom_path(aes(x=bias, y=crps_rel, color=model), alpha=0.6, linetype='dashed')+
  geom_point(aes(x=bias, y=crps_rel, color=model, size=horizon), alpha=0.6)+
  geom_hline(yintercept=1)+
  geom_vline(xintercept=0)+
  scale_y_continuous(trans='log2', name='')+
  scale_x_discrete(labels=NULL, name='')+
  scale_color_discrete(labels=c('CoMix contact data', 
                                'No contact data', 
                                'No interaction',
                                'Polymod contact data', 
                                'Exponential baseline', 
                                'Fixed value baseline'), 
                       name = 'Model')+
  facet_wrap(~periods, nrow=1)+
  ggtitle("A")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'none'
  )


outs_cas = score_forecasts(summary_preds_cas[date > as.Date('2020-10-01') & !(run %in% c(2,3) | run=='baseline_linex_lv')], 
                                    pandemic_periods, 
                                    suffix = '_cases', 
                                    age_groups = c('0-9',  '10-19',   '20-29',  '30-39',  '40-49',  '50-59',  '60-69', '70+'),
                                    labels = c('CoMix contact data', 
                                               'No contact data', 
                                               'No interaction',
                                               'Polymod contact data', 
                                               'Exponential baseline', 
                                               'Fixed value baseline'))
summary_scores_cas = outs_cas[[1]]


preds_plot_cas = outs_cas[[2]]
scores_plot_cas = outs_cas[[3]]
cov_plot_cas = outs_cas[[4]]
range_plot_cas = outs_cas[[5]]

overall_scores_cas = summary_scores_cas$score_overall[,c('model', 'horizon', 'crps', 'ae_median', 'bias')]

write.csv(dcast(overall_scores_cas, formula = horizon ~ model, value.var = 'crps'), file = 'outputs/overall_cases_table.csv')


age_scores_cas = summary_scores_cas$score_by_age[,c('model', 'age_group', 'horizon', 'crps_rel', 'ae_median_rel','crps_rellv', 'ae_median_rellv', 'bias')]



age_cas = ggplot(age_scores_cas[model != 'baseline_linex_lv',]) +
  #geom_segment(aes(x=horizon, xend=horizon, y=1, yend=crps_rel), alpha=0.3)+
  geom_point(aes(x=horizon, y=crps_rel, color=model), alpha=0.8)+
  geom_bump(aes(x=horizon, y=crps_rel, color=model), alpha=0.8)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(trans='log2', name='CRPS')+
  scale_x_discrete(labels=NULL, name='Horizon')+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)))+

  facet_wrap(~age_group, nrow=1)+
  ggtitle("C")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom'
  )

age_cas_lv = ggplot(age_scores_cas[model != 'baseline_linex_lv',]) +
  #geom_segment(aes(x=horizon, xend=horizon, y=1, yend=crps_rel), alpha=0.3)+
  geom_point(aes(x=horizon, y=crps_rellv, color=model), alpha=0.8)+
  geom_bump(aes(x=horizon, y=crps_rellv, color=model), alpha=0.8)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(trans='log2', name='CRPS')+
  scale_x_discrete(labels=NULL, name='Horizon')+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)))+
  
  facet_wrap(~age_group, nrow=1)+
  ggtitle("C")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom'
  )



period_scores_cas = summary_scores_cas$score_by_period[,c('model', 'periods', 'horizon', 'crps_rel', 'ae_median_rel', 'bias')]


period_cas = ggplot(period_scores_cas[model != 'baseline_linex_lv' & periods != 'Schools \nre-opened',]) +
  geom_bump(aes(x=horizon, y=crps_rel, color=model), alpha=0.6)+
  geom_point(aes(x=horizon, y=crps_rel, color=model), alpha=0.6)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(trans='log2', name='CRPS')+
  scale_x_discrete(name='Horizon')+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)))+

  facet_wrap(~periods, nrow=1)+
  ggtitle("B")+
  theme_minimal()+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom'
  )


period_inf/period_cas

ggsave('plots/periods.png',width = 7, height=5, units='in')

period_all_bias = ggplot(rbind(period_scores_cas[,type:='cases'], period_scores_inf[,type:='infections'])[model != 'baseline_linex_lv' & periods != 'Schools \nre-opened',]) +
  geom_path(aes(x=bias, y=crps_rel, color=model), alpha=0.6, linetype='dashed')+
  geom_point(aes(x=bias, y=crps_rel, color=model, size=horizon), alpha=0.6)+
  geom_hline(yintercept=1)+
  geom_vline(xintercept=0)+
  scale_y_continuous(trans='log2', name='CRPS')+
  scale_x_discrete(labels=NULL, name='Bias')+
  scale_size(breaks = c(0,7,14,21,28))+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)))+
  scale_fill_manual(labels=c('CoMix contact data', 
                             'No contact data', 
                             'No interaction',
                             'Polymod contact data', 
                             'Exponential baseline', 
                             'Fixed value baseline'), 
                    name = 'Model', values = as.vector(vibrant(6)))+

  facet_grid(type~periods, scale='free')+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom'
  )


ggsave('plots/periods_bias.png',width = 7, height=5, units='in')




overall_scores = rbind(overall_scores_inf[, type:='Infections'], overall_scores_cas[, type:='Cases'])

overall_scores = merge(overall_scores, overall_scores[model=='5'], by = c('type', 'horizon'), suffixes = c('','_baseline'))
overall_scores = merge(overall_scores, overall_scores[model=='5'], by = c('type', 'horizon'), suffixes = c('','_baselinelv'))

overall_scores[, crps_rel := crps / crps_baseline]
overall_scores[, crps_rellv := crps / crps_baselinelv]

overall_scores[, ae_median_rel := ae_median / ae_median_baseline]


overall_scores[, c('type' ,'horizon'  ,           'model' ,     'crps', 'ae_median',        'bias',  'crps_rel',  'ae_median_rel')]

write.csv(overall_scores[, c('type' ,'horizon'  ,           'model' ,     'crps', 'ae_median',         'bias',  'crps_rel' )]
, 'outputs/overall_scores_rel.csv')

age_scores = rbind(age_scores_inf[, -c('age_index')][, type:='Infections'], age_scores_cas[, type:='Cases'])


#age_scores[, c('type' ,'horizon'  ,           'model' ,  'age_group',  'crps', 'ae_median',         'bias',  'crps_rel' )]

write.csv(age_scores[, c('type' ,'horizon'  ,           'model' , 'age_group',    'crps_rel'   ,      'bias')]
          , 'outputs/age_scores_rel.csv')

period_scores = rbind(period_scores_inf[, type:='Infections'], period_scores_cas[, type:='Cases'])


#period_scores[, c('type' ,'horizon'  ,           'model' ,  'age_group',  'crps', 'ae_median',         'bias',  'crps_rel' )]

write.csv(period_scores[, c('type' ,'horizon'  ,           'model' , 'periods',   'crps_rel'   ,      'bias')]
          , 'outputs/period_scores_rel.csv')


overall_sp = 
  ggplot(overall_scores) + 
  #geom_curve(aes(x=model, xend=model, y=1, yend=crps_rel+horizon*0.00001, size=exp(horizon/7)), alpha=0.3, curvature = 0.15)+
  geom_path(aes(x=ae_median_rel, y=crps_rel, color=model), alpha=0.7)+
  geom_point(aes(x=ae_median_rel, y=crps_rel, color=model, size=horizon), alpha=0.7)+
  geom_abline(intercept=-1, slope=1, linetype=2)+
  geom_hline(yintercept = 1, alpha=0.5)+
  geom_vline(xintercept = 1, alpha=0.5)+
  facet_wrap(~type)+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)), guide='none')+
  scale_size(breaks = c(0,7,14,21,28))+
  scale_y_continuous(trans='log2', name='CRPS')+
  theme_minimal()+
  ggtitle('A')+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom'
  )

overall_sp_lv = 
  ggplot(overall_scores) + 
  #geom_curve(aes(x=model, xend=model, y=1, yend=crps_rel+horizon*0.00001, size=exp(horizon/7)), alpha=0.3, curvature = 0.15)+
  geom_path(aes(x=bias, y=crps_rellv, color=model), alpha=0.7)+
  geom_point(aes(x=bias, y=crps_rellv, color=model, size=horizon), alpha=0.7)+
  geom_hline(yintercept = 1, alpha=0.5)+
  geom_vline(xintercept = 0, alpha=0.5)+
  facet_wrap(~type)+
  scale_color_manual(labels=c('CoMix contact data', 
                              'No contact data', 
                              'No interaction',
                              'Polymod contact data', 
                              'Exponential baseline', 
                              'Fixed value baseline'), 
                     name = 'Model', values = as.vector(vibrant(6)), guide='none')+
  scale_size(breaks = c(0,7,14,21,28))+
  scale_y_continuous(trans='log2', name='CRPS')+
  theme_minimal()+
  ggtitle('A')+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom'
  )
layout_overall = 
"
AABB
AACC"
overageres = overall_sp + age_inf + age_cas + plot_layout(design = layout_overall)
overageres_lv = overall_sp_lv / age_inf_lv + age_cas_lv + plot_layout(design = layout_overall)


age_inf / age_cas

ggsave('plots/overall_results.pdf', plot = overageres, width = 12, height=7, units = 'in')
ggsave('plots/overall_results.png', plot = overageres, width = 12, height=7, units = 'in', dpi=300)

ggsave('plots/overall_results_lv.pdf', plot = overageres_lv, width = 12, height=7, units = 'in')
ggsave('plots/overall_results_lv.png', plot = overageres_lv, width = 12, height=7, units = 'in', dpi=300)






(((preds_plot_inf + scale_y_continuous(trans = 'log10', limits = c(1e2, 1e6),  labels=scales::comma) + ggtitle('A')) + (preds_plot_cas + scale_y_continuous(trans = 'log10', limits = c(1e2, 1e6),  labels=scales::comma) + ggtitle('B')))     /  
( ((scores_plot_inf + ggtitle('C')) + (scores_plot_cas + ggtitle('D'))))) / 
guide_area()  + plot_layout(guides = 'collect', heights=c(10,2,1))


ggsave('plots/preds_scores_both.png', width=15, height=20, units='in')

layout_covs = 
  "
AABB
CCDD
CCDD
CCDD
EEEE"

layout_covs = 
  "
AAAAACCCCDDDD
AAAAACCCCDDDD
AAAAACCCCDDDD
BBBBBCCCCDDDD
BBBBBCCCCDDDD
BBBBBCCCCDDDD
EEEEEEEEEEEEE"

cov_plot_inf = cov_plot_inf + plot_annotation(title='C')
cov_plot_cas = cov_plot_cas + plot_annotation(title='D')

(range_plot_inf  + xlab('Horizon (Weeks)') + ylab("Prop. of infection incidence \nin central range of projections") + theme(legend.position = 'none')) + 
  (range_plot_cas  + xlab('Horizon (Weeks)') + ylab("Prop. of case incidence \nin central range of projections") + theme(legend.position = 'none')) +  
  ((cov_plot_inf) & theme(legend.position = 'none')) + (cov_plot_cas & theme(legend.position = 'none'))  + 
  guide_area()  + plot_layout(guides = 'collect', design = layout_covs) + 
  plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')



ggsave('plots/coverage_plot.png', width=11, height=8, units='in')



inf_pars = readRDS('outputs/summary_pars_infweek.rds')
cas_pars = readRDS('outputs/summary_pars_cases.rds')


plot_parameters(inf_pars, dates, '_infections')
plot_parameters(cas_pars, dates, '_cases')


inf_pars_wide = dcast(unique(inf_pars[!is.na(index)]), formula = forecast_date + run + name ~ index , value.var = '50%')
cas_pars_wide = dcast(unique(cas_pars[!is.na(index)]), formula = forecast_date + run + name ~ index , value.var = '50%')

inf_pars_wide[, ad_rat := (`3`     +   `4`     +   `5`)*2 / (`1` + `2`)/3  ]
inf_pars_wide[, el_rat := (`6`     +   `7`  ) / (`1` + `2`) ]

cas_pars_wide[, ad_rat := (`3`     +   `4`     +   `5`)*2 / (`1` + `2`)/3  ]
cas_pars_wide[, el_rat := (`6`     +   `7`     +   `8`)*2 / (`1` + `2`) /3 ]


model_labs = c('CoMix', 'No data', 'No interaction', 'POLYMOD')
names(model_labs) = c(1,4,5,6)

par_labs = c('Infectiousness', 'Susceptibility')
names(par_labs) = c('inf_rate', 'susceptibility')

lab_dates

ggplot(inf_pars_wide) +
  geom_path(aes(x=ad_rat, y=el_rat, color=forecast_date), size=1)+
  scale_color_viridis(breaks = as.numeric(lab_dates), 
                                          labels = lab_dates, name='Infections')+
  new_scale_color()+
  
  geom_path(data=cas_pars_wide, aes(x=ad_rat, y=el_rat, color=forecast_date), size=1)+
  scale_color_viridis(breaks = as.numeric(lab_dates), 
                      labels = lab_dates, option="magma", name='Cases')+
  scale_y_continuous(trans='log2', name='Older adults relative to children', limits=c(1/3.25, 3))+
  scale_x_continuous(trans='log2', name='Younger adults relative to children', limits=c(1/3, 3))+
  geom_vline(xintercept=1, alpha=0.2)+
  geom_hline(yintercept=1, alpha=0.2)+
  facet_grid(name~run, scale='free', labeller=labeller(run=model_labs, name=par_labs))+
  theme_minimal()

ggsave('plots/paramplot_traj.png', width=10, height =6, units='in')
  

