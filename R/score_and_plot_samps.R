library(data.table)
library(ggplot2)
library(patchwork)
source('R/score_forecasts_samples.R')
source('R/plotting.R')


summary_preds_inf = readRDS('outputs/samples_preds_infweek.rds')
summary_preds_cas = readRDS('outputs/samples_preds_cases.rds')


age_labs_inf_reverse = 1:7
names(age_labs_inf_reverse) = c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' )

age_labs_inf_forward = names(age_labs_case_reverse)
names(age_labs_inf_forward) = age_labs_case_reverse


summary_scores_inf = score_forecasts(summary_preds_inf[date > as.Date('2020-10-01') & !(run %in% c(2,3) | run=='baseline_linex_lv')], 
                                 pandemic_periods,
                                 suffix = '_infections', 
                                 age_groups =  c('2-10', '11-15', '16-24', '25-34', '35-49', '50-69', '70+' ),
                                 labels = c('CoMix contact data', 
                                            'No contact data', 
                                            'No interaction',
                                            'Polymod contact data', 
                                            'Exponential baseline', 
                                            'Fixed value baseline'))


overall_scores_inf = summary_scores_inf$score_overall[,c('model', 'horizon', 'crps', 'ae_median', 'bias')]

#overall_scores_inf[, model := c('CoMix contact data', 
#                                'No contact data', 
#                                'No interaction',
#                                'Polymod contact data', 
#                                'Exponential baseline', 
#                                'Fixed value baseline')]
#
ggplot(overall_scores_inf) + 
  geom_point(aes(x=model, y=crps))


age_scores_inf = summary_scores_inf$score_by_age[,c('model', 'age_group', 'horizon', 'crps_rel', 'ae_median_rel', 'bias')]


age_scores_inf[, age_index := age_labs_inf_reverse[age_group]]

age_inf = ggplot(age_scores_inf[model != 'baseline_linex_lv',]) +
  #geom_segment(aes(x=horizon, xend=horizon, y=1, yend=crps_rel), alpha=0.3)+
  geom_point(aes(x=horizon, y=crps_rel, color=model))+
  geom_line(aes(x=horizon, y=crps_rel, color=model), linetype='dashed')+
  geom_hline(yintercept = 1)+
  scale_y_continuous(trans='log2', name='crps')+
  scale_x_discrete(labels=NULL, name='')+
  scale_color_discrete(labels=c('CoMix contact data', 
                                'No contact data', 
                                'No interaction',
                                'Polymod contact data', 
                                'Exponential baseline', 
                                'Fixed value baseline'), 
                       name = 'Model')+
  facet_wrap(~age_index, nrow=1, labeller = labeller(age_index = age_labs_inf_forward))+
  ggtitle("B")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'none'
  )


period_scores_inf = summary_scores_inf$score_by_period[,c('model', 'periods', 'horizon', 'crps_rel', 'ae_median_rel', 'bias')]



period_inf = ggplot(period_scores_inf[model != 'baseline_linex_lv',]) +
  geom_line(aes(x=horizon, y=crps_rel, color=model), alpha=0.6, linetype='dashed')+
  geom_point(aes(x=horizon, y=crps_rel, color=model), alpha=0.6)+
  geom_hline(yintercept = 1)+
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


period_inf_bias = ggplot(period_scores_inf[model != 'baseline_linex_lv',]) +
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


summary_scores_cas = score_forecasts(summary_preds_cas[date > as.Date('2020-10-01') & !(run %in% c(2,3) | run=='baseline_linex_lv')], 
                                    pandemic_periods, 
                                    suffix = '_cases', 
                                    age_groups = age_groups,
                                    labels = c('CoMix contact data', 
                                               'No contact data', 
                                               'No interaction',
                                               'Polymod contact data', 
                                               'Exponential baseline', 
                                               'Fixed value baseline'))


overall_scores_cas = summary_scores_cas$score_overall[,c('model', 'horizon', 'crps', 'ae_median', 'bias')]

#overall_scores_cas[, model := c('CoMix contact data', 
#                                'No contact data', 
#                                'No interaction',
#                                'Polymod contact data', 
#                                'Exponential baseline', 
#                                'Fixed value baseline')]
#

age_scores_cas = summary_scores_cas$score_by_age[,c('model', 'age_group', 'horizon', 'crps_rel', 'ae_median_rel', 'bias')]



age_cas = ggplot(age_scores_cas[model != 'baseline_linex_lv',]) +
  #geom_segment(aes(x=horizon, xend=horizon, y=1, yend=crps_rel), alpha=0.3)+
  geom_point(aes(x=horizon, y=crps_rel, color=model))+
  geom_line(aes(x=horizon, y=crps_rel, color=model), linetype='dashed')+
  geom_hline(yintercept = 1)+
  scale_y_continuous(trans='log2', name='crps')+
  scale_x_discrete(labels=NULL, name='Horizon')+
  scale_color_discrete(labels=c('CoMix contact data', 
                                'No contact data', 
                                'No interaction',
                                'Polymod contact data', 
                                'Exponential baseline', 
                                'Fixed value baseline'), 
                       name = 'Model')+
  facet_wrap(~age_group, nrow=1)+
  ggtitle("C")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom'
  )


period_scores_cas = summary_scores_cas$score_by_period[,c('model', 'periods', 'horizon', 'crps_rel', 'ae_median_rel', 'bias')]



period_cas = ggplot(period_scores_cas[model != 'baseline_linex_lv',]) +
  geom_line(aes(x=horizon, y=crps_rel, color=model), alpha=0.6, linetype='dashed')+
  geom_point(aes(x=horizon, y=crps_rel, color=model), alpha=0.6)+
  geom_hline(yintercept = 1)+
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
  ggtitle("B")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom', 
    strip.text.x = element_blank()
  )

period_inf/period_cas


period_all_bias = ggplot(rbind(period_scores_cas[,type:='cases'], period_scores_inf[,type:='infections'])[model != 'baseline_linex_lv',]) +
  geom_path(aes(x=bias, y=crps_rel, color=model), alpha=0.6, linetype='dashed')+
  geom_point(aes(x=bias, y=crps_rel, color=model, size=horizon), alpha=0.6)+
  geom_hline(yintercept=1)+
  geom_vline(xintercept=0)+
  scale_y_continuous(trans='log2', name='CRPS')+
  scale_x_discrete(labels=NULL, name='Bias')+
  scale_size(breaks = c(0,7,14,21,28))+
  scale_color_discrete(labels=c('CoMix contact data', 
                                'No contact data', 
                                'No interaction',
                                'Polymod contact data', 
                                'Exponential baseline', 
                                'Fixed value baseline'), 
                       name = 'Model')+
  facet_grid(type~periods, scale='free')+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom'
  )

ggsave('plots/periods.png',width = 10, height=5, units='in')


overall_scores = rbind(overall_scores_inf[, type:='Infections'], overall_scores_cas[, type:='Cases'])

overall_scores = merge(overall_scores, overall_scores[model=='baseline_expex_lv'], by = c('type', 'horizon'), suffixes = c('','_baseline'))
overall_scores[, crps_rel := crps / crps_baseline]

overall_sp = 
  ggplot(overall_scores) + 
  #geom_curve(aes(x=model, xend=model, y=1, yend=crps_rel+horizon*0.00001, size=exp(horizon/7)), alpha=0.3, curvature = 0.15)+
  geom_path(aes(x=bias, y=crps_rel, color=model), alpha=0.7)+
  geom_point(aes(x=bias, y=crps_rel, color=model, size=horizon), alpha=0.7)+
  geom_hline(yintercept = 1, alpha=0.5)+
  geom_vline(xintercept = 0, alpha=0.5)+
  facet_wrap(~type)+
  scale_color_discrete(labels=c('CoMix contact data', 
                                'No contact data', 
                                'No interaction',
                                'Polymod contact data', 
                                'Exponential baseline', 
                                'Fixed value baseline'), 
                       name = 'Model')+
  scale_y_continuous(trans='log2', name='crps')+
  theme_minimal()+
  ggtitle('A')+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'none'
  )

overageres = overall_sp / age_inf / age_cas + plot_layout(heights = c(2,1,1))

ggsave('plots/overall_results.pdf', plot = overageres, width = 7, height=7, units = 'in')
ggsave('plots/overall_results.png', plot = overageres, width = 7, height=7, units = 'in', dpi=300)




inf_pars =readRDS('outputs/summary_pars_infweek.rds')

inf_pars
