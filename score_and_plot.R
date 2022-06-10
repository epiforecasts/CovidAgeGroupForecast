library(data.table)
library(ggplot2)
library(patchwork)
source('R/score_forecasts.R')
source('R/plotting.R')


summary_preds_inf = readRDS('outputs/summary_preds_infweek.rds')
summary_preds_cas = readRDS('outputs/summary_preds_cases.rds')


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


overall_scores_inf = summary_scores_inf$score_overall[,c('model', 'interval_score', 'aem', 'bias')]

#overall_scores_inf[, model := c('CoMix contact data', 
#                                'No contact data', 
#                                'No interaction',
#                                'Polymod contact data', 
#                                'Exponential baseline', 
#                                'Fixed value baseline')]
#
ggplot(overall_scores_inf) + 
  geom_point(aes(x=model, y=interval_score))


age_scores_inf = summary_scores_inf$score_by_age[,c('model', 'age_group', 'interval_score_rel', 'aem_rel', 'bias')]



age_inf = ggplot(age_scores_inf[model != 'baseline_linex_lv',]) +
  geom_segment(aes(x=model, xend=model, y=1, yend=interval_score_rel), alpha=0.3)+
  geom_point(aes(x=model, y=interval_score_rel, color=model))+
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
  facet_wrap(~age_group, nrow=1)+
  ggtitle("B")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'none'
  )


period_scores_inf = summary_scores_inf$score_by_period[,c('model', 'periods', 'interval_score_rel', 'aem_rel', 'bias')]



period_inf = ggplot(period_scores_inf[model != 'baseline_linex_lv',]) +
  geom_segment(aes(x=model, xend=model, y=1, yend=interval_score_rel), alpha=0.3)+
  geom_point(aes(x=model, y=interval_score_rel, color=model))+
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


overall_scores_cas = summary_scores_cas$score_overall[,c('model', 'interval_score', 'aem', 'bias')]

#overall_scores_cas[, model := c('CoMix contact data', 
#                                'No contact data', 
#                                'No interaction',
#                                'Polymod contact data', 
#                                'Exponential baseline', 
#                                'Fixed value baseline')]
#

age_scores_cas = summary_scores_cas$score_by_age[,c('model', 'age_group', 'interval_score_rel', 'aem_rel', 'bias')]



age_cas = ggplot(age_scores_cas[model != 'baseline_linex_lv',]) +
  geom_segment(aes(x=model, xend=model, y=1, yend=interval_score_rel), alpha=0.3)+
  geom_point(aes(x=model, y=interval_score_rel, color=model))+
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
  facet_wrap(~age_group, nrow=1)+
  ggtitle("C")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'bottom'
  )


period_scores_cas = summary_scores_cas$score_by_period[,c('model', 'periods', 'interval_score_rel', 'aem_rel', 'bias')]



period_cas = ggplot(period_scores_cas[model != 'baseline_linex_lv',]) +
  geom_segment(aes(x=model, xend=model, y=1, yend=interval_score_rel), alpha=0.3)+
  geom_point(aes(x=model, y=interval_score_rel, color=model))+
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
    legend.position = 'none'
  )




overall_scores = rbind(overall_scores_inf[, type:='Infections'], overall_scores_cas[, type:='Cases'])

overall_scores = merge(overall_scores, overall_scores[model=='baseline_expex_lv'], by = c('type'), suffixes = c('','_baseline'))
overall_scores[, interval_score_rel := interval_score / interval_score_baseline]

overall_sp = 
  ggplot(overall_scores) + 
  geom_segment(aes(x=model, xend=model, y=1, yend=interval_score_rel), alpha=0.3)+
  geom_point(aes(x=model, y=interval_score_rel, color=model))+
  facet_wrap(~type)+
  scale_color_discrete(labels=c('CoMix contact data', 
                                'No contact data', 
                                'No interaction',
                                'Polymod contact data', 
                                'Exponential baseline', 
                                'Fixed value baseline'), 
                       name = 'Model')+
  scale_x_discrete(labels=NULL, name='')+
  scale_y_continuous(trans='log2', name='')+
  theme_minimal()+
  ggtitle('A')+
  theme(
    axis.text.x = element_text(angle=90), 
    legend.position = 'none'
  )

overageres = overall_sp / age_inf / age_cas + plot_layout(heights = c(2,1,1))

ggsave('plots/overall_results.pdf', plot = overageres, width = 10, height=7, units = 'in')
ggsave('plots/overall_results.png', plot = overageres, width = 7, height=7, units = 'in')


