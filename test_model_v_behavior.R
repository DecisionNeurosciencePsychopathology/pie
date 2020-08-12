library(rstan)
library(corrplot)
library(psych)
load("pie_data_processed.rdata")
fdffdf4=fdf[fdf$num_segments==4,]
data_use=fdf4
model_out=readRDS('models/RL_softmax_4_seg_1_points.rds') #'RL_pers_decay_softmax_4seg.rds')
params=summary(model_out,pars=c('alpha','beta'))$summary[,'50%'] #,'lambda','pers'
load('models/RL_softmax_4_seg_1_points_modelfit.RData') #pers_decay_

subjs=unique(data_use$ID)
num_subjs=length(model_out_LL_bysubj)
subj_perf=array(NA,dim=c(num_subjs,10))

for (s in 1:num_subjs) {
  subj_data=data_use[data_use$ID==subjs[s],]
  num_trials=dim(subj_data)[1]
  subj_perf[s,1]=as.numeric(subjs[s])
  subj_perf[s,2]=sum(subj_data$selected_prob>0.5)/num_trials
  subj_perf[s,3]=sum(subj_data$selected_prob>0.6)/num_trials
  subj_perf[s,4]=sum(subj_data$win)
  subj_data_nopts=subj_data[subj_data$show_points==0,]
  subj_perf[s,5]=sum(subj_data_nopts$selected_prob>0.5)/num_trials
  subj_perf[s,6]=sum(subj_data_nopts$selected_prob>0.6)/num_trials
  subj_perf[s,7]=sum(subj_data_nopts$win)
  subj_data_pts=subj_data[subj_data$show_points==1,]
  subj_perf[s,8]=sum(subj_data_pts$selected_prob>0.5)/num_trials
  subj_perf[s,9]=sum(subj_data_pts$selected_prob>0.6)/num_trials
  subj_perf[s,10]=sum(subj_data_pts$win)
}
subj_perf=data.frame(subj_perf)
names(subj_perf)=c('ID','all_prob50','all_best','all_win','nopts_prob50',
  'nopts_best','nopts_win','pts_prob50','pts_best','pts_win')
subj_perf$alpha=params[1:num_subjs]
subj_perf$beta=params[(num_subjs+1):(2*num_subjs)]
# subj_perf$lambda=params[(2*num_subjs+1):(3*num_subjs)]
# subj_perf$pers=params[(3*num_subjs+1):(4*num_subjs)]
subj_perf$AIC=model_out_V1_AIC
corr.test(subj_perf[,c(2,11:13)])
pairs(subj_perf[,c(2,11:13)],pch=19,col=rgb(.2,0,1,.5))
corrplot(cor(subj_perf[,c(2,11:13)]),diag=T,addCoef.col='black')


#test for bad subjects: using RL_softmax model w/pts shown only
dim(subj_perf[subj_perf$beta<0,])[1] #20
subj_perf$beta75=summary(model_out,pars='beta',probs=.75)$summary[,'75%']
subj_perf$beta25=summary(model_out,pars='beta',probs=.25)$summary[,'25%']
dim(subj_perf[subj_perf$beta75<0,])[1] #14
dim(subj_perf[subj_perf$beta25<0,])[1] #33
cor.test(subj_perf$beta75,subj_perf$pts_prob50) #r=0.71
subj_perf[subj_perf$beta75<0,]$ID # 6 11 14 15 23 48 62 65 66 88 89 90 95 40
subj_perf[subj_perf$beta<0,]$ID #6 10 11 14 15 19 23 26 28 48 62 65 66 74 88 89 90 95 97 40
# additional: 0010, 0019, 0026, 0028, 0074, 0097
subj_perf[subj_perf$beta25<0,]$ID # 2  6 10 11 14 15 16 19 21 22 23 26 28 34 39 
#41 48 52 62 63 65 66 67 74 81 87 88 89 90 93 95 97 40
# addiitional: 002,016,021,022,034,039,041,052,063,067,081,087,093
subj_perf[subj_perf$pts_prob50<0.25,]$ID
