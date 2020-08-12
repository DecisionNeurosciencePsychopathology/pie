library(tidyr)
library(ggplot2)
library(ggbeeswarm)
# all_models=c('Bayes_decay_softmax','Bayes_decay_Thompson_Qparam_softmax','Bayes_softmax',
#   'Bayes_pers_softmax','Bayes_pers_decay_softmax','RL_decay_softmax', 
#   'RL_decay_softmax_expbonus','RL_pers_decay_softmax','RL_pers_softmax','RL_softmax')

all_models=c('Bayes_decay_softmax','Bayes_decay_Thompson_Qparam_softmax','Bayes_softmax',
  'Bayes_pers_softmax','Bayes_pers_decay_softmax','RL_decay_softmax', 
  'RL_decay_softmax_expbonus','RL_pers_decay_softmax','RL_pers_softmax','RL_softmax',
  'RL_softmax_altAprior')
model_convergence=c(0,NA,1,1,0,0,NA,0,1,1,1) #all/14ex
# model_convergence=c(1,NA,1,1,1,0,NA,0,1,1) #20ex
inc_models=c(3,4,9,10,11) #c(1,3,4,5,6,8,9,10) #c(3,4,9,10)  #c(1,5,8) #

test_models=all_models[inc_models]
num_models=length(test_models)

load(paste0('models/',test_models[1],'_4_seg_0_points_modelfit.RData')) #_1_points_20_subjex
nS=length(model_out_LL_bysubj)
model_LL=array(data=NA,dim=c(nS,num_models*3))
model_AIC=array(data=NA,dim=c(nS,num_models*3))

model_LL[,1]=model_out_LL_bysubj
model_AIC[,1]=model_out_V1_AIC
model_LL[,(num_models*2+1)]=model_convergence[1]
model_AIC[,(num_models*2+1)]=model_convergence[1]

for (m in 2:num_models) {
  load(paste0('models/',test_models[m],'_4_seg_0_points_modelfit.RData')) #1_points_20_subjex_
  model_LL[,m]=model_out_LL_bysubj
  model_AIC[,m]=model_out_V1_AIC
  model_LL[,(num_models*2+m)]=model_convergence[m]
  model_AIC[,(num_models*2+m)]=model_convergence[m]
}

for (s in 1:nS) {
  model_LL[s,(num_models+1):(num_models*2)]=rank(model_LL[s,1:num_models])
  model_AIC[s,(num_models+1):(num_models*2)]=rank(model_AIC[s,1:num_models])
}

mean_AIC=colSums(model_AIC[,1:num_models])/nS 
median_AIC=apply(model_AIC[,1:num_models],2,FUN=median)#median AIC
best_AIC=colSums(model_AIC[,(num_models+1):(num_models*2)]==1)

model_AIC=data.frame(model_AIC)
names(model_AIC)=c(test_models,paste0('rank_',test_models),paste0('conv_',test_models))
model_AIC=model_AIC[order(model_AIC[,which(median_AIC==min(median_AIC))]),]
model_AIC$subj=rep(1:nS)
all_AIC_long=pivot_longer(model_AIC[c(1:num_models,(3*num_models+1))],test_models,names_to='model')

all_AIC_long$convergence=NA
for (l in 1:dim(all_AIC_long)[1]) {
  all_AIC_long$convergence[l]=model_convergence[which(all_models==all_AIC_long$model[l])]
}

ggplot(all_AIC_long,aes(x=model,y=value,color=subj,alpha=(convergence/3+0.6)))+ #fill=model
  geom_violin(draw_quantiles=0.5,color='gray90',fill='gray90',size=1.2)+
  geom_boxplot(fill=NA,width=0.2,fatten=3,color='gray30',outlier.size=0)+
  geom_beeswarm()+
  scale_color_gradient(low = "blue", high = "red")+
  theme_minimal()+
  scale_alpha_continuous(breaks=c(1/3+0.6,0.6),labels=c('yes','no'),name='model \nconvergence')+
  labs(color='ranked AIC \nfor best model',y='AIC per participant')

