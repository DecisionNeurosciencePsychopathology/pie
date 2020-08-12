library(rstan)
library(ggplot2)
# num_cores=3
# rstan_options(auto_write = TRUE)
# options(mc.cores = num_cores)
load("pie_data_processed.rdata")

segs=4
pts=0
excluded_IDs=NA
df_seg=df[df$num_segments==segs,]
df_seg_pts=df_seg[df_seg$show_points==pts,]
df_seg_pts_wex=df_seg_pts[!(df_seg_pts$ID %in% excluded_IDs),]
use_data=df_seg_pts_wex
seg_probs=c(0.36,0.43,0.56,0.65)

# make list of stan inputs ####
subj_IDs=unique(use_data$ID)
nS=length(subj_IDs)
nT=max(use_data$trial)*(dim(use_data)[1]/(nS*max(use_data$trial))) #max(use_data$block_num)/2

num_segments=array(0,c(nS,nT))
points_shown=array(0,c(nS,nT))
# choice=array(0,c(nS,nT))
# reward=array(0,c(nS,nT))
block_num=array(0,c(nS,nT))

for (s in 1:nS) {
  subj_data=use_data[use_data$ID==subj_IDs[s],]
  subj_blocks=unique(subj_data$block_num) #not all blocks if only 4 seg conditions
  num_segments[s,]=4 #subj_data$num_segments
  points_shown[s,]=as.numeric(subj_data$show_points)-1
  # choice[s,]=subj_data$selected_segment
  # reward[s,]=subj_data$win
  block_num[s,]=ifelse(subj_data$block_num==subj_blocks[1],1,
    ifelse(subj_data$block_num==subj_blocks[2],2,
      ifelse(subj_data$block_num==subj_blocks[3],3,4)))
}

data_toest=list(nS=nS,nT=nT,num_segments=num_segments,points_shown=points_shown,
                #choice=choice,reward=reward,
                block_num=block_num,seg_probs=seg_probs)

model_pp=stan_model(file='estimate_RL_decay_softmax_priorcheck.stan')
ppcheck_out=sampling(model_pp,data=data_toest,iter=1000,algorithm = "Fixed_param")
# ppcheck_out=sampling(model_pp,data=data_toest,iter=1,chains=1,algorithm = "Fixed_param")

#plot parameters
alphas=summary(ppcheck_out,pars='alpha')$summary[,'50%']
betas=summary(ppcheck_out,pars='beta')$summary[,'50%']
lambdas=summary(ppcheck_out,pars='lambda')$summary[,'50%']
hist(alphas)
hist(betas)
hist(lambdas)

#function to calculate mode
mode_func<-function(x) {
  tabresult <- tabulate(x)
    themode <- which(tabresult == max(tabresult))
    if(sum(tabresult == max(tabresult))>1) themode <- NA
    return(themode)
}
rand_choose<-function(x) {
  return(sample(x,1))
}
get_num<-function(x,num) {
  return(sum(x==num))
}


#add simulated choices to dataset
use_data$sim_choice_mode=NA
use_data$sim_choice_rand=NA
use_data$sim_reward_mode=NA
use_data$sim_reward_rand=NA
use_data$sim_choice_prop1=NA
use_data$sim_choice_prop2=NA
use_data$sim_choice_prop3=NA
use_data$sim_choice_prop4=NA
use_data$sim_reward_prop0=NA
use_data$sim_reward_prop1=NA
for (s in 1:nS) {
  choice_pars=paste0('choice[',s,',',seq(1,nT,by=1),']')
  # subj_sim_choices=summary(ppcheck_out,pars=choice_pars)$summary[,'50%']
  subj_sim_choices_all=extract(ppcheck_out,pars=choice_pars)
  subj_sim_choices_mode=as.numeric(apply(as.data.frame(subj_sim_choices_all),2,mode_func)) 
  subj_sim_choices_rand=as.numeric(apply(as.data.frame(subj_sim_choices_all),2,rand_choose))
  subj_sim_choice_prop1=as.numeric(apply(as.data.frame(subj_sim_choices_all),2,get_num,num=1))/2000
  subj_sim_choice_prop2=as.numeric(apply(as.data.frame(subj_sim_choices_all),2,get_num,num=2))/2000
  subj_sim_choice_prop3=as.numeric(apply(as.data.frame(subj_sim_choices_all),2,get_num,num=3))/2000
  subj_sim_choice_prop4=as.numeric(apply(as.data.frame(subj_sim_choices_all),2,get_num,num=4))/2000
  rew_pars=paste0('reward[',s,',',seq(1,nT,by=1),']')
  # subj_sim_rew=summary(ppcheck_out,pars=rew_pars)$summary[,'50%']
  subj_sim_rew_all=extract(ppcheck_out,pars=rew_pars)
  subj_sim_rew_mode=as.numeric(apply(as.data.frame(subj_sim_rew_all),2,mode_func)) 
  subj_sim_rew_rand=as.numeric(apply(as.data.frame(subj_sim_rew_all),2,rand_choose)) 
  subj_sim_reward_prop1=as.numeric(apply(as.data.frame(subj_sim_rew_all),2,get_num,num=1))/2000
  subj_sim_reward_prop0=as.numeric(apply(as.data.frame(subj_sim_rew_all),2,get_num,num=0))/2000
  subj_trials=which(use_data$ID==subj_IDs[s])
  use_data$sim_choice_mode[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_choices_mode
  use_data$sim_choice_rand[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_choices_rand
  use_data$sim_reward_mode[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_rew_mode
  use_data$sim_reward_rand[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_rew_rand
  use_data$sim_choice_prop1[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_choice_prop1
  use_data$sim_choice_prop2[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_choice_prop2
  use_data$sim_choice_prop3[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_choice_prop3
  use_data$sim_choice_prop4[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_choice_prop4
  use_data$sim_reward_prop0[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_reward_prop0
  use_data$sim_reward_prop1[subj_trials[1]:subj_trials[length(subj_trials)]]=subj_sim_reward_prop1
}

use_data_choice_long=pivot_longer(use_data,starts_with('sim_choice_prop'),
  names_to='choice',names_prefix='^sim_choice_prop',values_to='prop_choices')

#similarly convert real data
fdf_4seg=fdf[fdf$num_segments==4,]
choices_by_trial=as.data.frame(table(fdf_4seg$trial_adj,fdf_4seg$selected_prob))
names(choices_by_trial)=c('trial','choice','frequency')
choices_by_trial$prop=choices_by_trial$frequency/(nS*4)

fdf_4seg_nopts=fdf_4seg[fdf_4seg$show_points==0,]
choices_by_trial_nopts=as.data.frame(table(fdf_4seg_nopts$trial_adj,fdf_4seg_nopts$selected_prob))
names(choices_by_trial_nopts)=c('trial','choice','frequency')
choices_by_trial_nopts$prop=choices_by_trial_nopts$frequency/(nS*2)

# table(use_data$trial,use_data$sim_choice)

# use_data=use_data[!(use_data$sim_choice==2.5),]

#simulated data from prior
ggplot(use_data_choice_long[use_data_choice_long$trial>4,],aes(trial,prop_choices,color=choice))+
  geom_smooth(method='gam')+geom_point(alpha=0.1)+theme_classic()+coord_cartesian(ylim=c(0.1,0.4))
#real data
ggplot(choices_by_trial,aes(trial,prop,color=choice,fill=choice))+
  geom_line(aes(group=choice))+geom_point()+theme_classic()+coord_cartesian(ylim=c(0.1,0.4))
#real data: no points only
ggplot(choices_by_trial_nopts,aes(trial,prop,color=choice,fill=choice))+
  geom_line(aes(group=choice))+geom_point()+theme_classic()+coord_cartesian(ylim=c(0.1,0.4))


ggplot(use_data[use_data$trial>4,],aes(trial,sim_reward_rand))+ #sim_reward,
    # color=as.factor(sim_choice),fill=as.factor(sim_choice)))+
  geom_smooth(method='lm')+geom_count(aes(size = stat(prop)))+
  theme_classic()

ggplot(use_data[use_data$trial>4,],aes(trial,sim_reward_mode,
    color=as.factor(sim_choice_mode),fill=as.factor(sim_choice_mode)))+
  geom_smooth(method='lm')+geom_count(aes(size = stat(prop)))+
  theme_classic()
