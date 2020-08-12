source('sim_functions.R')
library(ggplot2)

#set up simulation parameters ####
num_trials=35
num_subjs=100
num_segments=4
if (num_segments==4) {
  seg_probs=c(0.36,0.43,0.56,0.65)
}

#distributions to sample parameter values
#learning rate
alpha_alpha=2
alpha_beta=2
#decay
decay_mean=0.8
decay_sd=0.05
#kappa (sutton bonus)
kappa_mean=0.25
kappa_sd=0.1
#epsilon (e-greedy)
epsilon_alpha=1
epsilon_beta=3
#beta (inverse temperature, gamma dist)
beta_alpha=6
beta_beta=2
#tau (perseveration)
tau_mean=0.25
tau_sd=0.2
#omega (UCB variance effect)
omega_mean=0.5
omega_sd=0.25

#list functions and note combination exclusions
value_functions=c('basic_RL_sim','RL_decay_sim','RL_SB_sim','ideal_bayes_sim',
                  'forgetting_bayes_sim')
value_nobayes=c(1,1,1,0,0)
value_decay=c(0,1,0,0,1)
value_SB=c(0,0,1,0,0)
choice_functions=c('greedy_sim','e_greedy_sim','softmax_sim','softmax_pers_sim',
'softmax_SB_sim','thompson_sim','thompson_softmax_sim','thompson_egreedy_sim',
'UCB_sim','UCB_egreedy_sim','UCB_softmax_sim','thompson_UCB_sim',
'thompson_UCB_egreedy_sim','thompson_UCB_softmax_sim')
choice_bayes=c(0,0,0,0,0,1,1,1,1,1,1,1,1,1)
choice_decay=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0)
choice_SB=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0)

plots=list()

# loop through value and choice function combinations and simulate
for (v in 1:length(value_functions)) {
  for (c in 1:length(choice_functions)) {
    
    #check that combo of functions is allowed (may be a better way to do this)
    if (max(c((value_nobayes[v]+choice_bayes[c]),(value_decay[v]+choice_decay[c]),
            (value_SB[v]+choice_SB[c])))<2) {
      
      #sample parameter values
      alpha_params=rbeta(num_subjs,alpha_alpha,alpha_beta)
      decay_params=rnorm(num_subjs,mean=decay_mean,sd=decay_sd)
      kappa_params=rnorm(num_subjs,mean=kappa_mean,sd=kappa_sd)
      epsilon_params=rbeta(num_subjs,epsilon_alpha,epsilon_beta)
      beta_params=rgamma(num_subjs,beta_alpha,rate=beta_beta)
      tau_params=rnorm(num_subjs,mean=tau_mean,sd=tau_sd)
      omega_params=rnorm(num_subjs,mean=omega_mean,sd=omega_sd)
      
      #set up variables
      all_alphas=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
      all_betas=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
      all_values=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
      seg_samplehx=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
      seg_samplehx_selected=array(data=NA,dim=c(num_trials,num_subjs))
      all_choices=array(data=NA,dim=c(num_trials,num_subjs))
      all_rews=array(data=NA,dim=c(num_trials,num_subjs))
      out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=13)
      
      for (s in 1:num_subjs) {
       all_alphas[1,,s]=1
       all_betas[1,,s]=1
       all_values[1,,s]=0.5
       seg_samplehx[1,,s]=0
        
        for (t in 1:num_trials) {
          
          #choice
          if (t<(num_segments+1)) { #forced sampling: assume even for now
            chosen_seg=t
          } else {
            if (value_nobayes[v]==0) {
              all_values[t,,s]=all_alphas[t,,s]/(all_alphas[t,,s]+all_betas[t,,s])
            }
            chosen_seg=do.call(choice_functions[c],list(in_values=all_values[t,,s],
              in_alphas=all_alphas[t,,s],in_betas=all_betas[t,,s],
              seg_samplehx=seg_samplehx[t,,s],last_chosen=all_choices[t-1,s],
              beta=beta_params[s],tau=tau_params[s],kappa=kappa_params[s],
              omega=omega_params[s],epsilon=epsilon_params[s]))
          }
          
          #update variables
          all_choices[t,s]=chosen_seg
          seg_samplehx[t+1,chosen_seg,s]=seg_samplehx[t,chosen_seg,s]+1
          seg_samplehx[t+1,-chosen_seg,s]=seg_samplehx[t,-chosen_seg,s]
          seg_samplehx_selected[t,s]=seg_samplehx[t+1,chosen_seg,s]
          
          #outcome
          rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
          all_rews[t,s]=rew

          #update values
          out_args=do.call(value_functions[v],list(in_values=all_values[t,,s],
            in_alpha=all_alphas[t,,s],in_beta=all_betas[t,,s],rew=rew,
            chosen_seg=chosen_seg,seg_samplehx=seg_samplehx[t,,s],
            alpha=alpha_params[s],decay=decay_params[s],kappa=kappa_params[s]))
          if (value_nobayes[v]==0) {
            all_alphas[t+1,,s]=out_args[1:num_segments]
            all_betas[t+1,,s]=out_args[(num_segments+1):(2*num_segments)]
          } else all_values[t+1,,s]=out_args
          
          out_info[(s-1)*num_trials+t,]=c(s,t,chosen_seg,
            seg_samplehx_selected[t,s],seg_probs[chosen_seg],rew,alpha_params[s],
            decay_params[s],kappa_params[s],epsilon_params[s],beta_params[s],
            tau_params[s],omega_params[s])
        }
      }
      
      #plot
      out_info_df=data.frame(out_info)
      names(out_info_df)=c('subj','trial','chosen_seg','samplehx_selected',
        'selected_prob','rew','alpha','decay','kappa','epsilon','beta','tau','omega')
      ggsave(paste('sim_plot_',c,v,'.png',sep='_'),ggplot(out_info_df,aes(trial,samplehx_selected,
        color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
        geom_smooth(method='gam')+geom_point(alpha=0.1)+theme_classic()+
        ylim(0,20)+labs(y='number of times option selected',
             title=paste(value_functions[v],choice_functions[c],sep=' with '))+
        scale_color_discrete(name='Probability of Reward')+
        scale_fill_discrete(name='Probability of Reward')+
        theme(legend.position="bottom"),
        width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
      
      # plots[(v-1)*length(choice_functions)+v]=ggplot(out_info_df,aes(trial,samplehx_selected,color=as.factor(selected_prob),
      #                        fill=as.factor(selected_prob)))+
      #   geom_smooth(method='gam')+geom_point(alpha=0.1)+theme_classic()+
      #   ylim(0,20)+labs(y='number of times option selected',
      #        title=paste(value_functions[v],choice_functions[c],sep=' with '))+
      #   scale_color_discrete(name='Probability of Reward')+
      #   scale_fill_discrete(name='Probability of Reward')
      
      #readline('press any key to move onto next combination of functions')
    }
  }
}
