source('sim_functions.R')
library(ggplot2)
library(tidyr)

#set up simulation parameters ####
num_trials=30
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
decay_mean=c(1.05,0.8,1)
decay_sd=0.05
#kappa (sutton bonus)
kappa_mean=c(1,-1,0)
kappa_sd=0.25
#epsilon (e-greedy)
epsilon_alpha=1
epsilon_beta=3
#beta (inverse temperature, gamma dist)
beta_alpha=c(2,20,8) #c(2,200,10)
beta_beta=2
#tau (perseveration)
tau_mean=c(0.75,-0.75,0)
tau_sd=0.25
#omega (UCB variance effect)
omega_mean=c(20,-20,0)
omega_sd=1

#list functions and note combination exclusions
value_functions=c('basic_RL_sim','RL_decay_sim','RL_SB_sim','ideal_bayes_sim',
                  'forgetting_bayes_sim')
value_nobayes=c(1,1,1,0,0)
value_decay=c(0,1,0,0,1)
value_SB=c(0,0,1,0,0)
choice_functions=c('greedy_sim','e_greedy_sim','softmax_sim','softmax_pers_sim',
'softmax_SB_sim','thompson_sim','thompson_softmax_sim','thompson_egreedy_sim',
'UCB_sim','UCB_egreedy_sim','UCB_softmax_sim','thompson_UCB_sim',
'thompson_UCB_egreedy_sim','thompson_UCB_softmax_sim','softmax_pers_SB_sim',
'softmax_pers_SBalt_sim','softmax_pers_fUCB_sim')
choice_bayes=c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,0)
choice_decay=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,1)
choice_SB=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,1)

v_RL=which(value_functions=='RL_decay_sim')
v_bayes=which(value_functions=='forgetting_bayes_sim')
c_RL=which(choice_functions=='softmax_pers_fUCB_sim')
c_bayes=which(choice_functions=='thompson_UCB_softmax_sim')
c_bayes_noT=which(choice_functions=='UCB_softmax_sim')

total_vars_alter=5 #decay, kappa, beta, tau, omega
var_levels=array(data=3,dim=total_vars_alter)
var_names=c('decay', 'kappa', 'beta', 'tau', 'omega')
level_names=c('high','low','med')
level_names_beta=c('low','high','med')

seg_dist_mat=cbind(c(0,1,2,1),c(1,0,1,2),c(2,1,0,1),c(1,2,1,0))

for (a in 1:total_vars_alter) { #loop through params to change
  for (level in 1:3) {
    var_levels[a]=level
    var_levels[-a]=3
    level_name_sel=ifelse(a==3,level_names_beta[level],level_names[level])
    
    #sample parameter values
    alpha_params=rbeta(num_subjs,alpha_alpha,alpha_beta)
    decay_params=rnorm(num_subjs,mean=decay_mean[var_levels[1]],sd=decay_sd)
    kappa_params=rnorm(num_subjs,mean=kappa_mean[var_levels[2]],sd=kappa_sd)
    epsilon_params=rbeta(num_subjs,epsilon_alpha,epsilon_beta)
    beta_params=rgamma(num_subjs,beta_alpha[var_levels[3]],rate=beta_beta)
    tau_params=rnorm(num_subjs,mean=tau_mean[var_levels[4]],sd=tau_sd)
    omega_params=rnorm(num_subjs,mean=omega_mean[var_levels[5]],sd=omega_sd)
    
    # #set up variables: RL model
    all_alphas=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    all_betas=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    all_values=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    seg_samplehx=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    seg_samplehx_selected=array(data=NA,dim=c(num_trials,num_subjs))
    all_choices=array(data=NA,dim=c(num_trials,num_subjs))
    all_rews=array(data=NA,dim=c(num_trials,num_subjs))
    seg_dist=array(data=NA,dim=c(num_trials,num_subjs))
    seg_dist_bin=array(data=NA,dim=c(num_trials,num_subjs))
    pers_exp_bins=array(data=NA,dim=c(num_trials,num_subjs))
    max_val=array(data=NA,dim=c(num_trials,num_subjs))
    stay_max=array(data=0,dim=c(num_trials,num_subjs))
    stay_nomax=array(data=0,dim=c(num_trials,num_subjs))
    switch_max=array(data=0,dim=c(num_trials,num_subjs))
    switch_nomax=array(data=0,dim=c(num_trials,num_subjs))
    postmax_val=array(data=NA,dim=c(num_trials,num_subjs))
    stay_postmax=array(data=0,dim=c(num_trials,num_subjs))
    stay_postnomax=array(data=0,dim=c(num_trials,num_subjs))
    switch_postmax=array(data=0,dim=c(num_trials,num_subjs))
    switch_postnomax=array(data=0,dim=c(num_trials,num_subjs))

    out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=29)

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
          if (value_nobayes[v_RL]==0) {
            all_values[t,,s]=all_alphas[t,,s]/(all_alphas[t,,s]+all_betas[t,,s])
          }
          chosen_seg=do.call(choice_functions[c_RL],list(in_values=all_values[t,,s],
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
        seg_dist[t,s]=ifelse(t==1,NA,seg_dist_mat[chosen_seg,all_choices[t-1,s]])
        seg_dist_bin[t,s]=ifelse(seg_dist[t,s]==0,0,1)
        max_val[t,s]=ifelse(chosen_seg %in% which(all_values[t,,s]==max(all_values[t,,s])),1,0)
        postmax_val[t,s]=ifelse(t==1,NA,
          ifelse(all_choices[t-1,s] %in% which(all_values[t-1,,s]==max(all_values[t-1,,s])),1,0))
        pers_exp_bins[t,s]=2*seg_dist_bin[t,s]+max_val[t,s]
        if (is.na(seg_dist_bin[t,s])) {
        } else if (max_val[t,s]==1&&seg_dist_bin[t,s]==0) {
          stay_max[t,s]=1
        } else if (max_val[t,s]==0&&seg_dist_bin[t,s]==0) {
          stay_nomax[t,s]=1
        } else if (max_val[t,s]==1&&seg_dist_bin[t,s]==1) {
          switch_max[t,s]=1
        } else if (max_val[t,s]==0&&seg_dist_bin[t,s]==1) {
          switch_nomax[t,s]=1
        }
        if (t==1||is.na(seg_dist_bin[t-1,s])) {
        } else if (postmax_val[t,s]==1&&seg_dist_bin[t-1,s]==0) {
          stay_postmax[t,s]=1
        } else if (postmax_val[t,s]==0&&seg_dist_bin[t-1,s]==0) {
          stay_postnomax[t,s]=1
        } else if (postmax_val[t,s]==1&&seg_dist_bin[t-1,s]==1) {
          switch_postmax[t,s]=1
        } else if (postmax_val[t,s]==0&&seg_dist_bin[t-1,s]==1) {
          switch_postnomax[t,s]=1
        }

        #outcome
        rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
        all_rews[t,s]=rew

        #update values
        out_args=do.call(value_functions[v_RL],list(in_values=all_values[t,,s],
          in_alpha=all_alphas[t,,s],in_beta=all_betas[t,,s],rew=rew,
          chosen_seg=chosen_seg,seg_samplehx=seg_samplehx[t,,s],
          alpha=alpha_params[s],decay=decay_params[s],kappa=kappa_params[s]))
        if (value_nobayes[v_RL]==0) {
          all_alphas[t+1,,s]=out_args[1:num_segments]
          all_betas[t+1,,s]=out_args[(num_segments+1):(2*num_segments)]
          all_values[t+1,,s]=all_alphas[t+1,,s]/(all_alphas[t+1,,s]+all_betas[t+1,,s])
        } else all_values[t+1,,s]=out_args

        out_info[(s-1)*num_trials+t,]=c(s,t,chosen_seg,seg_dist[t,s],
          seg_samplehx_selected[t,s],seg_probs[chosen_seg],rew,alpha_params[s],
          decay_params[s],kappa_params[s],epsilon_params[s],beta_params[s],
          tau_params[s],omega_params[s],seg_dist_bin[t,s],max_val[t,s],pers_exp_bins[t,s],
          stay_max[t,s],stay_nomax[t,s],switch_max[t,s],switch_nomax[t,s],
          stay_postmax[t,s],stay_postnomax[t,s],switch_postmax[t,s],switch_postnomax[t,s],
          all_values[t,,s])
      }
    }

    #create data for plotting
    out_info_df=data.frame(out_info)
    names(out_info_df)=c('subj','trial','chosen_seg','seg_dist','samplehx_selected',
      'selected_prob','rew','alpha','decay','kappa','epsilon','beta','tau','omega',
      'binary_seg_dist','max_value_chosen','pers_exp_bins','stay_max','stay_nomax',
      'switch_max','switch_nomax','stay_postmax','stay_postnomax','switch_postmax',
      'switch_postnomax','modelvalue_1','modelvalue_2','modelvalue_3','modelvalue_4')
    out_info_df$modelvaluediff_1=seg_probs[1]-out_info_df$modelvalue_1
    out_info_df$modelvaluediff_2=seg_probs[2]-out_info_df$modelvalue_2
    out_info_df$modelvaluediff_3=seg_probs[3]-out_info_df$modelvalue_3
    out_info_df$modelvaluediff_4=seg_probs[4]-out_info_df$modelvalue_4

    #make long format for with rows for each segment
    out_info_df$subj=factor(out_info_df$subj)
    out_info_long=pivot_longer(out_info_df,modelvalue_1:modelvaluediff_4,
      names_to=c('.value','segment'),names_sep="_")
    
    #plot
    ggsave(paste('plots/plots_04282020/sim_plot_RL_',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_df,aes(trial,samplehx_selected,
        color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
        geom_smooth(method='gam')+geom_point(alpha=0.1)+theme_classic()+
        ylim(0,20)+labs(y='number of times option selected',
          title=paste(value_functions[v_RL],choice_functions[c_RL],
                      var_names[a],level_name_sel,sep=' '))+
        scale_color_discrete(name='Probability of Reward')+
        scale_fill_discrete(name='Probability of Reward')+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
    
    ggsave(paste('plots/plots_04282020/sim_plot_RL_segvalues_',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_long,aes(trial,modelvaluediff,color=as.factor(segment),
        fill=as.factor(segment)))+geom_hline(yintercept=0)+
        geom_smooth(method='loess')+theme_classic()+#geom_point(alpha=0.1)+
        ylim(-0.2,0.2)+labs(y='difference in actual vs. learned segment probability',
          title=paste(value_functions[v_RL],choice_functions[c_RL],
                      var_names[a],level_name_sel,sep=' '))+
        scale_color_discrete(name='Segment Probabilities',labels=seg_probs)+
        scale_fill_discrete(name='Segment Probabilities',labels=seg_probs)+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')

    ggsave(paste('plots/plots_04282020/sim_plot_RL_segdist',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial,seg_dist,
          color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
          geom_smooth(method='gam')+theme_classic()+
          ylim(0,2)+labs(y='distance of selected segment from previous choice',
            title=paste(value_functions[v_RL],choice_functions[c_RL],
                        var_names[a],level_name_sel,sep=' '))+
          scale_color_discrete(name='Probability of Reward')+
          scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')

      ggsave(paste('plots/plots_04282020/sim_plot_RL_expchoices',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial))+
          geom_smooth(aes(y=stay_max),method='gam',color='maroon',linetype='solid')+
          geom_smooth(aes(y=stay_nomax),method='gam',color='maroon',linetype='dashed')+
          geom_smooth(aes(y=switch_max),method='gam',color='blue',linetype='solid')+
          geom_smooth(aes(y=switch_nomax),method='gam',color='blue',linetype='dashed')+
          theme_classic()+
          ylim(0,1)+labs(y='Proportion choices',
            title=paste(value_functions[v_RL],choice_functions[c_RL],
                        var_names[a],level_name_sel,sep=' '))+
          #scale_linetype_discrete(name='Choice',labels=c('stay','switch'))+
          #scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')

      ggsave(paste('plots/plots_04282020/sim_plot_RL_post_expchoices',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial))+
          geom_smooth(aes(y=stay_postmax),method='gam',color='maroon',linetype='solid')+
          geom_smooth(aes(y=stay_postnomax),method='gam',color='maroon',linetype='dashed')+
          geom_smooth(aes(y=switch_postmax),method='gam',color='blue',linetype='solid')+
          geom_smooth(aes(y=switch_postnomax),method='gam',color='blue',linetype='dashed')+
          theme_classic()+
          ylim(0,1)+labs(y='Proportion choices',
            title=paste(value_functions[v_RL],choice_functions[c_RL],
                        var_names[a],level_name_sel,sep=' '))+
          #scale_linetype_discrete(name='Choice',labels=c('stay','switch'))+
          #scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')

      ggsave(paste('plots/plots_04282020/sim_plot_RL_value_',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_df,aes(trial,rew))+
        geom_smooth(method='gam')+theme_classic()+
        coord_cartesian(ylim=c(0.4,0.6))+labs(y='Proportion Rewarded Choices',
          title=paste(value_functions[v_RL],choice_functions[c_RL],
                      var_names[a],level_name_sel,sep=' '))+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
    
    #set up variables: Bayesian model w/Thompson sampling
    all_alphas=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    all_betas=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    all_values=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    seg_samplehx=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    seg_samplehx_selected=array(data=NA,dim=c(num_trials,num_subjs))
    all_choices=array(data=NA,dim=c(num_trials,num_subjs))
    all_rews=array(data=NA,dim=c(num_trials,num_subjs))
    seg_dist=array(data=NA,dim=c(num_trials,num_subjs))
    seg_dist_bin=array(data=NA,dim=c(num_trials,num_subjs))
    pers_exp_bins=array(data=NA,dim=c(num_trials,num_subjs))
    max_val=array(data=NA,dim=c(num_trials,num_subjs))
    stay_max=array(data=0,dim=c(num_trials,num_subjs))
    stay_nomax=array(data=0,dim=c(num_trials,num_subjs))
    switch_max=array(data=0,dim=c(num_trials,num_subjs))
    switch_nomax=array(data=0,dim=c(num_trials,num_subjs))
    postmax_val=array(data=NA,dim=c(num_trials,num_subjs))
    stay_postmax=array(data=0,dim=c(num_trials,num_subjs))
    stay_postnomax=array(data=0,dim=c(num_trials,num_subjs))
    switch_postmax=array(data=0,dim=c(num_trials,num_subjs))
    switch_postnomax=array(data=0,dim=c(num_trials,num_subjs))
    out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=29)
    
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
          if (value_nobayes[v_bayes]==0) {
            all_values[t,,s]=all_alphas[t,,s]/(all_alphas[t,,s]+all_betas[t,,s])
          }
          chosen_seg=do.call(choice_functions[c_bayes],list(in_values=all_values[t,,s],
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
        seg_dist[t,s]=ifelse(t==1,NA,seg_dist_mat[chosen_seg,all_choices[t-1,s]])
        seg_dist_bin[t,s]=ifelse(seg_dist[t,s]==0,0,1)
        max_val[t,s]=ifelse(chosen_seg %in% which(all_values[t,,s]==max(all_values[t,,s])),1,0)
        postmax_val[t,s]=ifelse(t==1,NA,
          ifelse(all_choices[t-1,s] %in% which(all_values[t-1,,s]==max(all_values[t-1,,s])),1,0))
        pers_exp_bins[t,s]=2*seg_dist_bin[t,s]+max_val[t,s]
        if (is.na(seg_dist_bin[t,s])) {
        } else if (max_val[t,s]==1&&seg_dist_bin[t,s]==0) {
          stay_max[t,s]=1
        } else if (max_val[t,s]==0&&seg_dist_bin[t,s]==0) {
          stay_nomax[t,s]=1
        } else if (max_val[t,s]==1&&seg_dist_bin[t,s]==1) {
          switch_max[t,s]=1
        } else if (max_val[t,s]==0&&seg_dist_bin[t,s]==1) {
          switch_nomax[t,s]=1
        } 
        if (t==1||is.na(seg_dist_bin[t-1,s])) {
        } else if (postmax_val[t,s]==1&&seg_dist_bin[t-1,s]==0) {
          stay_postmax[t,s]=1
        } else if (postmax_val[t,s]==0&&seg_dist_bin[t-1,s]==0) {
          stay_postnomax[t,s]=1
        } else if (postmax_val[t,s]==1&&seg_dist_bin[t-1,s]==1) {
          switch_postmax[t,s]=1
        } else if (postmax_val[t,s]==0&&seg_dist_bin[t-1,s]==1) {
          switch_postnomax[t,s]=1
        } 
        
        #outcome
        rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
        all_rews[t,s]=rew
        
        #update values
        out_args=do.call(value_functions[v_bayes],list(in_values=all_values[t,,s],
          in_alpha=all_alphas[t,,s],in_beta=all_betas[t,,s],rew=rew,
          chosen_seg=chosen_seg,seg_samplehx=seg_samplehx[t,,s],
          alpha=alpha_params[s],decay=decay_params[s],kappa=kappa_params[s]))
        if (value_nobayes[v_bayes]==0) {
          all_alphas[t+1,,s]=out_args[1:num_segments]
          all_betas[t+1,,s]=out_args[(num_segments+1):(2*num_segments)]
          all_values[t+1,,s]=all_alphas[t+1,,s]/(all_alphas[t+1,,s]+all_betas[t+1,,s])
        } else all_values[t+1,,s]=out_args
        
        out_info[(s-1)*num_trials+t,]=c(s,t,chosen_seg,seg_dist[t,s],
          seg_samplehx_selected[t,s],seg_probs[chosen_seg],rew,alpha_params[s],
          decay_params[s],kappa_params[s],epsilon_params[s],beta_params[s],
          tau_params[s],omega_params[s],seg_dist_bin[t,s],max_val[t,s],pers_exp_bins[t,s],
          stay_max[t,s],stay_nomax[t,s],switch_max[t,s],switch_nomax[t,s],
          stay_postmax[t,s],stay_postnomax[t,s],switch_postmax[t,s],switch_postnomax[t,s],
          all_values[t,,s])
      }
    }
    
    #create data for plotting
    out_info_df=data.frame(out_info)
    names(out_info_df)=c('subj','trial','chosen_seg','seg_dist','samplehx_selected',
      'selected_prob','rew','alpha','decay','kappa','epsilon','beta','tau','omega',
      'binary_seg_dist','max_value_chosen','pers_exp_bins','stay_max','stay_nomax',
      'switch_max','switch_nomax','stay_postmax','stay_postnomax','switch_postmax',
      'switch_postnomax','modelvalue_1','modelvalue_2','modelvalue_3','modelvalue_4')
    out_info_df$modelvaluediff_1=seg_probs[1]-out_info_df$modelvalue_1
    out_info_df$modelvaluediff_2=seg_probs[2]-out_info_df$modelvalue_2
    out_info_df$modelvaluediff_3=seg_probs[3]-out_info_df$modelvalue_3
    out_info_df$modelvaluediff_4=seg_probs[4]-out_info_df$modelvalue_4

    #make long format for with rows for each segment
    out_info_df$subj=factor(out_info_df$subj)
    out_info_long=pivot_longer(out_info_df,modelvalue_1:modelvaluediff_4,
      names_to=c('.value','segment'),names_sep="_")
    
    #plot
    ggsave(paste('plots/plots_04282020/sim_plot_Bayes',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_df,aes(trial,samplehx_selected,
        color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
        geom_smooth(method='gam')+geom_point(alpha=0.1)+theme_classic()+
        ylim(0,20)+labs(y='number of times option selected',
          title=paste(value_functions[v_bayes],choice_functions[c_bayes],
                      var_names[a],level_name_sel,sep=' '))+
        scale_color_discrete(name='Probability of Reward')+
        scale_fill_discrete(name='Probability of Reward')+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
    
    ggsave(paste('plots/plots_04282020/sim_plot_Bayes_segvalues_',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_long,aes(trial,modelvaluediff,color=as.factor(segment),
        fill=as.factor(segment)))+geom_hline(yintercept=0)+
        geom_smooth(method='loess')+theme_classic()+#geom_point(alpha=0.1)+
        ylim(-0.2,0.2)+labs(y='difference in actual vs. learned segment probability',
          title=paste(value_functions[v_bayes],choice_functions[c_bayes],
                      var_names[a],level_name_sel,sep=' '))+
        scale_color_discrete(name='Segment Probabilities',labels=seg_probs)+
        scale_fill_discrete(name='Segment Probabilities',labels=seg_probs)+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
    
      ggsave(paste('plots/plots_04282020/sim_plot_Bayes_segdist',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial,seg_dist,
          color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
          geom_smooth(method='gam')+theme_classic()+
          ylim(0,2)+labs(y='distance of selected segment from previous choice',
            title=paste(value_functions[v_bayes],choice_functions[c_bayes],
                        var_names[a],level_name_sel,sep=' '))+
          scale_color_discrete(name='Probability of Reward')+
          scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
      
      ggsave(paste('plots/plots_04282020/sim_plot_Bayes_expchoices',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial))+
          geom_smooth(aes(y=stay_max),method='gam',color='maroon',linetype='solid')+
          geom_smooth(aes(y=stay_nomax),method='gam',color='maroon',linetype='dashed')+
          geom_smooth(aes(y=switch_max),method='gam',color='blue',linetype='solid')+
          geom_smooth(aes(y=switch_nomax),method='gam',color='blue',linetype='dashed')+
          theme_classic()+
          ylim(0,1)+labs(y='Proportion choices',
            title=paste(value_functions[v_bayes],choice_functions[c_bayes],
                        var_names[a],level_name_sel,sep=' '))+
          #scale_linetype_discrete(name='Choice',labels=c('stay','switch'))+
          #scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
      
      ggsave(paste('plots/plots_04282020/sim_plot_Bayes_post_expchoices',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial))+
          geom_smooth(aes(y=stay_postmax),method='gam',color='maroon',linetype='solid')+
          geom_smooth(aes(y=stay_postnomax),method='gam',color='maroon',linetype='dashed')+
          geom_smooth(aes(y=switch_postmax),method='gam',color='blue',linetype='solid')+
          geom_smooth(aes(y=switch_postnomax),method='gam',color='blue',linetype='dashed')+
          theme_classic()+
          ylim(0,1)+labs(y='Proportion choices',
            title=paste(value_functions[v_bayes],choice_functions[c_bayes],
                        var_names[a],level_name_sel,sep=' '))+
          #scale_linetype_discrete(name='Choice',labels=c('stay','switch'))+
          #scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
      
      ggsave(paste('plots/plots_04282020/sim_plot_Bayes_value_',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_df,aes(trial,rew))+
        geom_smooth(method='gam')+theme_classic()+
        coord_cartesian(ylim=c(0.4,0.6))+labs(y='Proportion Rewarded Choices',
          title=paste(value_functions[v_bayes],choice_functions[c_bayes],
                      var_names[a],level_name_sel,sep=' '))+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
      
    #set up variables: Bayesian model w/NO Thompson sampling
    all_alphas=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    all_betas=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    all_values=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    seg_samplehx=array(data=NA,dim=c(num_trials+1,num_segments,num_subjs))
    seg_samplehx_selected=array(data=NA,dim=c(num_trials,num_subjs))
    all_choices=array(data=NA,dim=c(num_trials,num_subjs))
    all_rews=array(data=NA,dim=c(num_trials,num_subjs))
    seg_dist=array(data=NA,dim=c(num_trials,num_subjs))
    seg_dist_bin=array(data=NA,dim=c(num_trials,num_subjs))
    pers_exp_bins=array(data=NA,dim=c(num_trials,num_subjs))
    max_val=array(data=NA,dim=c(num_trials,num_subjs))
    stay_max=array(data=0,dim=c(num_trials,num_subjs))
    stay_nomax=array(data=0,dim=c(num_trials,num_subjs))
    switch_max=array(data=0,dim=c(num_trials,num_subjs))
    switch_nomax=array(data=0,dim=c(num_trials,num_subjs))
    postmax_val=array(data=NA,dim=c(num_trials,num_subjs))
    stay_postmax=array(data=0,dim=c(num_trials,num_subjs))
    stay_postnomax=array(data=0,dim=c(num_trials,num_subjs))
    switch_postmax=array(data=0,dim=c(num_trials,num_subjs))
    switch_postnomax=array(data=0,dim=c(num_trials,num_subjs))
    out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=29)
    
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
          if (value_nobayes[v_bayes]==0) {
            all_values[t,,s]=all_alphas[t,,s]/(all_alphas[t,,s]+all_betas[t,,s])
          }
          chosen_seg=do.call(choice_functions[c_bayes_noT],list(in_values=all_values[t,,s],
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
        seg_dist[t,s]=ifelse(t==1,NA,seg_dist_mat[chosen_seg,all_choices[t-1,s]])
        seg_dist_bin[t,s]=ifelse(seg_dist[t,s]==0,0,1)
        max_val[t,s]=ifelse(chosen_seg %in% which(all_values[t,,s]==max(all_values[t,,s])),1,0)
        postmax_val[t,s]=ifelse(t==1,NA,
          ifelse(all_choices[t-1,s] %in% which(all_values[t-1,,s]==max(all_values[t-1,,s])),1,0))
        pers_exp_bins[t,s]=2*seg_dist_bin[t,s]+max_val[t,s]
        if (is.na(seg_dist_bin[t,s])) {
        } else if (max_val[t,s]==1&&seg_dist_bin[t,s]==0) {
          stay_max[t,s]=1
        } else if (max_val[t,s]==0&&seg_dist_bin[t,s]==0) {
          stay_nomax[t,s]=1
        } else if (max_val[t,s]==1&&seg_dist_bin[t,s]==1) {
          switch_max[t,s]=1
        } else if (max_val[t,s]==0&&seg_dist_bin[t,s]==1) {
          switch_nomax[t,s]=1
        } 
        if (t==1||is.na(seg_dist_bin[t-1,s])) {
        } else if (postmax_val[t,s]==1&&seg_dist_bin[t-1,s]==0) {
          stay_postmax[t,s]=1
        } else if (postmax_val[t,s]==0&&seg_dist_bin[t-1,s]==0) {
          stay_postnomax[t,s]=1
        } else if (postmax_val[t,s]==1&&seg_dist_bin[t-1,s]==1) {
          switch_postmax[t,s]=1
        } else if (postmax_val[t,s]==0&&seg_dist_bin[t-1,s]==1) {
          switch_postnomax[t,s]=1
        } 
        
        #outcome
        rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
        all_rews[t,s]=rew
        
        #update values
        out_args=do.call(value_functions[v_bayes],list(in_values=all_values[t,,s],
          in_alpha=all_alphas[t,,s],in_beta=all_betas[t,,s],rew=rew,
          chosen_seg=chosen_seg,seg_samplehx=seg_samplehx[t,,s],
          alpha=alpha_params[s],decay=decay_params[s],kappa=kappa_params[s]))
        if (value_nobayes[v_bayes]==0) {
          all_alphas[t+1,,s]=out_args[1:num_segments]
          all_betas[t+1,,s]=out_args[(num_segments+1):(2*num_segments)]
          all_values[t+1,,s]=all_alphas[t+1,,s]/(all_alphas[t+1,,s]+all_betas[t+1,,s])
        } else all_values[t+1,,s]=out_args
        
        out_info[(s-1)*num_trials+t,]=c(s,t,chosen_seg,seg_dist[t,s],
          seg_samplehx_selected[t,s],seg_probs[chosen_seg],rew,alpha_params[s],
          decay_params[s],kappa_params[s],epsilon_params[s],beta_params[s],
          tau_params[s],omega_params[s],seg_dist_bin[t,s],max_val[t,s],pers_exp_bins[t,s],
          stay_max[t,s],stay_nomax[t,s],switch_max[t,s],switch_nomax[t,s],
          stay_postmax[t,s],stay_postnomax[t,s],switch_postmax[t,s],switch_postnomax[t,s],
          all_values[t,,s])
      }
    }
    
    #create data for plotting
    out_info_df=data.frame(out_info)
    names(out_info_df)=c('subj','trial','chosen_seg','seg_dist','samplehx_selected',
      'selected_prob','rew','alpha','decay','kappa','epsilon','beta','tau','omega',
      'binary_seg_dist','max_value_chosen','pers_exp_bins','stay_max','stay_nomax',
      'switch_max','switch_nomax','stay_postmax','stay_postnomax','switch_postmax',
      'switch_postnomax','modelvalue_1','modelvalue_2','modelvalue_3','modelvalue_4')
    out_info_df$modelvaluediff_1=seg_probs[1]-out_info_df$modelvalue_1
    out_info_df$modelvaluediff_2=seg_probs[2]-out_info_df$modelvalue_2
    out_info_df$modelvaluediff_3=seg_probs[3]-out_info_df$modelvalue_3
    out_info_df$modelvaluediff_4=seg_probs[4]-out_info_df$modelvalue_4

    #make long format for with rows for each segment
    out_info_df$subj=factor(out_info_df$subj)
    out_info_long=pivot_longer(out_info_df,modelvalue_1:modelvaluediff_4,
      names_to=c('.value','segment'),names_sep="_")
    
    #plot
    ggsave(paste('plots/plots_04282020/sim_plot_Bayes_noT',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_df,aes(trial,samplehx_selected,
        color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
        geom_smooth(method='gam')+geom_point(alpha=0.1)+theme_classic()+
        ylim(0,20)+labs(y='number of times option selected',
          title=paste(value_functions[v_bayes],choice_functions[c_bayes_noT],
                      var_names[a],level_name_sel,sep=' '))+
        scale_color_discrete(name='Probability of Reward')+
        scale_fill_discrete(name='Probability of Reward')+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
    
    ggsave(paste('plots/plots_04282020/sim_plot_Bayes_noT_segvalues_',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_long,aes(trial,modelvaluediff,color=as.factor(segment),
        fill=as.factor(segment)))+geom_hline(yintercept=0)+
        geom_smooth(method='loess')+theme_classic()+#geom_point(alpha=0.1)+
        ylim(-0.2,0.2)+labs(y='difference in actual vs. learned segment probability',
          title=paste(value_functions[v_bayes],choice_functions[c_bayes_noT],
                      var_names[a],level_name_sel,sep=' '))+
        scale_color_discrete(name='Segment Probabilities',labels=seg_probs)+
        scale_fill_discrete(name='Segment Probabilities',labels=seg_probs)+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
    
      ggsave(paste('plots/plots_04282020/sim_plot_Bayes_noT_segdist',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial,seg_dist,
          color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
          geom_smooth(method='gam')+theme_classic()+
          ylim(0,2)+labs(y='distance of selected segment from previous choice',
            title=paste(value_functions[v_bayes],choice_functions[c_bayes_noT],
                        var_names[a],level_name_sel,sep=' '))+
          scale_color_discrete(name='Probability of Reward')+
          scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
      
      ggsave(paste('plots/plots_04282020/sim_plot_Bayes_noT_expchoices',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial))+
          geom_smooth(aes(y=stay_max),method='gam',color='maroon',linetype='solid')+
          geom_smooth(aes(y=stay_nomax),method='gam',color='maroon',linetype='dashed')+
          geom_smooth(aes(y=switch_max),method='gam',color='blue',linetype='solid')+
          geom_smooth(aes(y=switch_nomax),method='gam',color='blue',linetype='dashed')+
          theme_classic()+
          ylim(0,1)+labs(y='Proportion choices',
            title=paste(value_functions[v_bayes],choice_functions[c_bayes_noT],
                        var_names[a],level_name_sel,sep=' '))+
          #scale_linetype_discrete(name='Choice',labels=c('stay','switch'))+
          #scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
      
      ggsave(paste('plots/plots_04282020/sim_plot_Bayes_noT_post_expchoices',var_names[a],level,num_trials,'.png',sep='_'),
        ggplot(out_info_df,aes(trial))+
          geom_smooth(aes(y=stay_postmax),method='gam',color='maroon',linetype='solid')+
          geom_smooth(aes(y=stay_postnomax),method='gam',color='maroon',linetype='dashed')+
          geom_smooth(aes(y=switch_postmax),method='gam',color='blue',linetype='solid')+
          geom_smooth(aes(y=switch_postnomax),method='gam',color='blue',linetype='dashed')+
          theme_classic()+
          ylim(0,1)+labs(y='Proportion choices',
            title=paste(value_functions[v_bayes],choice_functions[c_bayes_noT],
                        var_names[a],level_name_sel,sep=' '))+
          #scale_linetype_discrete(name='Choice',labels=c('stay','switch'))+
          #scale_fill_discrete(name='Probability of Reward')+
          theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
      
      ggsave(paste('plots/plots_04282020/sim_plot_Bayes_noT_value_',var_names[a],level,num_trials,'.png',sep='_'),
      ggplot(out_info_df,aes(trial,rew))+
        geom_smooth(method='gam')+theme_classic()+
        coord_cartesian(ylim=c(0.4,0.6))+labs(y='Proportion Rewarded Choices',
          title=paste(value_functions[v_bayes],choice_functions[c_bayes_noT],
                      var_names[a],level_name_sel,sep=' '))+
        theme(legend.position="bottom"),
      width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
    
  }
}

