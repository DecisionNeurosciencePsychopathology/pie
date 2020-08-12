# simulate data, estimate parameters from simulated data, and assess parameter recovery

# set up ####
source('sim_functions.R')
library(ggplot2)
# library(rstan)
# library(shinystan)
# library(loo)
# options(mc.cores=2)
# rstan_options(auto_write=TRUE)

num_segments=8
num_trials_per_block=30+num_segments 
num_blocks=4
num_trials=num_trials_per_block*num_blocks
num_subjs=1000
if (num_segments==4) {
  seg_probs=c(0.36,0.43,0.56,0.65)
} else if (num_segments==8) {
  seg_probs=c(0.34,0.38,0.42,0.46,0.51,0.56,0.64,0.69)
}

out_compare=array(data=NA,dim=c(6*35,9))
count=1

# parameters: right now assume fixed across subjects
alpha=0.5
beta=5
kappa=1
chi=c(3,6,9)
# decay=0.95

seg_dist_mat4=cbind(c(0,1,2,1),c(1,0,1,2),c(2,1,0,1),c(1,2,1,0))
seg_dist_mat8=cbind(c(0,1,2,3,4,3,2,1),c(1,0,1,2,3,4,3,2),c(2,1,0,1,2,3,4,3),
                    c(3,2,1,0,1,2,3,4),c(4,3,2,1,0,1,2,3),c(3,4,3,2,1,0,1,2),
                    c(2,3,4,3,2,1,0,1),c(1,2,3,4,3,2,1,0))

num_segments_in=array(data=num_segments,dim=c(num_subjs,num_trials))
points_shown=array(data=0,dim=c(num_subjs,num_trials))
# block_num=array(data=1,dim=c(num_subjs,num_trials))

for (capacity in c(4,5,6,7,8)) {
  for (decay in c(1,0.95,0.9,0.8,0.7,0.6,0.5)) {
  
  print(paste('running models with capacity of',capacity,'and decay of',decay,sep=' '))
  
  # model: basic RL + softmax with exploration bonus ####
  out_compare[count,1:6]=c(10,capacity,decay,num_subjs,num_trials,num_blocks)
    
  all_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  all_values_pre=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  all_choices=array(data=NA,dim=c(num_subjs,num_trials))
  all_rews=array(data=NA,dim=c(num_subjs,num_trials))
  seg_samplehx=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  seg_dist=array(data=NA,dim=c(num_subjs,num_trials))
  seg_dist_bin=array(data=NA,dim=c(num_subjs,num_trials))
  max_val=array(data=NA,dim=c(num_subjs,num_trials))
  stay_max=array(data=0,dim=c(num_subjs,num_trials))
  stay_nomax=array(data=0,dim=c(num_subjs,num_trials))
  switch_max=array(data=0,dim=c(num_subjs,num_trials))
  switch_nomax=array(data=0,dim=c(num_subjs,num_trials))
  seg_samplehx_selected=array(data=NA,dim=c(num_subjs,num_trials))
  pers_exp_bins=array(data=NA,dim=c(num_subjs,num_trials))
  postmax_val=array(data=NA,dim=c(num_subjs,num_trials))
  stay_postmax=array(data=0,dim=c(num_subjs,num_trials))
  stay_postnomax=array(data=0,dim=c(num_subjs,num_trials))
  switch_postmax=array(data=0,dim=c(num_subjs,num_trials))
  switch_postnomax=array(data=0,dim=c(num_subjs,num_trials))
  block_num=array(data=NA,dim=c(num_subjs,num_trials))
  
  out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=20)
  
  for (s in 1:num_subjs) {
    all_values[s,1,]=0.5
    seg_samplehx[s,1,]=0
    
    block_tnum=0
    block=1
    # if (num_segments<4) seg_probs=block_seg_probs[block,]
    
    for (t in 1:num_trials) {
      if (block_tnum>num_trials_per_block) {
        block_tnum=1
        block=block+1
        # if (num_segments<4) seg_probs=block_seg_probs[block,]
        all_values[s,t,]=0.5
        seg_samplehx[s,t,]=0
      } else block_tnum=block_tnum+1
      block_num[s,t]=block
      
      #choice
      if (block_tnum<(num_segments+1)) { #forced sampling: assume even for now
        chosen_seg=block_tnum
      } else {
        # chosen_seg=softmax_sim(in_values=all_values[s,t,],beta=beta)
        chosen_seg=softmax_SBalt_sim(in_values=all_values[s,t,],
          seg_samplehx=seg_samplehx[s,t,],beta=beta,kappa=kappa)
      }
      
      #update variables
      all_choices[s,t]=chosen_seg
      seg_samplehx[s,t+1,chosen_seg]=seg_samplehx[s,t,chosen_seg]+1
      seg_samplehx[s,t+1,-chosen_seg]=seg_samplehx[s,t,-chosen_seg]
      seg_samplehx_selected[s,t]=seg_samplehx[s,t+1,chosen_seg]
      seg_dist[s,t]=ifelse(t==1,NA,seg_dist_mat8[chosen_seg,all_choices[s,t-1]])
      seg_dist_bin[s,t]=ifelse(seg_dist[s,t]==0,0,1)
      max_val[s,t]=ifelse(chosen_seg %in% which(all_values[s,t,]==max(all_values[s,t,])),1,0)
      
      if (is.na(seg_dist_bin[s,t])) {
      } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==0) {
        stay_max[s,t]=1
      } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==0) {
        stay_nomax[s,t]=1
      } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==1) {
        switch_max[s,t]=1
      } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==1) {
        switch_nomax[s,t]=1
      } 
      
      #outcome
      rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
      all_rews[s,t]=rew
      
      #update values
      all_values_pre[s,t+1,]=basic_RL_sim(in_values=all_values[s,t,],rew=rew,
        chosen_seg=chosen_seg,alpha=alpha)
      all_values[s,t+1,]=decay*all_values_pre[s,t+1,]
      
      out_info[(s-1)*num_trials+t,]=c(s,t,block,block_tnum,chosen_seg,seg_dist[s,t],
                                      seg_samplehx_selected[s,t],
        seg_probs[chosen_seg],rew,seg_dist_bin[s,t],max_val[s,t],pers_exp_bins[s,t],
        stay_max[s,t],stay_nomax[s,t],switch_max[s,t],switch_nomax[s,t],stay_postmax[s,t],
        stay_postnomax[s,t],switch_postmax[s,t],switch_postnomax[s,t])
    }
  }
  
  out_info_df=data.frame(out_info)
  names(out_info_df)=c('subj','trial','block','trial_in_block','chosen_seg',
                       'seg_dist','samplehx_selected',
    'selected_prob','rew','binary_seg_dist','max_value_chosen','pers_exp_bins',
    'stay_max','stay_nomax','switch_max','switch_nomax','stay_postmax',
    'stay_postnomax','switch_postmax','switch_postnomax')
  out_info_fc=out_info_df[out_info_df$trial_in_block>num_segments,]
  
  out_compare[count,7]=sum(out_info_fc$rew)/(num_trials*num_subjs)
  out_compare[count,8]=sum(out_info_fc$selected_prob>.5)/(num_trials*num_subjs)
  out_compare[count,9]=sum(out_info_fc$max_value_chosen)/(num_trials*num_subjs)
  count=count+1
  
  # model: capacity limited ####
  out_compare[count,1:6]=c(1,capacity,decay,num_subjs,num_trials,num_blocks)
    
  all_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  all_values_pre=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  all_choices=array(data=NA,dim=c(num_subjs,num_trials))
  all_rews=array(data=NA,dim=c(num_subjs,num_trials))
  seg_samplehx=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  seg_dist=array(data=NA,dim=c(num_subjs,num_trials))
  seg_dist_bin=array(data=NA,dim=c(num_subjs,num_trials))
  max_val=array(data=NA,dim=c(num_subjs,num_trials))
  stay_max=array(data=0,dim=c(num_subjs,num_trials))
  stay_nomax=array(data=0,dim=c(num_subjs,num_trials))
  switch_max=array(data=0,dim=c(num_subjs,num_trials))
  switch_nomax=array(data=0,dim=c(num_subjs,num_trials))
  seg_samplehx_selected=array(data=NA,dim=c(num_subjs,num_trials))
  pers_exp_bins=array(data=NA,dim=c(num_subjs,num_trials))
  postmax_val=array(data=NA,dim=c(num_subjs,num_trials))
  stay_postmax=array(data=0,dim=c(num_subjs,num_trials))
  stay_postnomax=array(data=0,dim=c(num_subjs,num_trials))
  switch_postmax=array(data=0,dim=c(num_subjs,num_trials))
  switch_postnomax=array(data=0,dim=c(num_subjs,num_trials))
  block_num=array(data=NA,dim=c(num_subjs,num_trials))
  
  out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=20)
  
  for (s in 1:num_subjs) {
    all_values[s,1,]=0.5
    seg_samplehx[s,1,]=0
    
    block_tnum=0
    block=1
    # if (num_segments<4) seg_probs=block_seg_probs[block,]
    
    for (t in 1:num_trials) {
      if (block_tnum>num_trials_per_block) {
        block_tnum=1
        block=block+1
        # if (num_segments<4) seg_probs=block_seg_probs[block,]
        all_values[s,t,]=0.5
        seg_samplehx[s,t,]=0
      } else block_tnum=block_tnum+1
      block_num[s,t]=block
      
      #choice
      if (block_tnum<(num_segments+1)) { #forced sampling: assume even for now
        chosen_seg=block_tnum
      } else {
        # chosen_seg=softmax_sim(in_values=all_values[s,t,],beta=beta)
        chosen_seg=softmax_SBalt_sim(in_values=all_values[s,t,],
          seg_samplehx=seg_samplehx[s,t,],beta=beta,kappa=kappa)
      }
      
      #update variables
      all_choices[s,t]=chosen_seg
      seg_samplehx[s,t+1,chosen_seg]=seg_samplehx[s,t,chosen_seg]+1
      seg_samplehx[s,t+1,-chosen_seg]=seg_samplehx[s,t,-chosen_seg]
      seg_samplehx_selected[s,t]=seg_samplehx[s,t+1,chosen_seg]
      seg_dist[s,t]=ifelse(t==1,NA,seg_dist_mat8[chosen_seg,all_choices[s,t-1]])
      seg_dist_bin[s,t]=ifelse(seg_dist[s,t]==0,0,1)
      max_val[s,t]=ifelse(chosen_seg %in% which(all_values[s,t,]==max(all_values[s,t,])),1,0)
      
      if (is.na(seg_dist_bin[s,t])) {
      } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==0) {
        stay_max[s,t]=1
      } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==0) {
        stay_nomax[s,t]=1
      } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==1) {
        switch_max[s,t]=1
      } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==1) {
        switch_nomax[s,t]=1
      } 
      
      #outcome
      rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
      all_rews[s,t]=rew
      
      #update values
      all_values_pre[s,t+1,]=basic_RL_sim(in_values=all_values[s,t,],rew=rew,
        chosen_seg=chosen_seg,alpha=alpha)
      all_values[s,t+1,]=decay*(capacity/num_segments)*all_values_pre[s,t+1,]
      
      out_info[(s-1)*num_trials+t,]=c(s,t,block,block_tnum,chosen_seg,seg_dist[s,t],
                                      seg_samplehx_selected[s,t],
        seg_probs[chosen_seg],rew,seg_dist_bin[s,t],max_val[s,t],pers_exp_bins[s,t],
        stay_max[s,t],stay_nomax[s,t],switch_max[s,t],switch_nomax[s,t],stay_postmax[s,t],
        stay_postnomax[s,t],switch_postmax[s,t],switch_postnomax[s,t])
    }
  }
  
  out_info_df=data.frame(out_info)
  names(out_info_df)=c('subj','trial','block','trial_in_block','chosen_seg',
                       'seg_dist','samplehx_selected',
    'selected_prob','rew','binary_seg_dist','max_value_chosen','pers_exp_bins',
    'stay_max','stay_nomax','switch_max','switch_nomax','stay_postmax',
    'stay_postnomax','switch_postmax','switch_postnomax')
  out_info_fc=out_info_df[out_info_df$trial_in_block>num_segments,]
  
  out_compare[count,7]=sum(out_info_fc$rew)/(num_trials*num_subjs)
  out_compare[count,8]=sum(out_info_fc$selected_prob>.5)/(num_trials*num_subjs)
  out_compare[count,9]=sum(out_info_fc$max_value_chosen)/(num_trials*num_subjs)
  count=count+1
  
  # additive limited capacity model ####
  out_compare[count,1:6]=c(2,capacity,decay,num_subjs,num_trials,num_blocks)
    
  all_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  all_values_pre=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  all_choices=array(data=NA,dim=c(num_subjs,num_trials))
  all_rews=array(data=NA,dim=c(num_subjs,num_trials))
  seg_samplehx=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  seg_dist=array(data=NA,dim=c(num_subjs,num_trials))
  seg_dist_bin=array(data=NA,dim=c(num_subjs,num_trials))
  max_val=array(data=NA,dim=c(num_subjs,num_trials))
  stay_max=array(data=0,dim=c(num_subjs,num_trials))
  stay_nomax=array(data=0,dim=c(num_subjs,num_trials))
  switch_max=array(data=0,dim=c(num_subjs,num_trials))
  switch_nomax=array(data=0,dim=c(num_subjs,num_trials))
  seg_samplehx_selected=array(data=NA,dim=c(num_subjs,num_trials))
  pers_exp_bins=array(data=NA,dim=c(num_subjs,num_trials))
  postmax_val=array(data=NA,dim=c(num_subjs,num_trials))
  stay_postmax=array(data=0,dim=c(num_subjs,num_trials))
  stay_postnomax=array(data=0,dim=c(num_subjs,num_trials))
  switch_postmax=array(data=0,dim=c(num_subjs,num_trials))
  switch_postnomax=array(data=0,dim=c(num_subjs,num_trials))
  block_num=array(data=NA,dim=c(num_subjs,num_trials))
  
  out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=20)
  
  for (s in 1:num_subjs) {
    all_values[s,1,]=0.5
    seg_samplehx[s,1,]=0
    
    block_tnum=0
    block=1
    # if (num_segments<4) seg_probs=block_seg_probs[block,]
    
    for (t in 1:num_trials) {
      if (block_tnum>num_trials_per_block) {
        block_tnum=1
        block=block+1
        # if (num_segments<4) seg_probs=block_seg_probs[block,]
        all_values[s,t,]=0.5
        seg_samplehx[s,t,]=0
      } else block_tnum=block_tnum+1
      block_num[s,t]=block
      
      #choice
      if (block_tnum<(num_segments+1)) { #forced sampling: assume even for now
        chosen_seg=block_tnum
      } else {
        # chosen_seg=softmax_sim(in_values=all_values[s,t,],beta=beta)
        chosen_seg=softmax_SBalt_sim(in_values=all_values[s,t,],
          seg_samplehx=seg_samplehx[s,t,],beta=beta,kappa=kappa)
      }
      
      #update variables
      all_choices[s,t]=chosen_seg
      seg_samplehx[s,t+1,chosen_seg]=seg_samplehx[s,t,chosen_seg]+1
      seg_samplehx[s,t+1,-chosen_seg]=seg_samplehx[s,t,-chosen_seg]
      seg_samplehx_selected[s,t]=seg_samplehx[s,t+1,chosen_seg]
      seg_dist[s,t]=ifelse(t==1,NA,seg_dist_mat8[chosen_seg,all_choices[s,t-1]])
      seg_dist_bin[s,t]=ifelse(seg_dist[s,t]==0,0,1)
      max_val[s,t]=ifelse(chosen_seg %in% which(all_values[s,t,]==max(all_values[s,t,])),1,0)
      
      if (is.na(seg_dist_bin[s,t])) {
      } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==0) {
        stay_max[s,t]=1
      } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==0) {
        stay_nomax[s,t]=1
      } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==1) {
        switch_max[s,t]=1
      } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==1) {
        switch_nomax[s,t]=1
      } 
      
      #outcome
      rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
      all_rews[s,t]=rew
      
      #update values
      all_values_pre[s,t+1,]=basic_RL_sim(in_values=all_values[s,t,],rew=rew,
        chosen_seg=chosen_seg,alpha=alpha)
      all_values[s,t+1,]=.5*decay*all_values_pre[s,t+1,]+
        .5*(capacity^decay/num_segments^decay)*all_values_pre[s,t+1,]
      
      out_info[(s-1)*num_trials+t,]=c(s,t,block,block_tnum,chosen_seg,seg_dist[s,t],
                                      seg_samplehx_selected[s,t],
        seg_probs[chosen_seg],rew,seg_dist_bin[s,t],max_val[s,t],pers_exp_bins[s,t],
        stay_max[s,t],stay_nomax[s,t],switch_max[s,t],switch_nomax[s,t],stay_postmax[s,t],
        stay_postnomax[s,t],switch_postmax[s,t],switch_postnomax[s,t])
    }
  }
  
  out_info_df=data.frame(out_info)
  names(out_info_df)=c('subj','trial','block','trial_in_block','chosen_seg',
                       'seg_dist','samplehx_selected',
    'selected_prob','rew','binary_seg_dist','max_value_chosen','pers_exp_bins',
    'stay_max','stay_nomax','switch_max','switch_nomax','stay_postmax',
    'stay_postnomax','switch_postmax','switch_postnomax')
  out_info_fc=out_info_df[out_info_df$trial_in_block>num_segments,]
  
  out_compare[count,7]=sum(out_info_fc$rew)/(num_trials*num_subjs)
  out_compare[count,8]=sum(out_info_fc$selected_prob>.5)/(num_trials*num_subjs)
  out_compare[count,9]=sum(out_info_fc$max_value_chosen)/(num_trials*num_subjs)
  count=count+1
  
  # compressed capacity model ####
  for (chi in c(3,6,9)) {
  out_compare[count,1:6]=c(chi,capacity,decay,num_subjs,num_trials,num_blocks)
  
  # if (chi==3) {
  #   beta_alt=5*beta
  # } else {
    beta_alt=beta
  # }
    
  all_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  all_values_pre=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  all_choices=array(data=NA,dim=c(num_subjs,num_trials))
  all_rews=array(data=NA,dim=c(num_subjs,num_trials))
  seg_samplehx=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
  seg_dist=array(data=NA,dim=c(num_subjs,num_trials))
  seg_dist_bin=array(data=NA,dim=c(num_subjs,num_trials))
  max_val=array(data=NA,dim=c(num_subjs,num_trials))
  stay_max=array(data=0,dim=c(num_subjs,num_trials))
  stay_nomax=array(data=0,dim=c(num_subjs,num_trials))
  switch_max=array(data=0,dim=c(num_subjs,num_trials))
  switch_nomax=array(data=0,dim=c(num_subjs,num_trials))
  seg_samplehx_selected=array(data=NA,dim=c(num_subjs,num_trials))
  pers_exp_bins=array(data=NA,dim=c(num_subjs,num_trials))
  postmax_val=array(data=NA,dim=c(num_subjs,num_trials))
  stay_postmax=array(data=0,dim=c(num_subjs,num_trials))
  stay_postnomax=array(data=0,dim=c(num_subjs,num_trials))
  switch_postmax=array(data=0,dim=c(num_subjs,num_trials))
  switch_postnomax=array(data=0,dim=c(num_subjs,num_trials))
  block_num=array(data=NA,dim=c(num_subjs,num_trials))
  
  out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=20)
  
  for (s in 1:num_subjs) {
    all_values[s,1,]=0.5
    seg_samplehx[s,1,]=0
    
    block_tnum=0
    block=1
    # if (num_segments<4) seg_probs=block_seg_probs[block,]
    
    for (t in 1:num_trials) {
      if (block_tnum>num_trials_per_block) {
        block_tnum=1
        block=block+1
        # if (num_segments<4) seg_probs=block_seg_probs[block,]
        all_values[s,t,]=0.5
        seg_samplehx[s,t,]=0
      } else block_tnum=block_tnum+1
      block_num[s,t]=block
      
      #choice
      if (block_tnum<(num_segments+1)) { #forced sampling: assume even for now
        chosen_seg=block_tnum
      } else {
        # chosen_seg=softmax_sim(in_values=all_values[s,t,],beta=beta)
        chosen_seg=softmax_SBalt_sim(in_values=all_values[s,t,],
          seg_samplehx=seg_samplehx[s,t,],beta=beta_alt,kappa=kappa)
      }
      
      #update variables
      all_choices[s,t]=chosen_seg
      seg_samplehx[s,t+1,chosen_seg]=seg_samplehx[s,t,chosen_seg]+1
      seg_samplehx[s,t+1,-chosen_seg]=seg_samplehx[s,t,-chosen_seg]
      seg_samplehx_selected[s,t]=seg_samplehx[s,t+1,chosen_seg]
      seg_dist[s,t]=ifelse(t==1,NA,seg_dist_mat8[chosen_seg,all_choices[s,t-1]])
      seg_dist_bin[s,t]=ifelse(seg_dist[s,t]==0,0,1)
      max_val[s,t]=ifelse(chosen_seg %in% which(all_values[s,t,]==max(all_values[s,t,])),1,0)
      
      if (is.na(seg_dist_bin[s,t])) {
      } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==0) {
        stay_max[s,t]=1
      } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==0) {
        stay_nomax[s,t]=1
      } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==1) {
        switch_max[s,t]=1
      } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==1) {
        switch_nomax[s,t]=1
      } 
      
      #outcome
      rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
      all_rews[s,t]=rew
      
      #update values
      all_values_pre[s,t+1,]=basic_RL_sim(in_values=all_values[s,t,],rew=rew,
        chosen_seg=chosen_seg,alpha=alpha)
      all_values[s,t+1,]=decay*(min(num_segments,(capacity*(chi^0.1)))/(num_segments))*
        (1/(1+exp(-1*chi*(all_values_pre[s,t+1,]-0.5))))
      
      out_info[(s-1)*num_trials+t,]=c(s,t,block,block_tnum,chosen_seg,seg_dist[s,t],
                                      seg_samplehx_selected[s,t],
        seg_probs[chosen_seg],rew,seg_dist_bin[s,t],max_val[s,t],pers_exp_bins[s,t],
        stay_max[s,t],stay_nomax[s,t],switch_max[s,t],switch_nomax[s,t],stay_postmax[s,t],
        stay_postnomax[s,t],switch_postmax[s,t],switch_postnomax[s,t])
    }
  }
  
  out_info_df=data.frame(out_info)
  names(out_info_df)=c('subj','trial','block','trial_in_block','chosen_seg',
                       'seg_dist','samplehx_selected',
    'selected_prob','rew','binary_seg_dist','max_value_chosen','pers_exp_bins',
    'stay_max','stay_nomax','switch_max','switch_nomax','stay_postmax',
    'stay_postnomax','switch_postmax','switch_postnomax')
  out_info_fc=out_info_df[out_info_df$trial_in_block>num_segments,]
  
  out_compare[count,7]=sum(out_info_fc$rew)/(num_trials*num_subjs)
  out_compare[count,8]=sum(out_info_fc$selected_prob>.5)/(num_trials*num_subjs)
  out_compare[count,9]=sum(out_info_fc$max_value_chosen)/(num_trials*num_subjs)
  count=count+1
  }
  
 # #plot
 #  ggsave(paste('sim_plot_RL_',capacity,'.png',sep=''),
 #    ggplot(out_info_fc,aes(trial_in_block,samplehx_selected,
 #      color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
 #      geom_smooth(method='gam')+geom_point(alpha=0.1)+theme_classic()+
 #      ylim(0,15)+labs(y='number of times option selected',
 #        title=paste0('segment capacity=',capacity,'.'))+
 #      scale_color_discrete(name='Probability of Reward')+
 #      scale_fill_discrete(name='Probability of Reward')+
 #      theme(legend.position="bottom"),
 #    width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
 #  
 #  ggsave(paste('sim_plot_RL_expchoices_',capacity,'.png',sep=''),
 #      ggplot(out_info_fc,aes(trial_in_block))+
 #        geom_smooth(aes(y=stay_max),method='gam',color='maroon',linetype='solid')+
 #        geom_smooth(aes(y=stay_nomax),method='gam',color='maroon',linetype='dashed')+
 #        geom_smooth(aes(y=switch_max),method='gam',color='blue',linetype='solid')+
 #        geom_smooth(aes(y=switch_nomax),method='gam',color='blue',linetype='dashed')+
 #        theme_classic()+
 #        ylim(0,1)+labs(y='Proportion choices',
 #          title=paste0('segment capacity=',capacity,'.'))+
 #        #scale_linetype_discrete(name='Choice',labels=c('stay','switch'))+
 #        #scale_fill_discrete(name='Probability of Reward')+
 #        theme(legend.position="bottom"),
 #    width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
 #    
    # ggsave(paste('sim_plot_RL_post_expchoices_',capacity,'.png',sep=''),
    #   ggplot(out_info_fc,aes(trial))+
    #     geom_smooth(aes(y=stay_postmax),method='gam',color='maroon',linetype='solid')+
    #     geom_smooth(aes(y=stay_postnomax),method='gam',color='maroon',linetype='dashed')+
    #     geom_smooth(aes(y=switch_postmax),method='gam',color='blue',linetype='solid')+
    #     geom_smooth(aes(y=switch_postnomax),method='gam',color='blue',linetype='dashed')+
    #     theme_classic()+
    #     ylim(0,1)+labs(y='Proportion choices',
    #       title=paste0('segment capacity=',capacity,'.'))+
    #     #scale_linetype_discrete(name='Choice',labels=c('stay','switch'))+
    #     #scale_fill_discrete(name='Probability of Reward')+
    #     theme(legend.position="bottom"),
    # width=2.5,height=2.5,dpi=300,scale=2,units="in",device='png')
  }
}

out_compare_df=data.frame(out_compare)
names(out_compare_df)=c('model','capacity','decay','subjects','trials','blocks',
                        'reward_per_trial','proportion_trials_prob_over_50',
                        'proportion_trials_max_prob')
out_compare_df$decayf=as.factor(out_compare_df$decay)
out_compare_df$model=as.factor(out_compare_df$model)
levels(out_compare_df$model)=c('capacity limited','adaptive capacity limited',
  'compressed capacity (3)','compressed capacity (6)','compressed capacity (9)',
  'standard RL')

ggplot(out_compare_df,aes(capacity,reward_per_trial,color=decayf,fill=decayf))+
  geom_point()+geom_smooth(method='gam',alpha=0.1)+
  facet_wrap(~model)+theme_bw()+labs(y='reward per trial',color='decay',fill='decay')

ggplot(out_compare_df,aes(capacity,proportion_trials_prob_over_50,color=decayf,fill=decayf))+
  geom_point()+geom_smooth(method='gam',alpha=0.1)+
  facet_wrap(~model)+theme_bw()+
  labs(y='proportion of choices with p(reward)>.5',color='decay',fill='decay')

ggplot(out_compare_df,aes(capacity,proportion_trials_max_prob,color=decayf,fill=decayf))+
  geom_point()+geom_smooth(method='gam',alpha=0.1)+
  facet_wrap(~model)+theme_bw()+
  labs(y='proportion of choices with max(prob)',color='decay',fill='decay')
