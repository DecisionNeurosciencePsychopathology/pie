# set up ####

source('sim_functions.R')
library(circular)
library(ggplot2)

num_trials_per_block=104
num_blocks=4
num_trials=num_trials_per_block*num_blocks
num_subjs=500
num_segments=8
if (num_segments==4) {
  seg_probs=c(0.36,0.43,0.56,0.65)
  # seg_probs=c(0,0.33,0.67,1.0)
} else if (num_segments<4) {
  seg_probs=c(0.2,0.8,0.56)[1:num_segments]
} else if (num_segments==8) {
  seg_probs=c(0.34,0.38,0.42,0.46,0.51,0.56,0.64,0.69)
  seg_means=c(0,45,90,135,180,225,270,315)
}

seg_dist_mat8=cbind(c(0,1,2,3,4,3,2,1),c(1,0,1,2,3,4,3,2),c(2,1,0,1,2,3,4,3),
                    c(3,2,1,0,1,2,3,4),c(4,3,2,1,0,1,2,3),c(3,4,3,2,1,0,1,2),
                    c(2,3,4,3,2,1,0,1),c(1,2,3,4,3,2,1,0))

# parameter values
alpha=0.5
beta=3
bf_sd=c(15,22.5,45,70)
num_bf=c(4,5,6,7,8)
rand_bf=0
value_bf=1
bf_start=0

all_cap=c(3,4,5,6,7,8)

out_compare=array(data=NA,dim=c(length(bf_sd)*length(num_bf)*length(all_cap),9))
count=1

for (c in 1:length(all_cap)) {
  for (sd in 1:length(bf_sd)) {
    for (n in 1:length(num_bf)) {
      
      print(paste('running models with capacity of',all_cap[c],',',num_bf[n],
                  'basis functions and sd of',bf_sd[sd],sep=' '))
      
      out_compare[count,1:6]=c(all_cap[c],num_bf[n],bf_sd[sd],num_subjs,num_trials,num_blocks)
      
      #specify fixed sampling - used to set basis function locations
      fixed_samp=sample(1:8,8,replace=F,prob=(rep(1/8,8)))
      
      #define locations of basis function means
      bf_means=array(NA,num_bf[n])
      bf_dist=360/num_bf[n]
      if (rand_bf==1) {
        bf_start=sample(seq(1,360,by=1),1)
        for (i in 1:num_bf[n]) bf_means[i]=bf_start+(i-1)*bf_dist
      }
      else if (value_bf==1) { #set means to first set of samples
        bf_segs=fixed_samp[1:num_bf[n]]
        
      } else {
        for (i in 1:num_bf[n]) bf_means[i]=bf_start+(i-1)*bf_dist
      }
      
      #define AUC of each basis function per segment - convert all to radians
      bf_sd_circ=exp(-.5*((bf_sd[sd]*pi)/180)^2) 
      auc_seg=array(NA,dim=c(num_segments,num_bf[n]))
      for (s in 1:num_segments) {
        seg_min=suppressWarnings(as.circular(seg_means[s]*pi/180-pi/num_segments,units='radians'))
        seg_max=suppressWarnings(as.circular(seg_means[s]*pi/180+pi/num_segments,units='radians'))
        bf_means_circ=suppressWarnings(as.circular(bf_means*pi/180,units='radians'))
        for (b in 1:num_bf[n]) { #something weird with probabilities here
          # p_min=ifelse(pwrappednormal(seg_min,mu=bf_means_circ[b],rho=bf_sd_circ)>0.99999,
          #          0,pwrappednormal(seg_min,mu=bf_means_circ[b],rho=bf_sd_circ))
          # p_max=ifelse(pwrappednormal(seg_max,mu=bf_means_circ[b],rho=bf_sd_circ)>0.99999,
          #              0,pwrappednormal(seg_max,mu=bf_means_circ[b],rho=bf_sd_circ))
          # auc_seg[s,b]=abs(p_max-p_min)
          auc_seg[s,b]=abs(pwrappednormal(seg_max,mu=bf_means_circ[b],rho=bf_sd_circ)-
            pwrappednormal(seg_min,mu=bf_means_circ[b],rho=bf_sd_circ))
          if (auc_seg[s,b]>0.999999) auc_seg[s,b]=0 #auc_seg[s,b]=auc_seg[s,b]+1
        }
      }
      
      # set up outputs
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
      bf_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_bf[n]))
      
      out_info=matrix(data=NA,nrow=(num_trials*num_subjs),ncol=20)
      
      # simulate data 
      for (s in 1:num_subjs) {
        all_values[s,1,]=0.5
        bf_values[s,1,]=0.5
        seg_samplehx[s,1,]=0
        
        block_tnum=1
        block=1
        # if (num_segments<4) seg_probs=block_seg_probs[block,]
        
        for (t in 1:num_trials) {
          if (block_tnum>num_trials_per_block) {
            block_tnum=1
            block=block+1
            # if (num_segments<4) seg_probs=block_seg_probs[block,]
            all_values[s,t,]=0.5
            bf_values[s,t,]=0.5
            seg_samplehx[s,t,]=0
          } else block_tnum=block_tnum+1
          block_num[s,t]=block
          
          #choice
          if (block_tnum<(num_segments+1)) { #forced sampling: assume even for now
            chosen_seg=block_tnum
          } else {
            chosen_seg=softmax_sim(in_values=all_values[s,t,],beta=beta)
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
          
          #filter reward through basis functions and update values of BFs 
          # w/capacity limitations
          rew_bf=rew*auc_seg[chosen_seg,]
          for (i in 1:num_bf[n]) {
            bf_values[s,t+1,i]=min(all_cap[c]/num_bf[n],1)*(bf_values[s,t,i]+
              alpha*(rew_bf[i]-bf_values[s,t,i]))
          }
          
          #filter updated BF values back into segment values for choices
          for (j in 1:num_segments) {
            all_values[s,t+1,j]=sum(bf_values[s,t+1,]*auc_seg[j,])
          }
          
          out_info[(s-1)*num_trials+t,]=c(s,t,block,block_tnum,chosen_seg,
            seg_dist[s,t],seg_samplehx_selected[s,t],seg_probs[chosen_seg],rew,
            seg_dist_bin[s,t],max_val[s,t],pers_exp_bins[s,t],stay_max[s,t],
            stay_nomax[s,t],switch_max[s,t],switch_nomax[s,t],stay_postmax[s,t],
            stay_postnomax[s,t],switch_postmax[s,t],switch_postnomax[s,t])
        }
      }
      
      out_info_df=data.frame(out_info)
      names(out_info_df)=c('subj','trial','block','trial_in_block','chosen_seg',
        'seg_dist','samplehx_selected','selected_prob','rew','binary_seg_dist',
        'max_value_chosen','pers_exp_bins','stay_max','stay_nomax','switch_max',
        'switch_nomax','stay_postmax','stay_postnomax','switch_postmax','switch_postnomax')
      out_info_fc=out_info_df[out_info_df$trial_in_block>num_segments,]
      
      out_compare[count,7]=sum(out_info_fc$rew)/(num_trials*num_subjs)
      out_compare[count,8]=sum(out_info_fc$selected_prob>.5)/(num_trials*num_subjs)
      out_compare[count,9]=sum(out_info_fc$max_value_chosen)/(num_trials*num_subjs)
      count=count+1
  }
  }
}

out_compare_df=data.frame(out_compare)
names(out_compare_df)=c('capacity','basis_functions','function_SD','subjects','trials','blocks',
                        'reward_per_trial','proportion_trials_prob_over_50',
                        'proportion_trials_max_prob')
out_compare_df$basis_functionsf=as.factor(out_compare_df$basis_functions)
out_compare_df$function_SDf=as.factor(out_compare_df$function_SD)

ggplot(out_compare_df,aes(capacity,reward_per_trial,color=basis_functionsf,
                          fill=basis_functionsf,group=basis_functionsf))+
  geom_point()+geom_line()+theme_bw()+facet_wrap(~function_SDf,nrow=1)+
  labs(x='capacity',y='reward per trial',color='# of basis functions',
       fill='# of basis functions')

ggplot(out_compare_df,aes(capacity,proportion_trials_prob_over_50,color=basis_functionsf,
                          fill=basis_functionsf,group=basis_functionsf))+
  geom_point()+geom_line()+theme_bw()+facet_wrap(~function_SDf,nrow=1)+
  labs(x='capacity',y='prop. trials with prob. > 0.5',color='# of basis functions',
       fill='# of basis functions')

ggplot(out_compare_df,aes(capacity,proportion_trials_max_prob,color=basis_functionsf,
                          fill=basis_functionsf,group=basis_functionsf))+
  geom_point()+geom_line()+theme_bw()+facet_wrap(~function_SDf,nrow=1)+
  labs(x='capacity',y='prop. trials with max. prob.',color='# of basis functions',
       fill='# of basis functions')

# out_compare_df_no15=out_compare_df[out_compare_df$function_SD>20,]
# ggplot(out_compare_df_no15,aes(capacity,reward_per_trial,color=basis_functionsf,
#                           fill=basis_functionsf,group=basis_functionsf))+
#   geom_point()+geom_line()+theme_bw()+facet_wrap(~function_SDf,nrow=1)+
#   labs(x='capacity',y='reward per trial',color='# of basis functions',
#        fill='# of basis functions')
# 
# ggplot(out_compare_df_no15,aes(capacity,proportion_trials_prob_over_50,color=basis_functionsf,
#                           fill=basis_functionsf,group=basis_functionsf))+
#   geom_point()+geom_line()+theme_bw()+facet_wrap(~function_SDf,nrow=1)+
#   labs(x='capacity',y='prop. trials with prob. > 0.5',color='# of basis functions',
#        fill='# of basis functions')
# 
# ggplot(out_compare_df_no15,aes(capacity,proportion_trials_max_prob,color=basis_functionsf,
#                           fill=basis_functionsf,group=basis_functionsf))+
#   geom_point()+geom_line()+theme_bw()+facet_wrap(~function_SDf,nrow=1)+
#   labs(x='capacity',y='prop. trials with max. prob.',color='# of basis functions',
#        fill='# of basis functions')
