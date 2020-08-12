# simulate data, estimate parameters from simulated data, and assess parameter recovery

# set up ####
source('sim_functions.R')
library(rstan)
library(shinystan)
library(loo)
options(mc.cores=2)
rstan_options(auto_write=TRUE)

num_loops=4
num_trials_per_block=34
num_blocks=4
num_trials=num_trials_per_block*num_blocks
num_subjs=50
num_segments=4
if (num_segments==4) {
  seg_probs=c(0.36,0.43,0.56,0.65)
  # seg_probs=c(0,0.33,0.67,1.0)
} else if (num_segments<4) {
  seg_probs=c(0.2,0.8,0.56)[1:num_segments]
} else if (num_segments==8) {
  seg_probs=c(0.34,0.38,0.42,0.46,0.51,0.56,0.64,0.69)
}

#distributions to sample parameter values
#learning rate
alpha_alpha=c(6,3,3)
alpha_beta=c(3,3,6)
# alpha_mean=c(log(1/3),log(1),log(3))
# alpha_sigma=0.33
#decay
decay_mean=c(0.9,0.95,1)
decay_sd=0.05
#kappa (sutton bonus)
kappa_mean=c(1,.5,0)
kappa_sd=0.25
#epsilon (e-greedy)
epsilon_alpha=1
epsilon_beta=3
#beta (inverse temperature, gamma dist)
beta_alpha=c(8,4,16) #c(3,6,12) #c(2,200,10)
beta_beta=2
beta_mean=c(2,4,8)
beta_sigma=1
#tau (perseveration)
tau_mean=c(-0.5,0.5,0)
tau_sd=0.25
#omega (UCB variance effect)
omega_mean=c(5,-5,0) #c(20,-20,0)
omega_sd=1
#lapse
epsilon_alpha=c(0.5,3,6)
epsilon_beta=6

seg_dist_mat=cbind(c(0,1,2,1),c(1,0,1,2),c(2,1,0,1),c(1,2,1,0))

num_segments_in=array(data=num_segments,dim=c(num_subjs,num_trials))
points_shown=array(data=0,dim=c(num_subjs,num_trials))
# block_num=array(data=1,dim=c(num_subjs,num_trials))

# model: RL  + softmax w/decay ####

for (a_levels in 2) { #1:3) {
  for (b_levels in 2) { #1:3) {
    for (l_levels in 2) { #1:3) {
      # if (a_levels>1||b_levels>1||e_levels>1) {

# for (x in 1:num_loops) {
# 
#   a_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
#   b_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
# 
#   l_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
#   p_levels=sample(c(1,2,3),1,prob=rep(1/3,3))

    #simulate data
    alpha_params=rbeta(num_subjs,alpha_alpha[a_levels],alpha_beta[a_levels])
    beta_params=rgamma(num_subjs,beta_alpha[b_levels],beta_beta)
    # eps_params=rbeta(num_subjs,epsilon_alpha[e_levels],epsilon_beta)
    # pers_params=rnorm(num_subjs,tau_mean[p_levels],tau_sd)
    lambda_params=rnorm(num_subjs,decay_mean[l_levels],decay_sd)
    # kappa_params=rnorm(num_subjs,kappa_mean[k_levels],kappa_sd)
    print(paste('running alpha value',alpha_alpha[a_levels]/(alpha_alpha[a_levels]+
      alpha_beta[a_levels]),', beta value',beta_alpha[b_levels]/beta_beta,
      # ', and lapse value',round(epsilon_alpha[e_levels]/(epsilon_alpha[e_levels]+epsilon_beta),3),
      # ', and perseveration',tau_mean[p_levels],
      ', and decay value',decay_mean[l_levels], #', and kappa value',kappa_mean[k_levels],
      '.',sep=" "))

    all_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
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
    block_num=array(data=NA,dim=c(num_subjs,num_trials))

    for (s in 1:num_subjs) {
      all_values[s,1,]=0.5
      seg_samplehx[s,1,]=0

      block_tnum=0
      block=1
      # if (num_segments<4) seg_probs=block_seg_probs[block,]

      for (t in 1:num_trials) {
        if (block_tnum>(num_trials_per_block-1)) {
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
          chosen_seg=softmax_sim(in_values=all_values[s,t,],last_chosen=
            all_choices[s,t-1],beta=beta_params[s]) #,epsilon=eps_params[s],tau=pers_params[s])
        }

        #update variables
        all_choices[s,t]=chosen_seg
        seg_samplehx[s,t+1,chosen_seg]=seg_samplehx[s,t,chosen_seg]+1
        seg_samplehx[s,t+1,-chosen_seg]=seg_samplehx[s,t,-chosen_seg]
        seg_dist[s,t]=ifelse(t==1,NA,seg_dist_mat[chosen_seg,all_choices[s,t-1]])
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
        all_values[s,t+1,]=RL_decay_sim(in_values=all_values[s,t,],rew=rew,
          chosen_seg=chosen_seg,alpha=alpha_params[s],decay=lambda_params[s])
      }
    }

    #estimate data
    data_toest=list(nS=num_subjs,nT=num_trials,num_segments=num_segments_in,
      points_shown=points_shown,choice=all_choices,reward=all_rews,block_num=block_num)
    model_out=stan(file='estimate_RL_decay_softmax.stan',data=data_toest,verbose=FALSE,
                   save_warmup=FALSE,iter=4000,chains=3,control=list(adapt_delta=0.9))
    # model_out=stan(file='estimate_RL_lapse_softmax.stan',data=data_toest,verbose=FALSE,
    #                save_warmup=FALSE,iter=1000,chains=1)

    pairs(model_out,pars=c('alpha_m','alpha_s','beta_m','beta_s',
                           # 'eps_m','eps_s'))
                           # 'pers_m','pers_s',
                          'lambda_m','lambda_s'))
    print(summary(model_out,pars=c('alpha_m','alpha_s','beta_m','beta_s',
                                   # 'eps_m','eps_s'))$summary)
                          # 'pers_m','pers_s',
                        'lambda_m','lambda_s'))$summary)
                          

    alpha_est_pre=summary(model_out,pars='alpha',probs=0.5)$summary
    alpha_est=alpha_est_pre[,4]
    beta_est_pre=summary(model_out,pars='beta',probs=0.5)$summary
    beta_est=beta_est_pre[,4]
    lambda_est_pre=summary(model_out,pars='lambda',probs=0.5)$summary
    lambda_est=lambda_est_pre[,4]
    # pers_est_pre=summary(model_out,pars='pers',probs=0.5)$summary
    # pers_est=pers_est_pre[,4]
    # epsilon_est_pre=summary(model_out,pars='eps',probs=0.5)$summary
    # epsilon_est=epsilon_est_pre[,4]

    par(mfrow=c(1,1))
    plot(alpha_params,alpha_est,pch=19,col='blue',ylim=c(0,1),xlim=c(0,1),
         main=paste('mean simulated alpha of',
           round(mean(alpha_params),3),', mean recovered alpha of',round(mean(alpha_est),3),'.',sep=" "))
    abline(a=0,b=1)
    plot(beta_params,beta_est,pch=19,col='maroon',ylim=c(0,15),xlim=c(0,15),
         main=paste('mean simulated beta of',
           round(mean(beta_params),3),', mean recovered beta of',round(mean(beta_est),3),'.',sep=" "))
    abline(a=0,b=1)
    plot(lambda_params,lambda_est,pch=19,col='forestgreen',ylim=c(0.8,1.2),xlim=c(0.8,1.2),
         main=paste('mean simulated lambda of',
           round(mean(lambda_params),3),', mean recovered lambda of',round(mean(lambda_est),3),'.',sep=" "))
    abline(a=0,b=1)
    # plot(pers_params,pers_est,pch=19,col='salmon',ylim=c(-1.5,1.5),xlim=c(-1.5,1.5),
    #      main=paste('mean simulated perseveration of',
    #        round(mean(pers_params),3),', mean recovered of',round(mean(pers_est),3),'.',sep=" "))
    # abline(a=0,b=1)
    # plot(eps_params,epsilon_est,pch=19,col='aquamarine',ylim=c(0,1),xlim=c(0,1),
    #      main=paste('mean simulated epsilon of',
    #        round(mean(eps_params),3),', mean recovered epsilon of',round(mean(epsilon_est),3),'.',sep=" "))
    # abline(a=0,b=1)

    par(mfrow=c(3,1))
    plot(density(alpha_params),col='blue',main='density: LR',xlim=c(0,1))
    par(new=TRUE)
    plot(density(alpha_est),col='blue',lty=2,main="",xlim=c(0,1))
    par(new=FALSE)

    plot(density(beta_params),col='maroon',main='density: IT',xlim=c(0,15))
    par(new=TRUE)
    plot(density(beta_est),col='maroon',lty=2,main="",xlim=c(0,15))
    par(new=FALSE)

    plot(density(lambda_params),col='forestgreen',main='density: decay',xlim=c(0.8,1.2))
    par(new=TRUE)
    plot(density(lambda_est),col='forestgreen',lty=2,main="",xlim=c(0.8,1.2))
    par(new=FALSE)
    # 
    # plot(density(pers_params),col='salmon',main='density: perseveration',xlim=c(-1.5,1.5))
    # par(new=TRUE)
    # plot(density(pers_est),col='salmon',lty=2,main="",xlim=c(-1.5,1.5))
    # par(new=FALSE)
    
    # plot(density(eps_params),col='aquamarine',main='density: lapse',xlim=c(0,1))
    # par(new=TRUE)
    # plot(density(epsilon_est),col='aquamarine',lty=2,main="",xlim=c(0,1))
    # par(new=FALSE)
      # }

    }
  }
}

# # model: forgetting Bayes + UCB + softmax ####
# 
# for (a_levels in 1:3) {
#   for (b_levels in 1:3) {
# # # l_levels=2
#     # for (l_levels in 2:3) {
#     
# # for (x in 1:num_loops) {
# # #   
# #   a_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
#   # b_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
#   # l_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
# #   k_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
# 
#     #simulate data
#     beta_params=rgamma(num_subjs,beta_alpha[b_levels],beta_beta)
#     # beta_params=exp(rnorm(num_subjs,log(beta_mean[b_levels]),
#     #                       log(beta_mean[b_levels]+beta_sigma)-log(beta_mean[b_levels])))
#     # lambda_params=rnorm(num_subjs,decay_mean[l_levels],decay_sd)
#     omega_params=rnorm(num_subjs,omega_mean[a_levels],omega_sd)
#     print(paste('running beta value',beta_alpha[b_levels]/beta_beta,
#       # ', decay value',decay_mean[l_levels],'.',sep=" ")) #,
#       ',and UCB omega value',omega_mean[a_levels],'.',sep=" "))
#     
#     all_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
#     all_alphas=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
#     all_betas=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
#     all_choices=array(data=NA,dim=c(num_subjs,num_trials))
#     all_rews=array(data=NA,dim=c(num_subjs,num_trials))
#     seg_samplehx=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
#     seg_dist=array(data=NA,dim=c(num_subjs,num_trials))
#     seg_dist_bin=array(data=NA,dim=c(num_subjs,num_trials))
#     max_val=array(data=NA,dim=c(num_subjs,num_trials))
#     stay_max=array(data=0,dim=c(num_subjs,num_trials))
#     stay_nomax=array(data=0,dim=c(num_subjs,num_trials))
#     switch_max=array(data=0,dim=c(num_subjs,num_trials))
#     switch_nomax=array(data=0,dim=c(num_subjs,num_trials))
#     block_num=array(data=NA,dim=c(num_subjs,num_trials))
#     
#     for (s in 1:num_subjs) {
#       all_alphas[s,1,]=1
#       all_betas[s,1,]=1
#       seg_samplehx[s,1,]=0
#       
#       block_tnum=0
#       block=1
#       # if (num_segments<4) seg_probs=block_seg_probs[block,]
#       
#       for (t in 1:num_trials) {
#         if (block_tnum>(num_trials_per_block-1)) {
#           block_tnum=1
#           block=block+1
#           # if (num_segments<4) seg_probs=block_seg_probs[block,]
#           all_values[s,t,]=0.5
#           seg_samplehx[s,t,]=0
#         } else block_tnum=block_tnum+1
#         block_num[s,t]=block
# 
#         #choice
#         if (block_tnum<(num_segments+1)) { #forced sampling: assume even for now
#           chosen_seg=block_tnum
#         } else {
#           all_values[s,t,]=all_alphas[s,t,]/(all_alphas[s,t,]+all_betas[s,t,])
#           chosen_seg=UCB_softmax_sim(in_values=all_values[s,t,],in_alphas=all_alphas[s,t,],
#             in_betas=all_betas[s,t,],beta=beta_params[s],omega=omega_params[s])
#         }
#         
#         #update variables
#         all_choices[s,t]=chosen_seg
#         seg_samplehx[s,t+1,chosen_seg]=seg_samplehx[s,t,chosen_seg]+1
#         seg_samplehx[s,t+1,-chosen_seg]=seg_samplehx[s,t,-chosen_seg]
#         seg_dist[s,t]=ifelse(t==1,NA,seg_dist_mat[chosen_seg,all_choices[s,t-1]])
#         seg_dist_bin[s,t]=ifelse(seg_dist[s,t]==0,0,1)
#         max_val[s,t]=ifelse(chosen_seg %in% which(all_values[s,t,]==max(all_values[s,t,])),1,0)
#         
#         if (is.na(seg_dist_bin[s,t])) {
#         } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==0) {
#           stay_max[s,t]=1
#         } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==0) {
#           stay_nomax[s,t]=1
#         } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==1) {
#           switch_max[s,t]=1
#         } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==1) {
#           switch_nomax[s,t]=1
#         } 
#         
#         #outcome
#         rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
#         all_rews[s,t]=rew
#         
#         #update values
#         # out_args=forgetting_bayes_sim(in_values=all_values[s,t,],
#         #     in_alpha=all_alphas[s,t,],in_beta=all_betas[s,t,],rew=rew,
#         #   chosen_seg=chosen_seg,decay=lambda_params[s])
#         out_args=ideal_bayes_sim(in_values=all_values[s,t,],
#             in_alpha=all_alphas[s,t,],in_beta=all_betas[s,t,],rew=rew,
#           chosen_seg=chosen_seg)
#         all_alphas[s,t+1,]=out_args[1:num_segments]
#         all_betas[s,t+1,]=out_args[(num_segments+1):(2*num_segments)]
#       }
#     }
#     
#     #estimate data
#     data_toest=list(nS=num_subjs,nT=num_trials,num_segments=num_segments_in,
#       points_shown=points_shown,choice=all_choices,reward=all_rews,block_num=block_num) #,rand_theta=rand_theta)
#     model_out=stan(file='estimate_Bayes_UCB_softmax.stan',data=data_toest,verbose=FALSE,
#                    save_warmup=FALSE,iter=4000,chains=3) #,control=list(adapt_delta=0.90))
#     # model_out=stan(file='estimate_Bayes_decay_UCB_softmax.stan',data=data_toest,verbose=FALSE,
#     #                save_warmup=FALSE,iter=1000,chains=2) #control=list(adapt_delta=0.99))
#     
#     # pairs(model_out,pars=c('beta_m','beta_s','lambda_m','lambda_s','omega_m','omega_s'))
#     # traceplot(model_out,pars=c('beta_m','beta_s','lambda_m','lambda_s','omega_m','omega_s'))
#     # summary(model_out,pars=c('beta_m','beta_s','lambda_m','lambda_s','omega_m','omega_s'))$summary
#     
#     pairs(model_out,pars=c('beta_m','beta_s','omega_m','omega_s'))
#     traceplot(model_out,pars=c('beta_m','beta_s','omega_m','omega_s'))
#     summary(model_out,pars=c('beta_m','beta_s','omega_m','omega_s'))$summary
# 
# 
#     beta_est_pre=summary(model_out,pars='beta',probs=0.5)$summary
#     beta_est=beta_est_pre[,4]
#     # lambda_est_pre=summary(model_out,pars='lambda',probs=0.5)$summary
#     # lambda_est=lambda_est_pre[,4]
#     omega_est_pre=summary(model_out,pars='omega',probs=0.5)$summary
#     omega_est=omega_est_pre[,4]
#     
#     par(mfrow=c(1,1))
#     plot(beta_params,beta_est,pch=19,col='maroon',ylim=c(0,15),xlim=c(0,15),
#          main=paste('mean sim. beta of',round(mean(beta_params),3),
#                     'and rec. beta of',round(mean(beta_est),3),'.',sep=" "))
#     abline(a=0,b=1)
#     # plot(lambda_params,lambda_est,pch=19,col='forestgreen',ylim=c(0.6,1.2),xlim=c(0.6,1.2),
#     #      main=paste('mean sim. lambda of',round(mean(lambda_params),3),
#     #       'and rec. lambda of',round(mean(lambda_est),3),'.',sep=" "))
#     # abline(a=0,b=1)
#     plot(omega_params,omega_est,pch=19,col='purple',ylim=c(-30,30),xlim=c(-30,30),
#          main=paste('mean sim. omega of',round(mean(omega_params),3),
#           'and rec. omega of',round(mean(omega_est),3),'.',sep=" "))
#     abline(a=0,b=1)
#     
#     par(mfrow=c(2,1))
#     plot(density(beta_params),col='maroon',main='density: IT',xlim=c(0,15))
#     par(new=TRUE)
#     plot(density(beta_est),col='maroon',lty=2,main="",xlim=c(0,15))
#     par(new=FALSE)
#     
#     # plot(density(lambda_params),col='forestgreen',main='density: decay',xlim=c(0.6,1.2))
#     # par(new=TRUE)
#     # plot(density(lambda_est),col='forestgreen',lty=2,main="",xlim=c(0.6,1.2))
#     # par(new=FALSE)
#     
#     plot(density(omega_params),col='purple',main='density: UCB omega',xlim=c(-30,30))
#     par(new=TRUE)
#     plot(density(omega_est),col='purple',lty=2,main="",xlim=c(-30,30))
#     par(new=FALSE)
#     
#     
#   #   }
#   }
# }

# # model: RL w/decay + softmax w/exploration bonus ####
# 
# # for (a_levels in 1:3) {
# #   for (b_levels in 1:3) {
# #     for (l_levels in 1:3) {
# 
# for (x in 1:num_loops) {
# 
#   a_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
#   b_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
#   l_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
#   k_levels=sample(c(1,2,3),1,prob=rep(1/3,3))
# 
#     #simulate data
#     alpha_params=rbeta(num_subjs,alpha_alpha[a_levels],alpha_beta[a_levels])
#     beta_params=rgamma(num_subjs,beta_alpha[b_levels],beta_beta)
#     lambda_params=rnorm(num_subjs,decay_mean[l_levels],decay_sd)
#     kappa_params=rnorm(num_subjs,kappa_mean[k_levels],kappa_sd)
#     print(paste('running alpha value',alpha_alpha[a_levels]/(alpha_alpha[a_levels]+
#       alpha_beta[a_levels]),', beta value',beta_alpha[b_levels]/beta_beta,
#       ', decay value',decay_mean[l_levels],', and kappa value',kappa_mean[k_levels],
#       '.',sep=" "))
# 
#     all_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
#     all_choices=array(data=NA,dim=c(num_subjs,num_trials))
#     all_rews=array(data=NA,dim=c(num_subjs,num_trials))
#     seg_samplehx=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
#     seg_dist=array(data=NA,dim=c(num_subjs,num_trials))
#     seg_dist_bin=array(data=NA,dim=c(num_subjs,num_trials))
#     max_val=array(data=NA,dim=c(num_subjs,num_trials))
#     stay_max=array(data=0,dim=c(num_subjs,num_trials))
#     stay_nomax=array(data=0,dim=c(num_subjs,num_trials))
#     switch_max=array(data=0,dim=c(num_subjs,num_trials))
#     switch_nomax=array(data=0,dim=c(num_subjs,num_trials))
#     block_num=array(data=NA,dim=c(num_subjs,num_trials))
# 
#     for (s in 1:num_subjs) {
#       all_values[s,1,]=0.5
#       seg_samplehx[s,1,]=0
# 
#       block_tnum=0
#       block=1
#       # if (num_segments<4) seg_probs=block_seg_probs[block,]
# 
#       for (t in 1:num_trials) {
#         if (block_tnum>(num_trials_per_block-1)) {
#           block_tnum=1
#           block=block+1
#           # if (num_segments<4) seg_probs=block_seg_probs[block,]
#           all_values[s,t,]=0.5
#           seg_samplehx[s,t,]=0
#         } else block_tnum=block_tnum+1
#         block_num[s,t]=block
# 
#         #choice
#         if (block_tnum<(num_segments+1)) { #forced sampling: assume even for now
#           chosen_seg=block_tnum
#         } else {
#           chosen_seg=softmax_SBalt_sim(in_values=all_values[s,t,],seg_samplehx=
#             seg_samplehx[s,t,],beta=beta_params[s],kappa=kappa_params[s])
#         }
# 
#         #update variables
#         all_choices[s,t]=chosen_seg
#         seg_samplehx[s,t+1,chosen_seg]=seg_samplehx[s,t,chosen_seg]+1
#         seg_samplehx[s,t+1,-chosen_seg]=seg_samplehx[s,t,-chosen_seg]
#         seg_dist[s,t]=ifelse(t==1,NA,seg_dist_mat[chosen_seg,all_choices[s,t-1]])
#         seg_dist_bin[s,t]=ifelse(seg_dist[s,t]==0,0,1)
#         max_val[s,t]=ifelse(chosen_seg %in% which(all_values[s,t,]==max(all_values[s,t,])),1,0)
# 
#         if (is.na(seg_dist_bin[s,t])) {
#         } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==0) {
#           stay_max[s,t]=1
#         } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==0) {
#           stay_nomax[s,t]=1
#         } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==1) {
#           switch_max[s,t]=1
#         } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==1) {
#           switch_nomax[s,t]=1
#         }
# 
#         #outcome
#         rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
#         all_rews[s,t]=rew
# 
#         #update values
#         all_values[s,t+1,]=RL_decay_sim(in_values=all_values[s,t,],rew=rew,
#           chosen_seg=chosen_seg,alpha=alpha_params[s],decay=lambda_params[s])
#       }
#     }
# 
#     #estimate data
#     data_toest=list(nS=num_subjs,nT=num_trials,num_segments=num_segments_in,
#       points_shown=points_shown,choice=all_choices,reward=all_rews,block_num=block_num)
#     model_out=stan(file='estimate_RL_decay_softmax_expbonus.stan',data=data_toest,verbose=FALSE,
#                    save_warmup=FALSE,iter=4000,control=list(adapt_delta=0.9))
#     # model_out=stan(file='estimate_RL_decay_softmax_expbonus.stan',data=data_toest,verbose=FALSE,
#     #                save_warmup=FALSE,iter=2000,chains=2)
# 
#     pairs(model_out,pars=c('alpha_m','alpha_s','beta_m','beta_s','lambda_m',
#                            'lambda_s','kappa_m','kappa_s'))
#     summary(model_out,pars=c('alpha_m','alpha_s','beta_m','beta_s','lambda_m',
#                              'lambda_s','kappa_m','kappa_s'))
# 
#     alpha_est_pre=summary(model_out,pars='alpha',probs=0.5)$summary
#     alpha_est=alpha_est_pre[,4]
#     beta_est_pre=summary(model_out,pars='beta',probs=0.5)$summary
#     beta_est=beta_est_pre[,4]
#     lambda_est_pre=summary(model_out,pars='lambda',probs=0.5)$summary
#     lambda_est=lambda_est_pre[,4]
#     kappa_est_pre=summary(model_out,pars='kappa',probs=0.5)$summary
#     kappa_est=kappa_est_pre[,4]
# 
#     par(mfrow=c(1,1))
#     plot(alpha_params,alpha_est,pch=19,col='blue',ylim=c(0,1),xlim=c(0,1),
#          main=paste('Simulated vs. recovered learning rates: mean simulated alpha of',
#            round(mean(alpha_params),3),'and mean recovered alpha of',round(mean(alpha_est),3),'.',sep=" "))
#     abline(a=0,b=1)
#     plot(beta_params,beta_est,pch=19,col='maroon',ylim=c(0,15),xlim=c(0,15),
#          main=paste('Simulated vs. recovered inverse temperatures: mean simulated beta of',
#            round(mean(beta_params),3),'and mean recovered beta of',round(mean(beta_est),3),'.',sep=" "))
#     abline(a=0,b=1)
#     plot(lambda_params,lambda_est,pch=19,col='forestgreen',ylim=c(0.5,1.5),xlim=c(0.5,1.5),
#          main=paste('Simulated vs. recovered decay: mean simulated lambda of',
#            round(mean(lambda_params),3),'and mean recovered lambda of',round(mean(lambda_est),3),'.',sep=" "))
#     abline(a=0,b=1)
#     plot(kappa_params,kappa_est,pch=19,col='goldenrod4',ylim=c(-0.5,1.5),xlim=c(-0.5,1.5),
#          main=paste('Simulated vs. recovered decay: mean simulated kappa of',
#            round(mean(kappa_params),3),'and mean recovered kappa of',round(mean(kappa_est),3),'.',sep=" "))
#     abline(a=0,b=1)
# 
#     par(mfrow=c(2,2))
#     plot(density(alpha_params),col='blue',main='density: LR',xlim=c(0,1))
#     par(new=TRUE)
#     plot(density(alpha_est),col='blue',lty=2,main="",xlim=c(0,1))
#     par(new=FALSE)
# 
#     plot(density(beta_params),col='maroon',main='density: IT',xlim=c(0,15))
#     par(new=TRUE)
#     plot(density(beta_est),col='maroon',lty=2,main="",xlim=c(0,15))
#     par(new=FALSE)
# 
#     plot(density(lambda_params),col='forestgreen',main='density: decay',xlim=c(0.5,1.5))
#     par(new=TRUE)
#     plot(density(lambda_est),col='forestgreen',lty=2,main="",xlim=c(0.5,1.5))
#     par(new=FALSE)
# 
#     plot(density(kappa_params),col='goldenrod4',main='density: exp. bonus',xlim=c(-0.5,1.5))
#     par(new=TRUE)
#     plot(density(kappa_est),col='goldenrod4',lty=2,main="",xlim=c(-0.5,1.5))
#     par(new=FALSE)
# 
#   #   }
#   # }
# }

# # model: basic RL + softmax ####
# 
# for (a_levels in 1:3) {
#   for (b_levels in 1:3) {
#     
#     #simulate data
#     alpha_params=rbeta(num_subjs,alpha_alpha[a_levels],alpha_beta[a_levels])
#     # alpha_params_pre=rnorm(num_subjs,alpha_mean[a_levels],alpha_sigma[a_levels])
#     # alpha_params=1/(1+exp(-1*alpha_params_pre))
#     beta_params=rgamma(num_subjs,beta_alpha[b_levels],beta_beta)
#     # beta_params=rnorm(num_subjs,beta_alpha[b_levels],beta_beta)
#     # for (b in 1:length(beta_params)) {
#     #   if (beta_params[b]<0) beta_params[b]=-1*beta_params[b]
#     # }
#     print(paste('running alpha value',alpha_alpha[a_levels]/(alpha_alpha[a_levels]+
#       alpha_beta[a_levels]),'and beta value',beta_alpha[b_levels]/beta_beta,'.',sep=" "))
#     
#     all_values=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
#     all_choices=array(data=NA,dim=c(num_subjs,num_trials))
#     all_rews=array(data=NA,dim=c(num_subjs,num_trials))
#     seg_samplehx=array(data=NA,dim=c(num_subjs,num_trials+1,num_segments))
#     seg_dist=array(data=NA,dim=c(num_subjs,num_trials))
#     seg_dist_bin=array(data=NA,dim=c(num_subjs,num_trials))
#     max_val=array(data=NA,dim=c(num_subjs,num_trials))
#     stay_max=array(data=0,dim=c(num_subjs,num_trials))
#     stay_nomax=array(data=0,dim=c(num_subjs,num_trials))
#     switch_max=array(data=0,dim=c(num_subjs,num_trials))
#     switch_nomax=array(data=0,dim=c(num_subjs,num_trials))
#     block_num=array(data=NA,dim=c(num_subjs,num_trials))
#     
#     for (s in 1:num_subjs) {
#       all_values[s,1,]=0.5
#       seg_samplehx[s,1,]=0
#       
#       block_tnum=1
#       block=1
#       # if (num_segments<4) seg_probs=block_seg_probs[block,]
#       
#       for (t in 1:num_trials) {
#         if (block_tnum>num_trials_per_block) {
#           block_tnum=1
#           block=block+1
#           # if (num_segments<4) seg_probs=block_seg_probs[block,]
#           all_values[s,t,]=0.5
#           seg_samplehx[s,t,]=0
#         } else block_tnum=block_tnum+1
#         block_num[s,t]=block
#         
#         #choice
#         if (block_tnum<(num_segments+1)) { #forced sampling: assume even for now
#           chosen_seg=block_tnum
#         } else {
#           chosen_seg=softmax_sim(in_values=all_values[s,t,],beta=beta_params[s])
#         }
#         
#         #update variables
#         all_choices[s,t]=chosen_seg
#         seg_samplehx[s,t+1,chosen_seg]=seg_samplehx[s,t,chosen_seg]+1
#         seg_samplehx[s,t+1,-chosen_seg]=seg_samplehx[s,t,-chosen_seg]
#         seg_dist[s,t]=ifelse(t==1,NA,seg_dist_mat[chosen_seg,all_choices[s,t-1]])
#         seg_dist_bin[s,t]=ifelse(seg_dist[s,t]==0,0,1)
#         max_val[s,t]=ifelse(chosen_seg %in% which(all_values[s,t,]==max(all_values[s,t,])),1,0)
#         
#         if (is.na(seg_dist_bin[s,t])) {
#         } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==0) {
#           stay_max[s,t]=1
#         } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==0) {
#           stay_nomax[s,t]=1
#         } else if (max_val[s,t]==1&&seg_dist_bin[s,t]==1) {
#           switch_max[s,t]=1
#         } else if (max_val[s,t]==0&&seg_dist_bin[s,t]==1) {
#           switch_nomax[s,t]=1
#         } 
#         
#         #outcome
#         rew=sample(c(0,1),1,prob=c(1-seg_probs[chosen_seg],seg_probs[chosen_seg]))
#         all_rews[s,t]=rew
#         
#         #update values
#         all_values[s,t+1,]=basic_RL_sim(in_values=all_values[s,t,],rew=rew,
#           chosen_seg=chosen_seg,alpha=alpha_params[s])
#       }
#     }
#     
#     #estimate data
#     data_toest=list(nS=num_subjs,nT=num_trials,num_segments=num_segments_in,
#       points_shown=points_shown,choice=all_choices,reward=all_rews,block_num=block_num)
#     model_out=stan(file='estimate_RL_softmax.stan',data=data_toest,verbose=FALSE,
#                    save_warmup=FALSE,iter=4000)
#     # model_out=stan(file='estimate_RL_softmax_wrescaling.stan',data=data_toest,verbose=FALSE,
#     #                save_warmup=FALSE,iter=4000)
#     
#     alpha_est_pre=summary(model_out,pars='alpha',probs=0.5)$summary
#     alpha_est=alpha_est_pre[,4]
#     beta_est_pre=summary(model_out,pars='beta',probs=0.5)$summary
#     beta_est=beta_est_pre[,4]
#     
#     par(mfrow=c(1,1))
#     plot(alpha_params,alpha_est,pch=19,col='blue',ylim=c(0,1),xlim=c(0,1),
#          main=paste('Simulated vs. recovered learning rates: mean simulated alpha of',
#            round(mean(alpha_params),3),'and mean recovered alpha of',round(mean(alpha_est),3),'.',sep=" "))
#     abline(a=0,b=1)
#     plot(beta_params,beta_est,pch=19,col='maroon',ylim=c(0,15),xlim=c(0,15),
#          main=paste('Simulated vs. recovered inverse temperatures: mean simulated beta of',
#            round(mean(beta_params),3),'and mean recovered beta of',round(mean(beta_est),3),'.',sep=" "))
#     abline(a=0,b=1)
#     
#     par(mfrow=c(2,1))
#     plot(density(alpha_params),col='blue',main='density: LR',xlim=c(0,1))
#     par(new=TRUE)
#     plot(density(alpha_est),col='blue',lty=2,main="",xlim=c(0,1))
#     par(new=FALSE)
#     
#     plot(density(beta_params),col='maroon',main='density: IT',xlim=c(0,15))
#     par(new=TRUE)
#     plot(density(beta_est),col='maroon',lty=2,main="",xlim=c(0,15))
#     par(new=FALSE)
#     
#     
#   }
# }
