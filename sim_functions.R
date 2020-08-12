# value functions ####
basic_RL_sim=function(in_values,rew,chosen_seg,alpha,...) {
  
  out_values=in_values
  out_values[chosen_seg]=in_values[chosen_seg]+alpha*(rew-in_values[chosen_seg])
  
  return(out_values)
}

RL_decay_sim=function(in_values,rew,chosen_seg,alpha,decay,...) {
  
  out_values=in_values
  out_values[chosen_seg]=in_values[chosen_seg]+alpha*(rew-in_values[chosen_seg])
  out_values[-chosen_seg]=decay*in_values[-chosen_seg]
  
  return(out_values)
}

RL_SB_sim=function(in_values,rew,chosen_seg,seg_samplehx,alpha,kappa,...) {
  
  out_values=in_values
  if (seg_samplehx[chosen_seg]==0) {
    num_chosen=kappa
  } else num_chosen=seg_samplehx[chosen_seg]
  out_values[chosen_seg]=in_values[chosen_seg]+
    alpha*(rew-in_values[chosen_seg])+
    kappa/num_chosen
  
  return(out_values)
}

ideal_bayes_sim=function(in_alpha,in_beta,rew,chosen_seg,...) {
  
  out_alpha=in_alpha
  out_alpha[chosen_seg]=in_alpha[chosen_seg]+rew
  out_beta=in_beta
  out_beta[chosen_seg]=in_beta[chosen_seg]+1-rew
  
  return(c(out_alpha,out_beta))
}

forgetting_bayes_sim=function(in_alpha,in_beta,rew,chosen_seg,decay,...) {
  
  out_alpha=in_alpha
  out_alpha[chosen_seg]=in_alpha[chosen_seg]+rew
  out_alpha[-chosen_seg]=decay*in_alpha[-chosen_seg]
  out_beta=in_beta
  out_beta[chosen_seg]=in_beta[chosen_seg]+1-rew
  out_beta[-chosen_seg]=decay*(in_beta[-chosen_seg])
  
  return(c(out_alpha,out_beta))
}

forgetting_all_bayes_sim=function(in_alpha,in_beta,rew,chosen_seg,decay,...) {
  
  out_alpha=in_alpha
  out_alpha[chosen_seg]=decay*(in_alpha[chosen_seg]+rew)
  out_alpha[-chosen_seg]=decay*in_alpha[-chosen_seg]
  out_beta=in_beta
  out_beta[chosen_seg]=decay*(in_beta[chosen_seg]+1-rew)
  out_beta[-chosen_seg]=decay*(in_beta[-chosen_seg])
  
  return(c(out_alpha,out_beta))
}

#choice functions ####

greedy_sim=function(in_values,...) {
  
  chosen_seg_all=which(in_values==max(in_values))
  #choose randomly among maximum valued segments if multiple
  if (length(chosen_seg_all)>1) {
    chosen_seg=sample(chosen_seg_all,1)
  } else chosen_seg=chosen_seg_all
  
  return(chosen_seg)
}

e_greedy_sim=function(in_values,epsilon,...) {
  
  num_segments=length(in_values)
  chosen_seg_all=which(in_values==max(in_values))
  prob_chosen=array(data=NA,dim=num_segments)
  #probability of non-max segments is epsilon/n 
  prob_chosen[-chosen_seg_all]=epsilon/num_segments
  #prob of max segments: (1-epsilon*(n-x)/n)*(1/x), where x = # of max segs
  prob_chosen[chosen_seg_all]=(1-epsilon*(num_segments-length(chosen_seg_all))/num_segments)*
    (1/length(chosen_seg_all))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

softmax_sim=function(in_values,beta,...) {
  
  prob_chosen=(exp(beta*in_values))/sum(exp(beta*in_values))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

egreedy_softmax_sim=function(in_values,epsilon,beta,...) {
  
  num_segments=length(in_values)
  prob_chosen=(1-epsilon)*(exp(beta*in_values))/sum(exp(beta*in_values)) + epsilon/num_segments
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

softmax_pers_sim=function(in_values,last_chosen,beta,tau,...) {
  
  #last_chosen here is the index of the last chosen segment
  pers=array(data=-1,dim=length(in_values))
  pers[last_chosen]=1
  
  prob_chosen=(exp(beta*in_values+tau*pers))/sum(exp(beta*in_values+tau*pers))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

egreedy_softmax_pers_sim=function(in_values,last_chosen,epsilon,beta,tau,...) {
  
  #last_chosen here is the index of the last chosen segment
  pers=array(data=-1,dim=length(in_values))
  pers[last_chosen]=1
  
  num_segments=length(in_values)
  prob_chosen=(1-epsilon)*(exp(beta*in_values+tau*pers))/sum(exp(beta*in_values+tau*pers)) + 
    epsilon/num_segments
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

softmax_SB_sim=function(in_values,seg_samplehx,beta,kappa,...) {
  
  prob_chosen=(exp(beta*in_values+kappa/seg_samplehx))/
    sum(exp(beta*in_values+kappa/seg_samplehx))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

softmax_SBalt_sim=function(in_values,seg_samplehx,beta,kappa,...) {
  
    prob_chosen=(exp(beta*in_values+kappa/(seg_samplehx^3)))/
    sum(exp(beta*in_values+kappa/(seg_samplehx^3)))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

softmax_pers_SB_sim=function(in_values,last_chosen,seg_samplehx,beta,kappa,tau,...) {
  
  #last_chosen here is the index of the last chosen segment
  pers=array(data=0,dim=length(in_values))
  pers[last_chosen]=1
  
  prob_chosen=(exp(beta*in_values+tau*pers+kappa/seg_samplehx))/
    sum(exp(beta*in_values+tau*pers+kappa/seg_samplehx))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

softmax_pers_SBalt_sim=function(in_values,last_chosen,seg_samplehx,beta,kappa,tau,...) {
  
  #last_chosen here is the index of the last chosen segment
  pers=array(data=0,dim=length(in_values))
  pers[last_chosen]=1
  
  prob_chosen=(exp(beta*in_values+tau*pers+kappa/(seg_samplehx^3)))/
    sum(exp(beta*in_values+tau*pers+kappa/(seg_samplehx^3)))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

softmax_pers_fUCB_sim=function(in_values,last_chosen,seg_samplehx,beta,kappa,tau,...) {
  
  #last_chosen here is the index of the last chosen segment
  pers=array(data=0,dim=length(in_values))
  pers[last_chosen]=1
  
  #UCB bonus (from Auer 2002)
  total_choices=sum(seg_samplehx)
  seg_samplehx[which(seg_samplehx==0)]=1
  UCB=sqrt(2*(log(seg_samplehx)/total_choices))
  
  prob_chosen=(exp(beta*in_values+tau*pers+kappa*UCB))/
    sum(exp(beta*in_values+tau*pers+kappa*UCB))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

thompson_sim=function(in_alphas,in_betas,...) {
  
  num_segments=length(in_alphas)
  sample_vals=array(NA,dim=num_segments)
  for (n in 1:num_segments) {
    sample_vals[n]=rbeta(1,in_alphas[n],in_betas[n])
  }
  
  chosen_seg_all=which(sample_vals==max(sample_vals))
  #choose randomly among maximum valued segments if multiple
  if (length(chosen_seg_all)>1) {
    chosen_seg=sample(chosen_seg_all,1)
  } else chosen_seg=chosen_seg_all
  
  return(chosen_seg)
}

thompson_softmax_sim=function(in_alphas,in_betas,beta,...) {
  
  num_segments=length(in_alphas)
  sample_vals=array(NA,dim=num_segments)
  for (n in 1:num_segments) {
    sample_vals[n]=rbeta(1,in_alphas[n],in_betas[n])
  }
  
  prob_chosen=(exp(beta*sample_vals))/sum(exp(beta*sample_vals))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

thompson_egreedy_sim=function(in_alphas,in_betas,epsilon,...) {
   
  num_segments=length(in_alphas)
  sample_vals=array(NA,dim=num_segments)
  for (n in 1:num_segments) {
    sample_vals[n]=rbeta(1,in_alphas[n],in_betas[n])
  }
  
  chosen_seg_all=which(sample_vals==max(sample_vals))
  prob_chosen=array(data=NA,dim=num_segments)
  #probability of non-max segments is epsilon/n 
  prob_chosen[-chosen_seg_all]=epsilon/num_segments
  #prob of max segments: (1-epsilon*(n-x)/n)*(1/x), where x = # of max segs
  prob_chosen[chosen_seg_all]=(1-epsilon*(num_segments-length(chosen_seg_all))/num_segments)*
    (1/length(chosen_seg_all))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

UCB_sim=function(in_alphas,in_betas,omega,...) {
  
  num_segments=length(in_alphas)
  mean_vals=in_alphas/(in_alphas+in_betas)
  var_vals=(in_alphas*in_betas)/
    ((in_alphas+in_betas)*(in_alphas+in_betas)*(in_alphas+in_betas+1))
  UCB_vals=mean_vals+omega*var_vals
  
  chosen_seg_all=which(UCB_vals==max(UCB_vals))
   #choose randomly among maximum valued segments if multiple
  if (length(chosen_seg_all)>1) {
    chosen_seg=sample(chosen_seg_all,1)
  } else chosen_seg=chosen_seg_all
  
  return(chosen_seg)
}

UCB_egreedy_sim=function(in_alphas,in_betas,omega,epsilon,...) {
  
  num_segments=length(in_alphas)
  mean_vals=in_alphas/(in_alphas+in_betas)
  var_vals=(in_alphas*in_betas)/
    ((in_alphas+in_betas)*(in_alphas+in_betas)*(in_alphas+in_betas+1))
  UCB_vals=mean_vals+omega*var_vals
  
  chosen_seg_all=which(UCB_vals==max(UCB_vals))
  prob_chosen=array(data=NA,dim=num_segments)
  
  #probability of non-max segments is epsilon/n 
  prob_chosen[-chosen_seg_all]=epsilon/num_segments
  #prob of max segments: (1-epsilon*(n-x)/n)*(1/x), where x = # of max segs
  prob_chosen[chosen_seg_all]=(1-epsilon*(num_segments-length(chosen_seg_all))/num_segments)*
    (1/length(chosen_seg_all))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

UCB_softmax_sim=function(in_alphas,in_betas,omega,beta,...) {
  
  num_segments=length(in_alphas)
  mean_vals=in_alphas/(in_alphas+in_betas)
  var_vals=(in_alphas*in_betas)/
    ((in_alphas+in_betas)*(in_alphas+in_betas)*(in_alphas+in_betas+1))
  UCB_vals=mean_vals+omega*var_vals
  
  prob_chosen=(exp(beta*UCB_vals))/sum(exp(beta*UCB_vals))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

thompson_UCB_sim=function(in_alphas,in_betas,omega,...) {
  
  num_segments=length(in_alphas)
    sample_vals=array(NA,dim=num_segments)
  for (n in 1:num_segments) {
    sample_vals[n]=rbeta(1,in_alphas[n],in_betas[n])
  }
  var_vals=(in_alphas*in_betas)/
    ((in_alphas+in_betas)*(in_alphas+in_betas)*(in_alphas+in_betas+1))
  UCB_vals=sample_vals+omega*var_vals
  
  chosen_seg_all=which(UCB_vals==max(UCB_vals))
   #choose randomly among maximum valued segments if multiple
  if (length(chosen_seg_all)>1) {
    chosen_seg=sample(chosen_seg_all,1)
  } else chosen_seg=chosen_seg_all
  
  return(chosen_seg)
}

thompson_UCB_egreedy_sim=function(in_alphas,in_betas,omega,epsilon,...) {
  
  num_segments=length(in_alphas)
    sample_vals=array(NA,dim=num_segments)
  for (n in 1:num_segments) {
    sample_vals[n]=rbeta(1,in_alphas[n],in_betas[n])
  }
  var_vals=(in_alphas*in_betas)/
    ((in_alphas+in_betas)*(in_alphas+in_betas)*(in_alphas+in_betas+1))
  UCB_vals=sample_vals+omega*var_vals
  
  chosen_seg_all=which(UCB_vals==max(UCB_vals))
  prob_chosen=array(data=NA,dim=num_segments)
  #probability of non-max segments is epsilon/n 
  prob_chosen[-chosen_seg_all]=epsilon/num_segments
  #prob of max segments: (1-epsilon*(n-x)/n)*(1/x), where x = # of max segs
  prob_chosen[chosen_seg_all]=(1-epsilon*(num_segments-length(chosen_seg_all))/num_segments)*
    (1/length(chosen_seg_all))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

thompson_UCB_softmax_sim=function(in_alphas,in_betas,omega,beta,...) {
  
  num_segments=length(in_alphas)
    sample_vals=array(NA,dim=num_segments)
  for (n in 1:num_segments) {
    sample_vals[n]=rbeta(1,in_alphas[n],in_betas[n])
  }
  var_vals=(in_alphas*in_betas)/
    ((in_alphas+in_betas)*(in_alphas+in_betas)*(in_alphas+in_betas+1))
  UCB_vals=sample_vals+omega*var_vals
  
  prob_chosen=(exp(beta*UCB_vals))/sum(exp(beta*UCB_vals))
  chosen_seg=sample(seq(1:num_segments),1,prob=prob_chosen)
  
  return(chosen_seg)
}

