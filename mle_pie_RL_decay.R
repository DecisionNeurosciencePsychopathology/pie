mle_pie_RL_decay = function (alpha, beta, lambda) {

  num_segments=numeric_segments
  choice=selected_segment
  reward=win
  block_num=block_num
  
  alpha_trans=1/(1+exp(-alpha))
  beta_trans=exp(beta)
  lambda_trans=1/(1+exp(-lambda))
  
  num_trials=length(choice)
  Q=matrix(NA,nrow=num_trials+1,ncol=8)
  P_chosen=array(NA,dim=num_trials)
  
  for (t in 1:num_trials) {
    if(t==1||(block_num[t]-block_num[t-1]>0)) { 
      for (i in 1:num_segments[t]) {
        Q[t,i]=0.5;
      }
      for (i in (num_segments[t]+1):8) {
        Q[t,i]=0;
      }
    }
    
    if(t>num_segments[t]&&(block_num[t]-block_num[t-num_segments[t]]==0)) {
      P_chosen[t]=exp(beta_trans*Q[t,choice[t]])/sum(exp(beta_trans*Q[t,1:num_segments[t]]))
    }
    #update values- done for free & forced choice trials
    for (j in 1:num_segments[t]) {
      Q[t+1,j]=ifelse(choice[t]==j,Q[t,j]+alpha_trans*(reward[t]-Q[t,j]),
                      lambda_trans*(Q[t,j]-0.5))
    }
  }
  
  LL=-1*(sum(log(P_chosen),na.rm=T))
  
  return(LL)

}