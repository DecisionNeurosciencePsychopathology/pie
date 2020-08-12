functions {
  real normal_lb_rng(real mu, real sigma, real lb) {
    real p = normal_cdf(lb, mu, sigma);  // cdf for bounds
    real u = uniform_rng(p, 1);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf for value
  }
}
data {
  int<lower=1> nS;
  int<lower=1> nT;
  int<lower=1> num_segments[nS,nT];
  int<lower=0,upper=1> points_shown[nS,nT];
  // int<lower=1,upper=8> choice[nS,nT]; //segment number of chosen option
  // int<lower=0,upper=1> reward[nS,nT];
  int<lower=1> block_num[nS,nT];
  vector[4] seg_probs;
  //int<lower=0,upper=1> missed_choice[nS,nT]; //are any trials missed? if so, 
  // we'll need this but it's left out for now
}

generated quantities {
  //draw RNGs from prior distribution and simulate choices
  
  //define variables
  int<lower=1,upper=8> choice[nS,nT]; //segment number of chosen option
  vector[8] Q;
  int<lower=0,upper=1> reward[nS,nT];
  
  //specify draws from priors
  real alpha_m = normal_rng(0,3);
  real beta_m = normal_lb_rng(0,5,0);
  real lambda_m = normal_rng(0,5);
  
  real alpha_s = student_t_rng(3,0,5); 
  real beta_s = student_t_rng(3,0,5);
  real lambda_s = student_t_rng(3,0,5);
  
  vector[nS] alpha_raw; 
  vector[nS] beta_raw;
  vector[nS] lambda_raw;
  
  vector<lower=0,upper=1>[nS] alpha;
  vector[nS] alpha_pre;
  vector[nS] beta;
  vector[nS] lambda;
  // vector[nS] lambda_pre;
  
  //re-draw variance parameters if negative (can't bound)
  while (alpha_s<0) {
    alpha_s = student_t_rng(3,0,5);
  }
  while (beta_s<0) {
    beta_s = student_t_rng(3,0,5);
  }
  while (lambda_s<0) {
    lambda_s = student_t_rng(3,0,5);
  }
  
  for (s in 1:nS) {
    alpha_raw[s] = normal_rng(0,1);
    beta_raw[s] = normal_rng(0,1);
    lambda_raw[s] = normal_rng(0,1);
  }
  
  //tranformed values
  alpha_pre=alpha_m + alpha_s*alpha_raw;
  alpha=inv_logit(alpha_pre);
  
  // lambda_pre=lambda_m + lambda_s*lambda_raw;
  // lambda=inv_logit(lambda_pre);
  
  beta=beta_m + beta_s*beta_raw;
  lambda=lambda_m + lambda_s*lambda_raw+1;
  // lambda=(lambda_m + lambda_s*lambda_raw)/10+1; //rescale so params are approx N(0,1)

  for (s in 1:nS) {
    for (t in 1:nT) {

      //new block: initialize Q values
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { //new block
        for (i in 1:num_segments[s,t]) {
          Q[i]=0; //0.5;
        }
        for (i in (num_segments[s,t]+1):8) {
          Q[i]=0;
        }
      }

      //predict only for free choice trials
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        choice[s,t] = categorical_logit_rng(beta[s]*Q[1:num_segments[s,t]]);
      } else {
        choice[s,t] = t-(block_num[s,t]-1)*34;
      }
      
      //reward based on given segment probabilities
      reward[s,t]=bernoulli_rng(seg_probs[choice[s,t]]);
      // print(t)
      // print(seg_probs[choice[s,t]])
      // print(choice[s,t])
      // print(reward[s,t])
      
      
      //update values- done for all trials
      for (j in 1:num_segments[s,t]) {
        Q[j] = choice[s,t]==j ? (Q[j] + alpha[s]*(reward[s,t]-Q[j])) : lambda[s]*Q[j];
      }
      // print(Q[1:4])
      
    }
  }
}
