data {
  int<lower=1> nS;
  int<lower=1> nT;
  int<lower=1> num_segments[nS,nT];
  int<lower=0,upper=1> points_shown[nS,nT];
  int<lower=1,upper=8> choice[nS,nT]; //segment number of chosen option
  int<lower=0,upper=1> reward[nS,nT];
  int<lower=1> block_num[nS,nT];
  //int<lower=0,upper=1> missed_choice[nS,nT]; //are any trials missed? if so, 
  // we'll need this but it's left out for now
}

parameters {
  //group-level means
  real alpha_m;
  real<lower=0> beta_m;
  real lambda_m;
  
  //group-level variances
  real<lower=0> alpha_s;
  real<lower=0> beta_s;
  real<lower=0> lambda_s;
  
  //subject-specific variances (for non-centered parameterization)
  vector[nS] alpha_raw;
  vector[nS] beta_raw;
  vector[nS] lambda_raw;
}

transformed parameters {
  vector<lower=0,upper=1>[nS] alpha;
  vector[nS] alpha_pre;
  vector<lower=0>[nS] beta;
  vector<lower=0,upper=1>[nS] lambda;
  vector[nS] lambda_pre;
  
  alpha_pre=alpha_m + alpha_s*alpha_raw;
  alpha=inv_logit(alpha_pre);
  
  lambda_pre=lambda_m + lambda_s*lambda_raw;
  lambda=inv_logit(lambda_pre);
  
  beta=beta_m + beta_s*beta_raw;
  // lambda=(lambda_m + lambda_s*lambda_raw)/10+1; //rescale so params are approx N(0,1)
}

model {
  //define variables needed for model estimation
  vector[8] Q;
  int unc; //indicator for unchosen options
  int count;

  //specify priors
  alpha_m~normal(0,3);
  beta_m~normal(0,10);
  lambda_m~normal(0,3);
  
  alpha_s~student_t(3,0,3); 
  beta_s~student_t(3,0,5);
  lambda_s~student_t(3,0,3);
  
  alpha_raw~std_normal();
  beta_raw~std_normal();
  lambda_raw~std_normal();
  
  for (s in 1:nS) {
    for (t in 1:nT) {
      
      //new block: initialize Q values at 0.5
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { 
        for (i in 1:num_segments[s,t]) {
            Q[i]=0.5;
        }
        for (i in (num_segments[s,t]+1):8) {
          Q[i]=0;
        }
      }
      
      //predict choice only for free choice trials 
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        choice[s,t] ~ categorical_logit(beta[s]*Q[1:num_segments[s,t]]);
      }
      //update values- done for free & forced choice trials
      for (j in 1:num_segments[s,t]) {
        Q[j] = choice[s,t]==j ? (Q[j] + alpha[s]*(reward[s,t]-Q[j])) : alpha[s]*lambda[s]*Q[j];
      }
      
    }
  }
}

generated quantities {
  //this section only computes what is estimated above- use for LL, posterior
  // checks, etc.
  //right now, this is only used to compute log likelihood- notice that LL is
  // computed based on choice given parameters & values, rather than predicting
  // choice as in model block above

  //define variables
  real log_lik[nS,nT];
  vector[8] Q;
  int unc[8]; //indicator for unchosen options

  for (s in 1:nS) {
    for (t in 1:nT) {

      //new block: initialize Q values
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { //new block
        for (i in 1:num_segments[s,t]) {
            Q[i]=0.5;
        }
        for (i in (num_segments[s,t]+1):8) {
          Q[i]=0;
        }
      }

      //calculate likelihood of choice only for free choice trials
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        log_lik[s,t] = categorical_logit_lpmf(choice[s,t]|beta[s]*Q[1:num_segments[s,t]]);
      } else {
        log_lik[s,t] = uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
      }
      
      //update values- done for all trials
      for (j in 1:num_segments[s,t]) {
        Q[j] = choice[s,t]==j ? (Q[j] + alpha[s]*(reward[s,t]-Q[j])) : alpha[s]*lambda[s]*Q[j];
      }
      
    }
  }
}
