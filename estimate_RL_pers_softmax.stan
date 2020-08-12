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
  real pers_m;
  
  //group-level variances
  real<lower=0> alpha_s;
  real<lower=0> beta_s;
  real<lower=0> pers_s;
  
  //subject-specific variances (for non-centered parameterization)
  vector[nS] alpha_raw;
  vector[nS] beta_raw;
  vector[nS] pers_raw;
}

transformed parameters {
  vector<lower=0,upper=1>[nS] alpha;
  vector[nS] alpha_pre;
  vector[nS] beta;
  vector[nS] pers;
  
  alpha_pre=alpha_m + alpha_s*alpha_raw;
  alpha=inv_logit(alpha_pre);
  
  beta=beta_m + beta_s*beta_raw;
  pers=pers_m + pers_s*pers_raw;
}

model {
  //define variables needed for model estimation
  vector[8] prev_choice;
  vector[8] Q;

  //specify priors
  alpha_m~normal(0,3);
  beta_m~normal(0,10);
  pers_m~normal(0,1);
  
  alpha_s~cauchy(0,3); 
  beta_s~cauchy(0,5);
  pers_s~cauchy(0,2);
  
  alpha_raw~normal(0,1);
  beta_raw~normal(0,1);
  pers_raw~normal(0,1);
  
  for (s in 1:nS) {
    for (t in 1:nT) {
      
      //new block: initialize Q values at 0.5
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { 
        for (i in 1:num_segments[s,t]) {
            Q[i]=0.5;
            prev_choice[i]=0;
        }
        for (i in (num_segments[s,t]):8) {
          Q[i]=0;
          prev_choice[i]=0;
        }
      }
      
      //predict choice only for free choice trials 
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]+1]==0)) {
        choice[s,t] ~ categorical_logit(beta[s]*Q[1:num_segments[s,t]] + pers[s]*prev_choice[1:num_segments[s,t]]);
          // print(t)
          // print(pers[s])
          // print(beta[s])
          // print(beta[s]*Q[1:num_segments[s,t]])
          // print(pers[s]*prev_choice[1:num_segments[s,t]])
          // for (j in 1:num_segments[s,t]) {
          // print(1/(1+exp(-1*(beta[s]*Q[j]+pers[s]*prev_choice[j]))))
          // }
      }
      
      //update values- done for free & forced choice trials
      Q[choice[s,t]] = Q[choice[s,t]] + alpha[s]*(reward[s,t]-Q[choice[s,t]]);
      for (j in 1:num_segments[s,t]) {
        prev_choice[j] = choice[s,t]==j ? 1 : -1;
      }
      // print(t)
      // print(choice[s,t])
      // print(prev_choice)
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
  vector[8] prev_choice;
  vector[8] Q;

  for (s in 1:nS) {
    for (t in 1:nT) {

      //new block: initialize Q values
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { //new block
        for (i in 1:num_segments[s,t]) {
            Q[i]=0.5;
            prev_choice[i]=0;
        }
        for (i in (num_segments[s,t]):8) {
          Q[i]=uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
          prev_choice[i]=uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings;
        }
      }

      //calculate likelihood of choice only for free choice trials
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]+1]==0)) {
        log_lik[s,t] = categorical_logit_lpmf(choice[s,t]|beta[s]*Q[1:num_segments[s,t]] + 
          pers[s]*prev_choice[1:num_segments[s,t]]);
      } else {
        log_lik[s,t] = uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
      }

      //update values- done for all trials
      Q[choice[s,t]] = Q[choice[s,t]] + alpha[s]*(reward[s,t]-Q[choice[s,t]]);
      for (j in 1:num_segments[s,t]) {
        prev_choice[j] = choice[s,t]==j ? 1 : -1;
      }
    }
  }


}
