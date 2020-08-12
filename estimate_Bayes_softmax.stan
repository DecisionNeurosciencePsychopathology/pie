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
  real<lower=0> beta_m;
  
  //group-level variances
  real<lower=0> beta_s;
  
  //subject-specific variances (for non-centered parameterization)
  vector[nS] beta_raw;
}

transformed parameters {
  vector[nS] beta;
  
  beta=beta_m + beta_s*beta_raw;
}

model {
  //define variables needed for model estimation
  vector[8] Q;
  vector[8] value_alpha;
  vector[8] value_beta;

  //specify priors
  beta_m~normal(0,10);
  
  beta_s~cauchy(0,5);
  
  beta_raw~normal(0,1);
  
  for (s in 1:nS) {
    for (t in 1:nT) {
      
      //new block: initialize Q values at 0.5
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { 
        for (i in 1:num_segments[s,t]) {
            value_alpha[i]=1;
            value_beta[i]=1;
        }
        for (i in (num_segments[s,t]+1):8) {
          value_alpha[i]=0;
          value_beta[i]=0;
        }
      }
      Q = value_alpha ./ (value_alpha+value_beta); //assume value is mean of dist.
      // print(t)
      // print(value_alpha)
      // print(value_beta)
      
      //predict choice only for free choice trials 
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        choice[s,t] ~ categorical_logit(beta[s]*Q[1:num_segments[s,t]]);
      }
      
      //update distributions per segment- done for free & forced choice trials
      value_alpha[choice[s,t]] = value_alpha[choice[s,t]]+reward[s,t];
      value_beta[choice[s,t]] = value_beta[choice[s,t]] - reward[s,t] + 1;

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
  vector[8] value_alpha;
  vector[8] value_beta;

  for (s in 1:nS) {
    for (t in 1:nT) {

      //new block: initialize Q values
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { //new block
        for (i in 1:num_segments[s,t]) {
          value_alpha[i]=1;
          value_beta[i]=1;
        }
        for (i in (num_segments[s,t]+1):8) {
          value_alpha[i]=0;
          value_beta[i]=0;
        }
      }
      Q = value_alpha ./ (value_alpha+value_beta); //assume value is mean of dist.
      for (i in (num_segments[s,t]+1):8) {
          Q[i]=0;
        }

      //calculate likelihood of choice only for free choice trials
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        log_lik[s,t] = categorical_logit_lpmf(choice[s,t]|beta[s]*Q[1:num_segments[s,t]]);
      } else {
        log_lik[s,t] = uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
      }

      //update distributions per segment- done for free & forced choice trials
      value_alpha[choice[s,t]] = value_alpha[choice[s,t]]+reward[s,t];
      value_beta[choice[s,t]] = value_beta[choice[s,t]] - reward[s,t] + 1;
    }
  }
}
