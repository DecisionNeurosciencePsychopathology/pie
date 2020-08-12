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

transformed data {
  matrix[4,4] dist4 = [[0,1,2,1],[1,0,1,2],[2,1,0,1],[1,2,1,0]];
  //matrix[8,8] dist8 = [[0,1,2,3,4,3,2,1],[1,0,1,2,3,4,3,2],[2,1,0,1,2,3,4,3],
  //  [3,2,1,0,1,2,3,4],[4,3,2,1,0,1,2,3],[3,4,3,2,1,0,1,2],[2,3,4,3,2,1,0,1].
  //  [1,2,3,4,3,2,1,0]];
}

parameters {
  //group-level means
  real<lower=0> beta_m;
  real<lower=0> lambda_m;
  
  //group-level variances
  real<lower=0> beta_s;
  real<lower=0> lambda_s;
  
  //subject-specific variances (for non-centered parameterization)
  vector[nS] beta_raw;
  vector[nS] lambda_raw;
}

transformed parameters {
  vector[nS] beta;
  vector[nS] lambda;
  
  beta=beta_m + beta_s*beta_raw;
  lambda=lambda_m + lambda_s*lambda_raw;
}

model {
  //define variables needed for model estimation
  vector[8] Q;
  vector[8] value_alpha;
  vector[8] value_beta;
  row_vector[8] seg_dist;
  
  //matrix[4,4] dist4 = [[0,1,2,1],[1,0,1,2],[2,1,0,1],[1,2,1,0]];
  //matrix[8,8] dist8 = [[0,1,2,3,4,3,2,1],[1,0,1,2,3,4,3,2],[2,1,0,1,2,3,4,3],
  //  [3,2,1,0,1,2,3,4],[4,3,2,1,0,1,2,3],[3,4,3,2,1,0,1,2],[2,3,4,3,2,1,0,1].
  //  [1,2,3,4,3,2,1,0]];

  //specify priors
  beta_m~normal(0,5);
  lambda_m~normal(0,2);
  
  beta_s~student_t(3,0,5);
  lambda_s~student_t(3,0,2);
  
  beta_raw~normal(0,1);
  lambda_raw~normal(0,1);
  
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
      
      //predict choice only for free choice trials 
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        choice[s,t] ~ categorical_logit(beta[s]*Q[1:num_segments[s,t]]);
      }
      
      //update distributions per segment- done for free & forced choice trials
      for (j in 1:num_segments[s,t]) {
        if (num_segments[s,t]==4) {
          seg_dist[1:4]=dist4[choice[s,t]];
        }
        value_alpha[j] = choice[s,t]==j ? value_alpha[j]+reward[s,t] : 
          lambda[s]*seg_dist[j]*value_alpha[j];
        value_beta[j] = choice[s,t]==j ? value_beta[j]-reward[s,t]+1 : 
          lambda[s]*seg_dist[j]*value_beta[j];
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
  vector[8] value_alpha;
  vector[8] value_beta;
  row_vector[8] seg_dist;

  for (s in 1:nS) {
    for (t in 1:nT) {

      //new block: initialize Q values
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { //new block
        for (i in 1:num_segments[s,t]) {
          value_alpha[i]=1;
          value_beta[i]=1;
        }
        for (i in (num_segments[s,t]+1):8) {
          value_alpha[i]=uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
          value_beta[i]=uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
        }
      }
      Q = value_alpha ./ (value_alpha+value_beta); //assume value is mean of dist.
      for (i in (num_segments[s,t]+1):8) {
          Q[i]=uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
        }

      //calculate likelihood of choice only for free choice trials
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        log_lik[s,t] = categorical_logit_lpmf(choice[s,t]|beta[s]*Q[1:num_segments[s,t]]);
      } else {
        log_lik[s,t] = uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
      }

      //update distributions per segment- done for free & forced choice trials
      for (j in 1:num_segments[s,t]) {
         if (num_segments[s,t]==4) {
          seg_dist[1:4]=dist4[choice[s,t]];
        }
        value_alpha[j] = choice[s,t]==j ? value_alpha[j]+reward[s,t] : 
          lambda[s]*seg_dist[j]*value_alpha[j];
        value_beta[j] = choice[s,t]==j ? value_beta[j]-reward[s,t]+1 : 
          lambda[s]*seg_dist[j]*value_beta[j];
      }
    }
  }
}
