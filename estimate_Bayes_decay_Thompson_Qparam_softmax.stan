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
  // real<lower=0,upper=1> rand_theta[nT,8]; //specify theta since stan can't compute with random 
  // draws in model block
}

parameters {
  //group-level means
  real<lower=0> beta_m;
  real lambda_m;
  
  //group-level variances
  real<lower=0> beta_s;
  real<lower=0> lambda_s;
  
  //subject-specific variances (for non-centered parameterization)
  vector[nS] beta_raw;
  vector[nS] lambda_raw;
  
  //Q values (as random draws from beta distributions of segments)
  real<lower=0,upper=1> Q[nS,nT,4];
}

transformed parameters {
  vector[nS] beta;
  vector[nS] lambda;
  
  beta=beta_m + beta_s*beta_raw;
  lambda=(lambda_m + lambda_s*lambda_raw)/10+1;
}

model {
  //define variables needed for model estimation
  // vector[4] Q;
  vector[4] value_alpha;
  vector[4] value_beta;
  vector[4] sm_values;

  //specify priors
  beta_m~normal(0,5);
  lambda_m~normal(0,1);
  
  beta_s~student_t(4,0,3);
  lambda_s~student_t(4,0,1);
  
  beta_raw~normal(0,1);
  lambda_raw~normal(0,1);
  
  for (s in 1:nS) {
    for (t in 1:nT) {
      
      //new block: initialize Q values at 0.5
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { 
        for (i in 1:4) { //num_segments[s,t]) {
            value_alpha[i]=1;
            value_beta[i]=1;
        }
        // for (i in (num_segments[s,t]+1):8) {
        //   value_alpha[i]=0;
        //   value_beta[i]=0;
        // }
      }
      
      //Thompsom sampling of segment values
      for (i in 1:num_segments[s,t]) {
        // print(value_alpha)
        // print(value_beta)
        Q[s,t,i] ~ beta(value_alpha[i],value_beta[i]);
      }
      // Q = value_alpha ./ (value_alpha+value_beta); //assume value is mean of dist.

      //predict choice only for free choice trials 
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        // print(block_num[s,t])
        // print(beta[s])
        // print(beta[s]*Q[1:num_segments[s,t]])
        for (i in 1:4) {
          sm_values[i]=beta[s]*Q[s,t,i];
        }
        // choice[s,t] ~ categorical_logit(beta[s]*Q[s,t,]); //[1:num_segments[s,t]]);
        choice[s,t] ~ categorical_logit(sm_values);
      }
      
      //update distributions per segment- done for free & forced choice trials
      for (j in 1:4) { //num_segments[s,t]) {
        value_alpha[j] = choice[s,t]==j ? value_alpha[j]+reward[s,t] : lambda[s]*value_alpha[j];
        value_beta[j] = choice[s,t]==j ? value_beta[j]-reward[s,t]+1 : lambda[s]*value_beta[j];
      }

    }
  }
}

// generated quantities {
//   //this section only computes what is estimated above- use for LL, posterior
//   // checks, etc.
//   //right now, this is only used to compute log likelihood- notice that LL is
//   // computed based on choice given parameters & values, rather than predicting
//   // choice as in model block above
// 
//   //define variables
//   real log_lik[nS,nT];
//   vector[4] Q;
//   vector[4] value_alpha;
//   vector[4] value_beta;
// 
//   for (s in 1:nS) {
//     for (t in 1:nT) {
// 
//       //new block: initialize Q values
//       if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { //new block
//         for (i in 1:4) { //num_segments[s,t]) {
//           value_alpha[i]=1;
//           value_beta[i]=1;
//         }
//         // for (i in (num_segments[s,t]+1):8) {
//         //   value_alpha[i]=uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
//         //   value_beta[i]=uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
//         // }
//       }
//       //Thompsom sampling of segment values
//       for (i in 1:num_segments[s,t]) {
//         Q[i] = beta_cdf(rand_theta[t,i],value_alpha[i],value_beta[i]);
//       }
//       // for (i in (num_segments[s,t]+1):8) {
//       //     Q[i]=uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
//       //   }
// 
//       //calculate likelihood of choice only for free choice trials
//       // (note: this assumes # of forced choice trials = # of segments)
//       if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
//         log_lik[s,t] = categorical_logit_lpmf(choice[s,t]|beta[s]*Q); //[1:num_segments[s,t]]);
//       } else {
//         log_lik[s,t] = uniform_rng(1e-16,1e-15); //0; //change from 0 to prevent Rhat warnings
//       }
// 
//       //update distributions per segment- done for free & forced choice trials
//       for (j in 1:4) { //num_segments[s,t]) {
//         value_alpha[j] = choice[s,t]==j ?  value_alpha[j]+reward[s,t] : lambda[s]*value_alpha[j];
//         value_beta[j] = choice[s,t]==j ? value_beta[j]-reward[s,t]+1 : lambda[s]*value_beta[j];
//       }
//     }
//   }
// }
