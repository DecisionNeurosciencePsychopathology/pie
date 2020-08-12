data {
  int<lower=1> nS;
  int<lower=1> nT;
  int<lower=1> num_segments[nS,nT];
  int<lower=0,upper=1> points_shown[nS,nT];
  int<lower=1,upper=4> choice[nS,nT]; //segment number of chosen option
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
  real kappa_m;
  
  //group-level variances
  real<lower=0> alpha_s;
  real<lower=0> beta_s;
  real<lower=0> lambda_s;
  real<lower=0> kappa_s;
  
  //subject-specific variances (for non-centered parameterization)
  vector[nS] alpha_raw;
  vector[nS] beta_raw;
  vector[nS] lambda_raw;
  vector[nS] kappa_raw;
}

transformed parameters {
  vector<lower=0,upper=1>[nS] alpha;
  vector[nS] alpha_pre;
  vector[nS] beta;
  vector[nS] lambda;
  vector[nS] kappa;
  
  alpha_pre=alpha_m + alpha_s*alpha_raw;
  alpha=inv_logit(alpha_pre);
  
  beta=beta_m + beta_s*beta_raw;
  lambda=(lambda_m + lambda_s*lambda_raw)/10+1; //rescale so params are approx N(0,1)
  kappa=kappa_m + kappa_s*kappa_raw;
}

model {
  //define variables needed for model estimation
  vector[4] Q;
  vector[4] samplehx;
  vector[4] trial_kappa;
  

  //specify priors
  alpha_m~normal(0,3);
  beta_m~normal(0,10);
  lambda_m~normal(0,1);
  kappa_m~normal(0,1);
  
  alpha_s~cauchy(0,3); 
  beta_s~cauchy(0,5);
  lambda_s~cauchy(0,2);
  kappa_s~cauchy(0,2);
  
  alpha_raw~normal(0,1);
  beta_raw~normal(0,1);
  lambda_raw~normal(0,1);
  kappa_raw~normal(0,1);
  
  for (s in 1:nS) {
    for (t in 1:nT) {
      
      //new block: initialize Q values at 0.5
      if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { 
        for (i in 1:4) {
            Q[i]=0.5;
            samplehx[i]=0;
        }
        // for (i in (num_segments[s,t]+1):8) {
        //   Q[i]=0;
        //   samplehx[i]=0;
        // }
      }
      
      //allocate kappa values to avoid int/real issues
      // for (i in 1:4) {
      //   trial_kappa[i]=samplehx[i]>0 ? kappa[s]*inv(samplehx[i]) : 0;
      // } 
      // for (i in (num_segments[s,t]+1):8) {
      //   trial_kappa[i]=0;
      // }
      // print(samplehx)
      // print(trial_kappa)
      // print(beta[s]*Q[1:num_segments[s,t]]+trial_kappa[1:num_segments[s,t]])
      
      //predict choice only for free choice trials 
      // (note: this assumes # of forced choice trials = # of segments)
      if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]]==0)) {
        print(t)
        print(samplehx)
        print(beta[s]*Q[1:num_segments[s,t]]+kappa[s]*inv(samplehx))
        choice[s,t] ~ categorical_logit(beta[s]*Q+kappa[s]*inv(samplehx));
      }
      
      //update sampling history
      samplehx[choice[s,t]]=samplehx[choice[s,t]]+1;
      
      //update values- done for free & forced choice trials
      for (j in 1:4) {
        Q[j] = choice[s,t]==j ? (Q[j] + alpha[s]*(reward[s,t]-Q[j])) : lambda[s]*Q[j];
      }
      
    }
  }
}
// 
// generated quantities {
//   //this section only computes what is estimated above- use for LL, posterior
//   // checks, etc.
//   //right now, this is only used to compute log likelihood- notice that LL is
//   // computed based on choice given parameters & values, rather than predicting
//   // choice as in model block above
// 
//   //define variables
//   real log_lik[nS,nT];
//   vector[8] Q;
//   vector[8] samplehx;
//   vector[8] trial_kappa;
// 
//   for (s in 1:nS) {
//     for (t in 1:nT) {
// 
//       //new block: initialize Q values
//       if(t==1||(block_num[s,t]-block_num[s,t-1]>0)) { //new block
//         for (i in 1:num_segments[s,t]) {
//             Q[i]=0.5;
//             samplehx[i]=0;
//         }
//         for (i in (num_segments[s,t]+1):8) {
//           Q[i]=0;
//           samplehx[i]=0;
//         }
//       }
//       
//       //allocate kappa values to avoid int/real issues
//       for (i in 1:num_segments[s,t]) {
//         trial_kappa[i]=samplehx[i]>0 ? kappa[s]*inv(samplehx[i]) : 0;
//       } 
//       for (i in (num_segments[s,t]+1):8) {
//         trial_kappa[i]=0;
//       }
//       
// 
//       //calculate likelihood of choice only for free choice trials
//       // (note: this assumes # of forced choice trials = # of segments)
//       if(t>num_segments[s,t]&&(block_num[s,t]-block_num[s,t-num_segments[s,t]+1]==0)) {
//         log_lik[s,t] = categorical_logit_lpmf(choice[s,t]|beta[s]*Q[1:num_segments[s,t]]+trial_kappa[1:num_segments[s,t]]);
//       } else {
//         log_lik[s,t] = 0;
//       }
//       samplehx[choice[s,t]]=samplehx[choice[s,t]]+1;
//       
//       //update values- done for all trials
//       for (j in 1:num_segments[s,t]) {
//         Q[j] = choice[s,t]==j ? (Q[j] + alpha[s]*(reward[s,t]-Q[j])) : lambda[s]*Q[j];
//       }
//       
//     }
//   }
// }
