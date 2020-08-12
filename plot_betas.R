library(rstan)
library(bayesplot)

model_out1=readRDS('models/RL_softmax_4_seg_1_points.rds')

beta_medians=summary(model_out1,pars="beta")$summary[,'50%']
beta_75s=summary(model_out1,pars="beta")$summary[,'75%']
color_betas=ifelse(beta_medians>0,1,ifelse(beta_75s<0,1.2,1.1))

mcmc_areas(model_out1,pars=vars(`beta[1]`:`beta[95]`),rhat=color_betas)
mcmc_areas_ridges(model_out1,pars=vars(`beta[1]`:`beta[95]`))

test=cbind(beta_medians,beta_75s,color_betas)
which(test[,2]<0)
which(test[,1]<0)
