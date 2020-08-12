library(rstan)
library(loo)

file_postfix='_4_seg_0_points' #_20_subjex'
nS=95 #81 #75

models=c('Bayes_decay_softmax','Bayes_decay_Thompson_Qparam_softmax','Bayes_softmax',
  'Bayes_pers_softmax','Bayes_pers_decay_softmax','RL_decay_softmax', 
  'RL_decay_softmax_expbonus','RL_pers_decay_softmax','RL_pers_softmax','RL_softmax')
model_pars=t(matrix(c(0,1,1,0,0, 0,1,1,0,0, 0,1,0,0,0, 
  0,1,0,0,1, 0,1,1,0,1, 1,1,1,0,0, 
  1,1,1,1,0, 1,1,1,0,1, 1,1,0,0,1, 1,1,0,0,0),ncol=length(models)))
par_names=c('alpha','beta','lambda','kappa','pers')

for (m in 1:length(models)) {
  if (file.exists(paste0('models/',models[m],file_postfix,'.rds'))) {
    cat('\n')
    print(paste0('NOW CHECKING MODEL:',models[m]))
    model_out=readRDS(paste0('models/',models[m],file_postfix,'.rds'))
    num_pars=sum(model_pars[m,])
    pars_out=array(NA,num_pars*2)
    count=1
    ind_pars='beta'
    for (n in 1:num_pars) {
      while (model_pars[m,][count]==0) count=count+1
      ind_pars_pre=c(ind_pars,par_names[count])
      ind_pars=unique(ind_pars_pre)
      pars_out[(n-1)*2+1]=paste0(ind_pars[n],'_m')
      pars_out[n*2]=paste0(ind_pars[n],'_s')
      count=count+1
    }
    num_divergent=sapply(get_sampler_params(model_out,inc_warmup=FALSE),
                        function(x) sum(x[,"divergent__"]))
    print(paste0('number of divergent transitions is:',sum(num_divergent)))
    pairs(model_out,pars=pars_out)
    traceplot(model_out,pars=pars_out)
    print(summary(model_out,pars=pars_out,probs=c(0.025,0.5,0.975))$summary)
    model_sum=summary(model_out,pars=c(pars_out,ind_pars))$summary
    model_sum_rhat=model_sum[model_sum[,10]>1.05,]
    if(dim(model_sum_rhat)[1]>5) {
      print(paste0('number of group and individual parameters with Rhat >1.05 is:',
                   dim(model_sum_rhat)[1],'.'))
    } else if (dim(model_sum_rhat)[1]>0) {
      print('parameters with Rhats >1.05:')
      print(model_sum_rhat[,10])
    } else {
      print('no Rhats greater than 1.05.')
    }
    
    params=array(NA,c(nS,num_pars))
    for (i in 1:num_pars) {
      params[,i]=as.numeric(summary(model_out,pars=ind_pars[i])$summary[,'50%'])
    }
    params=data.frame(params)
    names(params)=ind_pars
    if (num_pars>1) {
    pairs(params,pch=19)
    } else {
      hist(params[,1])
    }
  }
}

