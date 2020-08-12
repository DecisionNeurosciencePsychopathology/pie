# load data and set up ####

library(rstan)
library(loo)
num_cores=3
rstan_options(auto_write = TRUE)
options(mc.cores = num_cores)
load("pie_data_processed.rdata")

#specify data & models to use
segs=8
pts=0
all_model_num=1
excluded_IDs=c("010","016","023","034","041","048","052","061","063","066","067",
               "088","093") #avg prob on last 10 trials of 4 seg w/pts shown < 0.5
#c('006','011','014','015','023','048','062','065','066',
    # '088','089','090','095','40','010','019','026','028','074','097')
# 75% of beta > 0: 002,016,021,022,034,039,041,052,063,067,081,087,093
included_IDs=NA #c("004","008","009","013","014","017","045","049","053","054",
               #"056","059","064","071","073","075")

#select data and create model specifications
df_seg=df[df$num_segments==segs,]
df_seg_pts=df_seg[df_seg$show_points==pts,]
if (!is.na(included_IDs)) {
  df_seg_pts_wex=df_seg_pts[df_seg_pts$ID %in% included_IDs,]
} else {
  df_seg_pts_wex=df_seg_pts[!(df_seg_pts$ID %in% excluded_IDs),]
}
use_data=df_seg_pts_wex
models=c('Bayes_decay_softmax','Bayes_decay_Thompson_Qparam_softmax','Bayes_softmax',
  'Bayes_pers_softmax','Bayes_pers_decay_softmax','RL_decay_softmax', 
  'RL_decay_softmax_expbonus','RL_pers_decay_softmax','RL_pers_softmax','RL_softmax',
  'Bayes_spdecay_softmax')
#parameters in each model: learning rate, inverse temperature, decay, expbonus, perseveration
model_pars=t(matrix(c(0,1,1,0,0, 0,1,1,0,0, 0,1,0,0,0, 
  0,1,0,0,1, 0,1,1,0,1, 1,1,1,0,0, 
  1,1,1,1,0, 1,1,1,0,1, 1,1,0,0,1, 1,1,0,0,0,
  0,1,1,0,0),ncol=length(models)))
par_names=c('alpha','beta','lambda','kappa','pers')

# make list of stan inputs ####
subj_IDs=unique(use_data$ID)
nS=length(subj_IDs)
nT=max(use_data$trial)*max(use_data$block_num)/2

num_segments=array(0,c(nS,nT))
points_shown=array(0,c(nS,nT))
choice=array(0,c(nS,nT))
reward=array(0,c(nS,nT))
block_num=array(0,c(nS,nT))

for (s in 1:nS) {
  subj_data=use_data[use_data$ID==subj_IDs[s],]
  subj_blocks=unique(subj_data$block_num) #not all blocks if only 4 seg conditions
  num_segments[s,]=segs #subj_data$num_segments
  points_shown[s,]=as.numeric(subj_data$show_points)-1
  choice[s,]=subj_data$selected_segment
  reward[s,]=subj_data$win
  block_num[s,]=ifelse(subj_data$block_num==subj_blocks[1],1,
    ifelse(subj_data$block_num==subj_blocks[2],2,
      ifelse(subj_data$block_num==subj_blocks[3],3,4)))
}

data_toest=list(nS=nS,nT=nT,num_segments=num_segments,points_shown=points_shown,
                choice=choice,reward=reward,block_num=block_num)

# run stan model ####
for (m in 1:length(all_model_num)) {
  model_num=all_model_num[m]
  if (length(pts)>1) pt_name='all' else pt_name=pts
  if (length(segs)>1) seg_name='all' else seg_name=segs
  if (is.na(excluded_IDs)) {
    savename=paste(models[model_num],seg_name,'seg',pt_name,'points',sep='_')
  } else {
    savename=paste(models[model_num],seg_name,'seg',pt_name,'points',length(excluded_IDs),
                 'subjex',sep='_')
  }
model_out=stan(file=paste0('estimate_',models[model_num],'.stan'),data=data_toest,verbose=FALSE,
               save_warmup=FALSE,iter=4000,chains=3,control=list(adapt_delta=0.99))
# model_out=stan(file='estimate_RL_softmax.stan',data=data_toest,verbose=FALSE,
#                save_warmup=FALSE,iter=1000,chains=1)
saveRDS(model_out,file=paste0(savename,'.rds'))

# get posterior parameter values and LL ####
num_pars=sum(model_pars[model_num,])
pars_out=array(NA,num_pars*2)
count=1
ind_pars='beta'
for (n in 1:num_pars) {
  while (model_pars[model_num,][count]==0) count=count+1
  ind_pars_pre=c(ind_pars,par_names[count])
  ind_pars=unique(ind_pars_pre)
  pars_out[(n-1)*2+1]=paste0(ind_pars[n],'_m')
  pars_out[n*2]=paste0(ind_pars[n],'_s')
  count=count+1
}
# pars_out=c('alpha_m','alpha_s','beta_m','beta_s') #,'pers_m','pers_s','lambda_m','lambda_s')
# ind_pars=c('alpha','beta')

# unpack LL & calculate integrated AIC per person
model_out_LL=extract_log_lik(model_out,parameter_name="log_lik")
model_out_LL_mat=array(NA,dim=c(nS,nT))
count=1
for (t in 1:nT) {
  for (s in 1:nS) {
    model_out_LL_mat[s,t]=mean(model_out_LL[,count])
    count=count+1
  }
}
model_out_LL_bysubj=rowSums(model_out_LL_mat)
model_out_V1_AIC=2*num_pars-2*model_out_LL_bysubj
save(model_out_LL_bysubj,model_out_V1_AIC,file=paste0(savename,'_modelfit.RData'))

#plot parameters
pairs(model_out,pars=pars_out)
print(summary(model_out,pars=pars_out)$summary)

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