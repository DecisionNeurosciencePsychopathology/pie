library(bbmle)

load("pie_data_processed.rdata")

alpha_vals=c(0.25,0.5,0.75)
beta_vals=c(1,2,4,8)
lambda_vals=c(0.8,0.9,0.99)

#specify data
segs=4
pts=0
excluded_IDs=NA #c("010","016","023","034","041","048","052","061","063","066","067",
#"088","093") #avg prob on last 10 trials of 4 seg w/pts shown < 0.5
#c('006','011','014','015','023','048','062','065','066',
# '088','089','090','095','40','010','019','026','028','074','097')
# 75% of beta > 0: 002,016,021,022,034,039,041,052,063,067,081,087,093
included_IDs=NA 
#c("004","008","009","013","014","017","045","049","053","054",
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

# set up MLE output
num_subjs=length(unique(use_data$ID))
LL_matrix_mle=matrix(data=NA,nrow=num_subjs,ncol=6)

for (s in 1:num_subjs) {
  print(paste0('running subject #',s,' out of:',num_subjs,' total.'))
  LL_subj=matrix(data=NA,nrow=3*4*3,ncol=6)
  subj_data=use_data[use_data$ID==unique(use_data$ID)[s],]
  count=1 
  
  for (a_grid in 1:length(alpha_vals)) {
    for (b_grid in 1:length(beta_vals)) {
      for (l_grid in 1:length(lambda_vals)) {
        alpha=-log((1-alpha_vals[a_grid])/alpha_vals[a_grid])
        beta=log(beta_vals[b_grid])
        lambda=-log((1-lambda_vals[l_grid])/lambda_vals[l_grid])
        
        fit=bbmle::mle2(mle_pie_RL_decay,start=list(alpha=alpha,beta=beta,
            lambda=lambda),method='Nelder-Mead',data=subj_data,
          control=list(maxit=10000))
        
        #extract parameter estimates if model converged, otherwise NA
        if (exists('fit')&&fit@details$convergence==0) {
          iter_LL=fit@min
          iter_nt=dim(subj_data)[1]
          iter_alpha=1/(1+exp(as.numeric(fit@coef[1])))
          iter_beta=exp(as.numeric(fit@coef[2]))
          iter_lambda=1/(1+exp(as.numeric(fit@coef[3])))
          
          old_fit=fit
          rm(fit)
        } else {
          iter_LL=NA
          iter_nt=NA
          iter_alpha=NA
          iter_beta=NA
          iter_lambda=NA
        }
        
        #add to matrix
        LL_subj[count,]=c(as.numeric(subj_data$ID[1]),iter_LL,iter_nt,iter_alpha,
                          iter_beta,iter_lambda)
        count=count+1
      }
    }
  }
  
  #check if any iterations converged, if so pick lowest -LL and add to output
  if (sum(is.na(LL_subj[,2]))==(3*4*4)) {
    LL_matrix_mle[s,1]=subj_data$ID[1]
    LL_matrix_mle[s,2:(dim(LL_matrix_mle)[2])]=NA
  } else {
    best_LL=which(LL_subj[,2]==min(LL_subj[,2],na.rm=T))
    if (length(best_LL)>1) {
      best_LL=best_LL[1]
    }
    LL_matrix_mle[s,]=LL_subj[best_LL,]
  }
}
LL_matrix_mle=as.data.frame(LL_matrix_mle)
names(LL_matrix_mle)=c('Subj','LL','num_trials','alpha','beta','lambda')
pairs(LL_matrix_mle[,4:6],pch=19)
LL_matrix_mle_beta_lt1000=LL_matrix_mle[LL_matrix_mle$beta<1000,]
pairs(LL_matrix_mle_beta_lt1000[,4:6],pch=19)
cor(LL_matrix_mle$beta,LL_matrix_mle$lambda,method='spearman') #0.647
cor(LL_matrix_mle$alpha,LL_matrix_mle$lambda,method='spearman') #0.069
cor(LL_matrix_mle$alpha,LL_matrix_mle$beta,method='spearman') #0.485
saveRDS(LL_matrix_mle,file='RL_decay_MLE.rds')
hist(LL_matrix_mle$alpha)
hist(LL_matrix_mle$beta)
hist(LL_matrix_mle_beta_lt1000$beta)
hist(LL_matrix_mle$lambda)
