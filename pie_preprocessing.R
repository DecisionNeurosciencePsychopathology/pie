load("pie_data.rdata")
library(lme4)
library(ggplot2)
library(tidyverse)
source("pie_utility.R")

df <- as.tibble(pie_data_proc$df)
df = df %>% as_tibble %>% arrange(ID, block_num, trial)
df$num_segments <- as.factor(df$num_segments)
df$show_points <- as.factor(df$show_points)
df$even_uneven <- as.factor(df$even_uneven)
df$forced_sampling <- NA                #whether initial sampling is even or biased (uneven)
df$forced_sampling[df$even_uneven==0] <- 'uneven'
df$forced_sampling[df$even_uneven==1] <- 'even'

df$v_mean <- rowMeans(df[c('v_bayes1','v_bayes2','v_bayes3', 'v_bayes4','v_bayes5',
                           'v_bayes6','v_bayes7','v_bayes8')], na.rm =  T) 
df$v_diff <- df$v_bayes_selected - df$v_mean
df$mu_max <- apply(df[,c('dBetaMu1','dBetaMu2','dBetaMu3','dBetaMu4','dBetaMu5','dBetaMu6',
                         'dBetaMu7', 'dBetaMu8')], 1, function(x) {max(na.omit(x))})

df$n_unsampled <- apply(df[,c('samplehx1','samplehx2','samplehx3', 'samplehx4','samplehx5',
                              'samplehx6','samplehx7','samplehx8')],1, function(x) {sum(na.omit(x)==0)})
# value entropy
df$H <- apply(df[,c('dBetaMu1','dBetaMu2','dBetaMu3','dBetaMu4','dBetaMu5','dBetaMu6','dBetaMu7', 'dBetaMu8')], 
              1, function(x) {-sum(na.omit(x/norm(na.omit(as.matrix(x))))*log(na.omit(x/norm(na.omit(as.matrix(x))))))})
# need to mean-center entropy by # segments
df$Hscaled[df$num_segments=='8'] <- scale(df$H[df$num_segments=='8'])
df$Hscaled[df$num_segments=='4'] <- scale(df$H[df$num_segments=='4'])
# and by maintenance demand
df$Hscaled_show[df$show_points==1] <- scale(df$Hscaled[df$show_points==1])
df$Hscaled_show[df$show_points==0] <- scale(df$Hscaled[df$show_points==0])

# to eliminate (at least nominally) collinearity between trial and number of segments, adjust for condition
df$trial_adj <- df$trial 
df$trial_adj[df$num_segments=='8'] <- df$trial[df$num_segments=='8'] - 4

#RTs
df$badRT=ifelse(df$RT<.2,1,ifelse(df$RT>2,1,0))

#switch options
df$last_trial_block=c(ifelse(df[2:dim(df)[1],]$trial==df[1:dim(df)[1]-1,]$trial+1,0,1),1)
df$next_switch=c(ifelse(df[2:dim(df)[1],]$selected_segment==df[1:dim(df)[1]-1,]$selected_segment,0,1),NA)
df$next_switch=c(ifelse(df[2:dim(df)[1],]$trial==df[1:dim(df)[1]-1,]$trial+1,df$next_switch,NA),NA)
df$numeric_segments=ifelse(df$num_segments==4,4,8)
df$half=ifelse(df$trial-df$numeric_segments>15,2,1)
df$next_seg_dist=c(df[1:dim(df)[1]-1,]$selected_segment-df[2:dim(df)[1],]$selected_segment,NA)
df$next_seg_dist=ifelse(df$last_trial_block==1,NA,
  ifelse(df$next_seg_dist<0,-1*df$next_seg_dist,df$next_seg_dist))
df$next_seg_dist=
  ifelse(df$num_segments==4,
    ifelse(df$next_seg_dist<3,df$next_seg_dist,ifelse(df$next_seg_dist==3,1,NA)),
  ifelse(df$num_segments==8,
    ifelse(df$next_seg_dist<5,df$next_seg_dist,
      ifelse(df$next_seg_dist==5,3,
        ifelse(df$next_seg_dist==6,2,
          ifelse(df$next_seg_dist==7,1,NA)))),NA))

# note if segment was selected during initial forced choice in each block
df$forced_samplehx1=NA
df$forced_samplehx2=NA
df$forced_samplehx3=NA
df$forced_samplehx4=NA
df$forced_samplehx5=NA
df$forced_samplehx6=NA
df$forced_samplehx7=NA
df$forced_samplehx8=NA
df$forced_sample_sel=NA

for (t in 1:dim(df)[1]) {
  if (df$forced_choice[t]==1) {
    if (df$trial[t]==df$num_segments[t]) {
      df$forced_samplehx1[t]=df$samplehx1[t]
      df$forced_samplehx2[t]=df$samplehx2[t]
      df$forced_samplehx3[t]=df$samplehx3[t]
      df$forced_samplehx4[t]=df$samplehx4[t]
      df$forced_samplehx5[t]=df$samplehx5[t]
      df$forced_samplehx6[t]=df$samplehx6[t]
      df$forced_samplehx7[t]=df$samplehx7[t]
      df$forced_samplehx8[t]=df$samplehx8[t]
    }
  } else {
    df$forced_samplehx1[t]=df$forced_samplehx1[t-1]
    df$forced_samplehx2[t]=df$forced_samplehx2[t-1]
    df$forced_samplehx3[t]=df$forced_samplehx3[t-1]
    df$forced_samplehx4[t]=df$forced_samplehx4[t-1]
    df$forced_samplehx5[t]=df$forced_samplehx5[t-1]
    df$forced_samplehx6[t]=df$forced_samplehx6[t-1]
    df$forced_samplehx7[t]=df$forced_samplehx7[t-1]
    df$forced_samplehx8[t]=df$forced_samplehx8[t-1]
  }
  if (df$selected_segment[t]==1) {
    df$forced_sample_sel[t]=df$forced_samplehx1[t]
  } else if (df$selected_segment[t]==2) {
    df$forced_sample_sel[t]=df$forced_samplehx2[t]
  } else if (df$selected_segment[t]==3) {
    df$forced_sample_sel[t]=df$forced_samplehx3[t]
  } else if (df$selected_segment[t]==4) {
    df$forced_sample_sel[t]=df$forced_samplehx4[t]
  } else if (df$selected_segment[t]==5) {
    df$forced_sample_sel[t]=df$forced_samplehx5[t]
  } else if (df$selected_segment[t]==6) {
    df$forced_sample_sel[t]=df$forced_samplehx6[t]
  } else if (df$selected_segment[t]==7) {
    df$forced_sample_sel[t]=df$forced_samplehx7[t]
  } else if (df$selected_segment[t]==8) {
    df$forced_sample_sel[t]=df$forced_samplehx8[t]
  }
}

# only the free choices
fdf<-df[!as.logical(df$forced_choice),]

#check performance per person
num_subjs=length(unique(fdf$ID))
perf_check=matrix(data=NA,nrow=num_subjs,ncol=8)
for (n in 1:num_subjs) {
  subj_fdf=fdf[fdf$ID==unique(fdf$ID)[n],]
  perf_check[n,1]=unique(fdf$ID)[n]
  
  if (sum(subj_fdf$next_switch,na.rm=T)<dim(subj_fdf)[1]) {
  test_switch_after_wins=chisq.test(subj_fdf$next_switch,subj_fdf$win)
  perf_check[n,2]=test_switch_after_wins$p.value
  perf_check[n,3]=mean(subj_fdf[subj_fdf$win==0,]$next_switch)-mean(subj_fdf[subj_fdf$win==1,]$next_switch)
  perf_check[n,4]=ifelse(perf_check[n,2]<.05 && perf_check[n,3]<0,1,0)
  } 
  
  subj_half=subj_fdf[subj_fdf$half==2,]
  subj_half$higher_prob=ifelse(subj_half$selected_prob>.5,1,0)
  high_value_options_second_half=pbinom(dim(subj_half[subj_half$higher_prob==1,])[1],dim(subj_half)[1],.5)
  perf_check[n,5]=(dim(subj_half[subj_half$higher_prob==1,])[1])/(dim(subj_half)[1])
  perf_check[n,6]=high_value_options_second_half
  perf_check[n,7]=ifelse(high_value_options_second_half<.05,1,0)
  
  perf_check[n,8]=ifelse(perf_check[n,4]==1 && perf_check[n,7]==1,1,0)
}

# only the free choices
fdf<-df[!as.logical(df$forced_choice),]
fdf_r=fdf[fdf$badRT==0,]

#only non-free choices
nfdf<-df[as.logical(df$forced_choice),]
nfdf_r=nfdf[nfdf$badRT==0,]

# first free choice
ff <- as.tibble(df[(df$trial==5 & df$num_segments==4) | (df$trial==9 & df$num_segments==8),])
ff_r=ff[ff$badRT==0,]
#first free choice after unven sampling
uff <- ff[ff$forced_sampling=='uneven',]
uff_r=uff[uff$badRT==0,]

#long format with segments as separate rows, with value (v_bayes)
varyingvars<-names(df)[grep("[1-9]",names(df))]
ldf<-reshape2::melt(fdf, measure.vars = varyingvars)
ldf$type<-gsub("[0-9]*","",ldf$variable)
ldf <- ldf[ldf$type=='v_bayes',]
ldf <- ldf %>% arrange(ID,block_num,trial, type)
ldf_r<-reshape2::melt(fdf_r, measure.vars = varyingvars)
ldf_r$type<-gsub("[0-9]*","",ldf_r$variable)
ldf_r <- ldf_r[ldf_r$type=='v_bayes',]
ldf_r <- ldf_r %>% arrange(ID,block_num,trial, type)

#same for mean of beta dist 
mdf<-reshape2::melt(fdf, measure.vars = varyingvars)
mdf$type<-gsub("[0-9]*","",mdf$variable)
mdf <- mdf[mdf$type=='dBetaMu',]
mdf_r<-reshape2::melt(fdf_r, measure.vars = varyingvars)
mdf_r$type<-gsub("[0-9]*","",mdf_r$variable)
mdf_r <- mdf_r[mdf_r$type=='dBetaMu',]
# same for variance of beta dist
sdf<-reshape2::melt(fdf, measure.vars = varyingvars)
sdf$type<-gsub("[0-9]*","",sdf$variable)
sdf <- sdf[sdf$type=='dBetaSigmaSquare',]
sdf_r<-reshape2::melt(fdf_r, measure.vars = varyingvars)
sdf_r$type<-gsub("[0-9]*","",sdf_r$variable)
sdf_r <- sdf_r[sdf_r$type=='dBetaSigmaSquare',]

#create mean uncertainty & value variables
fdf$mean_u=apply(cbind(fdf$uncertainty1,fdf$uncertainty2,fdf$uncertainty3,fdf$uncertainty4,fdf$uncertainty5,fdf$uncertainty6,fdf$uncertainty7,fdf$uncertainty8),1,mean,na.rm=T)
fdf$median_u=apply(cbind(fdf$uncertainty1,fdf$uncertainty2,fdf$uncertainty3,fdf$uncertainty4,fdf$uncertainty5,fdf$uncertainty6,fdf$uncertainty7,fdf$uncertainty8),1,median,na.rm=T)
fdf$mean_mu=apply(cbind(fdf$dBetaMu1,fdf$dBetaMu2,fdf$dBetaMu3,fdf$dBetaMu4,fdf$dBetaMu5,fdf$dBetaMu6,fdf$dBetaMu7,fdf$dBetaMu8),1,mean,na.rm=T)
fdf$median_mu=apply(cbind(fdf$dBetaMu1,fdf$dBetaMu2,fdf$dBetaMu3,fdf$dBetaMu4,fdf$dBetaMu5,fdf$dBetaMu6,fdf$dBetaMu7,fdf$dBetaMu8),1,median,na.rm=T)
fdf$mean_sigma2=apply(cbind(fdf$dBetaSigmaSquare1,fdf$dBetaSigmaSquare2,fdf$dBetaSigmaSquare3,fdf$dBetaSigmaSquare4,fdf$dBetaSigmaSquare5,
                            fdf$dBetaSigmaSquare6,fdf$dBetaSigmaSquare7,fdf$dBetaSigmaSquare8),1,mean,na.rm=T)
fdf$median_sigma2=apply(cbind(fdf$dBetaSigmaSquare1,fdf$dBetaSigmaSquare2,fdf$dBetaSigmaSquare3,fdf$dBetaSigmaSquare4,fdf$dBetaSigmaSquare5,
                              fdf$dBetaSigmaSquare6,fdf$dBetaSigmaSquare7,fdf$dBetaSigmaSquare8),1,median,na.rm=T)

#convert data to long, with one row per segment per trial
fdf$trial_block_id=paste(fdf$ID,fdf$block_num,fdf$trial,sep="_")
varyingvars<-names(df)[grep("[1-9]",names(df))]
cldf=reshape(as.data.frame(fdf),varying=varyingvars,direction="long",idvar="trial_block_id",sep="",timevar="segment")
ldf=cldf[order(cldf$trial_block_id),]
ldf=ldf[!is.na(ldf$choice),]
ldf$seg_chosen=as.numeric(ldf$segment==ldf$selected_segment)

#add relative value & uncertainty- currently value/mean, but could use median or subtract instead of divide
ldf$rel_u=ldf$uncertainty/ldf$mean_u
ldf$rel_mu=ldf$dBetaMu/ldf$mean_mu
ldf$rel_sigma2=ldf$dBetaSigmaSquare/ldf$mean_sigma2

#rescale categorical variables to reduce collinearity
ldf$num_seg_resc=ifelse(ldf$num_segments==4,-.5,.5)
ldf$sampling_resc=ifelse(ldf$forced_sampling=="even",-.5,.5)
ldf$points_resc=ifelse(ldf$show_points==0,-.5,.5)

#normalize by trial number and then by number of segments
fdf$samplehx_selected_adj=fdf$samplehx_selected/fdf$trial
fdf$samplehx_selected_adj_ns=fdf$samplehx_selected/fdf$trial*as.numeric(fdf$numeric_segments)

fdf$pts_seg=as.factor(as.numeric(fdf$show_points)*3-3+as.numeric(fdf$num_segments))
levels(fdf$pts_seg)=c('4seg_nopts','8seg_nopts','4seg_pts','8seg_pts')

fdf$h_ic_selected=log2(1/fdf$dBetaMu_selected)


save.image('pie_data_processed.RData')
