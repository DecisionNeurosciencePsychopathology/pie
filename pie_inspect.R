##  inspect data from single subject

setwd("~/code/pie")
library(readr)
library(lme4)
# library(lmerTest)
library(ggplot2)
library(tidyverse)
library(readr)
library(multcompView)
library(stargazer)

#
load("pie_data.rdata")

df <- as.tibble(pie_data_proc$df)

df = df %>% as_tibble %>% arrange(ID, block_num, trial)

# value sampled on the first free choice as a function of even/uneven sampling -- should be lower in uneven

# inspect relative value and uncertainty signals
df$num_segments <- as.factor(df$num_segments)
df$show_points <- as.factor(df$show_points)
df$even_uneven <- as.factor(df$even_uneven)
df$forced_sampling <- NA
df$forced_sampling[df$even_uneven==0] <- 'uneven'
df$forced_sampling[df$even_uneven==1] <- 'even'

inx<-df[names(df)[grep("samplehx[0-9]",names(df))]]==0
df[names(df)[grep("v_bayes[0-9]",names(df))]][inx]<-NA

# calculate different value flavors

df$v_mean <- rowMeans(df[c('v_bayes1','v_bayes2','v_bayes3', 'v_bayes4',
                           'v_bayes5','v_bayes6','v_bayes7','v_bayes8')], na.rm =  T)
df$v_diff <- df$vbay_selected - df$v_mean

fdf<-df[!as.logical(df$forced_choice),]

ff <- as.tibble(df[(df$trial==5 & df$num_segments==4) | (df$trial==9 & df$num_segments==8),])


varyingvars<-names(df)[grep("[1-9]",names(df))]
ldf<-reshape2::melt(fdf, measure.vars = varyingvars)
ldf$type<-gsub("[0-9]*","",ldf$variable)
ldf <- ldf[ldf$type=='v_bayes',]

ggplot(ldf,aes(trial,value, color = variable)) + geom_smooth() + facet_wrap(~num_segments)

# their exploitation is helped by show_points in 8
ggplot(fdf,aes(trial,vbay_selected,color = num_segments, lty = show_points)) + geom_smooth(method = "loess") 

ggplot(fdf,aes(trial,v_diff,color = num_segments, lty = show_points)) + geom_smooth(method = "loess") 


# linear value-uncertainty relationship
ggplot(df,aes(vbay_selected,u,color = num_segments)) + geom_point() + facet_wrap(~ID)

# do they switch from exploration to exploitation
ggplot(df,aes(trial, selected_prob,color = num_segments, lty = show_points)) + geom_smooth() + facet_wrap(~ID)
ggplot(df,aes(trial, selected_prob,color = num_segments, lty = show_points)) + 
  geom_smooth(method = 'loess')

ggplot(df,aes(trial, v_bayes,color = num_segments, lty = forced_sampling)) + 
  geom_smooth(method = 'loess') + facet_wrap(~ID)

# formal look
m1 <- lmer(selected_prob ~ num_segments * show_points + trial + (1|ID), df)
summary(m1)
car::Anova(m1,'3')

# subjective value
m2 <- lmer(v_l ~ num_segments * show_points + trial + (1|ID), df)
summary(m2)
car::Anova(m2,'3')

m3 <- lmer(v_bayes ~ num_segments * show_points + trial + (1|ID), df)
summary(m3)
car::Anova(m3,'3')

m4 <- lmer(u ~ num_segments * show_points * trial + (1|ID), df)
summary(m4)
car::Anova(m4,'3')

ggplot(df,aes(trial, u,color = num_segments)) + geom_smooth(method = 'gam') + facet_wrap(~ID)
ggplot(df,aes(trial, 1-u,color = num_segments, lty = show_points)) + geom_smooth(method = 'gam')


######
# Find the Bob Wilson uncertainty-driven exploration effect
ggplot(ff[ff$vbay_selected==0,],aes(forced_sampling,samplehx_selected,color = vbay_selected)) + geom_jitter() + facet_wrap(show_points~num_segments)
ggplot(ff,aes(forced_sampling,vbay_selected,color = show_points)) + geom_jitter() + facet_wrap(~num_segments)

# run a logistic model
um1 <- glmer(samplehx_selected==0 ~ num_segments * show_points + (1|ID),ff[ff$forced_sampling=='uneven',],
             family = 'binomial')
summary(um1)

