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

# inx<-df[names(df)[grep("samplehx[0-9]",names(df))]]==0
# df[names(df)[grep("v_bayes[0-9]",names(df))]][inx]<-NA

# calculate different value flavors

df$v_mean <- rowMeans(df[c('v_bayes1','v_bayes2','v_bayes3', 'v_bayes4',
                           'v_bayes5','v_bayes6','v_bayes7','v_bayes8')], na.rm =  T)
df$v_diff <- df$vbay_selected - df$v_mean
df$v_max <- apply(df[,c('v_bayes1','v_bayes2','v_bayes3', 'v_bayes4',
                       'v_bayes5','v_bayes6','v_bayes7','v_bayes8')], 1, max, na.rm = T)

df$n_unsampled <- apply(df[,c('samplehx1','samplehx2','samplehx3', 'samplehx4',
                              'samplehx5','samplehx6','samplehx7','samplehx8')], 1, function(x) 
                        {sum(na.omit(x)==0)})
# value entropy
df$H <- apply(df[,c('dBetaMu1','dBetaMu2','dBetaMu3','dBetaMu4',
                    'dBetaMu5','dBetaMu6','dBetaMu7', 'dBetaMu8')], 1, function(x)
                      {-sum(na.omit(x)*log(na.omit(x)))})


fdf<-df[!as.logical(df$forced_choice),]

ff <- as.tibble(df[(df$trial==5 & df$num_segments==4) | (df$trial==9 & df$num_segments==8),])
uff <- ff[ff$forced_sampling=='uneven',]

# sanity check H timecourse plot -- large scaling difference between 4 and 8
ggplot(fdf, aes(trial,H, color = num_segments, lty = show_points)) + geom_smooth(method = "loess")


varyingvars<-names(df)[grep("[1-9]",names(df))]
ldf<-reshape2::melt(fdf, measure.vars = varyingvars)
ldf$type<-gsub("[0-9]*","",ldf$variable)
ldf <- ldf[ldf$type=='v_bayes',]

# how many remain unsampled
ggplot(fdf,aes(trial,n_unsampled, color = num_segments, lty = show_points)) + geom_smooth()
# subjective Bayesian probabilities by segment
ggplot(ldf,aes(trial,value, color = variable)) + geom_smooth() + facet_wrap(~num_segments)

# their exploitation is helped by show_points in 8
# selected value
ggplot(fdf,aes(trial,vbay_selected,color = num_segments, lty = show_points)) + geom_smooth(method = "loess") 
# initial values are inflated in even samplign by design
ggplot(fdf,aes(trial, vbay_selected,color = num_segments, lty = forced_sampling)) + 
  geom_smooth(method = 'loess') #+ facet_wrap(~ID)
# difference from mean value
ggplot(fdf,aes(trial,v_diff,color = num_segments, lty = show_points)) + geom_smooth(method = "loess") 
# objective value/probability
ggplot(fdf,aes(trial, selected_prob,color = num_segments, lty = show_points)) + geom_smooth() + facet_wrap(~ID)
ggplot(fdf,aes(trial, selected_prob,color = num_segments, lty = show_points)) + 
  geom_smooth(method = 'loess')


# value-uncertainty relationship
ggplot(fdf,aes(vbay_selected,u,color = num_segments, lty = show_points)) + geom_smooth(method = "loess") + facet_wrap(~ID)
ggplot(fdf,aes(vbay_selected,u,color = num_segments, lty = show_points)) + geom_smooth(method = "loess")
# full data 
ggplot(fdf[fdf$trial>10,],aes(vbay_selected,u, color = trial, shape = show_points)) + geom_point() + facet_grid(block_num~ID)
ggsave("uv_share.pdf", height = 20, width = 20)
# # right after forced sampling
# ggplot(ff,aes(vbay_selected,u,color = num_segments, lty = show_points)) + geom_smooth(method = "gam")
# ggplot(ff,aes(selected_prob,u,color = num_segments, lty = show_points)) + geom_smooth(method = "gam")

# do they switch from exploration to exploitation

########## 
# formal look at exploitation

# do they get design probabilities?
m1 <- lmer(selected_prob ~ num_segments * show_points + trial + (1|ID), fdf)
summary(m1)
car::Anova(m1,'3')

# ideal Bayesian observer value
m2 <- lmer(vbay_selected ~ num_segments * show_points + trial + (1|ID), fdf)
summary(m2)
car::Anova(m2,'3')

# value difference between chosen and best available
m3diff <- lmer(v_diff ~ num_segments * show_points + trial + (1|ID), fdf)
summary(m3diff)
car::Anova(m3diff,'3')

###########
# exploration
# crude measure of uncertainty: u = #samples_of_selected_segment/#trials(i.e. total # samples for normalization)
m4 <- lmer(u ~ num_segments * show_points * trial + (1|ID), fdf)
summary(m4)
car::Anova(m4,'3')
m4v <- lmer(u ~ v_max + num_segments * show_points * trial + (1|ID), fdf)
summary(m4v)
car::Anova(m4v,'3')
anova(m4,m4v)

# beta distribution uncertainty (mu reduces to the mean % reinforced, same as v_bayes)
# "s" stands for sigma^2
sm1 <- lmer(dBetaSigmaSquare_selected ~ num_segments * show_points * trial + (1|ID), fdf)
summary(sm1)
car::Anova(sm1,'3')
# this model is really circular -- the value of current choice should not predict its uncertainty, just for illustration
# interesting that covarying for mu does not change the predictors of sigma^2
sm2 <- lmer(dBetaSigmaSquare_selected ~ num_segments * show_points * trial + (1|ID), fdf)
summary(sm2)
car::Anova(sm2,'3')

# does value entropy predict exploration?
sm3 <- lmer(dBetaSigmaSquare_selected ~ H + num_segments * show_points * trial + (1|ID), fdf)
summary(sm3)
car::Anova(sm3,'3')

# entropy dynamics: entropy stays high in 8-show vs. 8-no-show (?selective maintenance)
hm1 <- lmer(H ~ num_segments * show_points * trial + (1|ID), fdf)
summary(hm1)
car::Anova(hm1,'3')


ggplot(fdf,aes(trial, u,color = num_segments, lty = show_points)) + geom_smooth(method = 'loess')
ggplot(fdf,aes(trial, dBetaSigmaSquare_selected,color = num_segments, lty = show_points)) + geom_smooth(method = 'loess')


######
# Find the Bob Wilson uncertainty-driven exploration effect
ggplot(ff[ff$vbay_selected==0,],aes(forced_sampling,samplehx_selected,color = vbay_selected)) + geom_jitter() + facet_wrap(show_points~num_segments)
ggplot(ff,aes(forced_sampling,u==1,color = vbay_selected)) + geom_jitter(width = .4, height = .03 ) + 
   facet_wrap(show_points~num_segments)

# run a logistic model
um1 <- glmer(samplehx_selected==0 ~ num_segments*show_points + (1|ID),uff[,],
             family = binomial(link = "logit"))
summary(um1)

# I don't think this is a valid model, just reminding myself that I did this
um2 <- lmer(u ~ forced_sampling * num_segments + forced_sampling * vbay_selected  
            + forced_sampling * show_points  + vbay_selected *  show_points + num_segments * show_points +  (1|ID),ff[,])
summary(um2)
car::Anova(um2,'3')

# the one I trust, not so much um2
um3 <- lmer(u ~ vbay_selected * num_segments + show_points * num_segments + (1|ID),uff)
summary(um3)
car::Anova(um3,'3')

# compare observed to expected exploration -- no clear prediction for expected because of value confound
u4plus <- sum(uff$u[uff$num_segments==4]==1)
u4minus <- sum(uff$u[uff$num_segments==4]<1)
observed = c(u4plus,u4minus)
expected = c(.25,.75)
chisq.test(x = observed, p = expected)

u8plus <- sum(uff$u[uff$num_segments==8]==1)
u8minus <- sum(uff$u[uff$num_segments==8]<1)
observed = c(u8plus,u8minus)
expected = c(.375,1-.375)
chisq.test(x = observed, p = expected)

