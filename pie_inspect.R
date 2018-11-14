##  inspect data from single subject

setwd("~/code/pie")
library(readr)
library(lme4)
# library(lmerTest)
library(ggplot2)
library(tidyverse)
library(multcompView)
library(stargazer)
source('~/code/R/vif.lme.R')
#
load("pie_data.rdata")

# more processing (can be offloaded to upstream scripts)
df <- as.tibble(pie_data_proc$df)
df = df %>% as_tibble %>% arrange(ID, block_num, trial)
df$num_segments <- as.factor(df$num_segments)
df$show_points <- as.factor(df$show_points)
df$even_uneven <- as.factor(df$even_uneven)
df$forced_sampling <- NA
df$forced_sampling[df$even_uneven==0] <- 'uneven'
df$forced_sampling[df$even_uneven==1] <- 'even'
df$v_mean <- rowMeans(df[c('v_bayes1','v_bayes2','v_bayes3', 'v_bayes4',
                           'v_bayes5','v_bayes6','v_bayes7','v_bayes8')], na.rm =  T)
df$v_diff <- df$v_bayes_selected - df$v_mean
df$mu_max <- apply(df[,c('dBetaMu1','dBetaMu2','dBetaMu3','dBetaMu4',
                         'dBetaMu5','dBetaMu6','dBetaMu7', 'dBetaMu8')], 1, function(x)
                  {max(na.omit(x))})

df$n_unsampled <- apply(df[,c('samplehx1','samplehx2','samplehx3', 'samplehx4',
                              'samplehx5','samplehx6','samplehx7','samplehx8')], 1, function(x) 
                        {sum(na.omit(x)==0)})
# value entropy
df$H <- apply(df[,c('dBetaMu1','dBetaMu2','dBetaMu3','dBetaMu4',
                    'dBetaMu5','dBetaMu6','dBetaMu7', 'dBetaMu8')], 1, function(x)
                      {-sum(na.omit(x)*log(na.omit(x)))})


# need to mean-center entropy by # segments
df$Hscaled[df$num_segments=='8'] <- scale(df$H[df$num_segments=='8'])
df$Hscaled[df$num_segments=='4'] <- scale(df$H[df$num_segments=='4'])

# and by maintenance demand
df$Hscaled_show[df$show_points==1] <- scale(df$Hscaled[df$show_points==1])
df$Hscaled_show[df$show_points==0] <- scale(df$Hscaled[df$show_points==0])



# to eliminate (at least nominally) collinearity between trial and number of segments, adjust for condition
df$trial_adj <- df$trial 
df$trial_adj[df$num_segments=='8'] <- df$trial[df$num_segments=='8'] - 4
# only the free choices
fdf<-df[!as.logical(df$forced_choice),]

ff <- as.tibble(df[(df$trial==5 & df$num_segments==4) | (df$trial==9 & df$num_segments==8),])
uff <- ff[ff$forced_sampling=='uneven',]

# sanity check H timecourse plot -- large scaling difference between 4 and 8
ggplot(fdf, aes(trial,H, color = num_segments, lty = show_points)) + geom_smooth(method = "loess")

ggplot(fdf, aes(trial_adj,Hscaled, color = num_segments, lty = show_points)) + geom_smooth(method = "loess")


varyingvars<-names(df)[grep("[1-9]",names(df))]
ldf<-reshape2::melt(fdf, measure.vars = varyingvars)
ldf$type<-gsub("[0-9]*","",ldf$variable)
ldf <- ldf[ldf$type=='v_bayes',]
ldf <- ldf %>% arrange(ID,block_num,trial, type)
# how many remain unsampled
# ggplot(fdf,aes(trial,n_unsampled, color = num_segments, lty = show_points)) + geom_smooth()

# beta mean
mdf<-reshape2::melt(fdf, measure.vars = varyingvars)
mdf$type<-gsub("[0-9]*","",mdf$variable)
mdf <- mdf[mdf$type=='dBetaMu',]
# beta variance
sdf<-reshape2::melt(fdf, measure.vars = varyingvars)
sdf$type<-gsub("[0-9]*","",sdf$variable)
sdf <- sdf[sdf$type=='dBetaSigmaSquare',]


# subjective Bayesian probabilities by segment
ggplot(ldf,aes(trial,value, color = variable)) + geom_smooth(method = "loess") + facet_wrap(~num_segments)
ggplot(mdf,aes(trial,value, color = variable)) + geom_smooth(method = "loess") + facet_wrap(~num_segments)
# ggplot(sdf,aes(trial,value, color = variable,size = num_segments, lty = show_points)) + geom_smooth() 

#########
# plots of exploitation
# their exploitation is helped by show_points in 8
# selected value
ggplot(fdf,aes(trial,dBetaMu_selected,color = num_segments, lty = show_points)) + geom_smooth(method = "loess") 
# initial values are inflated in even samplign by design
ggplot(fdf,aes(trial, dBetaMu_selected,color = num_segments, lty = show_points)) + 
  geom_smooth(method = 'loess') + facet_wrap(~forced_sampling)
# difference from mean value
ggplot(fdf,aes(trial,v_diff,color = num_segments, lty = show_points)) + geom_smooth(method = "loess") 
# objective value/probability
ggplot(fdf,aes(trial, selected_prob,color = num_segments, lty = show_points)) + geom_smooth() + facet_wrap(~ID)
ggplot(fdf,aes(trial, selected_prob,color = num_segments, lty = show_points)) + 
  geom_smooth(method = 'loess')

#########
# value-uncertainty relationship

# beta variance is by definition not epistemic uncertainty but risk
ggplot(fdf,aes(dBetaMu_selected,dBetaSigmaSquare_selected,color = num_segments, lty = show_points)) + 
  geom_point() + facet_wrap(~ID)

# u is truly uncorrelated with value and closer to epistemic uncertainty
ggplot(fdf,aes(dBetaMu_selected,u,color = num_segments, lty = show_points)) + 
  geom_point() + facet_wrap(~ID)

ggplot(fdf,aes(v_bayes_selected,u,color = num_segments, lty = show_points)) + 
  geom_smooth(method = "loess")
# full data 
ggplot(fdf[fdf$trial>10,],aes(v_bayes_selected,u, color = trial, shape = show_points)) + geom_point() + facet_grid(block_num~ID)
ggsave("uv_share.pdf", height = 20, width = 20)
# # right after forced sampling
# ggplot(ff,aes(v_bayes_selected,u,color = num_segments, lty = show_points)) + geom_smooth(method = "gam")
# ggplot(ff,aes(selected_prob,u,color = num_segments, lty = show_points)) + geom_smooth(method = "gam")

# do they switch from exploration to exploitation

########## 
# models of exploitation

# do they get design probabilities?
m1 <- lmer(selected_prob ~ num_segments * show_points * scale(trial) + 
             forced_sampling  +
             (1|ID), fdf)
summary(m1)
car::Anova(m1,'3')

# ideal Bayesian observer value
m2 <- lmer(dBetaMu_selected ~ num_segments * show_points + trial + (1|ID), fdf)
summary(m2)
car::Anova(m2,'3')

m2 <- lmer(dBetaMu_selected ~ num_segments * show_points + scale(trial_adj) + scale(mu_max) +
           Hscaled_show + (1|ID), fdf)
vif.lme(m2)
summary(m2)
car::Anova(m2,'3')


# does the forced sampling symmetry matter?
m2f <- lmer(dBetaMu_selected ~ num_segments * show_points * scale(trial) + 
              forced_sampling * num_segments   +
              (1|ID), fdf)
summary(m2f)
car::Anova(m2f,'3')
anova(m2,m2f)


# value difference between chosen and best available
m3diff <- lmer(v_diff ~ num_segments * show_points + trial + (1|ID), fdf)
summary(m3diff)
car::Anova(m3diff,'3')

###########
# exploration
# crude measure of uncertainty: u = #samples_of_selected_segment/#trials(i.e. total # samples for normalization)
# factors controlling choice uncertainty
m4 <- lmer(u ~ num_segments * show_points + show_points * scale(trial_adj) +  (1|ID), fdf)
vif.lme(m4)
summary(m4)
car::Anova(m4,'3')
m4v <- lmer(u ~ (num_segments + show_points + scale(mu_max) + scale(trial_adj) ) ^2 + (1|ID), fdf)
vif.lme(m4v)
summary(m4v)
car::Anova(m4v,'3')

m4h <- lmer(u ~ (num_segments + show_points +  scale(Hscaled_show) + scale(trial_adj)) ^2  + (1|ID), fdf)
vif.lme(m4h)
summary(m4h)

# stopped here, this model is interesting but runs into substantial collinearity problems (vif~5)
m4vh <- lmer(u ~ (num_segments + show_points + scale(mu_max) + scale(trial_adj) ) ^2 +
               (num_segments + show_points + scale(Hscaled) + scale(trial_adj) ) ^2 + (1|ID), fdf)
vif.lme(m4vh)
summary(m4vh)
car::Anova(m4vh,'3')

m4vhs <- lmer(u ~ (scale(mu_max) + scale(Hscaled_show) + scale(trial_adj)) ^2 + num_segments * show_points + (1|ID), fdf)
vif.lme(m4vhs)
summary(m4vhs)
car::Anova(m4vhs,'3')



anova(m4,m4v,m4h)

#? interaction with condition and trial
m5vh <- lmer(u ~ num_segments * show_points + show_points * trial * mu_max + num_segments * H + show_points * H + (1|ID), fdf)
summary(m5vh)
car::Anova(m5vh,'3')
plot(emmeans(m5vh, c("mu_max", "H"), by = c("num_segments", "show_points"), at = list(H = c())))

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

# when do they sample the most uncertain option? 'd' stands for directed
dm1 <- glmer(dBetaSigmaSquare_isSelectedMax ~ num_segments * show_points * scale(trial) + (1|ID), family = 'binomial', fdf)
summary(dm1)
car::Anova(dm1,'3')
# adjust for entropy
dm2 <- glmer(dBetaSigmaSquare_isSelectedMax ~ scale(H) * num_segments * scale(trial) + num_segments * show_points + (1|ID), family = 'binomial', 
             glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),data = fdf)
summary(dm2)
car::Anova(dm2,'3')

# and for max available value (temptation to exploit)
dm3 <- glmer(dBetaSigmaSquare_isSelectedMax ~ scale(H) * num_segments * scale(trial) + 
               scale(mu_max) * num_segments * scale(trial) + num_segments * show_points + (1|ID), family = 'binomial', 
             glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),data = fdf)
summary(dm3)
car::Anova(dm3,'3')


ggplot(fdf, aes(trial,as.integer(dBetaSigmaSquare_isSelectedMax), color = num_segments, lty = show_points)) +
       geom_smooth(method = "glm", method.args = list(family = "binomial"))

# entropy dynamics: entropy stays high in 8-show vs. 8-no-show (?selective maintenance)
hm1 <- lmer(H ~ num_segments * show_points * trial + (1|ID), fdf)
summary(hm1)
car::Anova(hm1,'3')


ggplot(fdf,aes(trial, u,color = num_segments, lty = show_points)) + geom_smooth(method = 'loess')
ggplot(fdf,aes(trial, dBetaSigmaSquare_selected,color = num_segments, lty = show_points)) + geom_smooth(method = 'loess')


######
# Find the Bob Wilson uncertainty-driven exploration effect
ggplot(ff[ff$v_bayes_selected==0,],aes(forced_sampling,samplehx_selected,color = v_bayes_selected)) + geom_jitter() + facet_wrap(show_points~num_segments)
ggplot(ff,aes(forced_sampling,u==1,color = v_bayes_selected)) + geom_jitter(width = .4, height = .03 ) + 
   facet_wrap(show_points~num_segments)

# look at beta distribution uncertainty and value statistics
ggplot(fdf,aes(trial, dBetaMu_selected, color = num_segments, lty = show_points)) + geom_smooth(method = "loess")
# NB: variance of the beta is not the same as epistemic uncertainty; it is closer to risk
ggplot(fdf,aes(dBetaMu_selected,u, color = num_segments, shape = show_points)) + geom_point()

ggplot(fdf,aes(trial, dBetaSigmaSquare_selected, color = num_segments, lty = show_points)) + geom_smooth(method = "loess")

sm1 <- lmer(dBetaSigmaSquare_selected ~ num_segments * show_points * trial + (1|ID), fdf)
summary(sm1)
car::Anova(sm1,'3')
m4v <- lmer(u ~ v_max * num_segments * show_points * trial +  (1|ID), fdf)
summary(m4v)
car::Anova(m4v,'3')
anova(m4,m4v)


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

