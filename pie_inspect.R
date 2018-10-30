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


df<-df[!as.logical(df$forced_choice),]
# value sampled on the first free choice as a function of even/uneven sampling -- should be lower in uneven

# inspect relative value and uncertainty signals
df$num_segments <- as.factor(df$num_segments)
df$show_points <- as.factor(df$show_points)
df$even_uneven <- as.factor(df$even_uneven)
df$forced_sampling <- NA
df$forced_sampling[df$even_uneven==0] <- 'uneven'
df$forced_sampling[df$even_uneven==1] <- 'even'


ff <- as.tibble(df[(df$trial==5 & df$num_segments==4) | (df$trial==9 & df$num_segments==8),])

# df$logVrel <- log(df$v_l)
# df$logUrel <- log(df$u_l)

# log-log value-uncertainty relationship
ggplot(df,aes(v_l,u_l,color = num_segments)) + geom_point() + facet_wrap(~ID)

# linear value-uncertainty relationship
ggplot(df,aes(v_l,u_l,color = num_segments)) + geom_point() + facet_wrap(~ID)

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

m4 <- lmer(u_l ~ num_segments * show_points * trial + (1|ID), df)
summary(m4)
car::Anova(m4,'3')

ggplot(df,aes(trial, u_l,color = num_segments)) + geom_smooth(method = 'gam') + facet_wrap(~ID)
ggplot(df,aes(trial, 1-u_l,color = num_segments, lty = show_points)) + geom_smooth(method = 'gam')

ggplot(ff,aes(forced_sampling,samphx_lag,color = show_points, shape = show_points)) + geom_jitter() + facet_wrap(~num_segments)

ggplot(ff,aes(even_uneven,rewhx_lag,color = show_points, shape = show_points)) + geom_jitter() + facet_wrap(~num_segments)

