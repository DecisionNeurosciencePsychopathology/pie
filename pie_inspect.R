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
ff <- as.tibble(pie_firstfree)


df<-df[!as.logical(df$forced_choice),]
# value sampled on the first free choice as a function of even/uneven sampling -- should be lower in uneven

# how often do they pick the never-sampled option in the uneven condition?
uff <- ff[ff$even_uneven==1,]

# inspect relative value and uncertainty signals
df$num_segments <- as.factor(df$num_segments)
df$show_points <- as.factor(df$show_points)


# df$logVrel <- log(df$v_l)
# df$logUrel <- log(df$u_l)

# log-log value-uncertainty relationship
ggplot(df,aes(logVrel,logUrel,color = num_segments)) + geom_point() + facet_wrap(~ID)

# linear value-uncertainty relationship
ggplot(df,aes(v_l,u_l,color = num_segments)) + geom_point() + facet_wrap(~ID)

# do they switch from exploration to exploitation
ggplot(df,aes(trial, selected_prob,color = num_segments, lty = show_points)) + geom_smooth() + facet_wrap(~ID)

# formal look
m1 <- lmer(selected_prob ~ num_segments + show_points + trial + (1|ID), df)
summary(m1)
car::Anova(m1,'3')

# subjective value
m2 <- lmer(v_l ~ num_segments * show_points + trial + (1|ID), df)
summary(m2)
car::Anova(m2,'3')



ggplot(df,aes(trial, logUrel,color = num_segments)) + geom_smooth() + facet_wrap(~ID)


ggplot(df,aes(trial,v_l,color = selected_segment)) + geom_smooth() + facet_wrap(~ID)



ggplot(df, aes(x = trial, y = selected_prob, color = as.factor(num_segments))) + geom_smooth(method = "gam") + facet_wrap(ID~show_points,ncol = 2)

# get lags
df = df %>% arrange(ID, block_num, trial) %>% group_by(ID, block_num) %>% 
  mutate(
    choice1 = selected_segment==1,
    choice2 = selected_segment==2,
    choice3 = selected_segment==3,
    choice4 = selected_segment==4,
    choice5 = selected_segment==5,
    choice6 = selected_segment==6,
    choice7 = selected_segment==7,
    choice8 = selected_segment==8
    # least_sampled = which(c(cum1,cum2,cum3,cum4))
      ) %>% ungroup()
df = df %>% arrange(ID, trial)

# ggplot(s2, aes(x = trial, y = selected_prob, color = as.factor(num_segments))) + geom_smooth(method = "gam") + facet_wrap(~show_points)
# ggplot(s4, aes(x = trial, y = selected_prob, color = as.factor(num_segments))) + geom_smooth(method = "gam") + facet_wrap(~show_points)

ggplot(df, aes(x = trial, y = selected_prob, color = as.factor(num_segments))) + geom_smooth(method = "gam") + facet_wrap(~show_points)


ggplot(s1[s1$forced_choice==0,], aes(x = selected_segment, color = as.factor(even_uneven))) + geom_histogram(position = "identity") + facet_wrap(~num_segments)
ggplot(s2[s2$forced_choice==0,], aes(x = selected_segment, color = as.factor(even_uneven))) + geom_histogram(position = "identity") + facet_wrap(~num_segments)
ggplot(s4[s4$forced_choice==0,], aes(x = selected_segment, color = as.factor(even_uneven))) + geom_histogram(position = "identity") + facet_wrap(~num_segments)

ggplot(df[df$forced_choice==0,], aes(x = selected_segment, color = as.factor(even_uneven))) + geom_histogram(position = "identity") + facet_wrap(~num_segments)
