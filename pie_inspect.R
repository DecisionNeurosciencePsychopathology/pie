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

ff$chose_unsampled <- ff$ch

# value sampled on the first free choice as a function of even/uneven sampling -- should be lower in uneven

# how often do they pick the never-sampled option in the uneven condition?
uff <- ff[ff$even_uneven==1,]


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
    choice8 = selected_segment==8,
    s1lag = lag(samplehx1),
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
