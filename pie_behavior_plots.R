library(ggplot2)
library(cowplot)
load("pie_data_processed.rdata")
subj_IDs=unique(fdf$ID)

fdf$show_points_alt=.3*(as.numeric(fdf$show_points)-1)+0.7
fdf$seg_dist=c(NA,fdf$next_seg_dist[1:(dim(fdf)[1]-1)])
fdf$switch=c(NA,fdf$next_switch[1:(dim(fdf)[1]-1)])

fdf$stay_max=ifelse(fdf$dBetaMu_isSelectedMax,(1-fdf$switch),0)
fdf$switch_max=ifelse(fdf$dBetaMu_isSelectedMax,fdf$switch,0)
fdf$stay_nomax=ifelse(!fdf$dBetaMu_isSelectedMax,(1-fdf$switch),0)
fdf$switch_nomax=ifelse(!fdf$dBetaMu_isSelectedMax,fdf$switch,0)

fdf$stay_postmax=ifelse(fdf$dBetaMu_isSelectedMax,(1-fdf$next_switch),0)
fdf$switch_postmax=ifelse(fdf$dBetaMu_isSelectedMax,fdf$next_switch,0)
fdf$stay_postnomax=ifelse(!fdf$dBetaMu_isSelectedMax,(1-fdf$next_switch),0)
fdf$switch_postnomax=ifelse(!fdf$dBetaMu_isSelectedMax,fdf$next_switch,0)

fdf4=fdf[fdf$num_segments==4,]
fdf8=fdf[fdf$num_segments==8,]

sample_subj=10

#plot probability of selected segment over time by show points & initial sampling, 
# separate plots for four and eight segments
subj_num=subj_IDs[sample_subj] #random subject
# ggplot(fdf[fdf$ID==subj_num,], aes(trial,selected_prob,color=as.factor(even_uneven)))+
#   geom_violin()+facet_grid(vars(num_segments),vars(show_points))+
#   theme_classic()

ggplot(fdf, aes(trial,selected_prob_alt,color=as.factor(show_points)))+
  geom_count(aes(size = stat(prop))) +
  scale_size_area(max_size = 5)+
  theme_classic()+facet_grid(rows=vars(num_segments))
# not much effect can be visualized

#plot number of times option was chosen over time
ggplot(fdf[fdf$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob)))+
  # geom_count(aes(size=stat(prop)))+scale_size_area(max_size=5)+
  geom_line()+geom_point()+
  theme_classic()+scale_color_manual(values=topo.colors(11))+
  facet_wrap(~num_segments*show_points*even_uneven)

ggplot(fdf,aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(method='lm')+#geom_count(aes(size = stat(prop)))+
  theme_classic()+#scale_color_manual(values=topo.colors(11))+scale_fill_manual(values=topo.colors(11))+
  facet_wrap(~num_segments*show_points*even_uneven)
ggplot(fdf,aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(method='lm')+#geom_count(aes(size = stat(prop)))+
  theme_classic()+#scale_color_manual(values=topo.colors(11))+scale_fill_manual(values=topo.colors(11))+
  facet_wrap(~num_segments*show_points)
ggplot(fdf,aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(method='lm')+#geom_count(aes(size = stat(prop)))+
  theme_classic()+#scale_color_manual(values=topo.colors(11))+scale_fill_manual(values=topo.colors(11))+
  facet_wrap(~num_segments*even_uneven)
ggplot(fdf,aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(method='lm')+#geom_count(aes(size = stat(prop)))+
  theme_classic()+#scale_color_manual(values=topo.colors(11))+scale_fill_manual(values=topo.colors(11))+
  facet_wrap(~num_segments)

# individual plots
ind4=ggplot(fdf4[fdf4$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(se=F,method='gam')+geom_point(aes(alpha=as.numeric(show_points),shape=as.factor(even_uneven)))+
  theme_classic()+scale_alpha(range = c(0.5, 1))+labs(y='number of times option selected',title=paste('Subject ID:',subj_num,sep=' '))
ind8=ggplot(fdf8[fdf8$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(se=F,method='gam')+geom_point(aes(alpha=as.numeric(show_points),shape=as.factor(even_uneven)))+
  theme_classic()+scale_alpha(range = c(0.5, 1))+ylab('number of times option selected')
plot_grid(ind4,ind8,ncol=1)

####use this
#loop through multiple subjects
subj_num=subj_IDs[sample(1:length(subj_IDs),1,replace=T)]
ind4=ggplot(fdf4[fdf4$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(se=F,method='gam')+geom_point(aes(alpha=as.numeric(show_points),shape=as.factor(even_uneven)))+
  theme_classic()+labs(y='number of times option selected',title=paste('Subject ID:',subj_num,sep=' '))+
  scale_alpha(name='Points Shown',breaks=c(1,2),labels=c('No points','Points'),range=c(0.5, 1))+
  scale_shape_manual(values=c(16,17),labels=c('Even','Uneven'),name='Sampling')+
  scale_color_discrete(name='Probability of Reward')+scale_fill_discrete(name='Probability of Reward')
ind8=ggplot(fdf8[fdf8$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(se=F,method='gam')+geom_point(aes(alpha=as.numeric(show_points),shape=as.factor(even_uneven)))+
  theme_classic()+ylab('number of times option selected')+
  scale_alpha(name='Points Shown',breaks=c(1,2),labels=c('No points','Points'),range=c(0.5, 1))+
  scale_shape_manual(values=c(16,17),labels=c('Even','Uneven'),name='Sampling')+
  scale_color_discrete(name='Probability of Reward')+scale_fill_discrete(name='Probability of Reward')
plot_grid(ind4,ind8,ncol=1)

#group segments with similar reward probabilities
fdf4$selected_prob_grp=ifelse(fdf4$selected_prob>.5,1,0)
fdf8$selected_prob_grp=ifelse(fdf8$selected_prob>.55,2,ifelse(fdf8$selected_prob>.45,1,0))

subj_num=subj_IDs[sample(1:length(subj_IDs),1,replace=T)]
ind4g=ggplot(fdf4[fdf4$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob_grp),fill=as.factor(selected_prob_grp)))+
  geom_smooth(se=F,method='gam')+geom_point(aes(alpha=as.numeric(show_points),shape=as.factor(even_uneven)))+
  theme_classic()+scale_alpha(range = c(0.5, 1))+labs(y='number of times option selected',title=paste('Subject ID:',subj_num,sep=' '))
ind8g=ggplot(fdf8[fdf8$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob_grp),fill=as.factor(selected_prob_grp)))+
  geom_smooth(se=F,method='gam')+geom_point(aes(alpha=as.numeric(show_points),shape=as.factor(even_uneven)))+
  theme_classic()+scale_alpha(range = c(0.5, 1))+ylab('number of times option selected')
plot_grid(ind4g,ind8g,ncol=1)

#separate out conditions
subj_num=subj_IDs[sample(1:length(subj_IDs),1,replace=T)]
ind4f=ggplot(fdf4[fdf4$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(se=F,method='gam')+geom_point()+
  theme_classic()+labs(y='number of times option selected',title=paste('Subject ID:',subj_num,sep=' '))+
  facet_wrap(~show_points*even_uneven)
ind8f=ggplot(fdf8[fdf8$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(se=F,method='gam')+geom_point()+
  theme_classic()+ylab('number of times option selected')+
  facet_wrap(~show_points*even_uneven)
plot_grid(ind4f,ind8f,ncol=1)

subj_num=subj_IDs[sample(1:length(subj_IDs),1,replace=T)]
ind4f=ggplot(fdf4[fdf4$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(se=F,method='gam')+geom_point()+
  theme_classic()+labs(y='number of times option selected',title=paste('Subject ID:',subj_num,sep=' '))+
  facet_wrap(~show_points)
ind8f=ggplot(fdf8[fdf8$ID==subj_num,],aes(trial,samplehx_selected,color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(se=F,method='gam')+geom_point()+
  theme_classic()+ylab('number of times option selected')+
  facet_wrap(~show_points)
plot_grid(ind4f,ind8f,ncol=1)


# other plots used to evaluate simulations ####
seg_dist4=ggplot(fdf4,aes(trial,seg_dist,
    color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(method='gam')+theme_classic()+
  ylim(0,2)+labs(y='distance of selected segment from previous choice')+
  scale_color_discrete(name='Probability of Reward')+
  scale_fill_discrete(name='Probability of Reward')+
  theme(legend.position="bottom")+
  facet_wrap(~show_points)
seg_dist8=ggplot(fdf8,aes(trial,seg_dist,
    color=as.factor(selected_prob),fill=as.factor(selected_prob)))+
  geom_smooth(method='gam',alpha=0.3)+theme_classic()+
  ylim(0,3)+labs(y='distance of selected segment from previous choice')+
  scale_color_discrete(name='Probability of Reward')+
  scale_fill_discrete(name='Probability of Reward')+
  theme(legend.position="bottom")+
  facet_wrap(~show_points)
plot_grid(seg_dist4,seg_dist8,ncol=1)

vmax_switch4=ggplot(fdf4,aes(trial))+
  geom_smooth(aes(y=stay_max),method='gam',color='maroon',linetype='solid')+
  geom_smooth(aes(y=stay_nomax),method='gam',color='maroon',linetype='dashed')+
  geom_smooth(aes(y=switch_max),method='gam',color='blue',linetype='solid')+
  geom_smooth(aes(y=switch_nomax),method='gam',color='blue',linetype='dashed')+
  theme_classic()+
  ylim(0,1)+labs(y='Proportion choices')+
  facet_wrap(~show_points)
vmax_switch8=ggplot(fdf8,aes(trial))+
  geom_smooth(aes(y=stay_max),method='gam',color='maroon',linetype='solid')+
  geom_smooth(aes(y=stay_nomax),method='gam',color='maroon',linetype='dashed')+
  geom_smooth(aes(y=switch_max),method='gam',color='blue',linetype='solid')+
  geom_smooth(aes(y=switch_nomax),method='gam',color='blue',linetype='dashed')+
  theme_classic()+
  ylim(0,1)+labs(y='Proportion choices')+
  facet_wrap(~show_points)
plot_grid(vmax_switch4,vmax_switch8,ncol=1)

samplehx4=ggplot(fdf4,aes(trial,samplehx_selected,color=as.factor(selected_prob),
    fill=as.factor(selected_prob)))+geom_smooth(method='lm')+
  theme_classic()+facet_wrap(~show_points)+ylim(0,20)
samplehx8=ggplot(fdf8,aes(trial,samplehx_selected,color=as.factor(selected_prob),
    fill=as.factor(selected_prob)))+geom_smooth(method='lm')+
  theme_classic()+facet_wrap(~show_points)+ylim(0,20)
plot_grid(samplehx4,samplehx8,ncol=1)

ggplot(fdf4,aes(trial))+
  geom_smooth(aes(y=stay_postmax),method='gam',color='maroon',linetype='solid')+
  geom_smooth(aes(y=stay_postnomax),method='gam',color='maroon',linetype='dashed')+
  geom_smooth(aes(y=switch_postmax),method='gam',color='blue',linetype='solid')+
  geom_smooth(aes(y=switch_postnomax),method='gam',color='blue',linetype='dashed')+
  theme_classic()+
  ylim(0,1)+facet_wrap(~show_points*even_uneven)

ggplot(fdf4,aes(trial,win))+
  geom_smooth(method='gam')+theme_classic()+
  coord_cartesian(ylim=c(0.4,0.6))+labs(y='Proportion Rewarded Choices')+
  facet_wrap(~show_points*even_uneven)
