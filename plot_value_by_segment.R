ggplot(fdf,aes(trial,v_bayes_selected,color=as.factor(selected_segment),
               fill=as.factor(selected_segment),linetype=as.factor(even_uneven)))+
  geom_smooth(method='lm',alpha=0.1)+#geom_count(aes(size = stat(prop)))+
  theme_classic()+#scale_color_manual(values=topo.colors(11))+scale_fill_manual(values=topo.colors(11))+
  facet_wrap(~num_segments*show_points) #*even_uneven)

library(tidyverse)
fdf_4seg=fdf[fdf$num_segments==4,]
fdf_8seg=fdf[fdf$num_segments==8,]
fdf_4seg_1=fdf_4seg_2=fdf_4seg_3=fdf_4seg
fdf_8seg_1=fdf_8seg_2=fdf_8seg_3=fdf_8seg_4=fdf_8seg_5=fdf_8seg_6=fdf_8seg_7=fdf_8seg
fdf_4seg$segment=fdf_8seg$segment=1
fdf_4seg_1$segment=fdf_8seg_1$segment=2
fdf_4seg_2$segment=fdf_8seg_2$segment=3
fdf_4seg_3$segment=fdf_8seg_3$segment=4
fdf_8seg_4$segment=5
fdf_8seg_5$segment=6
fdf_8seg_6$segment=7
fdf_8seg_7$segment=8
fdf_long=rbind(fdf_8seg,fdf_8seg_1,fdf_8seg_2,fdf_8seg_3,fdf_8seg_4,fdf_8seg_5,
               fdf_8seg_6,fdf_8seg_7,fdf_4seg,fdf_4seg_1,fdf_4seg_2,fdf_4seg_3)

fdf_long$seg_chosen=ifelse(fdf_long$selected_segment==fdf_long$segment,1,0)
ggplot(fdf_long,aes(trial,seg_chosen,color=as.factor(selected_segment),
               fill=as.factor(selected_segment),linetype=as.factor(even_uneven)))+
  geom_smooth(method='lm',alpha=0.1)+#geom_count(aes(size = stat(prop)))+
  theme_classic()+#scale_color_manual(values=topo.colors(11))+scale_fill_manual(values=topo.colors(11))+
  facet_wrap(~num_segments*show_points) #*even_uneven)

ggplot(fdf_long[fdf_long$num_segments==4,],aes(trial,seg_chosen,color=as.factor(selected_segment),
               fill=as.factor(selected_segment),linetype=as.factor(even_uneven)))+
  geom_smooth(method='lm',alpha=0.1)+#geom_count(aes(size = stat(prop)))+
  theme_classic()+#scale_color_manual(values=topo.colors(11))+scale_fill_manual(values=topo.colors(11))+
  facet_wrap(~show_points) #*even_uneven)
