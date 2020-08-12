load("pie_data_processed.rdata")

fdf4=fdf[fdf$num_segments==4,]
fdf8=fdf[fdf$num_segments==8,]
fdf4_nopts=fdf4[fdf4$show_points==0,]
fdf4_pts=fdf4[fdf4$show_points==1,]

#average probability of reward < 0.5 ####
nopts4_prob_bysubj=aggregate(fdf4_nopts$selected_prob,by=list(fdf4_nopts$ID),FUN=mean)
hist(nopts4_prob_bysubj$x)
dim(nopts4_prob_bysubj[nopts4_prob_bysubj$x<.5,])[1] #10
nopts4_prob_lt50=nopts4_prob_bysubj[nopts4_prob_bysubj$x<.5,]$Group.1

pts4_prob_bysubj=aggregate(fdf4_pts$selected_prob,by=list(fdf4_pts$ID),FUN=mean)
hist(pts4_prob_bysubj$x)
dim(pts4_prob_bysubj[pts4_prob_bysubj$x<.5,])[1] #7
pts4_prob_lt50=pts4_prob_bysubj[pts4_prob_bysubj$x<.5,]$Group.1

pts4_prob_lt50[pts4_prob_lt50 %in% nopts4_prob_lt50] #1 subject overlaps

a4_prob_bysubj=aggregate(fdf4$selected_prob,by=list(fdf4$ID),FUN=mean)
hist(a4_prob_bysubj$x)
dim(a4_prob_bysubj[a4_prob_bysubj$x<.5,])[1] #6
a4_prob_lt50=a4_prob_bysubj[a4_prob_bysubj$x<.5,]$Group.1

# probability of reward for last 10 trials in blocks
fdf4_nopts_last10=fdf4_nopts[fdf4_nopts$trial>24,]
fdf4_pts_last10=fdf4_pts[fdf4_pts$trial>24,]
fdf4_last10=fdf4[fdf4$trial>24,]

nopts4_prob_bysubj_last10=aggregate(fdf4_nopts_last10$selected_prob,
                                    by=list(fdf4_nopts_last10$ID),FUN=mean)
hist(nopts4_prob_bysubj_last10$x)
dim(nopts4_prob_bysubj_last10[nopts4_prob_bysubj_last10$x<.5,])[1] #17
nopts4_prob_lt50=nopts4_prob_bysubj_last10[nopts4_prob_bysubj_last10$x<.5,]$Group.1

pts4_prob_bysubj_last10=aggregate(fdf4_pts_last10$selected_prob,
                                  by=list(fdf4_pts_last10$ID),FUN=mean)
hist(pts4_prob_bysubj_last10$x)
dim(pts4_prob_bysubj_last10[pts4_prob_bysubj_last10$x<.5,])[1] #13
pts4_prob_lt50=pts4_prob_bysubj_last10[pts4_prob_bysubj_last10$x<.5,]$Group.1

pts4_prob_lt50[pts4_prob_lt50 %in% nopts4_prob_lt50] #4 subjects overlap

a4_prob_bysubj_last10=aggregate(fdf4_last10$selected_prob,
                                by=list(fdf4_last10$ID),FUN=mean)
hist(a4_prob_bysubj_last10$x)
dim(a4_prob_bysubj_last10[a4_prob_bysubj_last10$x<.5,])[1] #8
a4_prob_lt50=a4_prob_bysubj_last10[a4_prob_bysubj_last10$x<.5,]$Group.1

