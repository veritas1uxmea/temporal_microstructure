# https://stats.stackexchange.com/questions/190152/visualising-many-variables-in-one-plot/190328

## REVISED Figure 1. Individual Differences in Behavior Across Days*
#### Aggression
agg<-ggplot(dyaddf.summary1.AB %>% filter(state=="aggr"), 
            aes(x=day, y=total, group=idABx, color=idABx)) + 
  newggtheme +  
  facet_wrap(~dyadid, ncol=7) + 
  geom_path(lwd=1) + 
  scale_color_manual(values=c("firebrick", "dodgerblue")) +
  scale_y_continuous(breaks=c(0,200)) +
  ylab("Total duration per day")+
  xlab("Day")+
  ggtitle("a) Aggressive behaviors")+
  theme(plot.title = element_text(hjust=-0.09))


dyaddf.summary1.AB %>% 
  filter(state=="aggr") %>%
  group_by(day,idABx) %>% 
  summarize(mean_total=mean(total,na.rm = T)) %>% 
  mutate(dyadid = "Idontknow") -> mean_agg

dyaddf.summary1.AB %>% 
  filter(state=="aggr") %>% 
  ggplot(aes(x=day, y=total, group=dyadid))+
  geom_line(color="grey")+
  geom_line(data = mean_agg, aes(x=day, y=mean_total),color="red",size = 1.2)+
  newggtheme +  
  facet_wrap(~idABx)+
  ylab("Total duration per day")+
  xlab("Day")+
  ggtitle("a) Aggressive behaviors")+
  theme(plot.title = element_text(hjust=-0.09))+
  scale_y_continuous(breaks=c(0,50,100,150,200)) -> agg_new



#### Subordinate behavior 
sub<-ggplot(dyaddf.summary1.AB %>% filter(state=="sub"), 
            aes(x=day, y=total, group=idABx, color=idABx)) + 
  newggtheme +  
  facet_wrap(~dyadid, ncol=7) + 
  geom_path(lwd=1) + 
  scale_color_manual(values=c("firebrick", "dodgerblue")) +  
  scale_y_continuous(breaks=c(0,200)) +
  ylab("Total duration per day")+
  xlab("Day")+
  ggtitle("b) Subordinate behaviors")+
  theme(plot.title = element_text(hjust=-0.09))


dyaddf.summary1.AB %>% 
  filter(state=="sub") %>%
  group_by(day,idABx) %>% 
  summarize(mean_total=mean(total,na.rm = T)) %>% 
  mutate(dyadid = "Idontknow") -> mean_sub

dyaddf.summary1.AB %>% 
  filter(state=="sub") %>% 
  ggplot(aes(x=day, y=total, group=dyadid))+
  geom_line(color="grey")+
  geom_line(data = mean_sub, aes(x=day, y=mean_total),color="red",size = 1.2)+
  newggtheme +  
  facet_wrap(~idABx)+
  ylab("Total duration per day")+
  xlab("Day")+
  ggtitle("b) Subordinate behaviors")+
  theme(plot.title = element_text(hjust=-0.09))+
  scale_y_continuous(breaks=c(0,50,100,150,200,250)) -> sub_new

fig1_new<-arrangeGrob(agg_new,sub_new,ncol=1)

ggsave(fig1_new,file="img/figure1_revised.tiff",width=220,height=200,units="mm",dpi=300)
