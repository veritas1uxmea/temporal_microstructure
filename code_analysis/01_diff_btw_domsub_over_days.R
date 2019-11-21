# 01. Differences in behaviors of dominant and subordinate individuals over days


# States
dyaddf.summary.AB <- df %>% 
  group_by(dyadid,idAB,day,state) %>% 
  summarise(freq = n(), total = sum(duration)) %>% 
  as.data.frame()

#add zeros in
dyaddf.AB.allstate<-expand.grid(unique(dyaddf.summary.AB$dyadid), 
                                unique(dyaddf.summary.AB$idAB), 
                                unique(dyaddf.summary.AB$day), 
                                unique(dyaddf.summary.AB$state))

colnames(dyaddf.AB.allstate)<-c("dyadid","idAB", "day","state") 

dyaddf.AB.allstate$total <- dyaddf.AB.allstate$freq <- 0 

dyaddf.summary1.AB <- rbind(dyaddf.summary.AB,dyaddf.AB.allstate) %>% 
  as.data.frame %>%
  group_by(dyadid,idAB,day,state) %>% 
  summarize(freq=sum(freq), total=sum(total)) %>% #add zeros in
  mutate(idABx=ifelse(idAB=="A","Dominant","Subordinate")) %>% 
  mutate(subjectid=paste(dyadid,idABx,sep="-")) %>% 
  as.data.frame()



# Behaviors
dyaddfx.summary.AB <- df %>% 
  group_by(dyadid,idAB,day,behavior) %>% 
  summarise(freq = n(), total = sum(duration)) %>% 
  as.data.frame

#add zeros in
dyaddfx.AB.allbehav<-expand.grid(unique(dyaddfx.summary.AB$dyadid), 
                                 unique(dyaddfx.summary.AB$idAB), 
                                 unique(dyaddfx.summary.AB$day), 
                                 unique(dyaddfx.summary.AB$behavior))

colnames(dyaddfx.AB.allbehav)<-c("dyadid","idAB", "day","behavior") 
dyaddfx.AB.allbehav$total <- dyaddfx.AB.allbehav$freq <- 0 

dyaddfx.summary1.AB <- rbind(dyaddfx.summary.AB,dyaddfx.AB.allbehav) %>% 
  as.data.frame %>%
  group_by(dyadid,idAB,day,behavior) %>% 
  summarize(freq=sum(freq), total=sum(total)) %>% #add zeros in
  mutate(idABx=ifelse(idAB=="A","Dominant","Subordinate")) %>% 
  mutate(subjectid=paste(dyadid,idABx,sep="-"))

## aggressive behavior state: lognormal
aggrdf <- dyaddf.summary1.AB %>% filter(state=="aggr") %>% as.data.frame()
head(aggrdf)
aggrdf$total.x<-aggrdf$total+0.000001

# fit.agg.brms<-brm(total+0.000001~idABx+day +(1|dyadid/subjectid),
#                   family=gaussian(link="log"),data=aggrdf,
#                   control = list(adapt_delta =0.9))
# saveRDS(fit.agg.brms,"data/brms_results/fit.agg.brms.RDS")
fit.agg.brms<-readRDS("data/brms_results/fit.agg.brms.RDS")

summary(fit.agg.brms)


## subordinate behavior state: lognormal
subdf <- dyaddf.summary1.AB %>% filter(state=="sub") %>% as.data.frame()

# fit.sub.brms<-brm(total+0.000001~idABx+day+(1|dyadid/subjectid),
#                   family=gaussian(link="log"),
#                   control = list(adapt_delta =0.99),
#                   data=subdf)
# saveRDS(fit.sub.brms,"data/brms_results/fit.sub.brms.RDS")
fit.sub.brms <- readRDS("data/brms_results/fit.sub.brms.RDS")

summary(fit.sub.brms)



## Figure 1. Individual Differences in Behavior Across Days*
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


#### Subordinate
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

fig1<-arrangeGrob(agg,sub,ncol=1)

# ggsave(fig1,file="img/figure1.tiff",width=300,height=400,units="mm",dpi=300)


## Each behavior 

#### 1. Anogenital Sniffing: lognormal
agdf <- dyaddfx.summary1.AB %>% filter(behavior=="Sniff anogenital") %>% as.data.frame()

# fit.ag.brms<-brm(total+0.000001~idABx+day+(1|dyadid/subjectid),
#                     family=gaussian(link="log"),
#                     control = list(adapt_delta =0.9),
#                     data=agdf)
# saveRDS(fit.ag.brms,"data/brms_results/fit.ag.brms.RDS")
fit.ag.brms <- readRDS("data/brms_results/fit.ag.brms.RDS")

summary(fit.ag.brms)


 
#### 2. Sniff-follow: lognormal
sfdf <- dyaddfx.summary1.AB %>% filter(behavior=="Sniff follow (sniffing while following)") %>% as.data.frame()

# fit.sf.brms<-brm(total+0.000001~idABx+day+(1|dyadid/subjectid),
#                  family=gaussian(link="log"),
#                  control = list(adapt_delta =0.99),
#                  data=sfdf)
# saveRDS(fit.sf.brms,"data/brms_results/fit.sf.brms.RDS")
fit.sf.brms <- readRDS("data/brms_results/fit.sf.brms.RDS")

summary(fit.sf.brms)


#### 3. Sniff-head: Gaussian
 
shdf <- dyaddfx.summary1.AB %>% filter(behavior=="Sniff head") %>% as.data.frame()

# fit.sh.brms<-brm(total~idABx+day+(1|dyadid/subjectid),
#                  family=gaussian,
#                  data=shdf)
# saveRDS(fit.sh.brms,"data/brms_results/fit.sh.brms.RDS")
fit.sh.brms <- readRDS("data/brms_results/fit.sh.brms.RDS")

summary(fit.sh.brms)


#### 4. Sniff-body: Gamma distribution
sbdf <- dyaddfx.summary1.AB %>% filter(behavior=="Sniff body") %>% as.data.frame()

# fit.sb.brms<-brm(total~idABx+day+(1|dyadid/subjectid),
#                  family=gaussian,
#                  data=sbdf)
# saveRDS(fit.sb.brms,"data/brms_results/fit.sb.brms.RDS")
fit.sb.brms <- readRDS("data/brms_results/fit.sb.brms.RDS")

summary(fit.sb.brms)

#### 5. Allogroom: lognormal is the best fit
allodf <- dyaddfx.summary1.AB %>% filter(behavior=="Allogroom") %>% as.data.frame()

# fit.allo.brms<-brm(total+0.000001~idABx+day+(1|dyadid/subjectid),
#                    family=gaussian(),
#                    control = list(adapt_delta =0.99),
#                    data=allodf)
# saveRDS(fit.allo.brms,"data/brms_results/fit.allo.brms.RDS")
fit.allo.brms <- readRDS("data/brms_results/fit.allo.brms.RDS")

summary(fit.allo.brms)


#### 6. Side by side: Gamma  
sbsdf <- dyaddfx.summary1.AB %>% filter(behavior=="Contact side by side (without sniffing)") %>% as.data.frame()

# fit.sbs.brms<-brm(total~idABx+day+(1|dyadid/subjectid),
#                  family=gaussian,
#                  data=sbsdf)
# saveRDS(fit.sbs.brms,"data/brms_results/fit.sbs.brms.RDS")
fit.sbs.brms <- readRDS("data/brms_results/fit.sbs.brms.RDS")

summary(fit.sbs.brms)


## Suppl Figure S1. The durations of investigative and social behaviors exhibited by each mouse in each dyad across 5 days
sniffanogenital<-ggplot(dyaddfx.summary1.AB %>% filter(behavior=="Sniff anogenital"), 
                        aes(x=day, y=total, group=idABx, color=idABx)) + 
  newggtheme +  
  facet_wrap(~dyadid, ncol=7) + 
  geom_path(lwd=1) + 
  ylab("Total duration per day")+
  xlab("Day")+
  scale_color_manual(values=c("firebrick", "dodgerblue")) +  
  scale_y_continuous(breaks=c(0,200))+
  ggtitle("a) Anogenital sniffing")


snifffollow<-ggplot(dyaddfx.summary1.AB %>% filter(behavior=="Sniff follow (sniffing while following)"), 
                    aes(x=day, y=total, group=idABx, color=idABx)) + 
  newggtheme +  
  facet_wrap(~dyadid, ncol=7) + 
  geom_path(lwd=1) + 
  ylab("Total duration per day")+
  xlab("Day")+
  scale_color_manual(values=c("firebrick", "dodgerblue")) +  
  scale_y_continuous(breaks=c(0,150))+
  ggtitle("b) Sniff follow")

sniffhead<-ggplot(dyaddfx.summary1.AB %>% filter(behavior=="Sniff head"), 
                  aes(x=day, y=total, group=idABx, color=idABx)) + 
  newggtheme +  
  facet_wrap(~dyadid, ncol=7) + 
  geom_path(lwd=1) + 
  ylab("Total duration per day")+
  xlab("Day")+
  scale_color_manual(values=c("firebrick", "dodgerblue")) +  
  scale_y_continuous(breaks=c(0,400))+
  ggtitle("c) Sniff head")


sniffbody<-ggplot(dyaddfx.summary1.AB %>% filter(behavior=="Sniff body"), 
                  aes(x=day, y=total, group=idABx, color=idABx)) + 
  newggtheme +  
  facet_wrap(~dyadid, ncol=7) + 
  geom_path(lwd=1) + 
  ylab("Total duration per day")+
  xlab("Day")+
  scale_color_manual(values=c("firebrick", "dodgerblue")) +  
  scale_y_continuous(breaks=c(0,250))+
  ggtitle("d) Sniff body")


allogroom<-ggplot(dyaddfx.summary1.AB %>% filter(behavior=="Allogroom"), 
                  aes(x=day, y=total, group=idABx, color=idABx)) + 
  newggtheme +  
  facet_wrap(~dyadid, ncol=7) + 
  geom_path(lwd=1) + 
  ylab("Total duration per day")+
  xlab("Day")+
  scale_color_manual(values=c("firebrick", "dodgerblue")) +  
  scale_y_continuous(breaks=c(0,250)) +
  ggtitle("e) Allogroom")

sbs<-ggplot(dyaddfx.summary1.AB %>% filter(behavior=="Contact side by side (without sniffing)"), 
            aes(x=day, y=total, group=idABx, color=idABx)) + 
  newggtheme +  
  facet_wrap(~dyadid, ncol=7) + 
  geom_path(lwd=1) + 
  ylab("Total duration per day")+
  xlab("Day")+
  scale_color_manual(values=c("firebrick", "dodgerblue")) +  
  scale_y_continuous(breaks=c(0,300))+
  ggtitle("f) Side by side contact")

sfig1x<-arrangeGrob(sniffanogenital+
                      theme(plot.title = element_text(hjust=-0.07)),
                    snifffollow+
                      theme(plot.title = element_text(hjust=-0.07)),
                    sniffhead+
                      theme(plot.title = element_text(hjust=-0.07)))
sfig1y<-arrangeGrob(sniffbody+
                      theme(plot.title = element_text(hjust=-0.07)),
                    allogroom+
                      theme(plot.title = element_text(hjust=-0.07)),
                    sbs+
                      theme(plot.title = element_text(hjust=-0.07)))
# ggsave(sfig1x,file="img/sfig1x.tiff",width=300,height=600,units="mm",dpi=100)
# ggsave(sfig1y,file="img/sfig1y.tiff",width=300,height=600,units="mm",dpi=100)
