# 05. Behavioral changes across pre-reolution, middle, & post-resolution phases. 

# Setting up Times into DF.
# Time and Burst number of resolution for each method:

resolveddf <- cbind(rbind(resolvetimes, 
                          data.frame(dyadid="R",day=NA,time=NA)) %>% 
                      mutate(time1 = time + (1200*(day-1))) %>% 
                      arrange(dyadid), 
                    resolvetimesx %>% 
                      mutate(time1 = time + (1200*(day-1))) %>% 
                      arrange(dyadid)
)

resolveddf[,5]<-NULL #just deleting duplicate column 
colnames(resolveddf)<-c("dyadid", "dayphi", "timephi", "totaltimephi", "daydiff", "timediff", "totaltimediff")

resolveddf$burstphi <- outresx #phi method
resolveddf$burstdiff <- outres #diff method

## Add in the info regarding resolution period to pre- / middle- / post-

df0<-df
df0$end <- df0$time+df0$duration
df0$abstime <- df0$time + (1200 * (df0$day-1))
df0$totaltimephi <- resolveddf$totaltimephi[match(df0$dyadid, resolveddf$dyadid)]
df0$prepostphi <- ifelse(df0$abstime<df0$totaltimephi, "pre", "post")
df0$prepostphi[is.na(df0$prepostphi)] <- "pre" #R did not resolve
df0$totaltimediff <- resolveddf$totaltimediff[match(df0$dyadid, resolveddf$dyadid)]
df0$prepostdiff <- ifelse(df0$abstime<df0$totaltimediff, "pre", "post")
df0$phase <- ifelse(df0$prepostphi=="post", "post", 
                    ifelse(df0$prepostdiff=="pre", "pre",
                           "middle"))
write.csv(df0,"data/dominance_df.csv")

emptydf0z <- expand.grid(LETTERS[1:21], LETTERS[1:2], c("post", "middle", "pre"), unique(df0$behavior), 0,0)
colnames(emptydf0z) <- c("dyadid", "idAB", "phase", "behavior", "freq", "total")

df0sumz = df0 %>% group_by(dyadid,idAB,phase,behavior) %>% summarise(freq = n(), total = sum(duration)) %>% as.data.frame %>% rbind(emptydf0z) %>%
  group_by(dyadid,idAB,phase,behavior) %>% summarise(freq=sum(freq), total = sum(total))

df0sumz$phase <- factor(df0sumz$phase, levels=c("pre","middle", "post"))
df0sumz <- df0sumz %>% group_by(dyadid,phase,behavior) %>% mutate(pctfreq = freq / sum(freq))

df0sumz <- df0sumz %>% group_by(dyadid,phase,behavior) %>% mutate(pctfreq = freq / sum(freq), pctduration = total / sum(total))

df0sumz

df0sumz.meds <- df0sumz %>% 
  group_by(phase,behavior,idAB) %>% 
  summarise(medtotal = median(pctduration,na.rm=T), 
            uqr = quantile(pctduration, .75,na.rm=T),
            lqr = quantile(pctduration, .25,na.rm=T))

table(df0sumz.meds$phase)
sum(is.na(df0sumz.meds$phase))

#Order factors.
df0sumz.meds$behavior=factor(df0sumz.meds$behavior,levels=c(
  "Bite",   "Lunge",  "Tailrattle",   "Pursue (without sniffing)",  "Allogroom",  "Digging",                            
  "Flee" , "Defensive freeze",   "Subordinate posture", "Contact side by side (without sniffing)","Idle/nothing","Sniff follow (sniffing while following)","Sniff head","Sniff body" ,"Sniff anogenital","Self-grooming", "Moving", "Rearing",
  "Jumping"))

df0sumz.meds1 <- df0sumz.meds %>% filter(behavior!="Jumping")

### Stats - durations:
options(digits=2,scipen = 3)
df0sumz.list <- df0sumz %>% ungroup() %>% split(., list(.$behavior, .$phase))

df0sumz.list.res <- df0sumz.list %>%
  map(., ~select(., dyadid,idAB,pctduration)) %>%
  map(., ~spread(., idAB,pctduration)) %>%
  map(~ wilcox.test(.$A, .$B,paired=T,exact=F))

df0sumz.list.res.pvals <- df0sumz.list.res %>% map(~ round2(.[3][[1]],3)) %>% unlist

df0sumz.list.res.Wvals <- df0sumz.list.res %>% map(~ .[1][[1]]) %>% unlist

df0sumz.list.res.df <- data.frame(grp = names(df0sumz.list),  W= df0sumz.list.res.Wvals, pval = df0sumz.list.res.pvals)
rownames(df0sumz.list.res.df)<-1:nrow(df0sumz.list.res.df)

df0sumz.list.res.df

### Suppl Table S1. 
suppl_table_s3<-df0sumz.list.res.df %>% 
  mutate(grp=as.character(grp)) %>%
  mutate(phase=sub('.*\\.', '',grp),
         behavior=sub('\\..*', '', grp)) %>% 
  select(phase,behavior,W,pval) %>%
  mutate(sig=ifelse(pval<0.001,"***",ifelse(pval<0.01,"**",ifelse(pval<0.05,"*","")))) %>% 
  mutate(stats=paste(paste(W,pval,sep=" ("),sig,sep=")")) %>% 
  mutate(phase=factor(phase,levels=c("pre","middle","post"))) %>% 
  select(phase,behavior,stats) %>% 
  spread(phase,stats)

### Figure 5.
df0sumz.meds1x<-df0sumz.meds1 %>%
  ungroup %>%
  mutate(idABx=ifelse(idAB=="A","Dominant","Subordinate")) %>% 
  mutate(phase=ifelse(phase=="pre","Pre",ifelse(phase=="middle","Mid","Post"))) %>% 
  mutate(phase=factor(phase,levels=c("Pre","Mid","Post"))) %>% 
  mutate(behaviorx=ifelse(behavior=="Pursue (without sniffing)","Pursue",
                          ifelse(behavior=="Contact side by side (without sniffing)","Contact side by side",
                                 ifelse(behavior== "Sniff follow (sniffing while following)", "Sniff follow", 
                                        ifelse(behavior=="Tailrattle","Tail rattle", as.character(behavior))))))

unique(df0sumz.meds1x$behaviorx)
df0sumz.meds1x$behaviorx<-  factor(df0sumz.meds1x$behaviorx,
                                   levels=c(
  "Lunge","Bite","Tail rattle","Flee","Defensive freeze","Subordinate posture",
  "Sniff head","Sniff body","Sniff anogenital","Sniff follow" ,"Pursue","Allogroom",
  "Contact side by side","Digging" ,"Self-grooming","Moving","Rearing","Idle/nothing")
)

limits <- aes(ymax = uqr, ymin=lqr)

fig5<-ggplot(df0sumz.meds1x, aes(x=phase,y=medtotal, group=idABx,color=idABx)) + 
  geom_path() + 
  geom_point(size=3.5) +
  geom_pointrange(limits) +
  facet_wrap(~behaviorx, ncol=6) + 
  newggtheme +
  scale_color_manual(values=c("firebrick", "dodgerblue"))+
  ylab("Relative proportion of behavior")+
  xlab("Resolution Phase")


# ggsave(fig5,file="img/fig5.png",width=400,height=195,units="mm",dpi=600)



