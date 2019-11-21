# 07. Time-window cross-correlation analysis


# behavcodes=data.frame(behav=c("Lunge/Bite","Subordinate","Tailrattle","Pursuing","Allogroom","Side by side","Digging",
#                               "Sniffing"),
#                       code=c("A","B","C","D","E","F","G","H"))

Aggressive=c("Lunge","Bite")
names(Aggressive)<-"Lunge/Bite" #this is necessary for sttc_ppm_df function... trust me
Subordinate=c("Subordinate posture","Defensive freeze","Flee")
names(Subordinate)<-"Subordinate"
Tailrattle=c("Tailrattle")
names(Tailrattle)<-"Tail rattle"
Pursue=c("Pursue (without sniffing)")
names(Pursue)<-"Pursue"
Allogroom=c("Allogroom")
names(Allogroom)<-"Allogroom"
Sidebyside=c("Contact side by side (without sniffing)")
names(Sidebyside)<-"Side by side"
Digging=c("Digging")
names(Digging)<-"Digging"
Investigative=c("Sniff head","Sniff anogenital","Sniff body","Sniff follow (sniffing while following)")
names(Investigative)<-"Investigative"


behav_list<-list(Aggressive,Subordinate,Tailrattle,Pursue,Allogroom,Sidebyside,Digging,Investigative)

fsttc_initial_list<-list()
for(i in 1:8){
  for(j in 1:8){
    fsttc_initial_list[[((i-1)*8)+j]]<-sttc_ppm_df(df0,behav_list[[i]],behav_list[[j]])
  }
}

fsttc_df<-fsttc_initial_list %>% data.table::rbindlist() %>% select(-names) %>% unique()

# saveRDS(fsttc_df,"data/fsttc_df.RDS")
# 
# fsttc_df<-readRDS("data/fsttc_df.RDS")

table(fsttc_df$Contingency) #126 = 21*3phases*2directions

#### Stat tests - paired Wilcoxon signed rank test ==> Suppl Table S5. 
fsttc_list<-fsttc_df %>% 
  split(.,list(.$Contingency,.$phase)) %>% 
  map(~spread(.,direction,FSTTC)) 

fsttc_stat<-data.frame(phase=NA,Contingency=NA,V=NA,pval=NA,mean_DS=NA,sem_DS=NA,mean_SD=NA,sem_SD=NA)
for (i in 1:102){
  fsttc_stat[i,1]<-as.character(unique(fsttc_list[[i]]$phase))
  fsttc_stat[i,2]<-unique(fsttc_list[[i]]$Contingency)
  fsttc_stat[i,3]<-wilcox.test(fsttc_list[[i]]$'Dom->Sub',fsttc_list[[i]]$'Sub->Dom',exact=F,paired = T)[1]
  fsttc_stat[i,4]<-wilcox.test(fsttc_list[[i]]$'Dom->Sub',fsttc_list[[i]]$'Sub->Dom',exact=F,paired = T)[3]
  fsttc_stat[i,5]<-mean(fsttc_list[[i]]$`Dom->Sub`,na.rm = T)
  fsttc_stat[i,6]<-sem(fsttc_list[[i]]$`Dom->Sub`)
  fsttc_stat[i,7]<-mean(fsttc_list[[i]]$`Sub->Dom`,na.rm = T)
  fsttc_stat[i,8]<-sem(fsttc_list[[i]]$`Sub->Dom`)
}

for (i in 104:116){
  fsttc_stat[i,1]<-as.character(unique(fsttc_list[[i]]$phase))
  fsttc_stat[i,2]<-unique(fsttc_list[[i]]$Contingency)
  fsttc_stat[i,3]<-wilcox.test(fsttc_list[[i]]$'Dom->Sub',fsttc_list[[i]]$'Sub->Dom',exact=F,paired = T)[1]
  fsttc_stat[i,4]<-wilcox.test(fsttc_list[[i]]$'Dom->Sub',fsttc_list[[i]]$'Sub->Dom',exact=F,paired = T)[3]
  fsttc_stat[i,5]<-mean(fsttc_list[[i]]$`Dom->Sub`,na.rm = T)
  fsttc_stat[i,6]<-sem(fsttc_list[[i]]$`Dom->Sub`)
  fsttc_stat[i,7]<-mean(fsttc_list[[i]]$`Sub->Dom`,na.rm = T)
  fsttc_stat[i,8]<-sem(fsttc_list[[i]]$`Sub->Dom`)
}

for (i in 118:192){
  fsttc_stat[i,1]<-as.character(unique(fsttc_list[[i]]$phase))
  fsttc_stat[i,2]<-unique(fsttc_list[[i]]$Contingency)
  fsttc_stat[i,3]<-wilcox.test(fsttc_list[[i]]$'Dom->Sub',fsttc_list[[i]]$'Sub->Dom',exact=F,paired = T)[1]
  fsttc_stat[i,4]<-wilcox.test(fsttc_list[[i]]$'Dom->Sub',fsttc_list[[i]]$'Sub->Dom',exact=F,paired = T)[3]
  fsttc_stat[i,5]<-mean(fsttc_list[[i]]$`Dom->Sub`,na.rm = T)
  fsttc_stat[i,6]<-sem(fsttc_list[[i]]$`Dom->Sub`)
  fsttc_stat[i,7]<-mean(fsttc_list[[i]]$`Sub->Dom`,na.rm = T)
  fsttc_stat[i,8]<-sem(fsttc_list[[i]]$`Sub->Dom`)
}


for (i in c(103,117)){
  fsttc_stat[i,1]<-as.character(unique(fsttc_list[[i]]$phase))
  fsttc_stat[i,2]<-unique(fsttc_list[[i]]$Contingency)
  fsttc_stat[i,3]<-NA
  fsttc_stat[i,4]<-NA
  fsttc_stat[i,5]<-mean(fsttc_list[[i]]$`Dom->Sub`,na.rm = T)
  fsttc_stat[i,6]<-sem(fsttc_list[[i]]$`Dom->Sub`)
  fsttc_stat[i,7]<-mean(fsttc_list[[i]]$`Sub->Dom`,na.rm = T)
  fsttc_stat[i,8]<-sem(fsttc_list[[i]]$`Sub->Dom`)
}

fsttc_summary1<-fsttc_stat %>%filter(phase!="Mid") %>% 
  mutate(result=ifelse(mean_DS>mean_SD,"DS>SD",
                       ifelse(mean_DS<mean_SD,"DS<SD","-"))) %>% 
  mutate(result_pval=ifelse(is.nan(pval),"-",
                            ifelse(pval<0.001,paste(result,"***","(",V,")",sep=""),
                                   ifelse(pval<0.01,paste(result,"**","(",V,")",sep=""),
                                          ifelse(pval<0.05,paste(result,"*","(",V,")",sep="") ,"-"))))) %>%
  mutate(phase=factor(phase,levels=c("Pre","Post"))) %>% 
  select(phase,Contingency,result_pval) %>% 
  spread(phase,result_pval)

fsttc_summary2<-fsttc_stat %>%filter(phase!="Mid") %>% 
  mutate(numbers_ds=paste(round(mean_DS,2),"+/-",round(sem_DS,2))) %>% 
  mutate(phase=factor(phase,levels=c("Pre","Post"))) %>% 
  select(phase,Contingency,numbers_ds) %>% 
  spread(phase,numbers_ds)
colnames(fsttc_summary2)[2:3]<-c("Pre_DS","Post_DS")

fsttc_summary3<-fsttc_stat %>%filter(phase!="Mid") %>% 
  mutate(numbers_sd=paste(round(mean_SD,2),"+/-",round(sem_SD,2))) %>% 
  mutate(phase=factor(phase,levels=c("Pre","Post"))) %>% 
  select(phase,Contingency,numbers_sd) %>% 
  spread(phase,numbers_sd)
colnames(fsttc_summary3)[2:3]<-c("Pre_SD","Post_SD")

fsttc_summary<-left_join(fsttc_summary1,fsttc_summary2) %>% left_join(.,fsttc_summary3) %>% 
  select(Contingency,Pre_DS,Pre_SD,Post_DS,Post_SD,Pre,Post)


fsttc_summary -> suppl_table_s5

### Figure 7.
fig7_fsttc<-fsttc_df %>% 
  filter(phase!="Mid") %>% 
  group_by(phase,direction,Contingency) %>% 
  summarise(mean=mean(FSTTC,na.rm = T),stem=sem(FSTTC)) %>% 
  mutate(upr=mean+stem, lwr=mean-stem) %>%
  ungroup() %>% 
  mutate(phase=factor(phase,levels=c("Pre","Post"))) %>% 
  filter(Contingency=="Lunge/Bite->Lunge/Bite"|
           Contingency=="Lunge/Bite->Subordinate"|
           Contingency=="Lunge/Bite->Tail rattle"|
           Contingency=="Subordinate->Lunge/Bite"|
           Contingency=="Subordinate->Subordinate"|
           Contingency=="Subordinate->Tail rattle"|
           Contingency=="Tail rattle->Lunge/Bite"|
           Contingency=="Tail rattle->Subordinate"|
           Contingency=="Tail rattle->Tail rattle") %>% 
  mutate(Contingency=factor(Contingency,
                            levels = c( "Lunge/Bite->Subordinate", "Subordinate->Subordinate", "Tail rattle->Subordinate",
                                        "Lunge/Bite->Lunge/Bite", "Subordinate->Lunge/Bite", "Tail rattle->Lunge/Bite", 
                                        "Lunge/Bite->Tail rattle", "Subordinate->Tail rattle", "Tail rattle->Tail rattle" )))


tick<-c(0,0.25,0.5,0.75,1)


fig7<-ggplot(fig7_fsttc,aes(phase,mean,ymin=lwr,ymax=upr))+
  geom_path(aes(group=direction,color=direction))+
  geom_point(aes(color=direction),size=2.4)+
  geom_pointrange(aes(color=direction))+
  geom_hline(yintercept=0,linetype="dashed",color="grey")+
  scale_color_manual(values=c("firebrick","dodgerblue"))+
  scale_fill_manual(values=c("firebrick","dodgerblue"))+
  facet_wrap(~Contingency,ncol=3)+
  scale_y_continuous(breaks=tick)+
  ylab("")+xlab("")+
  theme(
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black"),
    plot.background = element_blank(), 
    text = element_text(color = "gray20", size = 5), 
    axis.text = element_text(size = rel(1)), 
    axis.text.x = element_text(color = "gray20", size = rel(1.8)), 
    axis.text.y = element_text(color = "gray20", size = rel(1.4)), 
    axis.title.x = element_text(size = rel(1.2), vjust = 0), 
    axis.title.y = element_text(size = rel(1.2), vjust = 1), 
    axis.ticks.y = element_blank(), 
    axis.ticks.x = element_blank(), 
    strip.text.x = element_text(size = rel(2)),
    legend.position = "bottom",
    legend.key=element_rect(fill=NA),
    legend.title = element_blank(),
    legend.text=element_text(size=rel(1.9))
  )


fig7
ggsave(fig7,file="img/fig7xx.tiff",width=130,height=150,units="mm",dpi=300)



#### Stat tests - paired Wilcoxon signed rank test ==> between pre- post-phases
fsttc_list2<-fsttc_df %>% 
  filter(phase!="Mid") %>% 
  filter(Contingency=="Lunge/Bite->Lunge/Bite"|
           Contingency=="Lunge/Bite->Subordinate"|
           Contingency=="Lunge/Bite->Tail rattle"|
           Contingency=="Subordinate->Lunge/Bite"|
           Contingency=="Subordinate->Subordinate"|
           Contingency=="Subordinate->Tail rattle"|
           Contingency=="Tail rattle->Lunge/Bite"|
           Contingency=="Tail rattle->Subordinate"|
           Contingency=="Tail rattle->Tail rattle") %>% 
  mutate(phase=factor(phase, levels = c("Pre","Post"))) %>% 
  split(.,list(.$Contingency,.$direction)) %>% 
  map(~spread(.,phase,FSTTC)) 
str(fsttc_list2)
summary(fsttc_list2)
length(fsttc_list2)
fsttc_stat2<-data.frame(phase=NA,Contingency=NA,V=NA,pval=NA,mean_pre=NA,sem_pre=NA,mean_post=NA,sem_post=NA)
for (i in 1:18){
  fsttc_stat2[i,1]<-as.character(unique(fsttc_list2[[i]]$direction))
  fsttc_stat2[i,2]<-unique(fsttc_list2[[i]]$Contingency)
  fsttc_stat2[i,3]<-wilcox.test(fsttc_list2[[i]]$'Pre',fsttc_list2[[i]]$'Post',exact=F,paired = T)[1]
  fsttc_stat2[i,4]<-wilcox.test(fsttc_list2[[i]]$'Pre',fsttc_list2[[i]]$'Post',exact=F,paired = T)[3]
  fsttc_stat2[i,5]<-mean(fsttc_list2[[i]]$Pre,na.rm = T)
  fsttc_stat2[i,6]<-sem(fsttc_list2[[i]]$Pre)
  fsttc_stat2[i,7]<-mean(fsttc_list2[[i]]$Post,na.rm = T)
  fsttc_stat2[i,8]<-sem(fsttc_list2[[i]]$Post)
}
 