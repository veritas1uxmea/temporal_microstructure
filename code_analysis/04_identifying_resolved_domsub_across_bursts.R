# 04. Identifying resolved dominant-subordinate relationships across bursts

detachAllPackages()
library(tidyverse)
library(scales)
source('functions/seqfunctions.R')

## 04.1 Phi coefficient method
dfburstA.sum <- getburstfreq0(dfburstANA) #- not including data outside of bursts


#### Phi Coefficient df - not including data outside of bursts
dfburstA.phi <- getphidf(dfburstA.sum, PVAL=.1)
dfphi <- dfburstA.phi %>% mutate(gamma='0.3')

#### Identify resolved bursts (the vertical lines on figure 3)
#first, add in number of behaviors in each burst.
dfburstA.phix <- dfburstA.phi %>% 
  left_join(dfburstA.sum %>% 
              group_by(dyadid,burstvalcume) %>% 
              summarise(totalbehav=sum(total), 
                        totalaggr=sum(total[state=="aggr"]),
                        totalsub=sum(total[state=="sub"]) )%>%     
              rename(burstno = burstvalcume))

head(dfburstA.phix)

resolvedA = findresolution_all(dfburstA.phix, P=.1, behavs=7) #only R not satisfied.
resolvedA  #K lands on an NA so should be +1 (11)
binom.test(1,21) #p-value = 2.098e-05 for only 1/21

resolvedA[11]<-11 #K lands on an NA so should be +1 (11)

maxburstsA <- dfburstA.phix %>% group_by(dyadid) %>% filter(burstno==max(burstno)) %>% .$burstno
names(maxburstsA)=names(resolvedA)
summary(round2(resolvedA/maxburstsA,2)) #average burst resolution occurrence

#### Find time of bursts resolved 
resolvedAx = resolvedA
resolvedAx[is.na(resolvedAx)]<-1000 #give R a ridiculously large number for not resolved relationship

dfburstANA$burstresolved = resolvedAx[match(dfburstANA$dyadid, names(resolvedAx))]
dfburstANA = dfburstANA %>% mutate(prepost = ifelse(burstvalcumeNA<burstresolved, 'pre', 'post'))

resolvetimes = dfburstANA %>% group_by(dyadid) %>% arrange(day,time) %>% filter(prepost=="post") %>%
  filter(row_number()==1) %>% select(dyadid, day, time)  %>% ungroup() %>% arrange(-day,-time) %>% data.frame()


#median and IQR of resolved times of 20 dyads (since one dyad never resolved according to phi coefficient method)
resolvetimes %>% mutate(time1 = time + (1200*(day-1))) %>% .$time1 %>% summary #median 3693=3d,33s;  lqr 2d,289s,  uqr=4day19s


#### Pre-switch versus Post-switch Matrices.
dfburstANApre = dfburstANA %>% filter(prepost=="pre")
dfburstANApost = dfburstANA %>% filter(prepost=="post")

preburstphi <- rbind(expand.grid(dfburstANApre$dyadid, c("A","B"), c("aggr", "sub")) %>% as.data.frame() %>% rename(dyadid=Var1, idAB=Var2, state=Var3) %>% mutate(total=0) ,
                     dfburstANApre %>% group_by(dyadid,idAB,state) %>% summarise(total=n()) %>% as.data.frame() ) %>% 
  group_by(dyadid,idAB,state) %>% summarise(total=sum(total)) %>% 
  split(., .$dyadid) %>%map(~ .$total) %>% map(~ matrix(., 2, 2, byrow=T) ) %>%  map(psych::phi) %>% unlist

postburstphi = rbind(expand.grid(dfburstANApost$dyadid, c("A","B"), c("aggr", "sub")) %>% as.data.frame() %>% rename(dyadid=Var1, idAB=Var2, state=Var3) %>% mutate(total=0) ,
                     dfburstANApost %>% group_by(dyadid,idAB,state) %>% summarise(total=n()) %>% as.data.frame() ) %>% group_by(dyadid,idAB,state) %>% summarise(total=sum(total)) %>% 
  split(., .$dyadid) %>%map(~ .$total) %>% map(~ matrix(., 2, 2, byrow=T) ) %>%  map(psych::phi) %>% unlist

wilcox.test(preburstphi[-17], postburstphi[-12],paired=T)


## Figure 3. Change Over Time in Phi Coefficient 
#dfburstA.phi$vline <- NA

# dfburstA.phi1<-rbind(dfburstA.phi,
#                      data.frame(val=NA,burstno=NA,dyadid=names(resolvedA),pval=NA,pvalx=NA,vline=resolvedA))

dfburstA.phi1<-left_join(dfburstA.phi,
                     data.frame(dyadid=names(resolvedA),vline=resolvedA))

dfburstA.phi1<- dfburstA.phi1 %>% 
  mutate(pvalx=ifelse(is.na(pvalx),NA,ifelse(pvalx=="NS","N.S.","Significant")))

dfburstA.phi1 %>% filter(is.na(vline))
dfburstA.phi1 %>% filter(dyadid=="T")

dfburstA.phi1 %>% filter(dyadid=="U")

library("scales")
integer_breaks <- function(n = 5, ...) {
  breaker <- pretty_breaks(n, ...)
  function(x) {
    breaks <- breaker(x)
    breaks[breaks == floor(breaks)]
  }
}

fig3<-ggplot(dfburstA.phi1 %>% filter(!is.na(pvalx)), aes(burstno, val, color=factor(pvalx))) + 
  geom_hline(yintercept=0, color="grey") +
  geom_vline(aes(xintercept=vline), color="black",size=1.2,linetype="dotted")+
  geom_point() + 
  scale_color_manual(values=c("black", "red"))+
  facet_wrap(~dyadid, scales="free_x", ncol=7) +
  scale_y_continuous(labels=c("-1","","0","","1"))+
  scale_x_continuous(breaks= integer_breaks())+
  newggtheme + xlab("Burst number") + ylab("Phi coefficient")

# ggsave(fig3,file="img/fig3.tiff",width=350,height=185,units="mm",dpi=600)

## 04.2 Absolute difference method 

outres<-NULL
for(i in 1:21){ outres[[i]] = findresolutionAS(dfburstA.sum, id=LETTERS[i])}
dev.off()
names(outres)=LETTERS[1:21]
outres # the number of the bursts when the relationship was resolved according to difference method

resolvedAx #by phi method, R is NA
outresx<-resolvedAx  # the number of the bursts when the relationship was resolved according to phi method
outresx["R"]=NA #by phi method, R is NA

#### resolved times
dfburstANA$burstresolvedx = outres[match(dfburstANA$dyadid, names(outres))]
dfburstANA = dfburstANA %>% mutate(prepostx = ifelse(burstvalcumeNA<burstresolvedx, 'pre', 'post'))

resolvetimesx = dfburstANA %>% group_by(dyadid) %>% arrange(day,time) %>% filter(prepostx=="post") %>%
  filter(row_number()==1) %>% select(dyadid, day, time)  %>% ungroup() %>% arrange(-day,-time) %>% data.frame()  

timesdf<-rbind(
  resolvetimes %>% mutate(time1 = time + (1200*(day-1))) %>% .$time1 %>% summary ,
  resolvetimesx %>% mutate(time1 = time + (1200*(day-1))) %>% .$time1 %>% summary) 

timesdf=as.data.frame.matrix(timesdf)
timesdf$method = c('Phi', "Difference")
colnames(timesdf)=c("min", "lqr", "median","mean","uqr","max","method")


## Figure 4.Differences in aggressive and subordinate behaviors between dominant and subordinate males by burst number within each dyad 
resolution<-data.frame(dyadid=LETTERS[1:21],resolved=c(13,25,9,10,6,12,12,8,1,5,2,1,5,16,1,8,1,9,1,5,1))
resolution$dyadid<-as.character(resolution$dyadid)

fig4df_agg<-dfburstA.sum  %>% filter(state=="aggr") %>% 
  group_by(dyadid,burstvalcume) %>% spread(idAB,total) %>% 
  mutate(lineval=B-A)
fig4df_sub<-dfburstA.sum  %>% filter(state=="sub") %>% 
  group_by(dyadid,burstvalcume) %>% spread(idAB,total) %>% 
  mutate(lineval=B-A)
fig4df<-rbind(fig4df_agg,fig4df_sub) %>% as.data.frame() %>% left_join(.,resolution)


fig4<-ggplot(fig4df,aes(burstvalcume,lineval,color=state))+
  geom_line()+
  geom_hline(yintercept = 0,linetype="dashed",color="grey")+
  facet_wrap(~dyadid,ncol=7,scales="free_x")+
  scale_color_manual(values=c("dodgerblue","red"),
                    labels=c("Aggressive behaviors(Dom-Sub)","Subordinate behaviors(Dom-Sub)") )+
  newggtheme+
  scale_x_continuous(breaks= pretty_breaks())+
  scale_y_continuous(limits = c(-40,40))+
  xlab("Burst number")+
  ylab("Difference")+
  geom_vline(aes(xintercept=resolved), color="black",size=1,linetype="dotted")
  

# ggsave(fig4,"img/fig4.png",width = 400,height = 195,units = "mm",dpi=600)  
