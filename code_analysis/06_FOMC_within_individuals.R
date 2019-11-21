# 06. First-order Markov Chains (Lag sequential) anlaysis within individuals

#### get list of behaviors as codes so can do sequence analysis:
behavcodes=data.frame(behavA=LETTERS[1:19],
                      behavB=letters[1:19],
                      behav=unique(df0$behavior))

df0$code <- ifelse(df0$idAB=="A",
                   LETTERS[1:19][match(df0$behavior,unique(df0$behavior))],
                   letters[1:19][match(df0$behavior,unique(df0$behavior))]
)


#### Make sequences for each dom/sub for pre-/post- resolution

df0 %>% filter(idAB=="A") %>% filter(phase=="pre") %>%
  group_by(dyadid,day) %>% arrange(day,abstime) %>%
  split(., list(.$dyadid)) %>%map(~ .$code) %>%
  map(~ seq_droprpt(.)) %>% map(~collapse(.))  -> domprecodes

df0 %>% filter(idAB=="A") %>% filter(phase=="post") %>%
  group_by(dyadid,day) %>% arrange(day,abstime) %>%
  split(., list(.$dyadid)) %>%map(~ .$code) %>%
  map(~ seq_droprpt(.)) %>% map(~collapse(.))-> dompostcodes

df0 %>% filter(idAB=="B") %>% filter(phase=="pre") %>%
  group_by(dyadid,day) %>% arrange(day,abstime) %>%
  split(., list(.$dyadid)) %>%map(~ .$code)%>%
  map(~ seq_droprpt(.)) %>% map(~ collapse(.)) -> subprecodes

df0 %>% filter(idAB=="B") %>% filter(phase=="post") %>%
  group_by(dyadid,day) %>% arrange(day,abstime) %>%
  split(., list(.$dyadid)) %>% map(~ .$code)%>%
  map(~  seq_droprpt(.)) %>% map(~ collapse(.)) -> subpostcodes

names(domprecodes) #all - but L has too short sequence
names(subpostcodes) #no R - as no R postcode.

#### Make observed frequency matrix
domprecodes.mats = lapply(domprecodes, function(x) obsmatrix(x,prob=F))
dompostcodes.mats = lapply(dom0postcodes, function(x) obsmatrix(x,prob=F))
subprecodes.mats = lapply(subprecodes, function(x) obsmatrix(x,prob=F))
subpostcodes.mats = lapply(subpostcodes, function(x) obsmatrix(x,prob=F))


#### Hierarchical G2 test for each with structural zeros along diagonal 
lrx_inc(domprecodes.mats[[3]]) #example

out1=list(); try(for(i in 1:21){out1[[i]] = lrx_inc(domprecodes.mats[[i]])}, TRUE)
out2=list(); try(for(i in 1:21){out2[[i]] = lrx_inc(dompostcodes.mats[[i]])}, TRUE)
out3=list(); try(for(i in 1:21){out3[[i]] = lrx_inc(subprecodes.mats[[i]])}, TRUE)
out4=list(); try(for(i in 1:21){out4[[i]] = lrx_inc(subpostcodes.mats[[i]])}, TRUE)

round2(lapply(out1, function(x) x$`p-value`) %>% unlist , 4) #1
round2(lapply(out2, function(x) x$`p-value`) %>% unlist , 4)
round2(lapply(out3, function(x) x$`p-value`) %>% unlist , 4) #12
round2(lapply(out4, function(x) x$`p-value`) %>% unlist , 4)

lapply(out1, function(x) x$`p-value`) %>% unlist %>% summary
lapply(out2, function(x) x$`p-value`) %>% unlist %>% summary
lapply(out3, function(x) x$`p-value`) %>% unlist %>% summary
lapply(out4, function(x) x$`p-value`) %>% unlist %>% summary

## Permutation Tests and analysis of the permuted data - ensuring no repeats in permuted sequences 
#### this is hashed out as the following takes a long time - to run, unhash.

# source('code_analysis/06_FOMC_within_individuals_permutation.R') 


### READ IN permutation data
domprecodespermed.vals.df1 <- readRDS("data/permutation/domprepermdf_new.RData")
subprecodespermed.vals.df1 <- readRDS("data/permutation/subprepermdf_new.RData")
dompostcodespermed.vals.df1 <- readRDS("data/permutation/dompostpermdf.RData")
subpostcodespermed.vals.df1 <- readRDS("data/permutation/subpostpermdf.RData")

domprecodespermed.vals.df2 = data.table::rbindlist(domprecodespermed.vals.df1)
subprecodespermed.vals.df2 = data.table::rbindlist(subprecodespermed.vals.df1)
dompostcodespermed.vals.df2 = data.table::rbindlist(dompostcodespermed.vals.df1)
subpostcodespermed.vals.df2 = data.table::rbindlist(subpostcodespermed.vals.df1)


# #if observed == NA should = 0 
domprecodespermed.vals.df2$Observed = ifelse(is.na(domprecodespermed.vals.df2$Observed), 0, domprecodespermed.vals.df2$Observed)
subprecodespermed.vals.df2$Observed = ifelse(is.na(subprecodespermed.vals.df2$Observed), 0, subprecodespermed.vals.df2$Observed)
dompostcodespermed.vals.df2$Observed = ifelse(is.na(dompostcodespermed.vals.df2$Observed), 0, dompostcodespermed.vals.df2$Observed)
subpostcodespermed.vals.df2$Observed = ifelse(is.na(subpostcodespermed.vals.df2$Observed), 0, subpostcodespermed.vals.df2$Observed)

#calculate pvalues
domprecodespermed.counts =  domprecodespermed.vals.df2 %>% group_by(Var1,Var2,dyadid) %>% summarise(higher = sum(value>=Observed)/1000, lower = sum(value<=Observed)/1000)
subprecodespermed.counts =  subprecodespermed.vals.df2 %>% group_by(Var1,Var2,dyadid) %>% summarise(higher = sum(value>=Observed)/1000, lower = sum(value<=Observed)/1000)
dompostcodespermed.counts =  dompostcodespermed.vals.df2 %>% group_by(Var1,Var2,dyadid) %>% summarise(higher = sum(value>=Observed)/1000, lower = sum(value<=Observed)/1000)
subpostcodespermed.counts =  subpostcodespermed.vals.df2 %>% group_by(Var1,Var2,dyadid) %>% summarise(higher = sum(value>=Observed)/1000, lower = sum(value<=Observed)/1000)


#median higher-p values
domprepvalsh = domprecodespermed.counts %>% ungroup() %>% group_by(Var1,Var2) %>%  summarize(higherp = median(higher), lowerp= median(lower)) %>% filter(higherp<.01) %>% data.frame   
subprepvalsh = subprecodespermed.counts %>% ungroup() %>% group_by(Var1,Var2) %>%  summarize(higherp = median(higher), lowerp= median(lower)) %>% filter(higherp<.01) %>% data.frame   
dompostpvalsh = dompostcodespermed.counts %>% ungroup() %>% group_by(Var1,Var2) %>%  summarize(higherp = median(higher), lowerp= median(lower)) %>% filter(higherp<.01) %>% data.frame   
subpostpvalsh = subpostcodespermed.counts %>% ungroup() %>% group_by(Var1,Var2) %>%  summarize(higherp = median(higher), lowerp= median(lower)) %>% filter(higherp<.01) %>% data.frame   

#median lower-p values
domprepvalsl = domprecodespermed.counts %>% ungroup() %>% group_by(Var1,Var2) %>%  summarize(higherp = median(higher), lowerp= median(lower)) %>% filter(lowerp<.01) %>% data.frame   
subprepvalsl = subprecodespermed.counts %>% ungroup() %>% group_by(Var1,Var2) %>%  summarize(higherp = median(higher), lowerp= median(lower)) %>% filter(lowerp<.01) %>% data.frame   
dompostpvalsl = dompostcodespermed.counts %>% ungroup() %>% group_by(Var1,Var2) %>%  summarize(higherp = median(higher), lowerp= median(lower)) %>% filter(lowerp<.01) %>% data.frame   
subpostpvalsl = subpostcodespermed.counts %>% ungroup() %>% group_by(Var1,Var2) %>%  summarize(higherp = median(higher), lowerp= median(lower)) %>% filter(lowerp<.01) %>% data.frame   

### Figure 6: Get df to plot kinetogram figure

dompostcodes.obsvals<-readRDS("data/permutation/dompostcodes.obsvals.RDS")
domprecodes.obsvals<-readRDS("data/permutation/domprecodes.obsvals.RDS")
subprecodes.obsvals<-readRDS("data/permutation/subprecodes.obsvals.RDS")
subpostcodes.obsvals<-readRDS("data/permutation/subpostcodes.obsvals.RDS")

### Make transition frequency matrix
domprecodes.tmats = lapply(domprecodes, function(x) obsmatrix(x,prob=T))
dompostcodes.tmats = lapply(dompostcodes, function(x) obsmatrix(x,prob=T))
subprecodes.tmats = lapply(subprecodes, function(x) obsmatrix(x,prob=T))
subpostcodes.tmats = lapply(subpostcodes, function(x) obsmatrix(x,prob=T))

#empty expanded df.
emptydf = as.data.frame(expand.grid(Var1=LETTERS[1:19],Var2=LETTERS[1:19]) %>% mutate(value=0))
emptydf <- zoo::coredata(emptydf)[rep(seq(nrow(emptydf)),21),]
emptydf$dyadid = rep(LETTERS[1:21],each=361)
emptydfx = emptydf[emptydf$dyadid!="R",]


emptydfx0 = emptydfx
emptydfx0[,1] <- tolower(emptydfx0[,1])
emptydfx0[,2] <- tolower(emptydfx0[,2])

emptydf0 = emptydf
emptydf0[,1] <- tolower(emptydf0[,1])
emptydf0[,2] <- tolower(emptydf0[,2])


dompre.trans<-rbind(do.call("rbind",Map(cbind, lapply(domprecodes.tmats, reshape2::melt), dyadid=LETTERS[1:21])),emptydf) %>%
  group_by(Var1,Var2,dyadid) %>% 
  summarize(value=sum(value,na.rm=T))
subpre.trans<-rbind(do.call("rbind",Map(cbind, lapply(subprecodes.tmats, reshape2::melt), dyadid=LETTERS[1:21])),emptydf0) %>%
  group_by(Var1,Var2,dyadid) %>% 
  summarize(value=sum(value))
dompost.trans<-rbind(do.call("rbind",Map(cbind, lapply(dompostcodes.tmats, reshape2::melt), dyadid=LETTERS[1:21][-18])),emptydfx) %>% 
  group_by(Var1,Var2,dyadid) %>% 
  summarize(value=sum(value))
subpost.trans<-rbind(do.call("rbind",Map(cbind, lapply(subpostcodes.tmats, reshape2::melt), dyadid=LETTERS[1:21][-18])),emptydfx0) %>% 
  group_by(Var1,Var2,dyadid) %>%
  summarize(value=sum(value))

#plot those that are p>.075 and show stars - edge thickness by transition value
dompre.trans.df = dompre.trans %>% group_by(Var1,Var2) %>% summarise(med = median(value)) %>% filter(med<0.075) %>% data.frame()
dompost.trans.df = dompost.trans %>% group_by(Var1,Var2) %>% summarise(med = median(value)) %>% filter(med>0.075) %>% data.frame()
subpre.trans.df = subpre.trans %>% group_by(Var1,Var2) %>% summarise(med = median(value)) %>% filter(med>0.075) %>% data.frame()
subpost.trans.df = subpost.trans %>% group_by(Var1,Var2) %>% summarise(med = median(value)) %>% filter(med>0.075) %>% data.frame()

#this will help plan the graph - 
g1=igraph::graph.edgelist(as.matrix(dompre.trans.df[,1:2]));plot(g1)
g2=igraph::graph.edgelist(as.matrix(dompost.trans.df[,1:2]));plot(g2)
g3=igraph::graph.edgelist(as.matrix(subpre.trans.df[,1:2]));plot(g3)
g4=igraph::graph.edgelist(as.matrix(subpost.trans.df[,1:2]));plot(g4)


#size of node by relative frequency of observed behavior (needs to be percentage of all behavior)
domprecodes.obsvals.df = do.call('rbind', Map(cbind,domprecodes.obsvals,dyadid=LETTERS[1:21])) %>% group_by(Var1,dyadid) %>% summarize(totalObserved = sum(Observed)) %>% ungroup() 
dompostcodes.obsvals.df = do.call('rbind', Map(cbind,dompostcodes.obsvals,dyadid=LETTERS[1:21][-18]))%>% group_by(Var1,dyadid) %>% summarize(totalObserved = sum(Observed)) %>% ungroup() 
subprecodes.obsvals.df = do.call('rbind', Map(cbind,subprecodes.obsvals,dyadid=LETTERS[1:21])) %>% group_by(Var1,dyadid) %>% summarize(totalObserved = sum(Observed)) %>% ungroup() 
subpostcodes.obsvals.df = do.call('rbind', Map(cbind,subpostcodes.obsvals,dyadid=LETTERS[1:21][-18])) %>% group_by(Var1,dyadid) %>% summarize(totalObserved = sum(Observed)) %>% ungroup() 

#add in emptydf
domprecodes.obsvals.df = rbind(domprecodes.obsvals.df, expand.grid(Var1=LETTERS[1:19], dyadid=LETTERS[1:21]) %>% mutate(totalObserved=0) %>% as.data.frame) %>% group_by(Var1,dyadid) %>% summarise(totalObserved=sum(totalObserved))
dompostcodes.obsvals.df = rbind(dompostcodes.obsvals.df, expand.grid(Var1=LETTERS[1:19], dyadid=LETTERS[1:21][-18]) %>% mutate(totalObserved=0) %>% as.data.frame) %>% group_by(Var1,dyadid) %>% summarise(totalObserved=sum(totalObserved))
subprecodes.obsvals.df = rbind(subprecodes.obsvals.df, expand.grid(Var1=LETTERS[1:19], dyadid=LETTERS[1:21]) %>% mutate(totalObserved=0) %>% as.data.frame) %>% group_by(Var1,dyadid) %>% summarise(totalObserved=sum(totalObserved))
subpostcodes.obsvals.df = rbind(subpostcodes.obsvals.df, expand.grid(Var1=LETTERS[1:19], dyadid=LETTERS[1:21][-18]) %>% mutate(totalObserved=0) %>% as.data.frame) %>% group_by(Var1,dyadid) %>% summarise(totalObserved=sum(totalObserved))

#get pct values
domprecodes.obsvals.df = domprecodes.obsvals.df  %>% group_by(dyadid) %>% arrange(dyadid,Var1) %>% mutate(totalObs = sum(totalObserved), pctObs=(100*totalObserved)/totalObs)
dompostcodes.obsvals.df = dompostcodes.obsvals.df  %>% group_by(dyadid) %>% arrange(dyadid,Var1) %>% mutate(totalObs = sum(totalObserved), pctObs=(100*totalObserved)/totalObs)
subprecodes.obsvals.df = subprecodes.obsvals.df  %>% group_by(dyadid) %>% arrange(dyadid,Var1) %>% mutate(totalObs = sum(totalObserved), pctObs=(100*totalObserved)/totalObs)
subpostcodes.obsvals.df = subpostcodes.obsvals.df  %>% group_by(dyadid) %>% arrange(dyadid,Var1) %>% mutate(totalObs = sum(totalObserved), pctObs=(100*totalObserved)/totalObs)

#get median
domprefreq = domprecodes.obsvals.df %>% group_by(Var1) %>% summarise(medpctObs = median(pctObs))
dompostfreq = dompostcodes.obsvals.df %>% group_by(Var1) %>% summarise(medpctObs = median(pctObs))
subprefreq = subprecodes.obsvals.df %>% group_by(Var1) %>% summarise(medpctObs = median(pctObs))
subpostfreq = subpostcodes.obsvals.df %>% group_by(Var1) %>% summarise(medpctObs = median(pctObs))

#nodesdf
domprefreq$medpctObslog =  round(5*log(domprefreq$medpctObs+1),0)
dompostfreq$medpctObslog =  round(5*log(dompostfreq$medpctObs+1),0)
subprefreq$medpctObslog =  round(5*log(subprefreq$medpctObs+1),0)
subpostfreq$medpctObslog =  round(5*log(subpostfreq$medpctObs+1),0)

#range for font sizes  0 - 17, 12px-62px
max(c(domprefreq$medpctObslog, dompostfreq$medpctObslog,subprefreq$medpctObslog,subpostfreq$medpctObslog))

fontsizes=data.frame(values=c(domprefreq$medpctObslog, dompostfreq$medpctObslog,subprefreq$medpctObslog,subpostfreq$medpctObslog),
                     fontsize = round(plotrix::rescale(c(domprefreq$medpctObslog, dompostfreq$medpctObslog,subprefreq$medpctObslog,subpostfreq$medpctObslog), c(12,62)),0)
) %>% unique %>% arrange(values)

domprefreq$fontsize <- fontsizes$fontsize[match(domprefreq$medpctObslog,fontsizes$values)]
dompostfreq$fontsize <- fontsizes$fontsize[match(dompostfreq$medpctObslog,fontsizes$values)]
subprefreq$fontsize <- fontsizes$fontsize[match(subprefreq$medpctObslog,fontsizes$values)]
subpostfreq$fontsize <- fontsizes$fontsize[match(subpostfreq$medpctObslog,fontsizes$values)]

dompre.trans.df$edgesize = dompre.trans.df$med*5
dompost.trans.df$edgesize = dompost.trans.df$med*5
subpre.trans.df$edgesize = subpre.trans.df$med*5
subpost.trans.df$edgesize = subpost.trans.df$med*5



### Supplementary table S2. 

dompre<-left_join(domprepvalsh,dompre.trans %>% group_by(Var1,Var2) %>% summarise(med = median(value))) %>% 
  mutate(sig=ifelse(higherp<0.001,paste(higherp,"***",sep=""),
                    ifelse(higherp<0.01,paste(higherp,"**",sep=""),
                           ifelse(higherp<0.05,paste(higherp,"*",sep=""),paste(higherp,"",sep=""))))) %>% 
  mutate(dompre=paste(format(round(med,3),nsmall=3),sig)) %>% 
  select(Var1,Var2,dompre) %>% 
  mutate(behavA=Var1) %>%  left_join(.,behavcodes) %>% mutate(Var1=behav) %>% select(Var1,Var2,dompre) %>% 
  mutate(behavA=Var2) %>% left_join(.,behavcodes) %>% mutate(Var2=behav) %>% select(Var1,Var2,dompre) %>% 
  mutate(contingency=paste(Var1,Var2,sep="-->")) %>% select(contingency,dompre)

subpre<-left_join(subprepvalsh,
                  subpre.trans %>% ungroup %>% mutate(Var1=toupper(Var1),Var2=toupper(Var2)) %>% 
                    group_by(Var1,Var2) %>% summarise(med = median(value))) %>% 
  mutate(sig=ifelse(higherp<0.001,paste(higherp,"***",sep=""),
                    ifelse(higherp<0.01,paste(higherp,"**",sep=""),
                           ifelse(higherp<0.05,paste(higherp,"*",sep=""),paste(higherp,"",sep=""))))) %>% 
  mutate(subpre=paste(format(round(med,3),nsmall=3),sig)) %>% 
  select(Var1,Var2,subpre) %>% 
  mutate(behavA=Var1) %>%  left_join(.,behavcodes) %>% mutate(Var1=behav) %>% select(Var1,Var2,subpre) %>% 
  mutate(behavA=Var2) %>% left_join(.,behavcodes) %>% mutate(Var2=behav) %>% select(Var1,Var2,subpre) %>% 
  mutate(contingency=paste(Var1,Var2,sep="-->")) %>% select(contingency,subpre)

dompost<-left_join(dompostpvalsh,dompost.trans %>% group_by(Var1,Var2) %>% summarise(med = median(value))) %>% 
  mutate(sig=ifelse(higherp<0.001,paste(higherp,"***",sep=""),
                    ifelse(higherp<0.01,paste(higherp,"**",sep=""),
                           ifelse(higherp<0.05,paste(higherp,"*",sep=""),paste(higherp,"",sep=""))))) %>% 
  mutate(dompost=paste(format(round(med,3),nsmall=3),sig)) %>% 
  select(Var1,Var2,dompost) %>% 
  mutate(behavA=Var1) %>%  left_join(.,behavcodes) %>% mutate(Var1=behav) %>% select(Var1,Var2,dompost) %>% 
  mutate(behavA=Var2) %>% left_join(.,behavcodes) %>% mutate(Var2=behav) %>% select(Var1,Var2,dompost) %>% 
  mutate(contingency=paste(Var1,Var2,sep="-->")) %>% select(contingency,dompost)

subpost<-left_join(subpostpvalsh,
                   subpost.trans %>% ungroup %>% 
                     mutate(Var1=toupper(Var1),Var2=toupper(Var2)) %>% 
                     group_by(Var1,Var2) %>% summarise(med = median(value)) ) %>% 
  mutate(sig=ifelse(higherp<0.001,paste(higherp,"***",sep=""),
                    ifelse(higherp<0.01,paste(higherp,"**",sep=""),
                           ifelse(higherp<0.05,paste(higherp,"*",sep=""),paste(higherp,"",sep=""))))) %>% 
  mutate(subpost=paste(format(round(med,3),nsmall=3),sig)) %>% 
  select(Var1,Var2,subpost) %>% 
  mutate(behavA=Var1) %>%  left_join(.,behavcodes) %>% mutate(Var1=behav) %>% select(Var1,Var2,subpost) %>% 
  mutate(behavA=Var2) %>% left_join(.,behavcodes) %>% mutate(Var2=behav) %>% select(Var1,Var2,subpost) %>% 
  mutate(contingency=paste(Var1,Var2,sep="-->")) %>% select(contingency,subpost)


suppl_table_s2<-full_join(dompre,subpre) %>% 
  full_join(.,dompost) %>% 
  full_join(.,subpost)

