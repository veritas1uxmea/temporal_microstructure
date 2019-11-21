
## Figure 2. Exemplar temporal pattern of aggressive and subrodinate behaviors by the eventual dominant and subordinate males on days  1-5 (Dyad C and Dyad G)
tmp_c=dfburstA %>% filter(dyadid=="C") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(time)) %>% filter(!is.na(burstval)) %>% data.frame() 

tmp_g=dfburstA %>% filter(dyadid=="G") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(time)) %>% filter(!is.na(burstval)) %>% data.frame() 

cc<-plotburstsX(df, tmp_c, ID="C",linebreak=T,legend=F,rel_size = 1.5)+ggtitle("Dyad C")+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))

gg<-plotburstsX(df, tmp_g, ID="G",linebreak=T,legend=F,rel_size = 1.5)+ggtitle("Dyad G")+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))


fig2 <- arrangeGrob(cc,gg,ncol=1)
#ggsave("figures/fig2_2.tiff",fig2,height=16,width=12)


## Suppl Figure S2. Temporal pattern of all dyads

tmp_a=dfburstA %>% filter(dyadid=="A") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_b=dfburstA %>% filter(dyadid=="B") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_c=dfburstA %>% filter(dyadid=="C") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_d=dfburstA %>% filter(dyadid=="D") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_e=dfburstA %>% filter(dyadid=="E") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_f=dfburstA %>% filter(dyadid=="F") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_g=dfburstA %>% filter(dyadid=="G") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_h=dfburstA %>% filter(dyadid=="H") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_i=dfburstA %>% filter(dyadid=="I") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_j=dfburstA %>% filter(dyadid=="J") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_k=dfburstA %>% filter(dyadid=="K") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_l=dfburstA %>% filter(dyadid=="L") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_m=dfburstA %>% filter(dyadid=="M") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_n=dfburstA %>% filter(dyadid=="N") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_o=dfburstA %>% filter(dyadid=="O") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_p=dfburstA %>% filter(dyadid=="P") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_q=dfburstA %>% filter(dyadid=="Q") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_r=dfburstA %>% filter(dyadid=="R") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_s=dfburstA %>% filter(dyadid=="S") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_t=dfburstA %>% filter(dyadid=="T") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 
tmp_u=dfburstA %>% filter(dyadid=="U") %>% group_by(day,burstval) %>% summarise(mint = min(time), maxt = max(end)) %>% filter(!is.na(burstval)) %>% data.frame() 



AA<-plotburstsX(df, tmp_a, ID="A",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad A")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
BB<-plotburstsX(df, tmp_b, ID="B",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad B")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
CC<-plotburstsX(df, tmp_c, ID="C",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad C")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
DD<-plotburstsX(df, tmp_d, ID="D",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad D")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
EE<-plotburstsX(df, tmp_e, ID="E",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad E")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
FF<-plotburstsX(df, tmp_f, ID="F",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad F")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
GG<-plotburstsX(df, tmp_g, ID="G",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad G")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
HH<-plotburstsX(df, tmp_h, ID="H",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad H")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
II<-plotburstsX(df, tmp_i, ID="I",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad I")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
JJ<-plotburstsX(df, tmp_j, ID="J",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad J")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
KK<-plotburstsX(df, tmp_k, ID="K",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad K")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
LL<-plotburstsX(df, tmp_l, ID="L",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad L")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
MM<-plotburstsX(df, tmp_m, ID="M",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad M")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
NN<-plotburstsX(df, tmp_n, ID="N",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad N")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
OO<-plotburstsX(df, tmp_o, ID="O",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad O")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
PP<-plotburstsX(df, tmp_p, ID="P",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad P")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
QQ<-plotburstsX(df, tmp_q, ID="Q",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad Q")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
RR<-plotburstsX(df, tmp_r, ID="R",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad R")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
SS<-plotburstsX(df, tmp_s, ID="S",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad S")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
TT<-plotburstsX(df, tmp_t, ID="T",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad T")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))
UU<-plotburstsX(df, tmp_u, ID="U",linebreak=T,legend=F,rel_size = 1.5)+
  ggtitle("Dyad U")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)))

sfig2x1 <- arrangeGrob(AA,BB,CC,DD,EE,FF,ncol=3)
sfig2x2 <- arrangeGrob(GG,HH,II,JJ,KK,LL,ncol=3)
sfig2x3 <- arrangeGrob(MM,NN,OO,PP,QQ,RR,ncol=3)
sfig2x4 <- arrangeGrob(SS,TT,UU,ncol=3)


# ggsave("img/sfig2x1.png",sfig2x1,height=450,width=1000,units = "mm")
# ggsave("img/sfig2x2.png",sfig2x2,height=450,width=1000,units = "mm")
# ggsave("img/sfig2x3.png",sfig2x3,height=450,width=1000,units = "mm")
# ggsave("img/sfig2x4.png",sfig2x4,height=225.5,width=1000,units = "mm")
