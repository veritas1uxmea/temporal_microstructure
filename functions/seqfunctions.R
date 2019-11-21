
collapse <-function(x){paste0(x,collapse="")}

uncollapse <-function(x,n=1){
  substring(x, seq(1,nchar(x),n), seq(n,nchar(x),n))
}

seq_droprpt <- function(x, coll=F) {
  if(coll==F){
    if(length(x) == 1) {return(x)}
    
    x <- as.vector(x)
    v <- c(x[1])
    for(i in 2:length(x)) {
      if(x[i] != v[length(v)]) {
        v <- c(v, x[i])
      }
    }
    return(v)
  }
  
  else
    if(coll==T){
      x <- uncollapse(x)
      
      if(length(x) == 1) {return(x)}
      
      x <- as.vector(x)
      v <- c(x[1])
      for(i in 2:length(x)) {
        if(x[i] != v[length(v)]) {
          v <- c(v, x[i])
        }
      }
      
      v1<-collapse(v)
      return(v1)
    }
}








sem <- function(x) sd(x, na.rm=T)/sqrt(length(x))

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

friedposthoc_sig <- function(obj)  {
  x<- obj$PostHoc.Test %>% as.data.frame %>% mutate(id = rownames(.)) %>% filter(V1<.05)
  colnames(x)[1]<-"pval"
  return(x)
}


getbursts <- function(df,ID, DAY=1, LEVEL=2, gamma=1){
  df1 <- df %>% filter(dyadid==ID) %>% filter(state=="sub"|state=="aggr") %>% arrange(day,time) %>% 
    select(dyadid,animalid,day,state,time,duration)
  df1empty <- data.frame(dyadid=ID, animalid=rep(c("S", "T"),each=10), day=1:5, state=rep(c("sub", "aggr"),10), time=NA, duration=NA)
  df2 <- rbind(df1,df1empty)
  
  df2time <- df2 %>% filter(day==DAY) %>% .$time
  df2timex <- ifelse(df2time == lag(df2time), df2time+0.0001, df2time) #don't exclude times, just add 0.0001
  df2timex1 <- c(df2time[1],df2timex[-1]) #make sure first time is included
  df2timex1 <- df2timex1[complete.cases(df2timex1)]
  
  if(length(df2timex1)==1) {
    df2timex1 <- NA #can't have bursts if only one time in a vector
    
  } 
  
  if(all(is.na(df2time))) {
    
    df2bursts.summary <-  data.frame(level=NA, start=NA, end=NA)
    warning("no events - empty dataframe produced")  
    
  } else
    
  {
    df2bursts <- bursts::kleinberg(df2timex1, gamma=gamma)
    df2bursts.summary <- df2bursts #%>% filter(level==LEVEL)
  }
  
  return(df2bursts.summary)
}


##########

plotbursts <- function(df,ID,legend=T,wd=NULL,linebreak=T){
  
  df1 <- df %>% filter(dyadid==ID) %>% filter(state=="sub"|state=="aggr") %>% arrange(day,time) %>% 
    mutate(ymin = ifelse(idAB=="B", 1, 2), ymax = ifelse(idAB=="A", 2, 3))
  
  #add empty df for days
  df1empty <- data.frame(dyadid=ID, animalid=rep(c("S", "T"),each=10), id=NA, day=1:5, 
                         behavior=rep(c('Flee','Bite','Defensive freeze','Tailrattle','Subordinate posture', 'Lunge'),each=10)
                         , state=rep(c("sub", "aggr"),10), time=NA, duration=NA,
                         idAB=rep(c("A", "B")), ymin=NA, ymax=NA)
  
  df1a<-df1[,colnames(df1empty)]
  
  df1x <- rbind(df1a,df1empty)
  df1x$behavior <- factor(df1x$behavior, levels=c('Flee','Subordinate posture','Defensive freeze','Tailrattle','Bite', 'Lunge'))
  
  
  P=ggplot(df1x, aes(x=time,y=ymin, color=behavior, fill=behavior))  + 
    scale_color_manual(values=c("cyan", "blue3", "dodgerblue", "firebrick", "red1", "indianred")) +
    scale_fill_manual(values=c("cyan", "blue3", "dodgerblue", "firebrick", "red1", "indianred")) +
    facet_wrap(~day, ncol=1) +
    theme_bw() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    ) + xlim(0,1200)+
    xlab("Time (s)") + ylab("")
  
  P=P+ annotate("text", label = "D", size = 4, x = -Inf, y = 2.1,hjust=-.45)+
       annotate("text", label = "S", size = 4, x = -Inf, y = 1,hjust=-.5)
  
  if(linebreak==T){ P = P + geom_hline(yintercept=1.5, lty=2) }
  
  
  if(is.null(wd) & legend==T){ P = P + geom_tile(aes(width=duration)) }  
  
  if(is.null(wd) & legend!=T){ P = P + geom_tile(aes(width=duration)) + theme(legend.position = "none") }
  
  if(!is.null(wd) & legend==T){ P = P + geom_tile(aes(width=wd))  }
  
  if(!is.null(wd) & legend!=T){ P = P + geom_tile(aes(width=wd)) + theme(legend.position = "none") }
  
  
  return(P)
}



####

plotburstsX <- function(df,linedf, ID,legend=T,wd=NULL,linebreak=T,rel_size=2){
  
  df1 <- df %>% filter(dyadid==ID) %>% filter(state=="sub"|state=="aggr") %>% arrange(day,time) %>% 
    mutate(ymin = ifelse(idAB=="B", 1, 2), ymax = ifelse(idAB=="A", 2, 3))
  
  #add empty df for days
  df1empty <- data.frame(dyadid=ID, animalid=rep(c("S", "T"),each=10), id=NA, day=1:5, 
                         behavior=rep(c('Flee','Bite','Defensive freeze','Tailrattle','Subordinate posture', 'Lunge'),each=10)
                         , state=rep(c("sub", "aggr"),10), time=NA, duration=NA,
                         idAB=rep(c("A", "B")), ymin=NA, ymax=NA)
  
  df1a<-df1[,colnames(df1empty)]
  
  df1x <- rbind(df1a,df1empty)
  
  df1x$behavior <- factor(df1x$behavior, levels=c('Flee','Subordinate posture','Defensive freeze','Tailrattle','Bite', 'Lunge'))
  
  df1x$mint=NA
  df1x$maxt=NA
  
  linedf = data.frame(dyadid=NA,animalid=NA,id=NA,day=linedf$day,behavior=NA,state=NA,time=NA,duration=NA,idAB=NA,ymin=NA, ymax=NA, mint=linedf$mint,maxt=linedf$maxt)
  
  df1x = rbind(df1x,linedf)
  
  P=ggplot(df1x, aes(x=time,y=ymin, color=behavior, fill=behavior))  + 
    scale_color_manual(values=c("cyan", "blue3", "dodgerblue", "firebrick", "red1", "indianred")) +
    scale_fill_manual(values=c("cyan", "blue3", "dodgerblue", "firebrick", "red1", "indianred")) +
    #geom_vline(aes(xintercept = mint), lty=1,lwd=1)+
    #geom_vline(aes(xintercept = maxt), lty=1,lwd=1)+
    facet_wrap(~day, ncol=1) +
    #theme_bw() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.text.y=element_blank(),
          axis.title.x = element_text(size=rel(rel_size)),
          axis.text.x=element_text(size=rel(rel_size)),
          strip.text = element_text(size=rel(rel_size)),
          axis.ticks.y=element_blank(),
          axis.line.x=element_line(color="white"),
          axis.line.y=element_line(color="white")
    ) + scale_x_continuous(limits = c(0, 1250))+
    xlab("Time (s)") + ylab("") +
    geom_rect(aes(xmin=mint, xmax=maxt,  ymin=2.7, ymax=3),fill = 'gold')
  
  P=P+ annotate("text", label = "D", size = rel_size*3.5, x = -Inf, y = 2.1,hjust=-.45)+
    annotate("text", label = "S", size = rel_size*3.5, x = -Inf, y = 1,hjust=-.5)
  
  if(linebreak==T){ P = P + geom_hline(yintercept=1.5, lty=2) +
    geom_hline(yintercept=2.6, lty=1)+
    geom_hline(yintercept=0.4, lty=1)
  }
  
  
  if(is.null(wd) & legend==T){ P = P + geom_tile(aes(width=duration)) }  
  
  if(is.null(wd) & legend!=T){ P = P + geom_tile(aes(width=duration)) + theme(legend.position = "none") }
  
  if(!is.null(wd) & legend==T){ P = P + geom_tile(aes(width=wd))  }
  
  if(!is.null(wd) & legend!=T){ P = P + geom_tile(aes(width=wd)) + theme(legend.position = "none") }
  
  
  return(P)
}




####

phicoeff = function(mat){
  phi =  ((mat[1]*mat[4])-(mat[2]*mat[3])) / sqrt(rowSums(mat)[1]*rowSums(mat)[2]*colSums(mat)[1]*colSums(mat)[2])
  N= sum(mat)
  chisq = N*(phi^2)
  pval=pchisq(chisq,df=1,lower.tail=FALSE) 
  return(list(matrix=mat, phi=phi, pval=pval, N=N))
}

#######

### Get Bursts df.

getburstsdf <- function(df, GM=0.3, LEVEL=2){
  
  dyadids <- unique(df$dyadid) #Get Dyad Ids as vector
  
  burstdfs<-vector("list", length(dyadids))
  for(i in seq_along(dyadids)){
    
    burstdfs[[i]] <- rbind(
      getbursts(df, ID=dyadids[i], DAY=1, gamma=GM) %>% mutate(day=1),
      getbursts(df, ID=dyadids[i], DAY=2, gamma=GM) %>% mutate(day=2),
      getbursts(df, ID=dyadids[i], DAY=3, gamma=GM) %>% mutate(day=3),
      getbursts(df, ID=dyadids[i], DAY=4, gamma=GM) %>% mutate(day=4),
      getbursts(df, ID=dyadids[i], DAY=5, gamma=GM) %>% mutate(day=5)
    )
  }
  
  names(burstdfs)<-dyadids
  
  burstdfs.A<- do.call('rbind', Map(cbind, burstdfs, dyadid=dyadids))
  burstdfs.A2 <- burstdfs.A %>% filter(level==LEVEL) %>% group_by(dyadid) %>% mutate(burstno = row_number())
  
  #burstdfs.A2sum <- burstdfs.A2 %>% group_by(dyadid,day) %>% summarize(total=n(), duration=mean(end-start)) %>% ungroup()
  
  return(burstdfs.A2)
}


#####

getrawburstsdf <- function(df, GM=0.3, LEVEL=2){
  
  dyadids <- unique(df$dyadid)
  burstdf <- getburstsdf(df,GM=GM,LEVEL=LEVEL)
  ##
  library(data.table)
  
  results.day1<-vector("list", length(dyadids))
  for(i in seq_along(dyadids)){
    tmpdf <- df %>% filter(dyadid==dyadids[i] & day==1 & (state=="aggr"|state=="sub")) %>% arrange(day,time)
    tmp <- burstdf %>% filter(dyadid==dyadids[i] & day==1 & level==2) %>% mutate(burstno = row_number())
    setDT(tmp)
    tmpdf$burstval <- tmp[.(start = tmpdf$time), on="start", roll=Inf][start > end, burstno := NA_integer_]$burstno  
    results.day1[[i]] <- tmpdf
  }
  ##
  results.day2<-vector("list", length(dyadids))
  for(i in seq_along(dyadids)){
    tmpdf <- df %>% filter(dyadid==dyadids[i] & day==2 & (state=="aggr"|state=="sub")) %>% arrange(day,time)
    tmp <- burstdf %>% filter(dyadid==dyadids[i] & day==2 & level==2) %>% mutate(burstno = row_number())
    setDT(tmp)
    tmpdf$burstval <- tmp[.(start = tmpdf$time), on="start", roll=Inf][start > end, burstno := NA_integer_]$burstno  
    results.day2[[i]] <- tmpdf
  }
  ##
  results.day3<-vector("list", length(dyadids))
  for(i in seq_along(dyadids)){
    tmpdf <- df %>% filter(dyadid==dyadids[i] & day==3 & (state=="aggr"|state=="sub")) %>% arrange(day,time)
    tmp <- burstdf %>% filter(dyadid==dyadids[i] & day==3 & level==2) %>% mutate(burstno = row_number())
    setDT(tmp)
    tmpdf$burstval <- tmp[.(start = tmpdf$time), on="start", roll=Inf][start > end, burstno := NA_integer_]$burstno  
    results.day3[[i]] <- tmpdf
  }
  ##
  results.day4<-vector("list", length(dyadids))
  for(i in seq_along(dyadids)){
    tmpdf <- df %>% filter(dyadid==dyadids[i] & day==4 & (state=="aggr"|state=="sub")) %>% arrange(day,time)
    tmp <- burstdf %>% filter(dyadid==dyadids[i] & day==4 & level==2) %>% mutate(burstno = row_number())
    setDT(tmp)
    tmpdf$burstval <- tmp[.(start = tmpdf$time), on="start", roll=Inf][start > end, burstno := NA_integer_]$burstno  
    results.day4[[i]] <- tmpdf
  }
  ##
  results.day5<-vector("list", length(dyadids))
  for(i in seq_along(dyadids)){
    tmpdf <- df %>% filter(dyadid==dyadids[i] & day==5 & (state=="aggr"|state=="sub")) %>% arrange(day,time)
    tmp <- burstdf %>% filter(dyadid==dyadids[i] & day==5 & level==2) %>% mutate(burstno = row_number())
    setDT(tmp)
    tmpdf$burstval <- tmp[.(start = tmpdf$time), on="start", roll=Inf][start > end, burstno := NA_integer_]$burstno  
    results.day5[[i]] <- tmpdf
  }
  ##
  names(results.day1)<-names(results.day2)<-names(results.day3)<-names(results.day4)<-names(results.day5)<-dyadids
  ##
  
  
  lowbursts<-
    rbind(
      do.call('rbind', results.day1),
      do.call('rbind', results.day2),
      do.call('rbind', results.day3),
      do.call('rbind', results.day4),
      do.call('rbind', results.day5)
    ) %>% as.data.frame()
  
  return(lowbursts)
}


####

addburstNAs <- function(dfburst){
  
  iter = max(dfburst$burstval,na.rm=T)
  
  if(iter > 26) stop("Max bursts for any dyad on any day is 26")
  
  lowbursts <- dfburst %>%
    ungroup() %>%
    arrange(dyadid,day,time) %>% 
    group_by(dyadid) %>% 
    mutate(burstvalcume = as.numeric(factor(ifelse(is.na(burstval), NA, 
                                                   as.numeric(factor(paste(day,LETTERS[burstval])))))) #using letters to unique id is ok here as max burstval =25         
    )
  
  #will remove ids with no bursts
  ids0 <- names(which(table(dfburst$dyadid,dfburst$burstval)[,1]==0) )
  
  lowbursts1 <- lowbursts %>% 
    arrange(dyadid,day,time) %>% 
    group_by(dyadid) %>% 
    filter(!dyadid %in% ids0) %>%
    mutate(burstvalNA = zoo::na.locf(zoo::na.locf(burstval, fromLast=T, na.rm=F)),
           burstvalcumeNA = zoo::na.locf(zoo::na.locf(burstvalcume, fromLast=T, na.rm=F)) ) %>% as.data.frame()
  
  return(lowbursts1)
}



#######


getburstfreq <- function(dfburst0){
  
  #1.  excluding NAs - need to make sure 0s are included for total counts.
  
  lowbursts1 <- dfburst0 %>% filter(!is.na(burstvalcume)) %>% ungroup() %>% as.data.frame()
  dyadidsx <- dfburst0$dyadid %>% unique
  
  
  lowbursts.totals <- rbind(
    lowbursts1 %>% as.data.frame()  %>% group_by(dyadid,burstvalcume,id,state) %>% summarise(total=n())
    ,  lowbursts1 %>% select(dyadid,burstvalcume)  %>% unique %>% mutate(id=paste0(dyadid,"S"), state="aggr", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcume)  %>% unique %>% mutate(id=paste0(dyadid,"S"), state="sub", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcume)  %>% unique %>% mutate(id=paste0(dyadid,"T"), state="aggr", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcume)  %>% unique %>% mutate(id=paste0(dyadid,"T"), state="sub", total=0)
  ) 
  
  lowbursts.totals <- lowbursts.totals %>% as.data.frame() %>%   group_by(dyadid, burstvalcume,id,state) %>%  summarise(total=sum(total)) 
  
  # get final phi coefficient for each [get final two]
  cmat_final=NULL
  for(i in seq_along(dyadidsx)){
    cmat_final[[i]]=lowbursts.totals %>% filter(dyadid==dyadidsx[[i]]) %>% ungroup() %>% 
      data.frame %>%
      filter(burstvalcume>=(max(burstvalcume)-1)) %>%
      group_by(id,state)    %>%
      summarise(total=sum(total)) %>%
      .$total %>% matrix(.,2,2,byrow=T) %>% phicoeff %>% .$phi %>% sign
    
  }
  
  dfsign<-data.frame(dyadid=dyadidsx, sign = cmat_final) %>% 
    mutate(A = ifelse(sign>0, paste0(dyadid,"S"), paste0(dyadid,"T")),
           B = ifelse(sign<0, paste0(dyadid,"S"), paste0(dyadid,"T")))
  
  dfsign0<-data.frame(id=c(dfsign$A,dfsign$B), idAB = rep(LETTERS[1:2],each=nrow(dfsign)))
  
  lowbursts.totalsZ <- lowbursts.totals %>% left_join(dfsign0) 
  
  return(lowbursts.totalsZ)
}

########


getburstfreq0 <- function(dfburst0){
  
  #1.  excluding NAs - need to make sure 0s are included for total counts.
  
  lowbursts1 <- dfburst0 %>% filter(!is.na(burstvalcume)) %>% ungroup() %>% as.data.frame()
  dyadidsx <- dfburst0$dyadid %>% unique
  
  
  lowbursts.totals <- rbind(
    lowbursts1 %>% as.data.frame()  %>% group_by(dyadid,burstvalcume,idAB,state) %>% summarise(total=n())%>% as.data.frame()
    ,  lowbursts1 %>% select(dyadid,burstvalcume)  %>% unique %>% mutate(idAB="A", state="aggr", total=0) %>% as.data.frame()
    ,  lowbursts1 %>% select(dyadid,burstvalcume)  %>% unique %>% mutate(idAB="A", state="sub", total=0)%>% as.data.frame()
    ,  lowbursts1 %>% select(dyadid,burstvalcume)  %>% unique %>% mutate(idAB="B", state="aggr", total=0)%>% as.data.frame()
    ,  lowbursts1 %>% select(dyadid,burstvalcume)  %>% unique %>% mutate(idAB="B", state="sub", total=0)%>% as.data.frame()
  ) 
  
  
  lowbursts.totals <- lowbursts.totals %>% as.data.frame() %>%   group_by(dyadid, burstvalcume,idAB,state) %>%  summarise(total=sum(total)) 
  
    return(lowbursts.totals)
}


#########

getphidf <- function(df,PVAL=.05){
  
  dfburst.sum <- df
  
  dyadidsx <- dfburst.sum$dyadid %>% unique
  
  resphiZ=NULL
  for(i in seq_along(dyadidsx)){
    resphiZ[[i]]=dfburst.sum %>% ungroup() %>% arrange(dyadid,burstvalcume,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcume) %>%
      map(~ psych::phi(matrix(.$total, 2, 2, byrow=T))) %>%    unlist
  }
  
  names(resphiZ)=dyadidsx
  
  resphiZdf = resphiZ %>%   
    map(~ data.frame(val=., burstno=names(.))) %>% 
    Map(cbind, ., dyadid=dyadidsx) %>%
    do.call('rbind', .)
  
  resphiZdf$burstno = as.numeric(as.character(resphiZdf$burstno))
  
  ## adding in p-value
  resphiZp=NULL
  for(i in seq_along(dyadidsx)){
    resphiZp[[i]]=dfburst.sum %>% ungroup() %>% arrange(dyadid,burstvalcume,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcume) %>%
      map(~ phicoeff(matrix(.$total, 2, 2, byrow=T))$pval) %>%    unlist
  }
  
  names(resphiZp)=dyadidsx
  
  resphiZdfp = resphiZp %>% 
    map(~ data.frame(pval=., burstno=names(.))) %>% 
    Map(cbind, ., dyadid=dyadidsx) %>%
    do.call('rbind', .)
  
  resphiZdfp$burstno = as.numeric(as.character(resphiZdfp$burstno))
  
  resphiZdfx = resphiZdf %>% left_join(resphiZdfp)
  resphiZdfx$pvalx = ifelse(resphiZdfx$pval<PVAL, "sig", "NS")
  
  return(resphiZdfx)
  
}



######

getfisherdf <- function(df){
  
  dfburst.sum <- df
  
  dyadidsx <- dfburst.sum$dyadid %>% unique
  
  resphiZ=NULL
  for(i in seq_along(dyadidsx)){
    resphiZ[[i]]=dfburst.sum %>% ungroup() %>% arrange(dyadid,burstvalcume,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcume) %>%
      map(~ fisher.test(matrix(.$total, 2, 2, byrow=T))[[3]][[1]]) %>%    unlist
  }
  
  names(resphiZ)=dyadidsx
  
  resphiZdf = resphiZ %>%   
    map(~ data.frame(val=., burstno=names(.))) %>% 
    Map(cbind, ., dyadid=dyadidsx) %>%
    do.call('rbind', .)
  
  resphiZdf$burstno = as.numeric(as.character(resphiZdf$burstno))
  
  ## adding in p-value
  resphiZp=NULL
  for(i in seq_along(dyadidsx)){
    resphiZp[[i]]=dfburst.sum %>% ungroup() %>% arrange(dyadid,burstvalcume,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcume) %>%
      map(~ fisher.test(matrix(.$total, 2, 2, byrow=T))[[1]]) %>%    unlist
  }
  
  names(resphiZp)=dyadidsx
  
  resphiZdfp = resphiZp %>% 
    map(~ data.frame(pval=., burstno=names(.))) %>% 
    Map(cbind, ., dyadid=dyadidsx) %>%
    do.call('rbind', .)
  
  resphiZdfp$burstno = as.numeric(as.character(resphiZdfp$burstno))
  
  resphiZdfx = resphiZdf %>% left_join(resphiZdfp)
  resphiZdfx$pvalx = ifelse(resphiZdfx$pval<0.05, "sig", "NS")
  colnames(resphiZdfx)<-c("Fval", "burstno", "dyadid", "Fpval", "Fpvalx")
  return(resphiZdfx)
  
}



######

getburstfreqNA <- function(dfburst0){
  
  #1.  excluding NAs - need to make sure 0s are included for total counts.
  
  lowbursts1 <- dfburst0  %>% ungroup() %>% as.data.frame()
  dyadidsx <- dfburst0$dyadid %>% unique
  
  
  ## 2. adding in NAs - need to make sure 0s are included for total counts.
  
  lowburstsNA.totals <- rbind(
    lowbursts1  %>% group_by(dyadid,burstvalcumeNA,id,state) %>% summarize(total=n())
    ,  lowbursts1 %>% select(dyadid,burstvalcumeNA) %>% unique %>% mutate(id=paste0(dyadid,"S"), state="aggr", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcumeNA) %>% unique %>% mutate(id=paste0(dyadid,"S"), state="sub", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcumeNA) %>% unique %>% mutate(id=paste0(dyadid,"T"), state="aggr", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcumeNA) %>% unique %>% mutate(id=paste0(dyadid,"T"), state="sub", total=0)
  ) %>% 
    group_by(dyadid, burstvalcumeNA,id,state) %>%
    summarize(total=sum(total)) 
  
  
  
  # get final phi coefficient for each [get final two]
  cmat_finalNA=NULL
  for(i in seq_along(dyadidsx)){
    cmat_finalNA[[i]]=lowburstsNA.totals %>% filter(dyadid==dyadidsx[[i]]) %>% ungroup() %>% 
      data.frame %>%
      filter(burstvalcumeNA>=(max(burstvalcumeNA)-1)) %>%
      group_by(id,state)    %>%
      summarize(total=sum(total)) %>%
      .$total %>% matrix(.,2,2,byrow=T) %>% phicoeff %>% .$phi %>% sign
    
  }
  
  dfsignNA<-data.frame(dyadid=dyadidsx, sign = cmat_finalNA) %>% 
    mutate(A = ifelse(sign>0, paste0(dyadid,"S"), paste0(dyadid,"T")),
           B = ifelse(sign<0, paste0(dyadid,"S"), paste0(dyadid,"T")))
  
  dfsignNA0<-data.frame(id=c(dfsignNA$A,dfsignNA$B), idAB = rep(LETTERS[1:2],each=nrow(dfsignNA)))
  
  
  lowburstsNA.totalsZ <- lowburstsNA.totals %>% left_join(dfsignNA0) 
  
  return(lowburstsNA.totalsZ)
}


#####
getburstfreqNA0 <- function(dfburst0){
  
  #1.  excluding NAs - need to make sure 0s are included for total counts.
  
  lowbursts1 <- dfburst0  %>% ungroup() %>% as.data.frame()
  dyadidsx <- dfburst0$dyadid %>% unique
  
  
  ## 2. adding in NAs - need to make sure 0s are included for total counts.
  
  lowburstsNA.totals <- rbind(
    lowbursts1  %>% group_by(dyadid,burstvalcumeNA,idAB,state) %>% summarize(total=n())
    ,  lowbursts1 %>% select(dyadid,burstvalcumeNA) %>% unique %>% mutate(idAB="A", state="aggr", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcumeNA) %>% unique %>% mutate(idAB="B", state="aggr", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcumeNA) %>% unique %>% mutate(idAB="A", state="sub", total=0)
    ,  lowbursts1 %>% select(dyadid,burstvalcumeNA) %>% unique %>% mutate(idAB="B", state="sub", total=0)
  ) %>% 
    group_by(dyadid, burstvalcumeNA,idAB,state) %>%
    summarise(total=sum(total)) 
  
  return(lowburstsNA.totals)
}





######


getphidfNA <- function(df){
  
  dfZ <- df
  
  dyadidsx <- dfburst.sum$dyadid %>% unique
  
  
  resNAphiZ=NULL
  for(i in seq_along(dyadidsx)){
    resNAphiZ[[i]]=dfZ %>% ungroup() %>% arrange(dyadid,burstvalcumeNA,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcumeNA) %>%
      map(~ psych::phi(matrix(.$total, 2, 2, byrow=T))) %>%    unlist
  }
  
  names(resNAphiZ)=dyadidsx
  
  resNAphiZdf = resNAphiZ %>%   
    map(~ data.frame(val=., burstno=names(.))) %>% 
    Map(cbind, ., dyadid=dyadidsx) %>%
    do.call('rbind', .)
  
  resNAphiZdf$burstno = as.numeric(as.character(resNAphiZdf$burstno))
  
  
  ## adding in p-value
  resNAphiZp=NULL
  for(i in seq_along(dyadidsx)){
    resNAphiZp[[i]]=dfZ %>% ungroup() %>% arrange(dyadid,burstvalcumeNA,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcumeNA) %>%
      map(~ phicoeff(matrix(.$total, 2, 2, byrow=T))$pval) %>%    unlist
  }
  
  names(resNAphiZp)=dyadidsx
  
  resNAphiZdfp = resNAphiZp %>% 
    map(~ data.frame(pval=., burstno=names(.))) %>% 
    Map(cbind, ., dyadid=dyadidsx) %>%
    do.call('rbind', .)
  
  resNAphiZdfp$burstno = as.numeric(as.character(resNAphiZdfp$burstno))
  
  resNAphiZdfx = resNAphiZdf %>% left_join(resNAphiZdfp)
  resNAphiZdfx$pvalx = ifelse(resNAphiZdfx$pval<0.05, "sig", "NS")
  
  
  
  return(resNAphiZdfx)
  
}


#####

getmats <- function(df){
  
  dfx <- df
  
  dyadidsx <- dfx$dyadid %>% unique
  
  
  #get matrices
  resphiZslide=NULL
  for(i in seq_along(dyadidsx)){
    resphiZslide[[i]]=dfx %>% ungroup() %>% 
      arrange(dyadid,burstvalcume,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcume) %>%
      map(~ .$total) %>% map(~ matrix(., 2, 2, byrow=T) ) 
  }
  
  names(resphiZslide)=dyadidsx
  
  return(resphiZslide)
}


###########

getphidf_slide <- function(df){
  
  dfx <- df
  
  dyadidsx <- dfx$dyadid %>% unique
  
  
  #get matrices
  resphiZslide=NULL
  for(i in seq_along(dyadidsx)){
    resphiZslide[[i]]=dfx %>% ungroup() %>% 
      arrange(dyadid,burstvalcume,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcume) %>%
      map(~ .$total) %>% map(~ matrix(., 2, 2, byrow=T) ) 
  }
  
  
  resphiZslideX <- lapply(resphiZslide, function(x) Map(`+`, head(x,-1), tail(x,-1))) #sliding function - every 2
  resphiZslideX_phi <- lapply(resphiZslideX, function(x) x %>% map(~ psych::phi(.)) %>%    unlist)
  resphiZslideX_phipval <- lapply(resphiZslideX, function(x) x %>% map(~ phicoeff(.)$pval) %>%    unlist)
  names(resphiZslideX_phi)=dyadidsx
  
  resphiZslideX_phidf <- resphiZslideX_phi %>%   
    Map(cbind, ., dyadid=dyadidsx) %>%
    do.call('rbind', .) %>% data.frame %>% group_by(dyadid) %>% mutate(burstno = row_number()) %>% data.frame %>%
    mutate(pval = resphiZslideX_phipval %>% unlist) 
  
  resphiZslideX_phidf$pval   <- as.numeric(resphiZslideX_phidf$pval )
  resphiZslideX_phidf$pvalx  <- ifelse(resphiZslideX_phidf$pval<0.05, "sig", "NS")
  
  return(resphiZslideX_phidf)
}

###########

getmats_slide <- function(df){
  
  dfx <- df
  
  dyadidsx <- dfx$dyadid %>% unique
  
  
  #get matrices
  resphiZslide=NULL
  for(i in seq_along(dyadidsx)){
    resphiZslide[[i]]=dfx %>% ungroup() %>% 
      arrange(dyadid,burstvalcume,idAB,state) %>%
      filter(dyadid==dyadidsx[[i]]) %>% split(., .$burstvalcume) %>%
      map(~ .$total) %>% map(~ matrix(., 2, 2, byrow=T) ) 
  }
  
  names(resphiZslide)=dyadidsx
  
  resphiZslideX <- lapply(resphiZslide, function(x) Map(`+`, head(x,-1), tail(x,-1))) #sliding function - every 2
  
  
  return(resphiZslideX)
}


#####


# Find relationship that is resolved phi coefficients.
# df = phicoeff results df, runs=consecutive runs conform to parameters, 
# consec - do they have to occur consecutivley, P=.05 pvalue to be under, 
# runs = "all" - whole length needs to be considered, not just next 3.
#NA are problems - not necessarily a negative result e.g. E burst 27

findresolution <- function(df, runs=3, consec=T, P=.05){
  
  dfburst.phi=df
  dfburst.phi$v <- ifelse(dfburst.phi$val>0 & dfburst.phi$pval<P, 1, 0) #could set p criteria to e.g. 0.1
  dfburst.phi$v[is.na(dfburst.phi$v)]<-0
  
  sp <- split(dfburst.phi$v, as.character(dfburst.phi$dyadid))
  spc <- lapply(sp, function(x) ifelse(x==0, 0, unlist(sapply(rle(x)$lengths, function(x) {return(1:x)}))))
  
  if(consec==T) {spw3 <- lapply(sp, function(x) which(cumsum(x)==runs)[1])}
  if(consec==F) {spw3 <- lapply(spc, function(x) which(x==runs)[1]) }
  
  spw4 <- lapply(spc, function(x) rev(which(x==1))[[1]])
  return(spw4)
}

#to ignore behaviors that have <5 or no phi-value-  these are usually values where only one animal behaves (is aggressive)
#however the run cannot start at these behaviors - as we do not have sufficient evidence to include them yet.
findresolution_all <- function(df, P=.05, behavs=5){

  newdf=df
  mynames=unique(df$dyadid)
  
  x <- ifelse(newdf$val>0 & newdf$pval<P, 1, 0) 
  
  sp <- split(x, as.character(newdf$dyadid))
  
  helper<-function(x,newdf){
    N<-which(x==1)[1]
    x[N:length(x)][is.na(x[N:length(x)])]<-1
    return(x)
  }
  
  sp <- lapply(sp, helper) #NA fill back in.
  
  x1=unlist(sp)
  x1[newdf$totalbehav < behavs  & newdf$val>0 & !is.na(x1)]<-1 # put 1's in for behavior but as long as NA not true
  x1[is.na(x1)]<-0 #remaining NAs should be zero?
  sp1 <- split(x1, as.character(newdf$dyadid))
  spc <- lapply(sp1, function(x) ifelse(x==0, 0, unlist(sapply(rle(x)$lengths, function(x) {return(1:x)}))))
  outp=unlist(lapply(spc, function(x)  ifelse(x[length(x)]==0, NA, ifelse(all(x)==T,  1, rev(which((x>0)==F))[[1]] + 1 )) ) )#technically this should add more if next value was an NA in sequence
  names(outp)=mynames
  return(outp)
  
  
}

##########

findresolutionAS = function(df, id="A", criterion="strict", buffer=3){
  BUFFER=buffer
  ID=id
  bursttotals=df
  xdf=subset(bursttotals,dyadid==ID)
  aggline=subset(xdf,idAB=="B" & state=="aggr")$total - subset(xdf,idAB=="A" & state=="aggr")$total
  subline=subset(xdf,idAB=="B" & state=="sub")$total - subset(xdf,idAB=="A" & state=="sub")$total
  
  if(criterion=="strict"){
    # STRICT CRITERION
    matplot(cbind(aggline,subline),ylim=c(-35,35),type="l",lty=c(1,1),col=c("dodgerblue","firebrick"),ylab="")
    title(ylab="Difference", line = 2)
    title(xlab="Burst Number", line = 2)
    title(ID, line = 0.5)
    abline(h=0,lty=2)
    cond1=(subline>=aggline)
    cond2=(subline>=0)
    resolved=(cond1&cond2)
    time=ifelse(any(resolved==F),max(which(resolved==FALSE))+1,1)
    abline(v=time,col="black")
    #cat("Resolved Time = ",time)
    return(time)
  }
  else
    if(criterion=="buffer"){
      # LOOSE CRITERION - only for sub diff > aggr diff.   Sub diff always has to be higher than 0. 
      BUFFER=3
      matplot(cbind(aggline,subline),ylim=c(-35,35),type="l",lty=c(1,1),col=c("dodgerblue","firebrick"),ylab="")
      title(ylab="Difference", line = 2)
      title(xlab="Burst Number", line = 2)
      title(ID, line = 0.5)
      abline(h=0,lty=2)
      cond3=(subline>=aggline-BUFFER)
      cond4=(subline>=0)
      resolved=(cond3&cond4)
      time=ifelse(any(resolved==F),max(which(resolved==FALSE))+1,1)
      abline(v=time,col="green")
     # cat("Resolved Time = ",time)
      return(time)
      }
}



#custom ggplot theme for these data
library(ggplot2)
newggtheme <- theme(
  plot.title = element_text(hjust = 0, vjust = 1, size = rel(2.3)), 
  panel.background = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black"),
  plot.background = element_blank(), 
  text = element_text(color = "gray20", size = 10), 
  axis.text = element_text(size = rel(1)), 
  axis.text.x = element_text(color = "gray20", size = rel(1.5)), 
  axis.text.y = element_text(color = "gray20", size = rel(1.5)), 
  axis.title.x = element_text(size = rel(1.5), vjust = 0), 
  axis.title.y = element_text(size = rel(1.5), vjust = 1), 
  axis.ticks.y = element_blank(), 
  axis.ticks.x = element_blank(), 
  strip.text.x = element_text(size = rel(2.3)),
  legend.position = "bottom",
  legend.key=element_rect(fill=NA),
  legend.title = element_blank(),
  legend.text=element_text(size=rel(1.7))
  )



######
plotbehav <- function(df,ID,legend=T,wd=NULL, behavA=NULL, behavB=NULL){
  
  df1 <- df %>% filter(dyadid==ID) %>% filter((behavior==behavA & idAB=="A") | (behavior==behavB & idAB=="B")) %>%
    arrange(day,time) %>% 
    mutate(ymin = ifelse(idAB=="B", 1, 2), ymax = ifelse(idAB=="A", 2, 3)) %>% select(dyadid,idAB,day,behavior,time,duration,ymin,ymax)
  
  #add empty df for days
  df1empty <- data.frame(dyadid=ID, 
                         idAB=rep(c("A", "B")),
                         day=1:5, 
                         behavior=rep(c(behavA,behavB),each=10),
                         time=NA, 
                         duration=NA,
                         ymin=NA, 
                         ymax=NA)
  
  
  
  df1x <- rbind(df1,df1empty)
  df1x$behavior <- factor(df1x$behavior)
  
  
  P=ggplot(df1x, aes(x=time,y=ymin, color=factor(behavior), fill=factor(behavior)))  + 
    scale_color_manual(values=c("red1", "blue")) +
    scale_fill_manual(values=c("red1", "blue")) +
    facet_wrap(~day, ncol=1) +
    theme_bw() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    ) + xlim(0,1200)+
    xlab("Time - seconds") + ylab("")
  
  P=P+ annotate("text", label = "A", size = 4, x = -Inf, y = 2,hjust=-.6)+
    annotate("text", label = "B", size = 4, x = -Inf, y = 1,hjust=-.5)
  
  if(is.null(wd) & legend==T){ P = P + geom_tile(aes(width=duration)) }  
  
  if(is.null(wd) & legend!=T){ P = P + geom_tile(aes(width=duration)) + theme(legend.position = "none") }
  
  if(!is.null(wd) & legend==T){ P = P + geom_tile(aes(width=wd))  }
  
  if(!is.null(wd) & legend!=T){ P = P + geom_tile(aes(width=wd)) + theme(legend.position = "none") }
  
  
  return(P)
}



#### Get STTC...

getsttc <- function(df,ID=NULL, behavA=NULL, behavB=NULL, dt=2, n=3){
  
  
  df$abstimeX <- ifelse(df$day==1, df$abstime, 
                        ifelse(df$day==2, df$abstime+n,
                               ifelse(df$day==3, df$abstime+(2*n),
                                      ifelse(df$day==4, df$abstime+(3*n),
                                             ifelse(df$day==5, df$abstime+(4*n), NA)))))
  
  
  ## preAB idA.behavA-->idB.behavB
  v1 = df %>% filter(dyadid==ID) %>% filter(prepostphi=="pre") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  v2 = df %>% filter(dyadid==ID) %>% filter(prepostphi=="pre") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  START = 0
  END = df %>% filter(dyadid==ID) %>% filter(prepostphi=="post") %>% filter(abstimeX == min(abstimeX)) %>% .$abstimeX
  DT = dt
  
  corrpreAB = .C("run_sttc2", as.integer(length(v1)), as.integer(length(v2)), 
                 as.double(DT), rec.time = range(c(START, END))
                 , coeff = double(1), as.double(v1), 
                 as.double(v2))
  
  preAB <- corrpreAB$coeff
  
  ## preBA idB.behavB-->idA.behavA
  corrpreBA = .C("run_sttc2", as.integer(length(v2)), as.integer(length(v1)), 
                 as.double(DT), rec.time = range(c(START, END))
                 , coeff = double(1), as.double(v2), 
                 as.double(v1))
  
  preBA <- corrpreBA$coeff
  

  
  
  ## postAB idA.behavA-->idB.behavB
  v1x = df %>% filter(dyadid==ID) %>% filter(prepostphi=="post") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  v2x = df %>% filter(dyadid==ID) %>% filter(prepostphi=="post") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  STARTx = df %>% filter(dyadid==ID) %>% filter(prepostphi=="post") %>% filter(abstimeX == min(abstimeX)) %>% .$abstimeX
  #ENDx = 6000
  ENDx = 6000 + (4*n)
  
  corrpostAB = .C("run_sttc2", as.integer(length(v1x)), as.integer(length(v2x)), 
                  as.double(DT), rec.time = range(c(STARTx, ENDx))
                  , coeff = double(1), as.double(v1x), 
                  as.double(v2x))
  
  postAB <- corrpostAB$coeff
  
  
  
  ## postBA idB.behavB-->idA.behavA
  
  corrpostBA = .C("run_sttc2", as.integer(length(v2x)), as.integer(length(v1x)), 
                  as.double(DT), rec.time = range(c(STARTx, ENDx))
                  , coeff = double(1), as.double(v2x), 
                  as.double(v1x))
  
  postBA <- corrpostBA$coeff
  
  
  
  
  ## preABrev idB.behavA-->idA.behavB
  z1 = df %>% filter(dyadid==ID) %>% filter(prepostphi=="pre") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  z2 = df %>% filter(dyadid==ID) %>% filter(prepostphi=="pre") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  corrpreABz = .C("run_sttc2", as.integer(length(z1)), as.integer(length(z2)), 
                 as.double(DT), rec.time = range(c(START, END))
                 , coeff = double(1), as.double(z1), 
                 as.double(z2))
  
  preABz <- corrpreABz$coeff
  
  ## preBA idA.behavB-->idB.behavA
  corrpreBAz = .C("run_sttc2", as.integer(length(z2)), as.integer(length(z1)), 
                 as.double(DT), rec.time = range(c(START, END))
                 , coeff = double(1), as.double(z2), 
                 as.double(z1))
  
  preBAz <- corrpreBAz$coeff
  
  
  
  
  ## postAB idB.behavA-->idA.behavB
  z1x = df %>% filter(dyadid==ID) %>% filter(prepostphi=="post") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  z2x = df %>% filter(dyadid==ID) %>% filter(prepostphi=="post") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  corrpostABz = .C("run_sttc2", as.integer(length(z1x)), as.integer(length(z2x)), 
                  as.double(DT), rec.time = range(c(STARTx, ENDx))
                  , coeff = double(1), as.double(z1x), 
                  as.double(z2x))
  
  postABz <- corrpostABz$coeff
  
  
  
  ## postBA idA.behavB-->idB.behavA
  
  corrpostBAz = .C("run_sttc2", as.integer(length(z2x)), as.integer(length(z1x)), 
                  as.double(DT), rec.time = range(c(STARTx, ENDx))
                  , coeff = double(1), as.double(z2x), 
                  as.double(z1x))
  
  postBAz <- corrpostBAz$coeff
  
  
  return(list("Pre: dom behavA --> sub behavB"=preAB, 
              "Pre: sub behavA --> dom behavB"=preABz, 
              "Pre: sub behavB --> dom behavA"=preBA, 
              "Pre: dom behavB --> sub behavA"=preBAz, 
              "Post: dom behavA --> sub behavB"=postAB, 
              "Post: sub behavA --> dom behavB"=postABz, 
              "Post: sub behavB --> dom behavA"=postBA, 
              "Post: dom behavB --> sub behavA"=postBAz)
         )
}







#### Get STTC with middle...

getsttcM <- function(df,ID=NULL, behavA=NULL, behavB=NULL, dt=2, n=3){
  
  
  df$abstimeX <- ifelse(df$day==1, df$abstime, 
                        ifelse(df$day==2, df$abstime+n,
                               ifelse(df$day==3, df$abstime+(2*n),
                                      ifelse(df$day==4, df$abstime+(3*n),
                                             ifelse(df$day==5, df$abstime+(4*n), NA)))))
  
  
  ## preAB idA.behavA-->idB.behavB
  v1 = df %>% filter(dyadid==ID) %>% filter(phase=="pre") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  v2 = df %>% filter(dyadid==ID) %>% filter(phase=="pre") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  START = 0
  END = df %>% filter(dyadid==ID) %>% filter(phase=="middle") %>% filter(abstimeX == min(abstimeX)) %>% .$abstimeX
  DT = dt
  
  corrpreAB = .C("run_sttc2", as.integer(length(v1)), as.integer(length(v2)), 
                 as.double(DT), rec.time = range(c(START, END))
                 , coeff = double(1), as.double(v1), 
                 as.double(v2))
  
  preAB <- corrpreAB$coeff
  
  ## preBA idB.behavB-->idA.behavA
  corrpreBA = .C("run_sttc2", as.integer(length(v2)), as.integer(length(v1)), 
                 as.double(DT), rec.time = range(c(START, END))
                 , coeff = double(1), as.double(v2), 
                 as.double(v1))
  
  preBA <- corrpreBA$coeff
  
  
  
  
  ## postAB idA.behavA-->idB.behavB
  v1x = df %>% filter(dyadid==ID) %>% filter(phase=="post") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  v2x = df %>% filter(dyadid==ID) %>% filter(phase=="post") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  STARTx = df %>% filter(dyadid==ID) %>% filter(phase=="post") %>% filter(abstimeX == min(abstimeX)) %>% .$abstimeX
  #ENDx = 6000
  ENDx = 6000 + (4*n)
  
  corrpostAB = .C("run_sttc2", as.integer(length(v1x)), as.integer(length(v2x)), 
                  as.double(DT), rec.time = range(c(STARTx, ENDx))
                  , coeff = double(1), as.double(v1x), 
                  as.double(v2x))
  
  postAB <- corrpostAB$coeff
  
  
  
  ## postBA idB.behavB-->idA.behavA
  
  corrpostBA = .C("run_sttc2", as.integer(length(v2x)), as.integer(length(v1x)), 
                  as.double(DT), rec.time = range(c(STARTx, ENDx))
                  , coeff = double(1), as.double(v2x), 
                  as.double(v1x))
  
  postBA <- corrpostBA$coeff
  
  
  
  
  ## preABrev idB.behavA-->idA.behavB
  z1 = df %>% filter(dyadid==ID) %>% filter(phase=="pre") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  z2 = df %>% filter(dyadid==ID) %>% filter(phase=="pre") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  corrpreABz = .C("run_sttc2", as.integer(length(z1)), as.integer(length(z2)), 
                  as.double(DT), rec.time = range(c(START, END))
                  , coeff = double(1), as.double(z1), 
                  as.double(z2))
  
  preABz <- corrpreABz$coeff
  
  ## preBA idA.behavB-->idB.behavA
  corrpreBAz = .C("run_sttc2", as.integer(length(z2)), as.integer(length(z1)), 
                  as.double(DT), rec.time = range(c(START, END))
                  , coeff = double(1), as.double(z2), 
                  as.double(z1))
  
  preBAz <- corrpreBAz$coeff
  
  
  
  
  ## postAB idB.behavA-->idA.behavB
  z1x = df %>% filter(dyadid==ID) %>% filter(phase=="post") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  z2x = df %>% filter(dyadid==ID) %>% filter(phase=="post") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  corrpostABz = .C("run_sttc2", as.integer(length(z1x)), as.integer(length(z2x)), 
                   as.double(DT), rec.time = range(c(STARTx, ENDx))
                   , coeff = double(1), as.double(z1x), 
                   as.double(z2x))
  
  postABz <- corrpostABz$coeff
  
  
  
  ## postBA idA.behavB-->idB.behavA
  
  corrpostBAz = .C("run_sttc2", as.integer(length(z2x)), as.integer(length(z1x)), 
                   as.double(DT), rec.time = range(c(STARTx, ENDx))
                   , coeff = double(1), as.double(z2x), 
                   as.double(z1x))
  
  postBAz <- corrpostBAz$coeff
  
  
  ### Middle phases:
  
  ## midAB idA.behavA-->idB.behavB
  v1m = df %>% filter(dyadid==ID) %>% filter(phase=="middle") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  v2m = df %>% filter(dyadid==ID) %>% filter(phase=="middle") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  STARTm = df %>% filter(dyadid==ID) %>% filter(phase=="middle") %>% filter(abstimeX == min(abstimeX)) %>% .$abstimeX
  ENDm = df %>% filter(dyadid==ID) %>% filter(phase=="post") %>% filter(abstimeX == min(abstimeX)) %>% .$abstimeX
  
  corrmidAB = .C("run_sttc2", as.integer(length(v1m)), as.integer(length(v2m)), 
                 as.double(DT), rec.time = range(c(STARTm, ENDm))
                 , coeff = double(1), as.double(v1m), 
                 as.double(v2m))
  
  midAB <- corrmidAB$coeff
  
  ## preBA idB.behavB-->idA.behavA
  corrmidBA = .C("run_sttc2", as.integer(length(v2m)), as.integer(length(v1m)), 
                 as.double(DT), rec.time = range(c(STARTm, ENDm))
                 , coeff = double(1), as.double(v2m), 
                 as.double(v1m))
  
  midBA <- corrmidBA$coeff
  
  
  ## midAB idB.behavA-->idA.behavB
  z1m = df %>% filter(dyadid==ID) %>% filter(phase=="middle") %>% filter(idAB=="B") %>%
    filter(behavior %in% behavA) %>% .$abstimeX %>% sort
  
  z2m = df %>% filter(dyadid==ID) %>% filter(phase=="middle") %>% filter(idAB=="A") %>%
    filter(behavior %in% behavB) %>% .$abstimeX %>% sort
  
  corrmidABz = .C("run_sttc2", as.integer(length(z1m)), as.integer(length(z2m)), 
                  as.double(DT), rec.time = range(c(STARTm, ENDm))
                  , coeff = double(1), as.double(z1m), 
                  as.double(z2m))
  
  midABz <- corrmidABz$coeff
  
  ## midBA idA.behavB-->idB.behavA
  corrmidBAz = .C("run_sttc2", as.integer(length(z2m)), as.integer(length(z1m)), 
                  as.double(DT), rec.time = range(c(STARTm, ENDm))
                  , coeff = double(1), as.double(z2m), 
                  as.double(z1m))
  
  midBAz <- corrmidBAz$coeff
  
  

  return(list("Pre: dom behavA --> sub behavB"=preAB, 
              "Pre: sub behavA --> dom behavB"=preABz, 
              "Pre: sub behavB --> dom behavA"=preBA, 
              "Pre: dom behavB --> sub behavA"=preBAz, 
              "Post: dom behavA --> sub behavB"=postAB, 
              "Post: sub behavA --> dom behavB"=postABz, 
              "Post: sub behavB --> dom behavA"=postBA, 
              "Post: dom behavB --> sub behavA"=postBAz,
         "Mid: dom behavA --> sub behavB"=midAB, 
         "Mid: sub behavA --> dom behavB"=midABz, 
         "Mid: sub behavB --> dom behavA"=midBA, 
         "Mid: dom behavB --> sub behavA"=midBAz)
  )
}


#needs output of for loop for calculating sttc values of 21 groups

plotsttc_means <- function(resoutM){
  
  #DOM aggr --> SUB sub behavior  == 1-pre, 5-post
  
  ppM <- data.frame(
    pre = lapply(resoutM, function(x) x[[1]]) %>% unlist ,
    mid = lapply(resoutM, function(x) x[[9]]) %>% unlist ,
    post = lapply(resoutM, function(x) x[[5]]) %>% unlist
  ) 
  
  
  #add blue ribbon
  rbind(
    ppM %>% summarise_each(funs(mean(., na.rm = TRUE))),
    ppM %>% summarise_each(funs(sem))
  ) %>% t %>% data.frame %>%
    select(mean=1,sem=2) %>%
    mutate(upr=mean+sem, lwr=mean-sem) -> rib
  
  rib$condition = factor(c("pre", "mid", "post"), levels=c("pre", "mid", "post"))
  
  ppMdf <- reshape2::melt(ppM) %>% mutate(id=rep(LETTERS[1:21],3))
  
  
  ### Compare to SUB aggr to DOM sub...
  
  ppMx <- data.frame(
    pre = lapply(resoutM, function(x) x[[2]]) %>% unlist ,
    mid = lapply(resoutM, function(x) x[[10]]) %>% unlist ,
    post = lapply(resoutM, function(x) x[[6]]) %>% unlist
  ) 
  
  
  #add blue ribbon
  rbind(
    ppMx %>% summarise_each(funs(mean(., na.rm = TRUE))),
    ppMx %>% summarise_each(funs(sem))
  ) %>% t %>% data.frame %>%
    select(mean=1,sem=2) %>%
    mutate(upr=mean+sem, lwr=mean-sem) -> ribx
  
  ribx$condition = factor(c("pre", "mid", "post"), levels=c("pre", "mid", "post"))
  
  ppMxdf <- reshape2::melt(ppMx) %>% mutate(id=rep(LETTERS[1:21],3))
  
  
  ppMboth <- rbind(ppM,ppMx)
  ppMboth$group <- rep(c("DOM-SUB", "SUB-DOM"),each=21)
  
  ribboth <- rbind(rib, ribx)
  ribboth$group <- rep(c("DOM-SUB", "SUB-DOM"),each=3)
  
  gp = ggplot() + 
    geom_ribbon(data=ribboth, aes(x=condition, ymin=lwr, ymax=upr, group=group, fill=group), alpha=0.2)+
    geom_line(data=ribboth, aes(x=condition, y = mean, group=group, color = group), lwd=1) + newggtheme + 
    scale_color_manual(values=c("red", "dodgerblue")) +
    scale_fill_manual(values=c("red", "dodgerblue")) 
  
  return(gp)
  
}


#### Get Results of STTC

sttc_ppm <- function(df, behavA = aggrbehav, behavB= subbehav,dt=2,n=2.01){
  resoutM=NULL
  for(i in 1:21){
    resoutM[[i]] <- getsttcM(df, ID=LETTERS[i], behavA, behavB,dt,n)
  }
  
  return(resoutM)
}






### Lag Seq----


#Get matrix - prob=F is observed frequencies, prob=T = transitional freq
obsmatrix <- function(x, prob=T) {
  X<-t(data.frame(uncollapse(x)))
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  if(prob) tt <- tt / rowSums(tt)
tt
}


## THIS IS FOR STRUCTURAL ZEROS ALONG THE DIAGONAL
lrx_inc <- function(obs){
  sz=length(diag(obs))
  expectm <- expected_inc(obs)
  aa3<-log(obs/expectm)
  aa3[aa3=='-Inf']<-0    #if get 0 observed value then get -Inf for log(obs/expected)
  LRX <- 2*(sum(obs*(aa3),na.rm=T))  #Likelihood Ratio
  dff <-  ((nrow(obs)-1)^2) - sz
  pvallrx<-1-pchisq(LRX,dff)
  return(list('LRX-sq'= LRX, 'df'=dff, 'p-value'= pvallrx))
}


expected_inc <- function(x){
  
  NN<-nrow(x)
  library(mipfp)
  seed.2d <- array(1,dim=c(NN,NN)) # desired targets (margins) : V1 and V2
  diag(seed.2d)<-0
  target.row <- rowSums(x)
  target.col <- colSums(x)
  tgt.data.2d <- list(target.row, target.col) # storing the margins in a list
  tgt.list.2d <- list(1,2) # list of dimensions of each marginal constrain
  
  # calling the Ipfp function, requesting the covariance matrices
  res.2d <- Ipfp(seed.2d, tgt.list.2d, tgt.data.2d)
  out <- res.2d$x.hat
  rownames(out)<-colnames(out)<-rownames(x)
  return(out)
}



epoch = function(x) {
  
  if(max(table(x)/sum(table(x)))>=0.5) stop("One code is greater than 50%  of the sequence - no solution possible")
  
  
  checkleft = function(p, i) {
    (p == 1) || x[p-1] != x[i]
  }
  
  checkright = function(p, i) {
    (p == nelem) || x[p+1] != x[i]
  }
  
  nelem = length(x)
  for (i in seq_len(nelem)) {
    looking = TRUE
    while(looking) {
      ok = TRUE
      p = sample(nelem, 1)
      ok = checkleft(p, i)  && checkright(p, i) &&
        checkleft(i, p)  && checkright(i, p)
      if (ok) {
        ## constraints satisfied, so swap elements i and p
        looking = FALSE
        old = x[i]
        x[i] = x[p]
        x[p] = old
      }
    }
  }
  x
}

permseq_norpt = function(x,Nperms=1000){
  
  vals <- table(uncollapse(x))
  freq = vals
  nletters = length(freq)
  
  
  x1 = rep(LETTERS[1:nletters], times=freq)
  
  permedseqs=list()
  anyleft<-list()
  equals <- list()
  
  for(i in 1:Nperms){
    
    x1x = sample(x1)
    x2 = epoch(x1x)
    permedseqs[[i]]=x2
    
    len = length(x2)
    anyleft[[i]] = sum(x2[1:(len-1)] == x2[2:len])
    equals[[i]] = all.equal(as.vector(table(x1x)), as.vector(table(x2)))
    
  }
  return(list('permedseqs'=permedseqs, 'anyleft'=anyleft, 'equals'=equals))
}



##### for 2 sequences synched

permseq_norpt2 = function(x,Nperms=1000){
  
  vals <- table(uncollapse(x))
  
  frequpper = vals[grepl("^[[:upper:]]+$", names(vals))]
  freqlower = vals[!grepl("^[[:upper:]]+$", names(vals))]
  
  nlettersupper = length(frequpper)
  nletterslower = length(freqlower)
  
#  freq = vals
#  nletters = length(freq)
  
#  x1 = rep(LETTERS[1:nletters], times=freq)
 
  x1 = c(rep(LETTERS[1:nlettersupper], times=frequpper),rep(letters[1:nletterslower], times=freqlower))
  
  permedseqs=list()
  anyleft<-list()
  equals <- list()
  
  for(i in 1:Nperms){
    
    x1x = sample(x1)
    x2 = epoch(x1x)
    permedseqs[[i]]=x2
    
    len = length(x2)
    anyleft[[i]] = sum(x2[1:(len-1)] == x2[2:len])
    equals[[i]] = all.equal(as.vector(table(x1x)), as.vector(table(x2)))
    
  }
  return(list('permedseqs'=permedseqs, 'anyleft'=anyleft, 'equals'=equals))
}




#### Multiple Post-Hoc Tests for Friedman's Test
# http://www.r-statistics.com/2010/02/post-hoc-analysis-for-friedmans-test-r-code/

friedman.test.with.post.hoc <- function(formu, data, to.print.friedman = T, to.post.hoc.if.signif = T,  to.plot.parallel = T, to.plot.boxplot = T, signif.P = .05, color.blocks.in.cor.plot = T, jitter.Y.in.cor.plot =F)
{
  # formu is a formula of the shape: 	Y ~ X | block
  # data is a long data.frame with three columns:    [[ Y (numeric), X (factor), block (factor) ]]
  
  # Note: This function doesn't handle NA's! In case of NA in Y in one of the blocks, then that entire block should be removed.
  
  
  # Loading needed packages
  if(!require(coin))
  {
    print("You are missing the package 'coin', we will now try to install it...")
    install.packages("coin")
    library(coin)
  }
  
  if(!require(multcomp))
  {
    print("You are missing the package 'multcomp', we will now try to install it...")
    install.packages("multcomp")
    library(multcomp)
  }
  
  if(!require(colorspace))
  {
    print("You are missing the package 'colorspace', we will now try to install it...")
    install.packages("colorspace")
    library(colorspace)
  }
  
  
  # get the names out of the formula
  formu.names <- all.vars(formu)
  Y.name <- formu.names[1]
  X.name <- formu.names[2]
  block.name <- formu.names[3]
  
  if(dim(data)[2] >3) data <- data[,c(Y.name,X.name,block.name)]	# In case we have a "data" data frame with more then the three columns we need. This code will clean it from them...
  
  # Note: the function doesn't handle NA's. In case of NA in one of the block T outcomes, that entire block should be removed.
  
  # stopping in case there is NA in the Y vector
  if(sum(is.na(data[,Y.name])) > 0) stop("Function stopped: This function doesn't handle NA's. In case of NA in Y in one of the blocks, then that entire block should be removed.")
  
  # make sure that the number of factors goes with the actual values present in the data:
  data[,X.name ] <- factor(data[,X.name ])
  data[,block.name ] <- factor(data[,block.name ])
  number.of.X.levels <- length(levels(data[,X.name ]))
  if(number.of.X.levels == 2) { warning(paste("'",X.name,"'", "has only two levels. Consider using paired wilcox.test instead of friedman test"))}
  
  # making the object that will hold the friedman test and the other.
  the.sym.test <- symmetry_test(formu, data = data,	### all pairwise comparisons
                                teststat = "max",
                                xtrafo = function(Y.data) { trafo( Y.data, factor_trafo = function(x) { model.matrix(~ x - 1) %*% t(contrMat(table(x), "Tukey")) } ) },
                                ytrafo = function(Y.data){ trafo(Y.data, numeric_trafo = rank, block = data[,block.name] ) }
  )
  # if(to.print.friedman) { print(the.sym.test) }
  
  
  if(to.post.hoc.if.signif)
  {
    if(pvalue(the.sym.test) < signif.P)
    {
      # the post hoc test
      The.post.hoc.P.values <- pvalue(the.sym.test, method = "single-step")	# this is the post hoc of the friedman test
      
      
      # plotting
      if(to.plot.parallel & to.plot.boxplot)	par(mfrow = c(1,2)) # if we are plotting two plots, let's make sure we'll be able to see both
      
      if(to.plot.parallel)
      {
        X.names <- levels(data[, X.name])
        X.for.plot <- seq_along(X.names)
        plot.xlim <- c(.7 , length(X.for.plot)+.3)	# adding some spacing from both sides of the plot
        
        if(color.blocks.in.cor.plot)
        {
          blocks.col <- rainbow_hcl(length(levels(data[,block.name])))
        } else {
          blocks.col <- 1 # black
        }
        
        data2 <- data
        if(jitter.Y.in.cor.plot) {
          data2[,Y.name] <- jitter(data2[,Y.name])
          par.cor.plot.text <- "Parallel coordinates plot (with Jitter)"
        } else {
          par.cor.plot.text <- "Parallel coordinates plot"
        }
        
        # adding a Parallel coordinates plot
        matplot(as.matrix(reshape(data2,  idvar=X.name, timevar=block.name,
                                  direction="wide")[,-1])  ,
                type = "l",  lty = 1, axes = FALSE, ylab = Y.name,
                xlim = plot.xlim,
                col = blocks.col,
                main = par.cor.plot.text)
        axis(1, at = X.for.plot , labels = X.names) # plot X axis
        axis(2) # plot Y axis
        points(tapply(data[,Y.name], data[,X.name], median) ~ X.for.plot, col = "red",pch = 4, cex = 2, lwd = 5)
      }
      
      if(to.plot.boxplot)
      {
        # first we create a function to create a new Y, by substracting different combinations of X levels from each other.
        subtract.a.from.b <- function(a.b , the.data)
        {
          the.data[,a.b[2]] - the.data[,a.b[1]]
        }
        
        temp.wide <- reshape(data,  idvar=X.name, timevar=block.name,
                             direction="wide") 	#[,-1]
        wide.data <- as.matrix(t(temp.wide[,-1]))
        colnames(wide.data) <- temp.wide[,1]
        
        Y.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 2, subtract.a.from.b, the.data =wide.data)
        names.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 2, function(a.b) {paste(a.b[2],a.b[1],sep=" - ")})
        
        the.ylim <- range(Y.b.minus.a.combos)
        the.ylim[2] <- the.ylim[2] + max(sd(Y.b.minus.a.combos))	# adding some space for the labels
        is.signif.color <- ifelse(The.post.hoc.P.values < .05 , "green", "grey")
        
        boxplot(Y.b.minus.a.combos,
                names = names.b.minus.a.combos ,
                col = is.signif.color,
                main = "Boxplots (of the differences)",
                ylim = the.ylim
        )
        legend("topright", legend = paste(names.b.minus.a.combos, rep(" ; PostHoc P.value:", number.of.X.levels),round(The.post.hoc.P.values,5)) , fill =  is.signif.color )
        abline(h = 0, col = "blue")
        
      }
      
      list.to.return <- list(Friedman.Test = the.sym.test, PostHoc.Test = The.post.hoc.P.values)
      if(to.print.friedman) {print(list.to.return)}
      return(list.to.return)
      
    }	else {
      print("The results where not significant, There is no need for a post hoc test")
      return(the.sym.test)
    }
  }
  
  # Original credit (for linking online, to the package that performs the post hoc test) goes to "David Winsemius", see:
  # http://tolstoy.newcastle.edu.au/R/e8/help/09/10/1416.html
}

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#get substring from right
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}





OR<-function(x){
  
  if(sum(rowSums(x)<5)>0 | sum(colSums(x)<5)>0) warning ("Some cells contain extremely low counts - proceed with caution")
  #warning/advice given on page 116 of Bakeman and Quera 2011
  
  a<-x[1,1]
  b<-x[1,2]
  c<-x[2,1]
  d<-x[2,2]
  
  
  
  if(a!=0 & b!=0 & c!=0 & d!=0){
    or<- (a*d) / (b*c)
    lnor<-log(or)
    se.ln.or<-sqrt((1/a) + (1/c) + (1/b) + (1/d))
    uci <- exp(log(or) + (1.96*se.ln.or))
    lci <- exp(log(or) - (1.96*se.ln.or))
    
    yq<- ((a*d) - (b*c)) / ((a*d) + (b*c))
    
    return(list('OR'= or,'95%CIs' = c(lci,uci),  'ln-OR' = lnor, 'YulesQ' = yq))
  }
  
  else
    if(b==0 | c==0){
      yq<- ((a*d) - (b*c)) / ((a*d) + (b*c))      
      return(list('Cannot calculate OR or lnOR - cells b or c equal to 0', 'YulesQ' = yq))
    }
  
  else
    if(a==0 | d==0){
      or<- (a*d) / (b*c)
      yq<- ((a*d) - (b*c)) / ((a*d) + (b*c))      
      return(list('Cannot calculate lnOR or CIs - cells a or d equal to 0', 
                  'OR' = or, 'YulesQ' = yq))
    }
  
}




####  FUNCTIONS for two sequences synched 


### TURN VECTOR OF CONTINGENCIES + GIVEN/TARGET CODES TO 2X2 MATRIX:
tbtmatrix <- function(behav1,behav2,letter1,letter2){
  behav12 <- paste0(behav1,behav2)
  given_yes_sum <- sum(behav1==letter1)   #27
  target_yes_sum <- sum(behav2==letter2)  #30
  sum_total <- length(behav1) #172
  giventarget <- paste0(letter1,letter2)
  cellA <- sum(behav12==giventarget)   #19
  cellC <- target_yes_sum - cellA  #11
  cellB <- given_yes_sum - cellA #8
  cellD <- sum_total - cellA - cellB - cellC
  tbtmat <- matrix(c(cellA,cellC,cellB,cellD),2,2)
  return(tbtmat)
}

### Make as function (sub given --> dom target):
subdom_tbtmat <- function(subdombehav1,subdombehav2){
  tbtsubdomresults<-vector("list",9)
  for(i in 1:9){
    for(j in 1:9){
      tbtsubdomresults[[i]][[j]] <- tbtmatrix(subdombehav1,subdombehav2,letters[i],LETTERS[j])
    }
  }
  return(tbtsubdomresults)
}

subdom_tbtmat_state <- function(subdombehav1,subdombehav2){
  tbtsubdomresults<-vector("list",7)
  for(i in 1:7){
    for(j in 1:7){
      tbtsubdomresults[[i]][[j]] <- tbtmatrix(subdombehav1,subdombehav2,letters[i+19],LETTERS[j+19])
    }
  }
  return(tbtsubdomresults)
}



### Make as function (dom given --> sub target):
domsub_tbtmat <- function(domsubbehav1,domsubbehav2){
  tbtdomsubresults<-vector("list",9)
  for(i in 1:9){
    for(j in 1:9){
      tbtdomsubresults[[i]][[j]] <- tbtmatrix(domsubbehav1,domsubbehav2,LETTERS[i],letters[j])
    }
  }
  return(tbtdomsubresults)
}

domsub_tbtmat_state <- function(domsubbehav1,domsubbehav2){
  tbtdomsubresults<-vector("list",7)
  for(i in 1:7){
    for(j in 1:7){
      tbtdomsubresults[[i]][[j]] <- tbtmatrix(domsubbehav1,domsubbehav2,LETTERS[i+19],letters[j+19])
    }
  }
  return(tbtdomsubresults)
}


## Function to test proportion of unique dyads showing behavior pre-resolution:
preprop<-function(behav1,behav2){
  #1:22 mid; 22-42 post, 43-63 pre
  subdompren <- unlist(lapply(seq2all_subdom[43:63], function(x) x[[behav1]][[behav2]][1]))  #
  domsubpren <- unlist(lapply(seq2all_domsub[43:63], function(x) x[[behav1]][[behav2]][1])) #
  val1<-sum(subdompren >0)
  val2<-sum(domsubpren >0)
  val3 <-length(subdompren)
  
  #matrix of proportions - 
  contig<-paste(LETTERS[behav1],LETTERS[behav2],sep="-")
  mat1 <- matrix(c(val1,val3-val1,val2,val3-val2),2,2) 
  if(val1==0 & val2==0) { prop1 <- NA }
  else  if(val1==val3 & val2==val3) { prop1 <- NA }
  else {
    prop1 <- prop.test(mat1,correct = F)
  } #
  rownames(mat1)<-c("Y","N")
  colnames(mat1)<-c("subdom","domsub")
  return(list(contig,mat1,prop1))
}


## Function to test proportion of unique dyads showing behavior post-resolution:
postprop<-function(behav1,behav2){
  #1:22 mid; 22-42 post, 43-63 pre
  subdommidn <- unlist(lapply(seq2all_subdom[1:21], function(x) x[[behav1]][[behav2]][1]))   #1-21 mid, 18=R, 3=C, 1=cellA   ## 0dyads
  subdompostn <- unlist(lapply(seq2all_subdom[22:42], function(x) x[[behav1]][[behav2]][1]))  #22-42 post, 18=R, 3=C, 1=cellA   ## 2dyads

  domsubmidn <- unlist(lapply(seq2all_domsub[1:21], function(x) x[[behav1]][[behav2]][1])) #1-21 post, 18=R, 3=C, 1=cellA   ## 6dyads
  domsubpostn <- unlist(lapply(seq2all_domsub[22:42], function(x) x[[behav1]][[behav2]][1])) #22-42 post, 18=R, 3=C, 1=cellA   ## 8dyads
  
  val1<-sum(subdommidn+subdompostn >0)
  val2<-sum(domsubmidn+domsubpostn >0)
  val3 <-length(subdompostn)
  
  #matrix of proportions - 
  contig<-paste(LETTERS[behav1],LETTERS[behav2],sep="-")
  mat1 <- matrix(c(val1,val3-val1,val2,val3-val2),2,2) 
  if(val1==0 & val2==0) { prop1 <- NA }
  else  if(val1==val3 & val2==val3) { prop1 <- NA }
  else {
    prop1 <- prop.test(mat1,correct = F) #
  }
  
  rownames(mat1)<-c("Y","N")
  colnames(mat1)<-c("subdom","domsub")
  return(list(contig,mat1,prop1))
}

#function to turn sttc results into df
sttcdf <- function(sttcres){
  
  tmpdf=data.frame(val=unlist(sttcres), names=names(unlist(sttcres)), 
                   phase = rep(c("Pre","Post","Mid"),each=4),
                   idA = rep(c("dom","sub","sub","dom")),
                   idB = rep(c("sub","dom","dom","sub")),
                   behav1 = rep(c("A","A","B","B")),
                   behav2 = rep(c("B","B","A","A")),
                   dyadid = rep(LETTERS[1:21],each=12)
  )  %>% mutate(direction = paste(idA,idB,sep="-"),
                tofrom = paste(behav1,behav2,sep="-"))
  tmpdf$phase <- factor(tmpdf$phase, levels=c("Pre","Mid", "Post"))
  return(tmpdf)
}




## fucntion to get frequency fo between individual transition in each dyad 
domsub_pre_freq<-function(dyadn=1){
  result=data.frame(behav1=NA,behav2=NA,freq=NA)
  for(x in 1:19){
    for(y in 1:19){
      result[c(((x-1)*19+1):((x)*19)),1]<-x
      result[(x-1)*19+y,2]<-y
      result[(x-1)*19+y,3]<-seq2all_domsub_pre[[i]][[x]][[y]][1]
    }
  }
  result1<-result %>% 
    mutate(dyad=dyadn,
           sum=sum(seq2all_domsub_pre[[dyadn]][[1]][[1]])) %>% 
    select(dyad,behav1,behav2,freq,sum) 
  return(result1)
}

subdom_pre_freq<-function(dyadn=1){
  result=data.frame(behav1=NA,behav2=NA,freq=NA)
  for(x in 1:19){
    for(y in 1:19){
      result[c(((x-1)*19+1):((x)*19)),1]<-x
      result[(x-1)*19+y,2]<-y
      result[(x-1)*19+y,3]<-seq2all_subdom_pre[[i]][[x]][[y]][1]
    }
  }
  result1<-result %>% 
    mutate(dyad=dyadn,
           sum=sum(seq2all_subdom_pre[[dyadn]][[1]][[1]])) %>% 
    select(dyad,behav1,behav2,freq,sum)
  
  return(result1)
}


domsub_post_freq<-function(dyadn=1){
  result=data.frame(behav1=NA,behav2=NA,freq=NA)
  for(x in 1:19){
    for(y in 1:19){
      result[c(((x-1)*19+1):((x)*19)),1]<-x
      result[(x-1)*19+y,2]<-y
      result[(x-1)*19+y,3]<-seq2all_domsub_post[[i]][[x]][[y]][1]
    }
  }
  result1<-result %>% 
    mutate(dyad=dyadn,
           sum=sum(seq2all_domsub_post[[dyadn]][[1]][[1]])) %>% 
    select(dyad,behav1,behav2,freq,sum)

  return(result1)
}

subdom_post_freq<-function(dyadn=1){
  result=data.frame(behav1=NA,behav2=NA,freq=NA)
  for(x in 1:19){
    for(y in 1:19){
      result[c(((x-1)*19+1):((x)*19)),1]<-x
      result[(x-1)*19+y,2]<-y
      result[(x-1)*19+y,3]<-seq2all_subdom_post[[i]][[x]][[y]][1]
    }
  }
  result1<-result %>% 
    mutate(dyad=dyadn,
           sum=sum(seq2all_subdom_post[[dyadn]][[1]][[1]])) %>% 
    select(dyad,behav1,behav2,freq,sum)
  
  return(result1)
}


#for state
domsub_state_pre_freq<-function(dyadn=1){
  result=data.frame(state1=NA,state2=NA,freq=NA)
  for(x in 1:7){
    for(y in 1:7){
      result[c(((x-1)*7+1):((x)*7)),1]<-x
      result[(x-1)*7+y,2]<-y
      result[(x-1)*7+y,3]<-seq2all_domsub_state_pre[[i]][[x]][[y]][1]
    }
  }
  result1<-result %>% 
    mutate(dyad=dyadn,
           sum=sum(seq2all_domsub_state_pre[[dyadn]][[1]][[1]])) %>% 
    select(dyad,state1,state2,freq,sum) 
  return(result1)
}

subdom_state_pre_freq<-function(dyadn=1){
  result=data.frame(state1=NA,state2=NA,freq=NA)
  for(x in 1:7){
    for(y in 1:7){
      result[c(((x-1)*7+1):((x)*7)),1]<-x
      result[(x-1)*7+y,2]<-y
      result[(x-1)*7+y,3]<-seq2all_subdom_state_pre[[i]][[x]][[y]][1]
    }
  }
  result1<-result %>% 
    mutate(dyad=dyadn,
           sum=sum(seq2all_subdom_state_pre[[dyadn]][[1]][[1]])) %>% 
    select(dyad,state1,state2,freq,sum)
  
  return(result1)
}


domsub_state_post_freq<-function(dyadn=1){
  result=data.frame(state1=NA,state2=NA,freq=NA)
  for(x in 1:7){
    for(y in 1:7){
      result[c(((x-1)*7+1):((x)*7)),1]<-x
      result[(x-1)*7+y,2]<-y
      result[(x-1)*7+y,3]<-seq2all_domsub_state_post[[i]][[x]][[y]][1]
    }
  }
  result1<-result %>% 
    mutate(dyad=dyadn,
           sum=sum(seq2all_domsub_state_post[[dyadn]][[1]][[1]])) %>% 
    select(dyad,state1,state2,freq,sum)
  
  return(result1)
}

subdom_state_post_freq<-function(dyadn=1){
  result=data.frame(state1=NA,state2=NA,freq=NA)
  for(x in 1:7){
    for(y in 1:7){
      result[c(((x-1)*7+1):((x)*7)),1]<-x
      result[(x-1)*7+y,2]<-y
      result[(x-1)*7+y,3]<-seq2all_subdom_state_post[[i]][[x]][[y]][1]
    }
  }
  result1<-result %>% 
    mutate(dyad=dyadn,
           sum=sum(seq2all_subdom_state_post[[dyadn]][[1]][[1]])) %>% 
    select(dyad,state1,state2,freq,sum)
  
  return(result1)
}


#### to make sttc dataframe at once

sttc_ppm_df<-function(df0, behav1,behav2){
  
temp<-sttc_ppm(df0, behavA = behav1, behavB = behav2,dt=2.0,n=2.01)

temp_df<-data.frame(FSTTC=unlist(temp), names=names(unlist(temp)), 
           phase = rep(c("Pre","Post","Mid"),each=4),
           idA = rep(c("Dom","Sub","Sub","Dom")),
           idB = rep(c("Sub","Dom","Dom","Sub")),
           given = rep(c(names(behav1)[1],names(behav1)[1],names(behav2)[1],names(behav2)[1])),
           target = rep(c(names(behav2)[1],names(behav2)[1],names(behav1)[1],names(behav1)[1])),
           dyadid = rep(LETTERS[1:21],each=12)
)  %>% 
  mutate(direction = paste(idA,idB,sep="->"),Contingency = paste(given,target,sep="->")) %>% 
  select(FSTTC,dyadid,phase,direction,Contingency,names) %>% 
  as.data.frame() %>% 
  unique(.)
return(temp_df)
}


####
get_z_value<-function(domsubmatrix,subdommatrix){ 
  a1<-domsubmatrix[1,1]
  c1<-domsubmatrix[1,2]
  b1<-domsubmatrix[2,1]
  d1<-domsubmatrix[2,2]
  
  a2<-subdommatrix[1,1]
  c2<-subdommatrix[1,2]
  b2<-subdommatrix[2,1]
  d2<-subdommatrix[2,2]
  if(b1==0|c1==0|b2==0|c2==0){
    z=NA
  }
  else{
  z=( (log(OR(domsubmatrix)$OR)-log(OR(subdommatrix)$OR))/sqrt((1/a1) + (1/c1) + (1/b1) + (1/d1)+(1/a2) + (1/c2) + (1/b2) + (1/d2)))}
  return(z)
}



