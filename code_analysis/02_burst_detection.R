# 03. Burst detection

### Calculate Bursts with Kleinberg's algorithm 
dfburstA <-  getrawburstsdf(df) #gamma=0.3
dfburstANA <- addburstNAs(dfburstA) # Add cumulative burst number for each and add in NAs


#frequency of bursts
dfbursts.freq = rbind(
  dfburstANA %>% 
    filter(!is.na(burstval)) %>% 
    select(dyadid, day, burstval) %>% 
    unique %>% 
    mutate(gamma='0.3')) %>% 
  group_by(dyadid,day,gamma) %>% 
  filter(burstval==max(burstval)
         )

# add zeros in
dfbursts.freq3 <- dfbursts.freq %>% filter(gamma==0.3)
emptydf.freq3<-cbind(expand.grid(LETTERS[1:21],1:5),0,"0.3") %>% as.data.frame.matrix
colnames(emptydf.freq3) <- c("dyadid", "day", "burstval", "gamma")

dfbursts.freq3 <- rbind(dfbursts.freq3 %>% 
                          as.data.frame ,
                        emptydf.freq3 %>% 
                          as.data.frame) %>% 
  group_by(dyadid,day) %>% 
  summarise(burstval= sum(burstval))

dfbursts.freq3 %>% 
  group_by(dyadid) %>% 
  summarize(total=sum(burstval))  %>% 
  summarize(median = median(total),
            lqr = quantile(total,0.25), 
            uqr = quantile(total,0.75)) %>% 
  mutate(day="All") -> burst.freq_all

dfbursts.freq3 %>% 
  group_by(day) %>% 
  summarize(median = median(burstval),
            lqr = quantile(burstval,0.25), 
            uqr = quantile(burstval,0.75)) -> burst.freq_days

#durations - df for mixed effects model
dfbursts.dur<-rbind(
  dfburstANA %>% filter(!is.na(burstval)) %>% 
    group_by(dyadid,day,burstval) %>% 
    summarize(starttime = min(time), endtime = max(time)) %>% 
    mutate(burstduration = endtime-starttime, gamma='0.3')
  )


dfbursts.dur %>% 
  group_by(dyadid) %>% 
  summarise(bd =median(burstduration)) %>% 
  ungroup() %>% 
  summarise(median  = median(bd), 
            lqr = quantile(bd,.25), 
            uqr=quantile(bd,.75)) %>% 
  mutate(day="All") -> burst.dur_all

dfbursts.dur %>% 
  group_by(day,dyadid) %>% 
  summarise(bd=median(burstduration)) %>% 
  group_by(day) %>%
  summarise(median  = median(bd), 
            lqr = quantile(bd,.25), 
            uqr=quantile(bd,.75)) -> burst.dur_days



### Proportion of agg/sub behaviors in/out bursts.
propburst = dfburstANA %>% 
  group_by(dyadid,day,state) %>% 
  summarise(outburst = sum(is.na(burstval)), 
            total = n() ) %>% 
  mutate(inburst = total-outburst, propin = inburst/total)


propburst  %>%
  group_by(state) %>%  
  summarize(median = median(propin), 
            lqr=quantile(propin,.25), 
            uqr=quantile(propin,.75)) %>% 
  data.frame %>% 
  mutate(day="All") -> burst.prop_all

propburst %>% 
  group_by(day,state) %>%  
  summarize(median = median(propin), 
            lqr=quantile(propin,.25), 
            uqr=quantile(propin,.75)) %>% 
  data.frame -> burst.prop_days


## Table 2. Descriptive statistics of aggressive and subordinate behavior by day
table2<-left_join(
rbind(burst.freq_days,burst.freq_all) %>% 
  mutate('Burst frequency'=glue("{round(median,0)} [{round(lqr,0)}, {round(uqr,0)}]")) %>% 
  select(day,'Burst frequency'),

rbind(burst.dur_days,burst.dur_all) %>% 
  mutate('Burst duration (s)'=glue("{round(median,1)} [{round(lqr,1)}, {round(uqr,1)}]")) %>% 
  select(day,'Burst duration (s)'), 

by="day") %>% left_join(.,
                        
rbind(burst.prop_days,burst.prop_all) %>% 
  mutate(value=glue("{round(median*100,1)} [{round(lqr*100,1)}, {round(uqr*100,1)}]")) %>% 
  select(day,state,value) %>% 
  spread(state,value) %>% 
  mutate('Proportion of aggressive behaviors in bursts (%)'=aggr,
         'Proportion of subordinate behaviors in bursts (%)'=sub) %>% 
  select(day,'Proportion of aggressive behaviors in bursts (%)','Proportion of subordinate behaviors in bursts (%)'), 
by="day")


### Proportion of Duration of agg/sub behaviors in/out bursts.
propburstd=dfburstANA %>% 
  group_by(dyadid,day,state) %>% 
  summarise(outburst = sum(duration[is.na(burstval)]), 
            total = sum(duration) ) %>% 
  mutate(inburst = total-outburst, propin = inburst/total)

propburstd %>% 
  group_by(day,state) %>% 
  summarize(median = median(propin), 
            lqr=quantile(propin,.25), 
            uqr=quantile(propin,.75)) %>% 
  data.frame


dfburstANA %>% 
  group_by(dyadid,state) %>% 
  summarise(outburst = sum(duration[is.na(burstval)]), 
            total = sum(duration) ) %>% 
  mutate(inburst = total-outburst, 
         propin = inburst/total) %>% 
  group_by(state) %>% 
  summarize(median = median(propin), 
            lqr=quantile(propin,.25), 
            uqr=quantile(propin,.75)) %>% 
  data.frame



#### Burst frequency over days: zero-inflated poisson

# hist(dfbursts.freq3$burstval,breaks=20) #zero-inflated? poisson 


dfbursts.freq3 <-dfbursts.freq3 %>% mutate(burstvalx=as.integer(burstval),burstval=as.integer(burstval),
                                           dayx=as.integer(day)) %>%   as.data.frame()

# fit.burstfreq.brms <- brm(burstvalx ~ dayx + (1|dyadid),
#                           data=dfbursts.freq3,
#                           family=zero_inflated_poisson())
#                                          
# saveRDS(fit.burstfreq.brms,"data/brms_results/fit.burstfreq.brms.RDS")
fit.burstfreq.brms <- readRDS("data/brms_results/fit.burstfreq.brms.RDS")

summary(fit.burstfreq.brms)



#### Effect of Day on Burst Duration
# hist(dfbursts.dur$burstduration,breaks = 100)

# fit.burstdur.brms<-brm(burstduration  ~ day + (1|dyadid),
#                           family=gaussian(link="log"),
#                           data=dfbursts.dur)
# saveRDS(fit.burstdur.brms,"data/brms_results/fit.burstdur.brms.RDS")
fit.burstdur.brms <- readRDS("data/brms_results/fit.burstdur.brms.RDS")

summary(fit.burstdur.brms)


#### The effect of days on the proportion of aggressive and subordinate behaviors transpired within bursts

# fit.burstin.brms <- brm(inburst|trials(total) ~ day +  (1 | dyadid), family = binomial, data = propburst)
# saveRDS(fit.burstin.brms,"data/brms_results/fit.burstin.brms.RDS")
fit.burstin.brms <- readRDS("data/brms_results/fit.burstin.brms.RDS")

summary(fit.burstin.brms)

