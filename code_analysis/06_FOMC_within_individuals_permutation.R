#### Permutation Tests - ensuring no repeats in permuted sequences ----

domprecodes1 <- permseq_norpt(domprecodes[[1]],1000)
domprecodes2 <- permseq_norpt(domprecodes[[2]],1000)
domprecodes3 <- permseq_norpt(domprecodes[[3]],1000)
domprecodes4 <- permseq_norpt(domprecodes[[4]],1000)
domprecodes5 <- permseq_norpt(domprecodes[[5]],1000)
domprecodes6 <- permseq_norpt(domprecodes[[6]],1000)
domprecodes7 <- permseq_norpt(domprecodes[[7]],1000)
domprecodes8 <- permseq_norpt(domprecodes[[8]],1000)
domprecodes9 <- permseq_norpt(domprecodes[[9]],1000)
domprecodes10 <- permseq_norpt(domprecodes[[10]],1000)
domprecodes11 <- permseq_norpt(domprecodes[[11]],1000)
##domprecodes12 <- permseq_norpt(domprecodes[[12]],1000) #not 1000 possibles.
domprecodes13 <- permseq_norpt(domprecodes[[13]],1000)
domprecodes14 <- permseq_norpt(domprecodes[[14]],1000)
domprecodes15 <- permseq_norpt(domprecodes[[15]],1000)
domprecodes16 <- permseq_norpt(domprecodes[[16]],1000)
domprecodes17 <- permseq_norpt(domprecodes[[17]],1000)
domprecodes18 <- permseq_norpt(domprecodes[[18]],1000)
domprecodes19 <- permseq_norpt(domprecodes[[19]],1000)
domprecodes20 <- permseq_norpt(domprecodes[[20]],1000)
domprecodes21 <- permseq_norpt(domprecodes[[21]],1000)
domprecodes.permed <- list(domprecodes1,domprecodes2,domprecodes3,domprecodes4,domprecodes5,
                           domprecodes6,domprecodes7,domprecodes8,domprecodes9,domprecodes10,
                           domprecodes11,#domprecodes12,
                           domprecodes13,domprecodes14,domprecodes15,
                           domprecodes16,domprecodes17,domprecodes18,domprecodes19,domprecodes20,
                           domprecodes21)

saveRDS(domprecodes.permed, file="data/permutation/domprecodespermed1.RData")

 
subpostcodes3 <- permseq_norpt(subpostcodes[[3]],1000)
subpostcodes4 <- permseq_norpt(subpostcodes[[4]],1000)
subpostcodes5 <- permseq_norpt(subpostcodes[[5]],1000)
subpostcodes6 <- permseq_norpt(subpostcodes[[6]],1000)
subpostcodes7 <- permseq_norpt(subpostcodes[[7]],1000)
subpostcodes8 <- permseq_norpt(subpostcodes[[8]],1000)
subpostcodes9 <- permseq_norpt(subpostcodes[[9]],1000)
subpostcodes10 <- permseq_norpt(subpostcodes[[10]],1000)
subpostcodes11 <- permseq_norpt(subpostcodes[[11]],1000)
subpostcodes12 <- permseq_norpt(subpostcodes[[12]],1000)
subpostcodes13 <- permseq_norpt(subpostcodes[[13]],1000)
subpostcodes14 <- permseq_norpt(subpostcodes[[14]],1000)
subpostcodes15 <- permseq_norpt(subpostcodes[[15]],1000)
subpostcodes16 <- permseq_norpt(subpostcodes[[16]],1000)
subpostcodes17 <- permseq_norpt(subpostcodes[[17]],1000)
subpostcodes18 <- permseq_norpt(subpostcodes[[18]],1000)
subpostcodes19 <- permseq_norpt(subpostcodes[[19]],1000)
subpostcodes20 <- permseq_norpt(subpostcodes[[20]],1000)
subpostcodes.permed <- list(subpostcodes1,subpostcodes2,subpostcodes3,subpostcodes4,subpostcodes5,
                            subpostcodes6,subpostcodes7,subpostcodes8,subpostcodes9,subpostcodes10,
                            subpostcodes11,subpostcodes12,subpostcodes13,subpostcodes14,subpostcodes15,
                            subpostcodes16,subpostcodes17,subpostcodes18,subpostcodes19,subpostcodes20)

saveRDS(subpostcodes.permed, file="data/permutation/subpostcodespermed.RData")




#### Analyze permuted sequences - if don't want to run above, can bring in data from RData files:

domprecodespermed <- readRDS("data/permutation/domprecodespermed1.RData")
dompostcodespermed <- readRDS("data/permutation/dompostcodespermed.RData")
subprecodespermed <- readRDS("data/permutation/subprecodespermed1.RData")
subpostcodespermed <- readRDS("data/permutation/subpostcodespermed.RData")

# no dyad 18 R in POST--
# no dyad 12 L in PRE--

#check perms worked:
sum(unlist(lapply(domprecodespermed, function(z) z[[3]])))   #20000
sum(unlist(lapply(domprecodespermed, function(z) z[[2]]==0))) #20000
sum(unlist(lapply(subprecodespermed, function(z) z[[3]])))   #20000
sum(unlist(lapply(subprecodespermed, function(z) z[[2]]==0))) #20000

sum(unlist(lapply(dompostcodespermed, function(z) z[[3]])))   #20000
sum(unlist(lapply(dompostcodespermed, function(z) z[[2]]==0))) #20000
sum(unlist(lapply(subpostcodespermed, function(z) z[[3]])))   #20000
sum(unlist(lapply(subpostcodespermed, function(z) z[[2]]==0))) #20000


## calculate pvalues
domprecodes.obsvals <- lapply(domprecodes.mats, reshape2::melt)
subprecodes.obsvals <- lapply(subprecodes.mats, reshape2::melt)
dompostcodes.obsvals <- lapply(dompostcodes.mats, reshape2::melt)
subpostcodes.obsvals <- lapply(subpostcodes.mats, reshape2::melt)

domprecodes.obsvals <- lapply(domprecodes.obsvals, setNames, c("Var1","Var2","Observed"))
subprecodes.obsvals <- lapply(subprecodes.obsvals, setNames, c("Var1","Var2","Observed"))
dompostcodes.obsvals <- lapply(dompostcodes.obsvals, setNames, c("Var1","Var2","Observed"))
subpostcodes.obsvals <- lapply(subpostcodes.obsvals, setNames, c("Var1","Var2","Observed"))

goupper<-function(df){
  df$Var1<-toupper(df$Var1)
  df$Var2<-toupper(df$Var2)
  return(df)
}

subprecodes.obsvals <- lapply(subprecodes.obsvals, goupper)
subpostcodes.obsvals <- lapply(subpostcodes.obsvals, goupper)

saveRDS(domprecodes.obsvals,"data/permutation/domprecodes.obsvals.RDS")
saveRDS(dompostcodes.obsvals,"data/permutation/dompostcodes.obsvals.RDS")
saveRDS(subprecodes.obsvals,"data/permutation/subprecodes.obsvals.RDS")
saveRDS(subpostcodes.obsvals,"data/permutation/subpostcodes.obsvals.RDS")


library(doParallel)
detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

## 1:20 as no L pre, no R post.
domprecodespermed.vals = list(); foreach(i=1:20)  %do%  {domprecodespermed.vals[[i]] = lapply(domprecodespermed[[i]][[1]], function(x) reshape2::melt(obsmatrix( collapse(x),prob=F)))}
subprecodespermed.vals = list(); foreach(i=1:20)  %do%  {subprecodespermed.vals[[i]] = lapply(subprecodespermed[[i]][[1]], function(x) reshape2::melt(obsmatrix( collapse(x),prob=F)))}
dompostcodespermed.vals = list(); foreach(i=1:20)  %do%  {dompostcodespermed.vals[[i]] = lapply(dompostcodespermed[[i]][[1]], function(x) reshape2::melt(obsmatrix( collapse(x),prob=F)))}
subpostcodespermed.vals = list();  foreach(i=1:20)  %do%  {subpostcodespermed.vals[[i]] = lapply(subpostcodespermed[[i]][[1]], function(x) reshape2::melt(obsmatrix( collapse(x),prob=F)))}

stopCluster(cl)

# Make permutations into list of data.frames
domprecodespermed.vals.df = lapply(domprecodespermed.vals, function(x) do.call('rbind',Map(cbind, x, PermNo = 1:1000)))
subprecodespermed.vals.df = lapply(subprecodespermed.vals, function(x) do.call('rbind',Map(cbind, x, PermNo = 1:1000)))
dompostcodespermed.vals.df = lapply(dompostcodespermed.vals, function(x) do.call('rbind',Map(cbind, x, PermNo = 1:1000)))
subpostcodespermed.vals.df = lapply(subpostcodespermed.vals, function(x) do.call('rbind',Map(cbind, x, PermNo = 1:1000)))


# add in observed value - "observed" from original matrices

domprecodespermed.vals.df1 <- map2(domprecodespermed.vals.df, domprecodes.obsvals[-12], left_join)
subprecodespermed.vals.df1 <- map2(subprecodespermed.vals.df, subprecodes.obsvals[-12], left_join)
dompostcodespermed.vals.df1 <- map2(dompostcodespermed.vals.df, dompostcodes.obsvals, left_join)
subpostcodespermed.vals.df1 <- map2(subpostcodespermed.vals.df, subpostcodes.obsvals, left_join)


domprecodespermed.vals.df1 <- Map(cbind,domprecodespermed.vals.df1,dyadid=LETTERS[1:21][-12]) #NO L
subprecodespermed.vals.df1 <- Map(cbind,subprecodespermed.vals.df1,dyadid=LETTERS[1:21][-12]) #NO L
dompostcodespermed.vals.df1 <- Map(cbind,dompostcodespermed.vals.df1,dyadid=LETTERS[1:21][-18]) #NO R
subpostcodespermed.vals.df1 <- Map(cbind,subpostcodespermed.vals.df1,dyadid=LETTERS[1:21][-18]) #NO R

saveRDS(domprecodespermed.vals.df1,"data/permutation/domprepermdf_new.RData")
saveRDS(subprecodespermed.vals.df1,"data/permutation/subprepermdf_new.RData")
saveRDS(dompostcodespermed.vals.df1,"data/permutation/dompostpermdf.RData")
saveRDS(subpostcodespermed.vals.df1,"data/permutation/subpostpermdf.RData")
