library(psych)
library(magrittr)
library(MASS)
library(data.table)
library(bursts)
library(tidyverse) #make sure select function is not masked by MASS package
library(gridExtra)
library(brms)
library(viridis)
library(glue)
library(mipfp)


dyn.load('functions/run_sttc2_64.dll') #C++ code to run STTC
source('functions/seqfunctions.R')


df<-read.csv("data/rawdata/dominance.csv", stringsAsFactors=F)


source('code_analysis/01_diff_btw_domsub_over_days.R')
source('code_analysis/02_burst_detection.R')
source('code_analysis/03_burst_detection_figures.R')
source('code_analysis/04_identifying_resolved_domsub_across_bursts.R')
source('code_analysis/05_changes_across_pre_mid_post_resolution.R')
source('code_analysis/06_FOMC_within_individuals.R')
source('code_analysis/07_timed-window_cross-correlation_FSTTC.R')
source('code_analysis/08_save_figures.R')
