##### setup #####
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(vegan)
# setwd("C:/Users/Casey/Desktop/research") for laptop
# outpatient2020_01_25 <- read.csv("C:/Users/Casey/Desktop/research/20200128plots/outpatient2020_01_25.csv")
outpatient_main = outpatient2020_01_25 %>% 
  select(sid, perc_total, last, visit, lq_all,pcr_neg, seq_count, name)%>%
  filter(!(sid == "JAR018" | sid == "BAR019" | sid == "NAO032" | sid == "EIL010" ))


# those with multiple visits will have more than 224 samples for a given visit
outpatient_2vv = outpatient_main %>% group_by(sid,visit) %>% filter(n() > 224)
# take the max lq_all for that visit
outpatient_2vv = outpatient_2vv %>% group_by(sid, visit) %>% filter(lq_all == max(lq_all))
# combine this with our other visit
tempdf = outpatient_main %>% group_by(sid,visit) %>% filter(n() <= 224)
outpatient_2vv = rbind(tempdf, outpatient_2vv)
rm(tempdf)
outpatient_2vv = outpatient_2vv %>%  filter(!(
  sid == "BAM014"|sid=="BSN024"|sid=="ESP008"|sid=="MEB022"|
    sid=="NDR020"|sid=="SKH012"|sid=="VCR009"|sid=="VAC001"|sid=="A_MO25"))

# for use in correlations dataset
# setwd("C:/Users/Casey/Desktop/WorkData/Microbiome Brandie 2019-2020")# for desktop
write.csv(outpatient_2vv, "outpatient_2vis.csv")
# ============================================== #


