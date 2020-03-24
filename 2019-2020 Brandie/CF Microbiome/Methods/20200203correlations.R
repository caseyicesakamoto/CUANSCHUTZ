setwd("C:/Repositories/CUANSCHUTZ/2019-2020 Brandie/CF Microbiome")
# correlations
# import T1 data and MH data
##### setup ######
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(vegan)
library(readxl)
cult3 <- read_excel("C:/Users/Casey/Desktop/research/cult3.xlsx")
patient_t1_df <- read.csv("C:/Users/Casey/Desktop/research/patient_t1_df.csv")
MH_subj <- read.csv("C:/Users/Casey/Desktop/research/MH_subj.csv")
outpatient_2vis <- read.csv("C:/Users/Casey/Desktop/research/outpatient_2vis.csv")
sublist = c("BAM014", "BSN024", "ESP008", "MEB022", "NDR020", "SKH012", 
            "VAC001", "VCR009", "A_MO25", "JAR018", "BAR019", "NAO032", "EIL010")
outpatient_2vis = outpatient_2vis %>% select(lq_all, sid, visit) %>%
  filter(!(sid %in% sublist))%>% unique.data.frame()
cult3 = cult3 %>% filter(!(SID %in% sublist)) %>% mutate(sid = SID) %>% select(-SID)
cult3 = cult3 %>% mutate(changeCFU = after - before) %>% select(sid, changeCFU)
patient_load = outpatient_2vis %>% spread(visit, lq_all) %>%
  mutate(load_dif = `2`-`1`) %>% select(-`1`,-`2`)
MH_subj = MH_subj%>% select(-X)

# data frame
patient_df = full_join(MH_subj, patient_t1_df, by="sid") %>% select(-c("X", "visit.x"))
patient_df = full_join(patient_load, patient_df, by="sid")
patient_df = full_join(patient_df, cult3, by="sid")
rm(MH_subj,patient_t1_df,patient_load, cult3)
#####---------------#######
###### Change in FEV vs change in MH #####
FEV_BD_cor = ggplot(data = patient_df) + 
  geom_point(aes(MH_BD,change_fev_perc_v1v2)) + theme_classic() + 
  coord_cartesian(xlim = c(0.45,1)) + labs(title = "Percent Change in FEV1 (V1-V2) vs MH B-Diversity",
                                           x = "Morisita-Horn Beta-Diversity", y = "%-Change in FEV1")
print(FEV_BD_cor)
# 3 missing a v2
# corr = 0.006
cor(patient_df$MH_BD, patient_df$change_fev_perc_v1v2, use = "complete.obs", method = "spearman")
#####----------------#####
##### Change in FEV vs change in load #####
FEV_load_cor = ggplot(data = patient_df) + 
  geom_point(aes(load_dif,change_fev_perc_v1v2)) + theme_classic() +
  labs(title = "Percent Change in FEV1 (V1-V2) vs Change in Load",
       x = "Change in Log-Copies/Reaction", y = "%-Change in FEV1")
print(FEV_load_cor)
# 3 missing a v2
# corr = -.047
cor(patient_df$load_dif, patient_df$change_fev_perc_v1v2, use = "complete.obs", method = "spearman")

##### change in dominant organism ####
FEV_cfu_cor = ggplot(data = patient_df) + 
  geom_point(aes(changeCFU,change_fev_perc_v1v2)) + theme_classic() +
  labs(title = "Percent Change in FEV1 (V1-V2) vs Change in Dominant Organism",
       x = "Change in Log CFU of Dominant Bacteria", y = "%-Change in FEV1")
print(FEV_cfu_cor)
# 3 missing a v2
# corr = -0.232
# 9 subs missing an obs (10 total points)
cor(patient_df$changeCFU, patient_df$change_fev_perc_v1v2, use = "complete.obs", method = "spearman")

####### save plots
ggsave("20200204FEV_BD_corr.tiff", plot = FEV_BD_cor, height = 12, width = 12)
ggsave("20200204FEV_load_corr.tiff", plot = FEV_load_cor, height = 12, width = 12)
ggsave("20200204FEV_cfu_corr.tiff", plot = FEV_cfu_cor, height = 12, width = 12)



##### alpha diversity #####
# used shannons H
outpatient2020_01_25 <- read.csv("C:/Users/Casey/Desktop/research/20200128plots/outpatient2020_01_25.csv")
shm = outpatient2020_01_25 %>% select(visit,sid, ShannonH_Median, ShannonH_Mean) %>%
  filter(!(sid %in% sublist)&visit==1)%>% unique.data.frame()%>%na.omit()
# similar enough
fivenum(shm$ShannonH_Median);fivenum(shm$ShannonH_Mean)
