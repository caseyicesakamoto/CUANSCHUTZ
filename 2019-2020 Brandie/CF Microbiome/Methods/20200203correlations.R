setwd("C:/Repositories/CUANSCHUTZ/2019-2020 Brandie/CF Microbiome")
# correlations
# import T1 data and MH data
##### setup ######
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(vegan)
library(readxl)
cult3 <- read_excel("C:/Users/Casey/Desktop/research/cult3.xlsx") # laptop
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

###### corrs for inflammatory markers added 4/2020####
# did this on desktop
outpatient_2vis <- read.csv("C:/Users/Casey/Desktop/WorkData/Microbiome Brandie 2019-2020/outpatient_2vis.csv")
sublist = c("BAM014", "BSN024", "ESP008", "MEB022", "NDR020", "SKH012", 
            "VAC001", "VCR009", "A_MO25", "JAR018", "BAR019", "NAO032", "EIL010")
outpatient_2vis = outpatient_2vis %>% select(lq_all, sid, visit, elas, IL1B, IL8, HMGB_1,
                                             Cell_Counts_x_104,X_Neutrophils___) %>%
  filter(!(sid %in% sublist))%>% unique.data.frame()

# load difference visit 2 - 1
patient_load = outpatient_2vis %>% select(visit, lq_all, sid) %>%spread(visit, lq_all) %>%
  mutate(load_dif = `2`-`1`) %>% select(-`1`,-`2`)

# inflamm difference
patient_elas  = outpatient_2vis %>% select(visit, elas, sid) %>%spread(visit, elas) %>%
  mutate(elas_dif = `2`-`1`) %>% select(-`1`,-`2`)

patient_IL1B = outpatient_2vis %>% select(visit, IL1B, sid) %>%spread(visit, IL1B) %>%
  mutate(IL1B_dif = `2`-`1`) %>% select(-`1`,-`2`)

patient_IL8 = outpatient_2vis %>% select(visit, IL8, sid) %>%spread(visit, IL8) %>%
  mutate(IL8_dif = `2`-`1`) %>% select(-`1`,-`2`)

patient_HMGB_1 = outpatient_2vis %>% select(visit, HMGB_1, sid) %>%spread(visit,HMGB_1) %>%
  mutate(HMGB_1_dif = `2`-`1`) %>% select(-`1`,-`2`)

patient_cellct = outpatient_2vis %>% select(visit, Cell_Counts_x_104, sid) %>%spread(visit,Cell_Counts_x_104) %>%
  mutate(cellct_dif = `2`-`1`) %>% select(-`1`,-`2`)

patient_neutrophil = outpatient_2vis %>% select(visit, X_Neutrophils___, sid) %>%spread(visit,X_Neutrophils___) %>%
  mutate(neutrophils_dif = `2`-`1`) %>% select(-`1`,-`2`)

# df
# quite a few missing a visit for neutrophils and white cell ct
# subjects typically missing all inflam data or neutrophils or white cell cct
patient_inflam = full_join(patient_load, patient_elas, by = "sid")
patient_inflam = full_join(patient_inflam, patient_IL1B, by = "sid")
patient_inflam = full_join(patient_inflam, patient_IL8, by = "sid")
patient_inflam = full_join(patient_inflam, patient_HMGB_1, by = "sid")
patient_inflam = full_join(patient_inflam, patient_cellct, by = "sid")
patient_inflam = full_join(patient_inflam, patient_neutrophil, by = "sid")

par(mfrow = c(3,2))
load_elas_cor = ggplot(data = patient_inflam) + 
  geom_point(aes(load_dif,elas_dif)) + theme_classic() +
  labs(title = "Change in Log-Copies/Reaction (V2-V1) vs Change in Elastase",
       x = "Change in Log-Copies/Reaction", y = "Change in Elastase")
print(load_elas_cor) # aa007 big outlier

load_IL1B_cor = ggplot(data = patient_inflam) + 
  geom_point(aes(load_dif,IL1B_dif)) + theme_classic() +
  labs(title = "Change in Log-Copies/Reaction (V2-V1) vs Change in IL-1B",
       x = "Change in Log-Copies/Reaction", y = "Change in IL-1B")
print(load_IL1B_cor)# lme026 big outlier

load_IL8_cor = ggplot(data = patient_inflam) + 
  geom_point(aes(load_dif,IL8_dif)) + theme_classic() +
  labs(title = "Change in Log-Copies/Reaction (V2-V1) vs Change in IL-8",
       x = "Change in Log-Copies/Reaction", y = "Change in IL-8")
print(load_IL8_cor) # this one is strange

load_HMGB_cor = ggplot(data = patient_inflam) + 
  geom_point(aes(load_dif,HMGB_1_dif)) + theme_classic() +
  labs(title = "Change in Log-Copies/Reaction (V2-V1) vs Change in HMGB-1",
       x = "Change in Log-Copies/Reaction", y = "Change in HMGB-1")
print(load_HMGB_cor) # JWS CAM outliers

load_cellct_cor = ggplot(data = patient_inflam) + 
  geom_point(aes(load_dif,cellct_dif)) + theme_classic() +
  labs(title = "Change in Log-Copies/Reaction (V2-V1) vs Change in White Blood Cell Ct",
       x = "Change in Log-Copies/Reaction", y = "Change in White Blood Cell ct")
print(load_cellct_cor) # lots of missing here

load_neutro_cor = ggplot(data = patient_inflam) + 
  geom_point(aes(load_dif,neutrophils_dif)) + theme_classic() +
  labs(title = "Change in Log-Copies/Reaction (V2-V1) vs Change in % Neutrophil",
       x = "Change in Log-Copies/Reaction", y = "Change in % Neutrophil")
print(load_neutro_cor) # lots of missing data

# save plots
setwd("D:/Repositories/CUANSCHUTZ/2019-2020 Brandie/CF Microbiome/Output")
ggsave("20200420load_elas_cor.tiff", plot = load_elas_cor, height = 12, width = 12)
ggsave("20200420load_il1b_cor.tiff", plot = load_IL1B_cor, height = 12, width = 12)
ggsave("20200420load_il8_cor.tiff", plot = load_IL8_cor, height = 12, width = 12)
ggsave("20200420load_hmgb_cor.tiff", plot = load_HMGB_cor, height = 12, width = 12)
ggsave("20200420load_cellct_cor.tiff", plot = load_cellct_cor, height = 12, width = 12)
ggsave("20200420load_neutro_cor.tiff", plot = load_neutro_cor, height = 12, width = 12)


### correlations inflam markers
cor(patient_inflam$load_dif,patient_inflam$elas_dif, use = "complete.obs", method = "spearman") # -0.214
cor(patient_inflam$load_dif,patient_inflam$IL1B_dif, use = "complete.obs", method = "spearman") # 0.0330
cor(patient_inflam$load_dif,patient_inflam$IL8_dif, use = "complete.obs", method = "spearman") # 0.2
cor(patient_inflam$load_dif,patient_inflam$HMGB_1_dif, use = "complete.obs", method = "spearman") #0.213
cor(patient_inflam$load_dif,patient_inflam$cellct_dif, use = "complete.obs", method = "spearman")# -0.0714
cor(patient_inflam$load_dif,patient_inflam$neutrophils_dif, use = "complete.obs", method = "spearman") # -0.429

