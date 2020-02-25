# Setup
# This table has been updated as of 1/30/2020 to reflect the 38 samples
setwd("C:/Users/Casey/Desktop/research")
library(tidyverse)
library(tableone)
outpatient2020_01_25 <- read.csv("C:/Users/Casey/Desktop/research/20200128plots/outpatient2020_01_25.csv")
# age, gender, base fev1% (m r), fev1 pred(mr)
# change pred fev (mr), v1 pes score, days between visit (mr)
# ? genotyp, cf resp culture v1, positie viral pcr
# antibiotic class prescribed
outpatient_main = outpatient2020_01_25 %>% 
  select(sid,molecularid, perc_total, last, visit, lq_all, fev1_perc_6m_prior,sex , pes_total, age_yrs,
         fev1_perc_pre, fev1_perc_post,Antibiotic_class, visit_dif,pcr_neg, genotype1,genotype2 )%>%
  filter(!(sid == "JAR018" | sid == "BAR019" | sid == "NAO032" | sid == "EIL010" ))

outpatient_2vv = outpatient_main %>% group_by(sid,visit) %>% filter(n() > 224)
# take the max lq_all
outpatient_2vv = outpatient_2vv %>% group_by(sid, visit) %>% filter(lq_all == max(lq_all))

# combine this with our other patient data
tempdf = outpatient_main %>% group_by(sid,visit) %>% filter(n() <= 224)
outpatient_2vv = rbind(tempdf, outpatient_2vv)
rm(tempdf)
# subjects not in our 19
sublist = c("BAM014", "BSN024", "ESP008", "MEB022", "NDR020", "SKH012", "VAC001", "VCR009", "A_MO25")
outpatient_2vv = outpatient_2vv %>%  filter(!(sid %in% sublist))

summary(outpatient_2vv)

outpatient_2vv = outpatient_2vv %>% mutate(
  genotype = case_when(genotype1 == "F508del" & genotype2 == "F508del" ~ "F508del/F508del",
                       genotype1 == "F508del" & genotype2 != "F508del" ~ "F508del/Other",
                       visit==1 & sid =="SLP006"~"Other"))
# observations for table 1
unique_df = outpatient_2vv %>% select(-c(last, perc_total)) %>% unique.data.frame()
# sid and mid
molecular_id = unique_df[,1:2]
write.csv(molecular_id, "molecularid.csv")
#
df_v1 = unique_df %>% filter(visit==1) %>% select(-visit_dif, -visit)
df_v2 = unique_df %>% filter(visit==2) 
df_v2 = df_v2 %>%  select(sid, visit_dif)

patient_t1_df = full_join(df_v1, df_v2, by = 'sid')
patient_t1_df = patient_t1_df %>% select(-visit.x, -visit.y)
rm(df_v1, df_v2, unique_df)
patient_t1_df = patient_t1_df %>% mutate(change_fev_perc = fev1_perc_pre - fev1_perc_6m_prior,
                                         change_fev_perc_v1v2 = fev1_perc_post - fev1_perc_pre)
patient_t1_df$age_yrs = floor(patient_t1_df$age_yrs)
# create table 1
all_var = c("age_yrs", "sex","genotype", "fev1_perc_6m_prior", "fev1_perc_pre", "change_fev_perc", "pes_total", "Antibiotic_class", "visit_dif", "pcr_neg")
med_iqr_vars = c("age_yrs", "fev1_perc_6m_prior","fev1_perc_pre", "change_fev_perc", "pes_total", "visit_dif")
cat_vars = c("sex","genotype", "antibiotic_class","pcr_neg")
t1 = CreateTableOne(data = patient_t1_df, vars = all_var, factorVars = cat_vars)
final = print(t1, nonnormal = med_iqr_vars, showAllLevels = TRUE, minMax = TRUE)
write.csv(final, "table1.csv")


# looking at lq_all
unique_df = outpatient_2vv %>% select(-c(last, perc_total)) %>% unique.data.frame()
df_v1 = unique_df %>% filter(visit==1) %>% select(-visit_dif, -visit)
df_v2 = unique_df %>% filter(visit==2) 

df_v1 = df_v1 %>% select(sid, lq_all)
df_v2 = df_v2 %>% select(sid, lq_all)
load_df = full_join(df_v1,df_v2,by='sid') # x is visit 1, y is v2
load_df = load_df %>% select(sid, lq_all.x, lq_all.y) %>% mutate(dif = lq_all.y - lq_all.x)
mean(load_df$dif); sd(load_df$dif)
hist(df_v1$lq_all)
hist(df_v2$lq_all)

# wilcoxon signed rank test
wilcox.test(df_v1$lq_all, df_v2$lq_all, paired = TRUE)
# 

rm(df_v1,df_v2,unique_df)

write.csv(patient_t1_df, "patient_t1_df.csv")
