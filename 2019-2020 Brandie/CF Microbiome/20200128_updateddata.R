##### setup #####
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(vegan)
setwd("C:/Users/Casey/Desktop/research")
outpatient2020_01_25 <- read.csv("C:/Users/Casey/Desktop/research/20200128plots/outpatient2020_01_25.csv")
outpatient_main = outpatient2020_01_25 %>% 
  select(sid, perc_total, last, visit, lq_all,pcr_neg, seq_count)%>%
  filter(!(sid == "JAR018" | sid == "BAR019" | sid == "NAO032" | sid == "EIL010" ))


# those with multiple visits will have more than 224 samples for a given visit
outpatient_2vv = outpatient_main %>% group_by(sid,visit) %>% filter(n() > 224)
# take the max lq_all for that visit
outpatient_2vv = outpatient_2vv %>% group_by(sid, visit) %>% filter(lq_all == max(lq_all))
# combine this with our other visit
tempdf = outpatient_main %>% group_by(sid,visit) %>% filter(n() <= 224)
outpatient_2vv = rbind(tempdf, outpatient_2vv)
rm(tempdf)
write.csv(outpatient_2vv, "outpatient_2vis.csv")
####################################################

##### spaghetti plots #####
library(data.table)
sublist = c("BAM014", "BSN024", "ESP008", "MEB022", "NDR020", "SKH012", "VAC001", "VCR009", "A_MO25")
long_op = outpatient_2vv %>% select(c("sid", "visit", "lq_all")) %>% filter(!(sid %in% sublist))
wide_op = long_op %>% dcast(sid~visit, fun=mean) 
long_op = wide_op %>% gather(visit, lq_all, `1`:`2`)
lq_all_noodle = ggplot(data = long_op, 
                       aes(x = as.numeric(visit), y = lq_all, group = sid, label = sid)) + 
  geom_line() + 
  stat_summary(aes(group = 1), 
               geom = "line", fun.y = mean, color = "blue") 
#geom_text(aes(label = sid), size = 2)
lq_all_noodle = lq_all_noodle + labs(title = "Change in Log Copies/Reaction Within Subjects",
                                     x = "Visit", 
                                     y = "Log Copies per Reaction")

lq_all_noodle = lq_all_noodle  + geom_label_repel(point.padding = NA,
                                                #direction = "y",
                                                  size = 4)
#print(lq_all_noodle)
ggsave("20200128_spaghetti_plot.png", plot = lq_all_noodle, height = 12, width = 12)
# ==============================================================================#

# ==============================================#

##### mh barchart #####
# load in cleaned data
outpatient2v = read.csv("C:/Users/Casey/Desktop/research/outpatient_2v")
outpatient2v = outpatient2v %>% select(-X)
# load in MH diversity matrix
Outpatient_MH<- read.delim("C:/Users/Casey/Desktop/research/OutpatientAbx_Morisita-Horn_Morisita-Horn.txt")

# -------------------------------------------------------------------------
##### Data Cleaning #####
# get a reference for "name" and "sid"
sid_name_key = outpatient2v %>%select(sid, visit, name)%>% group_by(sid) %>% distinct()

# set the matrix rownames for easier data filtering
rownames(Outpatient_MH) = colnames(Outpatient_MH)

# get the final matrix
MH_final = Outpatient_MH %>% select(levels(sid_name_key$name)) %>%filter(rownames(Outpatient_MH)
                                                                         %in% levels(sid_name_key$name))
rownames(MH_final)=colnames(MH_final)

# create a dictionary-like object
sid_name = as.data.frame(sid_name_key) %>% mutate(sid_vis = paste0(sid,"-",visit)) %>% mutate(sid=NULL, visit=NULL)
sid_name$name = as.character(sid_name$name)
dictionary = sid_name$name
names(dictionary) = sid_name$sid_vis

############## Janky code but it works! 
# take the subid mapping in dictionary
# match it to visits 1 & 2
# find corresponding MH cell in matrix
# return as named list

# still need to delete duplicates so n = 19
y=NULL
for(i in names(dictionary)){
  x = NULL
  a = i
  b = str_split_fixed(a, "-", n=2)
  c = ifelse(grepl(b[1], names(dictionary)),1,0 )
  for(j in 1:38){
    if(c[j]==1){
      p = dictionary[j]
      x = c(x,p)
    } 
  }
  l = MH_final %>% select(x[1]) %>% filter(rownames(MH_final)%in%x[2])
  y = c(y,l)
}
q = as.data.frame(y)

# long format for boxplot
MH_subj = gather(q, key = "Subject", value = "MH-BD",1:38)
MH_subj = MH_subj %>% arrange(Subject) %>% slice(seq(1,37, 2)) 
# get original sid back
firstpart = function(string_vector){
  a = NULL
  for (i in 1:length(string_vector)){
    a[i] = unlist(strsplit(string_vector[i],'.',fixed = TRUE))[1]
  }
  return(a)
}
# final dataframe for boxplot
MH_subj = MH_subj %>% mutate(sid = firstpart(Subject)) %>% rename("MH_BD"=`MH-BD`)%>% select(sid,MH_BD)
write.csv(MH_subj, "MH_subj.csv")
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#### PLOTS #####
MH_subj%>% arrange(MH_subj$MH_BD)

MH_box = ggplot(data = MH_subj, aes(x = "", y = MH_BD, label= sid)) + geom_boxplot(outlier.shape = NA) + geom_point()
MH_box = MH_box + labs(title = "Beta Diversity Within Subjects",
                       subtitle = "Visit 1 - Visit 2",
                       x = "", 
                       y = "Morisita-Horn Beta Diversity")

# outside middle 50%
outer_quart = MH_box + geom_label_repel(data = subset(MH_subj,!between(MH_BD,quantile(MH_BD,.25),
                                                                       quantile(MH_BD, .75))),
                                        direction = "x")

# bottom 25%
lower_quart = MH_box + geom_label_repel(data = subset(MH_subj, MH_BD < quantile(MH_BD, .25)),
                                        direction = "x")
# save plots
ggsave("MH_outer_qtls.png", plot = outer_quart, width = 5, height = 7)
ggsave("MH_lower_qtl.png", plot = lower_quart, width = 5, height = 7)

##### stacked barcharts #####
# dropping those without 2 visits
outpatient_2vv = outpatient_2vv %>%  filter(!(
  sid == "BAM014"|sid=="BSN024"|sid=="ESP008"|sid=="MEB022"|
    sid=="NDR020"|sid=="SKH012"|sid=="VCR009"|sid=="VAC001"|sid=="A_MO25"))

# new sid (MH, viral) MH_subj
name_df = outpatient_2vv %>% filter(visit==1)%>%select(sid, pcr_neg) %>% unique.data.frame()
name_df = full_join(name_df, MH_subj, by = "sid")
name_df = name_df%>% mutate(sid_new = paste0(MH_BD,"_",pcr_neg,"_", sid))
name_df = name_df %>% select(sid, sid_new)
outpatient_2vv = full_join(outpatient_2vv, name_df, by= "sid")
outpatient_cutoff5 = outpatient_2vv %>% filter(perc_total > 5)
#describe(outpatient_cutoff5$last)

outpatient_cutoff5 = outpatient_cutoff5 %>% group_by(last) %>% mutate(pos = cumsum(perc_total) - (0.5)* perc_total)

# stack labels
taxa_sbc_c5 = ggplot(data = outpatient_cutoff5,
                     aes(x = visit.x, y = perc_total, fill = last, label = last)) + 
  geom_bar(stat = "identity") + facet_wrap(~sid_new) +
  theme(legend.position = "none")+ 
  geom_text(size = 3 , position = position_stack(vjust = 0.5)) 

taxa_sbc_c5 = taxa_sbc_c5 + labs(title = "Change in Community Composition Within Subjects",
                                 x = "Visit", 
                                 y = "Total Percent of Community")
ggsave("20200128_sbc_c5.png", plot = taxa_sbc_c5, height = 12, width = 12)
