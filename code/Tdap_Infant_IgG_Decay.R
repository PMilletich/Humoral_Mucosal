#---- Script Metadata #----
# Title: Maternal_Infant_Breastmilk_5.1.R
# Author: Trish Milletich, PhD
# Date: 2024-10-10
#--------------------------

setwd("C:/Users/PMilletich/OneDrive - University of Maryland School of Medicine/Desktop/PATTI/HuMo/")
#setwd("~/Downloads/")

library(ggplot2)
library(ggpubr)
library(tidyr)
library(DescTools)
library(plyr)


Antigen_list = c("PT_Titer","FH_Titer","PRN_Titer", "Tetanus_Titer", "Diphtheria_Titer")
Metadata_list = c("Treatment_group", "Antibody_Type", "Sample_Type", "PARYST", 
                  "USUBJID", "Planned_Time", "RSUBJID",
                  "BMI", "AGE", "AGEU", "IWGTKG", "PREGNN")


OG = read.csv("../Mother_Infant_All_3.csv")
OG = OG[,c(Antigen_list, Metadata_list)]


Infant = subset(OG, OG$AGEU != "YEARS")
Maternal = subset(OG, OG$AGEU == "YEARS")

#####################################
#Longitudinal consecutive change
#####################################
BM = subset(Maternal, Maternal$Sample_Type != "Serum")
All_Antibodies = data.frame()
current_titer = "PT_Titer"
for(current_titer in c("PT_Titer","FH_Titer","PRN_Titer", 
                       "Tetanus_Titer", "Diptheria_Titer")) {
  All_long_Serum = data.frame()
  current_patient = unique(BM$USUBJID)[1]
  Only_1 = c()
  for (current_patient in unique(BM$USUBJID)) {
    ID_subset = subset(BM, BM$USUBJID == current_patient)
    ID_subset = ID_subset[,c("Days", "Antibody_Type", "Sample_Type", "Treatment_group", current_titer)]
    ID_subset_wide = spread(data = ID_subset, 
                            key = "Antibody_Type", 
                            value = current_titer)
    ID_subset_wide$ID = current_patient
    if (0 %in% ID_subset_wide$Days) {
      IgA_Base = ID_subset_wide[ID_subset_wide$Days == 0,"IgA"]
      IgG_Base = ID_subset_wide[ID_subset_wide$Days == 0,"IgG"]
      
      ID_subset_wide$IgA_FR = ID_subset_wide$IgA/IgA_Base
      ID_subset_wide$IgG_FR = ID_subset_wide$IgG/IgG_Base
      All_long_Serum = rbind(All_long_Serum, ID_subset_wide)
    } else {
      Only_1 = c(Only_1, current_patient)
    }
  }
  
  GMT = data.frame()
  current_time = unique(All_long_Serum$Days)[1]
  for (current_time in unique(All_long_Serum$Days)) {
    current_time_subset = subset(All_long_Serum, All_long_Serum$Days == current_time)
    
    Breastmilk = subset(current_time_subset, current_time_subset$Sample_Type != "Serum")
    
    
    BM_GMT_IgA = suppressWarnings(Gmean(Breastmilk$IgA[Breastmilk$IgA>0],
                                        method = c("classic"), conf.level = 0.95,
                                        sides = c("two.sided"), na.rm = T))
    BM_GMT_IgG = suppressWarnings(Gmean(Breastmilk$IgG[Breastmilk$IgG>0],
                                        method = c("classic"), conf.level = 0.95,
                                        sides = c("two.sided"), na.rm = T))
    
    current_row = data.frame("Days" = current_time, 
                             "GMT_Mean" = c(BM_GMT_IgA[1], BM_GMT_IgG[1]),
                             "GMT_Lower" = c(BM_GMT_IgA[2], BM_GMT_IgG[2]),
                             "GMT_Upper" = c(BM_GMT_IgA[3], BM_GMT_IgG[3]), 
                             "Group" = c("Breastmilk:IgA", "Breastmilk:IgG"))
    GMT = rbind(GMT, current_row)
  }
  
  
  GMT$Days = gsub("Birth", "Birth\n", GMT$Days)
  unique(GMT$Days)
  
  
  GMT$Days = factor(GMT$Days, levels = c("Prior to\nPrenatal\nVaccination", 
                                         "1 month post\nPrenatal\nVaccination",
                                         0, 42, 70, 130, 180))
  
  BM_Box = BM
  BM_Box$Titer = BM_Box[,current_titer]
  
  compare_list = list(c("0","42"), c("42","70"), c("70", "130"), c("130", "180"))
  
  Titer_plot =  ggplot(data = BM_Box, aes(x = Days, y = log10(Titer),
                                          fill = Antibody_Type))+
    theme_bw() + facet_wrap(~Treatment_group, nrow = 2) + 
    scale_fill_manual( name = 'Titer',  breaks = c("IgA", "IgG"),
                       values = c("darkblue", "orange")) +
    geom_boxplot(position = position_dodge(0.75), width = 0.75, color= "black") + 
    
    stat_compare_means(data = BM_Box[BM_Box$Antibody_Type == "IgA",],
                       comparisons = compare_list, tip.length = 0, 
                       step.increase = 0.05,size = 3, 
                       color = "darkblue", 
                       method = "wilcox.test") + 
    
    stat_compare_means(data = BM_Box[BM_Box$Antibody_Type == "IgG",],
                       comparisons = compare_list, tip.length = 0, 
                       step.increase = 0.05, vjust = 2, size = 3, 
                       color = "darkorange1", 
                       method = "wilcox.test") + 
    
    xlab("Days after Delivery(0)") + ylab("log10(Titer)") + 
    ggtitle(paste(gsub("_Titer", "", current_titer), "(IU/mL)")); Titer_plot
  
  GMT$Measure = gsub("_Titer", "", current_titer)
  All_Antibodies = rbind(All_Antibodies, GMT)
  assign(paste(current_titer, "_plot", sep = ""), Titer_plot)
  print(paste(current_titer, "_plot", sep = ""))
}


All_Proportions = ggarrange(PT_Titer_plot, FH_Titer_plot, PRN_Titer_plot, 
                            Tetanus_Titer_plot, Diptheria_Titer_plot,
                            ncol = 5, common.legend = T, align = "hv", 
                            legend = "bottom"); All_Proportions

jpeg("~/Desktop/BM_GMT_Stat_5.jpeg", res = 400, height = 3500, width = 8000)
All_Proportions
dev.off()


#####################################
#Boostrix vs TD
#####################################
BM = subset(Maternal, Maternal$Sample_Type != "Serum")
All_Antibodies = data.frame()
current_titer = "PT_Titer"
for(current_titer in c("PT_Titer","FH_Titer","PRN_Titer", 
                       "Tetanus_Titer", "Diptheria_Titer")) {
  All_long_Serum = data.frame()
  current_patient = unique(BM$USUBJID)[1]
  Only_1 = c()
  for (current_patient in unique(BM$USUBJID)) {
    ID_subset = subset(BM, BM$USUBJID == current_patient)
    ID_subset = ID_subset[,c("Days", "Antibody_Type", "Sample_Type", "Treatment_group", current_titer)]
    ID_subset_wide = spread(data = ID_subset, 
                            key = "Antibody_Type", 
                            value = current_titer)
    ID_subset_wide$ID = current_patient
    if (0 %in% ID_subset_wide$Days) {
      IgA_Base = ID_subset_wide[ID_subset_wide$Days == 0,"IgA"]
      IgG_Base = ID_subset_wide[ID_subset_wide$Days == 0,"IgG"]
      
      ID_subset_wide$IgA_FR = ID_subset_wide$IgA/IgA_Base
      ID_subset_wide$IgG_FR = ID_subset_wide$IgG/IgG_Base
      All_long_Serum = rbind(All_long_Serum, ID_subset_wide)
    } else {
      Only_1 = c(Only_1, current_patient)
    }
  }
  
  GMT = data.frame()
  current_time = unique(All_long_Serum$Days)[1]
  for (current_time in unique(All_long_Serum$Days)) {
    current_time_subset = subset(All_long_Serum, All_long_Serum$Days == current_time)
    
    Breastmilk = subset(current_time_subset, current_time_subset$Sample_Type != "Serum")
    
    
    BM_GMT_IgA = suppressWarnings(Gmean(Breastmilk$IgA[Breastmilk$IgA>0],
                                        method = c("classic"), conf.level = 0.95,
                                        sides = c("two.sided"), na.rm = T))
    BM_GMT_IgG = suppressWarnings(Gmean(Breastmilk$IgG[Breastmilk$IgG>0],
                                        method = c("classic"), conf.level = 0.95,
                                        sides = c("two.sided"), na.rm = T))
    
    current_row = data.frame("Days" = current_time, 
                             "GMT_Mean" = c(BM_GMT_IgA[1], BM_GMT_IgG[1]),
                             "GMT_Lower" = c(BM_GMT_IgA[2], BM_GMT_IgG[2]),
                             "GMT_Upper" = c(BM_GMT_IgA[3], BM_GMT_IgG[3]), 
                             "Group" = c("Breastmilk:IgA", "Breastmilk:IgG"))
    GMT = rbind(GMT, current_row)
  }
  
  
  GMT$Days = gsub("Birth", "Birth\n", GMT$Days)
  unique(GMT$Days)
  
  
  GMT$Days = factor(GMT$Days, levels = c("Prior to\nPrenatal\nVaccination", 
                                         "1 month post\nPrenatal\nVaccination",
                                         0, 42, 70, 130, 180))
  
  BM_Box = BM
  BM_Box$Titer = BM_Box[,current_titer]
  compare_list = list(c("0","42"), c("42","70"), c("70", "130"), c("130", "180"))
  
  Titer_plot =  ggplot(data = BM_Box, aes(x = Days, y = log10(Titer),
                                          fill = Treatment_group))+
    theme_bw() + facet_wrap(~Antibody_Type, nrow = 2) + 
    scale_fill_manual( name = 'Titer',  breaks = c("BOOSTRIX", "Td"),
                       values = c("magenta", "darkgreen")) +
    geom_boxplot(position = position_dodge(0.75), width = 0.75, color= "black") + 
    
    stat_compare_means(aes(label = paste0(after_stat(p.format))), 
                       #tip.length = 0, step.increase = 0.05,
                       size = 3,method = "wilcox.test") + 
    xlab("Days after Delivery(0)") + ylab("log10(Titer)") + 
    ggtitle(paste(gsub("_Titer", "", current_titer), "(IU/mL)")); Titer_plot
  
  GMT$Measure = gsub("_Titer", "", current_titer)
  All_Antibodies = rbind(All_Antibodies, GMT)
  assign(paste(current_titer, "_plot", sep = ""), Titer_plot)
  print(paste(current_titer, "_plot", sep = ""))
}


All_Proportions = ggarrange(PT_Titer_plot, FH_Titer_plot, PRN_Titer_plot,
                            ncol = 3, common.legend = T, align = "hv", 
                            legend = "bottom")

All_Proportions

jpeg("~/Desktop/Boostrix_TD_5.jpeg", res = 400, height = 3500, width = 6000)
All_Proportions
dev.off()
