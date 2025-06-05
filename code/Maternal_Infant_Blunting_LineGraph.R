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

##########################
#Convert Wide to Long

Difference_data = OG

All_long = data.frame()
current_patient = unique(Difference_data$USUBJID)[1]

for (current_patient in unique(Difference_data$USUBJID)) {
  ID_subset = subset(Difference_data, Difference_data$USUBJID == current_patient)
  ID_subset_meta = ID_subset[,c(Metadata_list)]
  ID_subset = ID_subset[,c("Planned_Time", "Sample_Type", "Antibody_Type",  Antigen_list)]
  ID_subset_long = reshape2::melt(data = ID_subset, 
                                  id.vars = c("Planned_Time", "Sample_Type", "Antibody_Type"), 
                                  value = "Antigen")
  ID_subset_meta = merge(ID_subset_meta, ID_subset_long)
  
  All_long = rbind.fill(All_long, ID_subset_meta)
}

head(All_long)


####################################################################################
#Serum 
####################################################################################
All_long_Serum = All_long #subset(All_long, All_long$Sample_Type == "Serum")


All_Antibodies_OG = data.frame()
current_titer = "PT_Titer"; current_treatment = "Tdap"; current_time = "Delivery/Birth"
for(current_titer in c("PT_Titer","FH_Titer","PRN_Titer",
                       "Tetanus_Titer","Diphtheria_Titer")) {
  GMT = data.frame()
  for (current_treatment in unique(All_long_Serum$Treatment_group)  ) {
    for (current_time in unique(All_long_Serum$Planned_Time)) {
      current_time_subset = subset(All_long_Serum, All_long_Serum$Planned_Time == current_time &
                                     All_long_Serum$Treatment_group == current_treatment)
      current_time_subset = subset(current_time_subset, current_time_subset$variable == current_titer)
      
      Serum_M = subset(current_time_subset, current_time_subset$Sample_Type == "Serum" & 
                       current_time_subset$AGEU == 'YEARS')
      Serum_I = subset(current_time_subset, current_time_subset$Sample_Type == "Serum" & 
                         current_time_subset$AGEU != 'YEARS')
      
      Breastmilk_IgA = subset(current_time_subset, current_time_subset$Sample_Type != "Serum" &
                                current_time_subset$Antibody_Type == "IgA")
      Breastmilk_IgG = subset(current_time_subset, current_time_subset$Sample_Type != "Serum" &
                                current_time_subset$Antibody_Type == "IgG")
      
      if (nrow(Serum_M) > 0) {
        SER_GMT_IgG = suppressWarnings(Gmean(Serum_M$value[Serum_M$value>0],
                                             method = c("classic"), conf.level = 0.95,
                                             sides = c("two.sided"), na.rm = T))
      } else {
        SER_GMT_IgG = c(NaN, NaN, NaN)
      }
      
      if (nrow(Breastmilk_IgA) > 0) {
        Maternal_GMT_IgA = suppressWarnings(Gmean(Breastmilk_IgA$value[Breastmilk_IgA$value>0],
                                                  method = c("classic"), conf.level = 0.95,
                                                  sides = c("two.sided"), na.rm = T))
        Maternal_GMT_IgG = suppressWarnings(Gmean(Breastmilk_IgG$value[Breastmilk_IgG$value>0],
                                                  method = c("classic"), conf.level = 0.95,
                                                  sides = c("two.sided"), na.rm = T))
      } else {
        Maternal_GMT_IgA = c(NaN, NaN, NaN)
        Maternal_GMT_IgG = c(NaN, NaN, NaN)
        
      }
      
      
      if (nrow(Serum_I) > 0) {
        SER_GMT_IgG_Infant = suppressWarnings(Gmean(Serum_I$value[Serum_I$value>0],
                                             method = c("classic"), conf.level = 0.95,
                                             sides = c("two.sided"), na.rm = T))
      } else {
        SER_GMT_IgG_Infant = c(NaN, NaN, NaN)
      }
      
      Maternal_GMT_IgA = unname(Maternal_GMT_IgA)
      Maternal_GMT_IgG = unname(Maternal_GMT_IgG)
      SER_GMT_IgG = unname(SER_GMT_IgG)
      SER_GMT_IgG_Infant = unname(SER_GMT_IgG_Infant)
      
      current_row = data.frame("Treatment" = current_treatment, 
                               "Days" = current_time, 
                               "GMT_Mean" = c(Maternal_GMT_IgA[1], Maternal_GMT_IgG[1], 
                                              SER_GMT_IgG[1], SER_GMT_IgG_Infant[1]),
                               "GMT_Lower" = c(Maternal_GMT_IgA[2], Maternal_GMT_IgG[2], 
                                               SER_GMT_IgG[2], SER_GMT_IgG_Infant[2] ),
                               "GMT_Upper" = c(Maternal_GMT_IgA[3], Maternal_GMT_IgG[3], 
                                               SER_GMT_IgG[3], SER_GMT_IgG_Infant[3] ), 
                               "Group" = c("Breastmilk\nIgA", "Breastmilk\nIgG",
                                           "Maternal Serum\nIgG", "Infant Serum\nIgG"))
      GMT = rbind(GMT, current_row)
    }
  }
  unique(GMT$Days)
  GMT = GMT[complete.cases(GMT),]
  
  # GMT$Days = factor(GMT$Days, levels = c("Prior to\nPrenatal\nVaccination", 
  #                                        "1 month post\nPrenatal\nVaccination",
  #                                        0, 42, 70, 130, 180))
  GMT$Measure = gsub("_Titer", "", current_titer)
  All_Antibodies_OG = rbind(All_Antibodies_OG, GMT)
}

All_Antibodies = All_Antibodies_OG
All_Antibodies$Grouping = paste(All_Antibodies$Group, All_Antibodies$Measure, All_Antibodies$Treatment, sep = "_")

Proper_Labels = data.frame("Measure" = c("PT","FH","PRN","Tetanus", "Diphtheria"), 
                           "variable" = c("PT_Titer","FH_Titer","PRN_Titer","Tetanus_Titer", "Diphtheria_Titer"), 
                           "Measure_2" = c("Pertussis Toxin",
                                           "Filamentous Hemagluttinin","Pertactin",
                                           "Tetanus Toxin", "Diptheria Toxin"))
All_Antibodies = merge(All_Antibodies, Proper_Labels)
All_long = merge(All_long, Proper_Labels)

All_Antibodies$Measure_2 = factor(All_Antibodies$Measure_2, 
                                    levels =c("Pertussis Toxin","Filamentous Hemagluttinin","Pertactin",
                                              "Tetanus Toxin", "Diptheria Toxin") )



Days_Labels = data.frame("Days" = c("Day 1","Day 31",
                                       "Delivery/Birth", "Birth+42d","Birth+70d",
                                       "Birth+130d", "Birth+6m"), 
                         "Planned_Time" = c("Day 1","Day 31",
                                    "Delivery/Birth", "Birth+42d","Birth+70d",
                                    "Birth+130d", "Birth+6m"), 
                           "Days.Since.Birth" = c(-90, -60, 0, 42, 70, 130, 180))

All_Antibodies = merge(All_Antibodies, Days_Labels)
All_long = merge(All_long, Days_Labels)
All_long$Grouping = paste(All_long$AGEU, All_long$Days.Since.Birth, All_long$Treatment, sep = "_")





All_Antibodies_TP = All_Antibodies[grepl("Serum", All_Antibodies$Grouping),]

All_Antibodies_TP[All_Antibodies_TP$GMT_Upper == max(All_Antibodies_TP$GMT_Upper), "GMT_Upper"] = 120


Titer_plot_P = ggplot(data = All_Antibodies_TP, 
                    aes(x = Days.Since.Birth, y = GMT_Mean, group = Grouping, 
                        color = Group, linetype = Treatment)) + 
  geom_vline(xintercept = -90, color = "#6F0202", alpha = 0.35, linewidth = 1.6) + 
  geom_vline(xintercept = 42, color = "#9C6EFD", alpha = 0.35, linewidth = 1.6) + 
  geom_vline(xintercept = 70, color = "#9C6EFD", alpha = 0.35, linewidth = 1.6) + 
  geom_vline(xintercept = 100, color = "#9C6EFD", alpha = 0.35, linewidth = 1.6) + 
  geom_hline(data = subset(All_Antibodies_TP, All_Antibodies_TP$Measure %in% c("PT")),
             aes(yintercept = 10), linewidth = 1)+ 
  geom_hline(data = subset(All_Antibodies_TP,  All_Antibodies_TP$Measure %in% c("Tetanus", "Diphtheria")),
             aes(yintercept = 0.1), linewidth = 1)+ 
  geom_line(linewidth = 0.85) + geom_point() + 
  scale_color_manual(name = "Antibodies",
                     breaks = c("Breastmilk\nIgA", "Breastmilk\nIgG", 
                                "Maternal Serum\nIgG", "Infant Serum\nIgG"),
                     values = c("#A3DEFF", "#FB7676", "#6F0202", "#9C6EFD")) +
  scale_linetype_manual(name = "Maternal\nVaccination", 
                        breaks = c("Td", "Tdap"), 
                        values = c("11", "solid")) + 
  facet_wrap(~Measure_2,  scale = "free", nrow = 2) + 
  scale_y_continuous(transform=scales::pseudo_log_trans(base = 10),
                     breaks=c(1, 10, 100, 1000),
                     labels = expression(1, 10, 10^2, 10^3),
                     expand = c(0, 0), 
                     limits = c(-1, 1000)) +
  # geom_errorbar(data = All_Antibodies_TP,
  #               aes(x = Days.Since.Birth, y = GMT_Mean,
  #                   ymin = GMT_Lower, ymax =GMT_Upper,  group = Group),
  #               linetype = "solid", width = 0.05, show.legend = F, linewidth = 0.75) +
  theme_bw() + 
  theme(#legend.text = element_text(size = 12),
        #legend.title = element_text(size = 13), 
        legend.position = "bottom", 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 10,
                                  margin = margin(1,1,1,1, "mm"))) + 
  xlab("Days after Delivery(0)")  + ylab("GMT"); Titer_plot_P

jpeg("./Figures/Infant_Mother_Serum_Longitudinal.jpeg", res = 400, 
     height = 2000, width = 3000)
Titer_plot_P
dev.off()

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
