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
library(ggh4x)



Antigen_list = c("PT_Titer","FH_Titer","PRN_Titer", "Tetanus_Titer", "Diphtheria_Titer")
Metadata_list = c("Treatment_group", "Antibody_Type", "Sample_Type", "PARYST", 
                  "USUBJID", "Planned_Time", "RSUBJID",
                  "BMI", "AGE", "AGEU", "IWGTKG", "PREGNN")


OG = read.csv("../Mother_Infant_All_4.csv")
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

All_Antibodies_Updated = All_Antibodies_OG

All_Antibodies_OG$Grouping = paste(All_Antibodies_OG$Group, 
                                   All_Antibodies_OG$Measure, 
                                   All_Antibodies_OG$Treatment, sep = "_")

Proper_Labels = data.frame("Measure" = c("PT","FH","PRN","Tetanus", "Diphtheria"), 
                           "variable" = c("PT_Titer","FH_Titer","PRN_Titer","Tetanus_Titer", "Diphtheria_Titer"), 
                           "Measure_2" = c("Pertussis\nToxin",
                                           "Filamentous\nHemagglutinin","Pertactin",
                                           "Tetanus", "Diphtheria"))
All_Antibodies_OG = merge(All_Antibodies_OG, Proper_Labels)

All_Antibodies_OG$Measure_2 = factor(All_Antibodies_OG$Measure_2, 
                                  levels =c("Filamentous\nHemagglutinin","Pertactin","Pertussis\nToxin",
                                            "Diphtheria", "Tetanus") )

Days_Labels = data.frame("Days" = c("Day 1","Day 31",
                                    "Delivery/Birth", "Birth+42d","Birth+70d",
                                    "Birth+130d", "Birth+6m"), 
                         "Planned_Time" = c("Day 1","Day 31",
                                            "Delivery/Birth", "Birth+42d","Birth+70d",
                                            "Birth+130d", "Birth+6m"), 
                         "Days.Since.Birth" = c(-90, -60, 0, 42, 70, 130, 180))

All_Antibodies_OG = merge(All_Antibodies_OG, Days_Labels)

# All_long = merge(All_long, Proper_Labels)
# All_long = merge(All_long, Days_Labels)
# All_long$Grouping = paste(All_long$AGEU, All_long$Days.Since.Birth, All_long$Treatment, sep = "_")

# All_Antibodies = subset(All_Antibodies_OG, grepl("Infant", All_Antibodies_OG$Group ) == F)
# All_Antibodies = subset(All_Antibodies, All_Antibodies$Group != "Breastmilk\nIgA")
# All_Antibodies$Group= gsub("Breastmilk\nIgG", "HBM IgG", All_Antibodies$Group)
# All_Antibodies$Group= gsub("\n", " ", All_Antibodies$Group)

# All_long_1 = subset(All_long, All_long$Antibody_Type == "IgG" & All_long$AGEU == "YEARS" & 
#                     All_long$Days.Since.Birth >= 0) 
# 
# All_long_1$Group = ifelse(All_long_1$Sample_Type == "Serum", "Maternal Serum IgG", "HBM IgG")
# 
# All_long_1_180 = data.frame(table(All_long_1$USUBJID))
# All_long_1_180 = subset(All_long_1_180, All_long_1_180$Freq == 30)
# table(All_long_1_180$Freq)
# All_long_1 = subset(All_long_1, All_long_1$USUBJID %in% All_long_1_180$Var1)

write.csv(All_Antibodies_OG, "Longitudinal_GMTs_SampleTypes.csv", row.names = F)


Pertussis = ggplot(data = All_Antibodies_OG[grepl("Infant", All_Antibodies_OG$Group) == F,], 
                      aes(x = Days.Since.Birth, y = GMT_Mean, group = Grouping, 
                          color = Group, linetype = Treatment)) + 
  geom_hline(data = All_Antibodies_OG[grepl("Infant", All_Antibodies_OG$Group) == F &
                                 All_Antibodies_OG$Measure %in% c("PT"),], 
             aes(yintercept = 10), linewidth = 1, color = "grey") + 
  geom_hline(data = All_Antibodies_OG[grepl("Infant", All_Antibodies_OG$Group) == F &
                                        All_Antibodies_OG$Measure %in% c("Tetanus", "Diphtheria"),], 
             aes(yintercept = 0.10), linewidth = 1, color = "grey") + 
  
  geom_errorbar(data = All_Antibodies_OG[grepl("Infant", All_Antibodies_OG$Group) == F,],
                aes(x = Days.Since.Birth, y = GMT_Mean,
                    ymin = GMT_Lower, ymax =GMT_Upper,  group = Group),
                linetype = "solid", width = 0.05, show.legend = F, linewidth = 0.75) +
  geom_line(linewidth = 1) + 
  geom_point() + 
  scale_color_manual(name = "Antibodies",
                     breaks = c("Maternal Serum\nIgG", "Breastmilk\nIgG", "Breastmilk\nIgA"),
                     values = c("darkgreen", "#FB7676", "#81C2E6")) +
  scale_linetype_manual(name = "Maternal\nVaccination", 
                        breaks = c("Td", "Tdap"), 
                        values = c("11", "solid")) + 
  scale_y_continuous(transform=scales::pseudo_log_trans(base = 10),
                     breaks=c(1, 10, 100, 1000),
                     labels = expression(1, 10, 10^2, 10^3),
                     expand = c(0, 0), 
                     limits = c(-1, 1000)) +
  scale_x_continuous(breaks = c(-90, -60, 0, 42,  70, 100, 130, 180 ),
                     label = c("MV", "PMV", 0, 42,  70, 100, 130, 180)) +
  facet_wrap(~Measure_2, ncol = 5) + 
  theme_bw() + 
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        axis.text.x = element_text(size = 8),
        legend.position = "bottom", 
        legend.margin=margin(c(0,0,0,0)), 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 12,
                                  margin = margin(1,1,1,1, "mm"))) + 
  xlab("Days after Delivery(0)")  + ylab("GMT"); Pertussis


jpeg("./Figures/Longitudinal_Line_HBM_Serum_Clean.jpeg", res = 400, 
     height = 2000, width = 5000)
Pertussis
dev.off()




