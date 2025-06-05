#---- Script Metadata #----
# Title: IgG_Decay.R
# Author: Trish Milletich, PhD
# Date: 2025-05-05
#--------------------------
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

setwd("C:/Users/PMilletich/OneDrive - University of Maryland School of Medicine/Desktop/PATTI/HuMo/")

OG = read.csv("../Mother_Infant_All_3.csv")
OG = OG[,c(Antigen_list, Metadata_list)]
#####################################################################################
# Decay Comparison 
#####################################################################################
Maternal = subset(OG, OG$AGEU == "YEARS" & 
                    OG$Planned_Time %in% c("Delivery/Birth", "Birth+6m") &
                    OG$Antibody_Type == "IgG")

Fold_Difference = data.frame()
current_ID = unique(Maternal$USUBJID)[1]
current_ID = "CPA.05788"
for (current_ID in unique(Maternal$USUBJID)) {
  ID_subset = subset(Maternal, Maternal$USUBJID == current_ID)
  for (current_Sample in unique(ID_subset$Sample_Type)) {
    Sample_subset = subset(ID_subset, ID_subset$Sample_Type == current_Sample)
    if (nrow(Sample_subset) == 1 ) {next}
    Sample_subset_FR = Sample_subset[Sample_subset$Planned_Time != "Delivery/Birth", Antigen_list]/
      Sample_subset[Sample_subset$Planned_Time == "Delivery/Birth", Antigen_list]
    
    
    Sample_subset_FR_T = data.frame(t(Sample_subset_FR)); colnames(Sample_subset_FR_T) = "m6.Delivery"
    Sample_subset_FR_T$Ag = rownames(Sample_subset_FR_T)
    Sample_subset_FR_T$Sample = current_Sample 
    Sample_subset_FR_T$ID = current_ID
    Sample_subset_FR_T$Treatment = unique(ID_subset$Treatment_group)
    Fold_Difference = rbind(Fold_Difference,Sample_subset_FR_T)
  }
}


# for (current_ID in unique(Maternal$USUBJID)) {
#   ID_subset = subset(Maternal, Maternal$USUBJID == current_ID)
#   for (current_time in unique(ID_subset$Planned_Time)) {
#     Time_subset = subset(ID_subset, ID_subset$Planned_Time == current_time)
#     if (nrow(Time_subset) == 1 ) {next}
#     Time_subset_FR = Time_subset[Time_subset$Sample_Type != "Serum", Antigen_list]/
#       Time_subset[Time_subset$Sample_Type == "Serum", Antigen_list]
#     Time_subset_FR_T = data.frame(t(Time_subset_FR)); colnames(Time_subset_FR_T) = "HBM.Serum"
#     Time_subset_FR_T$Ag = rownames(Time_subset_FR_T)
#     Time_subset_FR_T$Time = current_time 
#     Time_subset_FR_T$ID = current_ID
#     Time_subset_FR_T$Treatment = unique(ID_subset$Treatment_group)
#     Fold_Difference = rbind(Fold_Difference,Time_subset_FR_T)
#   }
# }
Fold_Difference$Sample = ifelse(Fold_Difference$Sample == "Serum", "Serum", "HBM")

Fold_Difference$Ag = gsub("_Titer", "", Fold_Difference$Ag)
Fold_Difference$Ag = gsub("PT", "Pertussis Toxin", Fold_Difference$Ag)
Fold_Difference$Ag = gsub("PRN", "Pertactin", Fold_Difference$Ag)
Fold_Difference$Ag = gsub("FH", "Filamentous\nHemagglutinin", Fold_Difference$Ag)
Fold_Difference$Ag = factor(Fold_Difference$Ag, 
                                levels = c("Filamentous\nHemagglutinin", 
                                           "Pertactin", "Pertussis Toxin", 
                                           "Diphtheria", "Tetanus"))
boxplot_decay = ggplot(Fold_Difference, 
       aes(x = Ag, y = m6.Delivery, fill = Sample)) + 
  geom_hline(yintercept = 1) + 
  geom_boxplot(color = "black") + 
  facet_wrap(~Treatment, ncol = 2, strip.position = "top") + 
  theme_bw() + 
  ylab("6m/Birth") + 
  coord_cartesian(ylim = c(0, 3.5)) + 
  stat_compare_means(aes(label = paste0("p=", after_stat(p.format))),
                     label.y = 3.24) + 
  scale_fill_manual(breaks = c("HBM", "Serum"), 
                    values = c("#FB7676", "#6F0202")) + 
  theme(axis.title.x = element_blank(), 
    strip.background =element_rect(fill="white"),
    strip.text = element_text(size = 12,
                              margin = margin(1,1,1,1, "mm")),
    legend.position = "bottom",
    legend.margin=margin(t = -0.15, unit='cm'))


####################################################################################
#Line Graphs
####################################################################################
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
                                           "Filamentous\nHemagluttinin","Pertactin",
                                           "Tetanus", "Diptheria"))
All_Antibodies_OG = merge(All_Antibodies_OG, Proper_Labels)
All_long = merge(All_long, Proper_Labels)

All_Antibodies_OG$Measure_2 = factor(All_Antibodies_OG$Measure_2, 
                                     levels =c("Filamentous\nHemagluttinin","Pertactin","Pertussis\nToxin",
                                               "Diptheria", "Tetanus") )



Days_Labels = data.frame("Days" = c("Day 1","Day 31",
                                    "Delivery/Birth", "Birth+42d","Birth+70d",
                                    "Birth+130d", "Birth+6m"), 
                         "Planned_Time" = c("Day 1","Day 31",
                                            "Delivery/Birth", "Birth+42d","Birth+70d",
                                            "Birth+130d", "Birth+6m"), 
                         "Days.Since.Birth" = c(-90, -60, 0, 42, 70, 130, 180))

All_Antibodies_OG = merge(All_Antibodies_OG, Days_Labels)
All_long = merge(All_long, Days_Labels)
All_long$Grouping = paste(All_long$AGEU, All_long$Days.Since.Birth, All_long$Treatment, sep = "_")

All_Antibodies = subset(All_Antibodies_OG, grepl("Infant", All_Antibodies_OG$Group ) == F)
All_Antibodies = subset(All_Antibodies, All_Antibodies$Group != "Breastmilk\nIgA")

All_Antibodies$Group= gsub("Breastmilk\nIgG", "HBM IgG", All_Antibodies$Group)
All_Antibodies$Group= gsub("\n", " ", All_Antibodies$Group)


All_long_1 = subset(All_long, All_long$Antibody_Type == "IgG" & All_long$AGEU == "YEARS" & 
                      All_long$Days.Since.Birth >= 0) 

All_long_1$Group = ifelse(All_long_1$Sample_Type == "Serum", "Maternal Serum IgG", "HBM IgG")

All_long_1_180 = data.frame(table(All_long_1$USUBJID))
All_long_1_180 = subset(All_long_1_180, All_long_1_180$Freq == 30)
table(All_long_1_180$Freq)
All_long_1 = subset(All_long_1, All_long_1$USUBJID %in% All_long_1_180$Var1)

Titer_plot_P = ggplot(data = All_Antibodies[! All_Antibodies$Days.Since.Birth %in% c(-90, -60),], 
                      aes(x = Days.Since.Birth, y = (GMT_Mean), group = Grouping, 
                          color = Group, linetype = Treatment)) + 
  geom_line(linewidth = 1) + 
  geom_point() + 
  scale_color_manual(name = "Antibodies",
                     breaks = c("Maternal Serum IgG", "HBM IgG"),
                     values = c("#6F0202", "#FB7676")) +
  scale_linetype_manual(name = "Maternal\nVaccination", 
                        breaks = c("Td", "Tdap"), 
                        values = c("11", "solid")) + 
  ggh4x::facet_grid2(Measure_2~Group,  scales = "free_y", independent = "y") + 
  theme_bw() + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 13), 
        legend.position = "none", 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 12,
                                  margin = margin(1,1,1,1, "mm"))) + 
  xlab("Days after Delivery(0)")  + ylab("GMT [log10]"); Titer_plot_P

jpeg("./Figures/IgG_Decay_HBM_Serum.jpeg", res = 400,
     height = 3000, width = 7000)
ggarrange(Titer_plot_P, boxplot_decay, ncol = 2, common.legend = T, widths  = c(1, 2))
dev.off()

