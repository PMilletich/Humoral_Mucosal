#---- Script Metadata #----
# Title: Correlation_LineScatter.R
# Author: Trish Milletich, PhD
# Date: 2025-05-05
#--------------------------

setwd("C:/Users/PMilletich/OneDrive - University of Maryland School of Medicine/Desktop/PATTI/HuMo/")
#setwd("~/Downloads/")

library(ggplot2)
library(ggpubr)
library(tidyr)
library(DescTools)

OG = read.csv("../Mother_Infant_All_4.csv")
OG = subset(OG, OG$Treatment_group == "Tdap")
Maternal = subset(OG, OG$AGEU == "YEARS")

Fim = subset(Maternal, is.na(Maternal$Fimbraie_Titer) == F)



Maternal = subset(Maternal, Maternal$Planned_Time %in% c("Delivery/Birth", "Birth+6m"))
#Maternal = subset(Maternal, Maternal$Treatment_group == "Tdap")
#####################################
#Correlation of IgG serum vs BM
#####################################
current_titer = "PT_Titer"
for(current_titer in c("PT_Titer","FH_Titer","PRN_Titer")) {
  
  All_wide = data.frame()
  missing_data = c()
  current_patient = unique(Maternal$Subject_ID)[1]
  for (current_patient in unique(Maternal$Subject_ID)) {
    ID_subset = subset(Maternal, Maternal$Subject_ID == current_patient)
    if (nrow(ID_subset) > 1){
      ID_subset = ID_subset[,c("Planned_Time", "Sample_Type", "Antibody_Type", 
                               "Treatment_group", current_titer)]
      ID_subset$ST = paste(ID_subset$Sample_Type, ID_subset$Antibody_Type, sep = "_")
      ID_subset$ST = gsub("Breastmilk or colostrum", "Breastmilk", ID_subset$ST)
      ID_subset$Sample_Type = NULL; ID_subset$Antibody_Type = NULL
      
      
      ID_subset_wide = spread(data = ID_subset, 
                              key = "ST", 
                              value = current_titer )
      ID_subset_wide$ID = current_patient
      All_wide = rbind(All_wide, ID_subset_wide)
    } else {
      missing_data = c(missing_data, current_patient)
    }
  }
  
  current_pair = list(c("Serum_IgG", "Breastmilk_IgG"), 
                      c("Breastmilk_IgG", "Breastmilk_IgA"))[[1]]
  for (current_pair in list( c("Breastmilk_IgG", "Breastmilk_IgA"),
                             c("Serum_IgG", "Breastmilk_IgG")
                           )) {
    
    print(paste(current_titer, current_pair[2]))
    
    
    if (current_pair[2] == "Breastmilk_IgA") {
      xlim_label = -1.5; ylim_label = 2
      ylim_max = 2.5; ylim_min = -1.25
    } else {
      xlim_label = 0.25
      ylim_label = 0.5
      ylim_min = -3.5;  ylim_max = 1.5
    }
    
    Pair_data = All_wide
    
    Pair_data$Pair1 = Pair_data[,current_pair[1]]
    Pair_data$Pair2 = Pair_data[,current_pair[2]]
    Pair_data = subset(Pair_data, is.na(Pair_data$Pair1) == F)
    Pair_data = subset(Pair_data, is.na(Pair_data$Pair2) == F)
    #Delivery 
    Delivery = subset(Pair_data, Pair_data$Planned_Time == "Delivery/Birth")
    Delivery_corr = cor.test(x = log10(Delivery[,current_pair[1]]), 
                             y = log10(Delivery[,current_pair[2]]), 
                             method = "pearson", exact = F)
    Delivery_corr_rho = Delivery_corr$estimate
    Delivery_corr_pvalue = ifelse(Delivery_corr$p.value< 0.001, "<0.001", 
                                  paste("=", round(Delivery_corr$p.value, 3), sep = ""))
    
    #Month 6 
    Month6 = subset(Pair_data, Pair_data$Planned_Time == "Birth+6m")
    Month6_corr = cor.test(x = log10(Month6[,current_pair[1]]), 
                           y = log10(Month6[,current_pair[2]]), 
                           method = "pearson", exact = F)
    Month6_corr_rho = Month6_corr$estimate
    Month6_corr_pvalue = ifelse(Month6_corr$p.value< 0.001, "<0.001", 
                                paste("=", round(Month6_corr$p.value, 3), sep = ""))
    
    correlation_plot = ggplot(Pair_data, aes(x = log10(Pair1), 
                                             y = log10(Pair2),
                                             color = Planned_Time, 
                                             shape = Planned_Time)) + 
      coord_cartesian(ylim = c(ylim_min, ylim_max)) + 
      geom_point(size = 1) + theme_bw() + 
      geom_label(label = paste(gsub("_", " ", current_titer), 
                               "\nDelivery: rs = ", round(Delivery_corr_rho, 2), ", p",Delivery_corr_pvalue,
                               "\nBirth+6m: rs = ", round(Month6_corr_rho, 2), ", p",Month6_corr_pvalue, sep = ""), 
                 aes(x = xlim_label, y = ylim_label), color = "black") + 
      scale_color_manual(name = "Time", 
                         breaks = c("Delivery/Birth", "Birth+6m"), 
                         values = c("#4B0092", "#2a6439")) + 
      scale_shape_manual(name = "Time", 
                         breaks = c("Delivery/Birth", "Birth+6m"), 
                         values = c(16, 17)) + 
      xlab(current_pair[1]) + ylab(current_pair[2]) + 
      
      theme(legend.text = element_text(size=11)) + 
      guides(colour = guide_legend(override.aes = list(size=3))) + 
      
      geom_smooth(data = Delivery,  aes(x= log10(Pair1), y = log10(Pair2)),
                  color = "#4B0092", fill = "#4B0092", 
                  method = "lm", formula = "y~x") + 
      geom_smooth(data = Month6, aes(x= log10(Pair1), y = log10(Pair2)),
                  color = "#40B0A6", fill = "#40B0A6",
                  method = "lm", formula = "y~x"); correlation_plot
    
    current_pair_title = paste(current_pair, collapse = "_")
    
    assign(paste(current_titer, current_pair_title,"_correlation", sep = ""), correlation_plot)
    print(paste(current_titer,current_pair_title,  "_correlation", sep = ""))
    
  }
}


BM.IgA_BM.IgG = ggarrange(
  PT_TiterSerum_IgG_Breastmilk_IgG_correlation, 
  FH_TiterSerum_IgG_Breastmilk_IgG_correlation, 
  PRN_TiterSerum_IgG_Breastmilk_IgG_correlation, 
  PT_TiterBreastmilk_IgG_Breastmilk_IgA_correlation,
  FH_TiterBreastmilk_IgG_Breastmilk_IgA_correlation, 
  PRN_TiterBreastmilk_IgG_Breastmilk_IgA_correlation,
  ncol = 3,nrow = 2, align = "hv", common.legend = T);BM.IgA_BM.IgG

jpeg(paste("./Figures/Correlation_All_Pertussis_OnlyTdap.jpeg", sep = ""),
     res = 400, height = 4000, width = 4000)
plot(BM.IgA_BM.IgG)
dev.off()
# 
# BM.IgG_Serum.IgG = ggarrange(PT_TiterSerum_IgG_Breastmilk_IgG_correlation, 
#                              FH_TiterSerum_IgG_Breastmilk_IgG_correlation, 
#                              PRN_TiterSerum_IgG_Breastmilk_IgG_correlation, 
#                              Tetanus_TiterSerum_IgG_Breastmilk_IgG_correlation, 
#                              Diphtheria_TiterSerum_IgG_Breastmilk_IgG_correlation, 
#                              ncol = 3, nrow = 2, common.legend = T, align = "hv")
# jpeg(paste("./Figures/Correlation_SerumIgG_BMIgG.jpeg", sep = ""),
#      res = 400, height = 3000, width = 4000)
# plot(BM.IgG_Serum.IgG)
# dev.off()
# 
# BM.IgA_Serum.IgG = ggarrange(PT_TiterSerum_IgG_Breastmilk_IgA_correlation,
#                      FH_TiterSerum_IgG_Breastmilk_IgA_correlation, 
#                      PRN_TiterSerum_IgG_Breastmilk_IgA_correlation, 
#                      Tetanus_TiterSerum_IgG_Breastmilk_IgA_correlation, 
#                      Diphtheria_TiterSerum_IgG_Breastmilk_IgA_correlation, 
#                      
#                      ncol = 3, nrow = 2,  common.legend = T, align = "hv")
#   
# jpeg(paste("./Figures/Correlation_SerumIgG_BMIgA.jpeg", sep = ""),
#      res = 400, height = 3000, width = 4000)
# plot(BM.IgA_Serum.IgG)
# dev.off()
# 
#   
# jpeg(paste("./Figures/Correlation_2Times_BM_Serum.jpeg", sep = ""),
#      res = 400, height = 3000, width = 4000)
# plot(All_corr)
# dev.off()
# 
# 
# #####################################
# #Ratio of IgG serum vs BM
# #####################################
# current_titer = "PT_Titer"
# for(current_titer in c("PT_Titer","FH_Titer","PRN_Titer")) {
#   All_wide = data.frame()
#   missing_data = c()
#   current_patient = unique(Maternal$Subject_ID)[1]
#   for (current_patient in unique(Maternal$Subject_ID)) {
#     ID_subset = subset(Maternal, Maternal$Subject_ID == current_patient)
#     if (nrow(ID_subset) > 1){
#       ID_subset = ID_subset[,c("Planned_Time", "Sample_Type", "Treatment_group", current_titer)]
#       ID_subset_wide = spread(data = ID_subset, 
#                               key = "Sample_Type", 
#                               value = current_titer )
#       ID_subset_wide$ID = current_patient
#       ID_subset_wide$S.BM = ID_subset_wide$Serum/ID_subset_wide$`Breastmilk or colostrum`
#       
#       ID_subset_wide$BM.S = ID_subset_wide$`Breastmilk or colostrum`/ID_subset_wide$Serum
#       All_wide = rbind(All_wide, ID_subset_wide)
#     } else {
#       missing_data = c(missing_data, current_patient)
#     }
#   }
#   All_wide = subset(All_wide, All_wide$Serum >= 10)
#   Delivery = subset(All_wide, All_wide$Planned_Time == "Delivery/Birth")
#   Birth6mo = subset(All_wide, All_wide$Planned_Time == "Birth+6m")
#   
#   summary(Delivery$BM.S)
#   summary(Delivery$S.BM)
#   
#   
#   
#   
#   Delivery_spearman = cor.test(x = log10(Delivery$`Breastmilk or colostrum`), 
#                                y = log10(Delivery$Serum), 
#                                method = "spearman", exact = F)
#   Delivery_rho = Delivery_spearman$estimate
#   Delivery_pvalue = Delivery_spearman$p.value
#   Delivery_pvalue = ifelse(Delivery_pvalue< 0.001, "<0.001", paste("=", round(Delivery_pvalue, 3), sep = ""))
#   
#   
#   Birth6mo_spearman = cor.test(x = log10(Birth6mo$`Breastmilk or colostrum`), 
#                                y = log10(Birth6mo$Serum), 
#                                method = "spearman", exact = F)
#   Birth6mo_rho = Birth6mo_spearman$estimate
#   Birth6mo_pvalue = Birth6mo_spearman$p.value
#   Birth6mo_pvalue = ifelse(Birth6mo_pvalue< 0.001, "<0.001", paste("=", round(Birth6mo_pvalue, 3), sep = ""))
#   
#   
#   correlation_plot = ggplot(All_wide, aes(x = log10(`Breastmilk or colostrum`), 
#                                           y = log10(Serum),
#                                           color = Planned_Time)) + 
#     geom_point(alpha = 0.8) + theme_bw() + 
#     ggtitle(paste(gsub("_", " ", current_titer), 
#                   "\nDelivery: rs = ", round(Delivery_rho, 2), ", p",Delivery_pvalue,
#                   "\nBirth+6m: rs = ", round(Birth6mo_rho, 2), ", p",Birth6mo_pvalue, sep = "")) + 
#     scale_color_manual(breaks = c("Delivery/Birth", "Birth+6m"), 
#                        values = c("mediumaquamarine", "grey4")) + 
#     xlab("log10(Breast milk IgG)") + ylab("log10(Serum IgG)") + 
#     geom_smooth(data = Delivery, aes(x= log10(`Breastmilk or colostrum`), y = log10(Serum)),
#                 color = "mediumaquamarine", method = "lm") + 
#     theme(legend.text = element_text(size=11)) + 
#     guides(colour = guide_legend(override.aes = list(size=3))) + 
#     geom_smooth(data = Birth6mo,  aes(x= log10(`Breastmilk or colostrum`), y = log10(Serum)),
#                 color = "grey4", method = "lm"); correlation_plot
#   
#   assign(paste(current_titer, "_correlation", sep = ""), correlation_plot)
#   print(paste(current_titer, "_correlation", sep = ""))
# }
# 
# 
