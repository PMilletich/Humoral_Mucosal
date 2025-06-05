#---- Script Metadata #----
# Title: Maternal_Meta.R
# Author: Trish Milletich, PhD
# Date: 2024-12-20
#--------------------------

library(ggplot2)
library(ggpubr)
library(data.table)

setwd("C:/Users/PMilletich/OneDrive - University of Maryland School of Medicine/Desktop/PATTI/HuMo/")
Metadata = read.csv("../Mother_Infant_All_3.csv") 

##############################################################################
#
Metadata = subset(Metadata, Metadata$Planned_Time %in% c("Delivery/Birth", "Birth+6m"))
Metadata = subset(Metadata, Metadata$Antibody_Type == "IgG")
Maternal = subset(Metadata, Metadata$AGEU == "YEARS")


#########################################
All_Antibodies = data.frame()
current_titer = "PT_Titer"
for(current_titer in c("PT_Titer","FH_Titer","PRN_Titer",#"Fimbraie_Titer",
                       "Tetanus_Titer","Diphtheria_Titer")) {
  All_wide = data.frame()
  current_patient = "CPA.05788" #unique(Maternal$Subject_ID)[1]
  Only_1 = c()
  for (current_patient in unique(Maternal$Subject_ID)) {
    ID_subset = subset(Maternal, Maternal$Subject_ID == current_patient)
    
    if (nrow(ID_subset[ID_subset$Sample_Type == "Breastmilk or colostrum",]) == 2) {
      HBM_FC=ID_subset[ID_subset$Planned_Time == "Birth+6m" & 
                        ID_subset$Sample_Type == "Breastmilk or colostrum",current_titer]/
        ID_subset[ID_subset$Planned_Time == "Delivery/Birth" & 
                    ID_subset$Sample_Type == "Breastmilk or colostrum",current_titer]
    } else {
      HBM_FC= NA
    }
    
    
    if (nrow(ID_subset[ID_subset$Sample_Type == "Serum",]) == 2) {
      Serum_FC=ID_subset[ID_subset$Planned_Time == "Birth+6m" & 
                         ID_subset$Sample_Type == "Serum",current_titer]/
        ID_subset[ID_subset$Planned_Time == "Delivery/Birth" & 
                    ID_subset$Sample_Type == "Serum",current_titer]
    } else {
      Serum_FC= NA
    }
    
    
    
    All_wide = rbind(All_wide, 
                     data.frame("ID" = current_patient, 
                                "Times" = c("Birth", "Month 6"), 
                                "FC_Samples" = c(HBM_FC, Serum_FC),
                                "Samples" = c("HBM", "Serum"),
                                "Treatment" = unique(ID_subset$Treatment_group), 
                                "Antigen" = current_titer))
  }
  All_Antibodies = rbind(All_Antibodies, All_wide)
}

All_Antibodies = All_Antibodies[complete.cases(All_Antibodies),]
All_Antibodies$Sample = ifelse(All_Antibodies$Samples == "HBM", "Breastmilk or colostrum", "Serum")
All_Antibodies$Ag = All_Antibodies$Antigen
#colnames(All_Antibodies) = c("ID", "Times", "FC_Samples", "Samples", "Treatment",  "Antigen", "Sample", "Ag")
data_merge = merge(All_Antibodies, Fold_Difference, all = T)
mixup = subset(data_merge, data_merge$m6.Delivery != data_merge$FC_Samples)


All_Antibodies$Antigen = gsub("_Titer", "", All_Antibodies$Antigen)
All_Antibodies$Antigen = gsub("PT", "Pertussis Toxin", All_Antibodies$Antigen)
All_Antibodies$Antigen = gsub("PRN", "Pertactin", All_Antibodies$Antigen)
All_Antibodies$Antigen = gsub("FH", "Filamentous\nHemagglutinin", All_Antibodies$Antigen)

All_Antibodies$Antigen = factor(All_Antibodies$Antigen, 
                                levels = c("Filamentous\nHemagglutinin", 
                                           "Pertactin", "Pertussis Toxin", 
                                           "Diphtheria", "Tetanus"))


#jpeg("./Figures/m6.Delivery_IgG_Serum.HBM_All.jpeg", res = 400, height = 3000, width = 3000)
ggplot(All_Antibodies, 
       aes(x = Antigen, y = FC_Samples, fill = Samples)) + 
  geom_hline(yintercept = 1) + 
  geom_boxplot(color = "black") + 
  facet_wrap(~Treatment, nrow = 2, strip.position = "right") + 
  theme_bw() + 
  ylab("6m/Birth") + 
  coord_cartesian(ylim = c(0, 3.5)) + 
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     label.y = 3.24) + 
  scale_fill_manual(breaks = c("HBM", "Serum"), 
                    values = c("#FF2828", "#6F0202")) + 
  theme(
    strip.background =element_rect(fill="white"),
    strip.text = element_text(size = 12,
                              margin = margin(1,1,1,1, "mm")),
    legend.position = "bottom",
    legend.margin=margin(t = -0.15, unit='cm'))
#dev.off()





