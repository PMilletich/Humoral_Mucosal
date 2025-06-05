#---- Script Metadata #----
# Title: Maternal_Meta.R
# Author: Trish Milletich, PhD
# Date: 2024-12-20
#--------------------------

library(ggplot2)
library(ggpubr)
library(data.table)

setwd("C:/Users/PMilletich/OneDrive - University of Maryland School of Medicine/Desktop/PATTI/HuMo/")
Metadata = read.csv("../Mother_Infant_All_4.csv") 

Baseline = subset(Metadata, Metadata$Planned_Time == "Day 31")
table(Baseline[Baseline$PT_Titer < 10,"Treatment_group"])
table(Baseline$Treatment_group)


##############################################################################
#
unique(Metadata$Planned_Time)
Metadata = subset(Metadata, Metadata$Planned_Time %in% c("Delivery/Birth", "Birth+42d", "Birth+6m"))
#Metadata = subset(Metadata, Metadata$Antibody_Type == "IgG")
Maternal = subset(Metadata, Metadata$AGEU == "YEARS")


#########################################
All_Antibodies_OG = data.frame()
current_titer = "PT_Titer"
for(current_titer in c("PT_Titer","FH_Titer","PRN_Titer",#"Fimbraie_Titer",
                       "Tetanus_Titer","Diphtheria_Titer")) {
  All_wide = data.frame()
  current_patient = "CPA.05788" #unique(Maternal$Subject_ID)[1]
  Only_1 = c()
  for (current_patient in unique(Maternal$Subject_ID)) {
    ID_subset = subset(Maternal, Maternal$Subject_ID == current_patient)

    ID_subset_HBM_IgG = subset(ID_subset, ID_subset$Sample_Type != "Serum" & 
                             ID_subset$Antibody_Type == "IgG")
    ID_subset_HBM_IgA = subset(ID_subset, ID_subset$Sample_Type != "Serum" & 
                                 ID_subset$Antibody_Type == "IgA")
    
    
    #####################
    #BM IgG 6M/Day42 
    if (nrow(ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time != "Delivery/Birth",]) == 2) {
      HBM_FC_IgG =ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time == "Birth+6m",current_titer]/
        ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time == "Birth+42d",current_titer]
    } else { HBM_FC_IgG= NA }
    
    
    #####################
    #BM IgG Day42/Delivery 
    
    if (nrow(ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time != "Birth+6m",]) == 2) {
      HBM_FC_IgG_2=ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time == "Birth+42d",current_titer]/
        ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time == "Delivery/Birth",current_titer]
    } else { HBM_FC_IgG_2= NA }
    
    #####################
    #BM IgG 6M/Delivery 
    
    if (nrow(ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time != "Birth+42d",]) == 2) {
      HBM_FC_IgG_3 =ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time == "Birth+6m",current_titer]/
        ID_subset_HBM_IgG[ID_subset_HBM_IgG$Planned_Time == "Delivery/Birth",current_titer]
    } else { HBM_FC_IgG_3= NA }
    
    #####################
    #BM IgA 6M/Day42 
    
    if (nrow(ID_subset_HBM_IgA[ID_subset_HBM_IgA$Planned_Time != "Delivery/Birth",]) == 2) {
      HBM_FC_IgA =ID_subset_HBM_IgA[ID_subset_HBM_IgA$Planned_Time == "Birth+6m",current_titer]/
        ID_subset_HBM_IgA[ID_subset_HBM_IgA$Planned_Time == "Birth+42d",current_titer]
    } else {
      HBM_FC_IgA= NA
    }
    
    #####################
    #BM IgA Day42/Delivery 
    if (nrow(ID_subset_HBM_IgA[ID_subset_HBM_IgA$Planned_Time != "Birth+6m",]) == 2) {
      HBM_FC_IgA_2=ID_subset_HBM_IgA[ID_subset_HBM_IgA$Planned_Time == "Birth+42d",current_titer]/
        ID_subset_HBM_IgA[ID_subset_HBM_IgA$Planned_Time == "Delivery/Birth",current_titer]
    } else {
      HBM_FC_IgA_2= NA
    }
    
    #####################
    #Serum
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
                                "FC_Samples" = c(HBM_FC_IgG, HBM_FC_IgG_2, HBM_FC_IgG_3,
                                                 HBM_FC_IgA, HBM_FC_IgA_2, 
                                                 Serum_FC),
                                "Samples" = c("HBM IgG\n6m/Day42", "HBM IgG\nMature/Colostrum", "HBM IgG\n6m/Colostrum",
                                              "HBM IgA\n6m/Day42", "HBM IgA\nMature/Colostrum",
                                              "Serum IgG"),
                                "Color_Group" = c("HBM IgG", "HBM IgG", "HBM IgG",
                                                  "HBM IgA", "HBM IgA", 
                                                  "Serum IgG"),

                                "Treatment" = unique(ID_subset$Treatment_group), 
                                "Antigen" = current_titer))
  }
  All_Antibodies_OG = rbind(All_Antibodies_OG, All_wide)
}

Num_data = data.frame("X_Group" = c(1,2,4,5,7,8),
                      "Samples" = c( "HBM IgG\nMature/Colostrum","HBM IgA\nMature/Colostrum",
                                     "HBM IgG\n6m/Day42","HBM IgA\n6m/Day42",
                                     "HBM IgG\n6m/Colostrum", "Serum IgG"))



All_Antibodies = merge(All_Antibodies_OG,Num_data)#[complete.cases(All_Antibodies_OG),]
#All_Antibodies$Sample = ifelse(All_Antibodies$Samples == "HBM", "Breastmilk or colostrum", "Serum")
All_Antibodies$Ag = All_Antibodies$Antigen
#colnames(All_Antibodies) = c("ID", "Times", "FC_Samples", "Samples", "Treatment",  "Antigen", "Sample", "Ag")
# data_merge = merge(All_Antibodies, Fold_Difference, all = T)
# mixup = subset(data_merge, data_merge$m6.Delivery != data_merge$FC_Samples)


All_Antibodies$Antigen = gsub("_Titer", "", All_Antibodies$Antigen)
All_Antibodies$Antigen = gsub("PT", "Pertussis Toxin", All_Antibodies$Antigen)
All_Antibodies$Antigen = gsub("PRN", "Pertactin", All_Antibodies$Antigen)
All_Antibodies$Antigen = gsub("FH", "Filamentous Hemagglutinin", All_Antibodies$Antigen)

All_Antibodies$Antigen = factor(All_Antibodies$Antigen, 
                                levels = c("Filamentous Hemagglutinin", 
                                           "Pertactin", "Pertussis Toxin", 
                                           "Diphtheria", "Tetanus"))
All_Antibodies$Samples = factor(All_Antibodies$Samples , 
                                levels = c("HBM IgG\nMature/Colostrum", "HBM IgA\nMature/Colostrum",
                                           "HBM IgG\n6m/Day42","HBM IgA\n6m/Day42",
                                           "HBM IgG\n6m/Colostrum", "Serum IgG"))


#All_Antibodies$X_Group = factor(All_Antibodies$X_Group , 
#                                levels = c(2,5,1,4,3,6))



All_Antibodies = subset(All_Antibodies, All_Antibodies$ID != "CPA.06883")


All_Antibodies = subset(All_Antibodies, All_Antibodies$X_Group != 7)

Pertussis_1 = All_Antibodies[All_Antibodies$Treatment == "Tdap" & 
                               All_Antibodies$Antigen %in% c("Filamentous Hemagglutinin", 
                                                           "Pertactin", "Pertussis Toxin"),]

#Pertussis_1 = Pertussis_1[order(Pertussis_1$X_Group),]


jpeg("./Figures/m6.Delivery_IgG_Serum.HBM_All.jpeg", res = 400, height = 1200, width = 3000)
Tdap_Pert = ggplot(Pertussis_1, 
       aes(x = X_Group, y = FC_Samples, fill = Color_Group , group = Samples)) + 
  geom_hline(yintercept = 1) + 
  #ggtitle(paste("Tdap (n=", length(unique(Pertussis_1$ID)), ")", sep = "")) + 
  geom_boxplot(color = "black") + 
  facet_grid(~Antigen) + 
  theme_bw() + ylab("FC") + 
  coord_cartesian(ylim = c(0, 3.0)) + 
  stat_compare_means(comparisons = list(
    c(1,2), c(1,4),
    c(4,5), c(2,5),
    c(1,8),c(4,8)),
    #paired = T,
    size = 3,
    label.y = c(0, 0.3, 0, 0.5, 0.7, 0.9),
    tip.length = 0.0002, 
    quiet = T) + 
  scale_x_continuous(breaks = c(1.5, 4.5, 8), 
                     labels = c("Mature Milk/\nColostrum", "6m/\nMature Milk", "6m/\nDelivery")) + 
  scale_fill_manual(breaks = c("HBM IgG", "HBM IgA", "Serum IgG"),
                    values = c( "#FB7676","#81C2E6", "darkgreen")) +
  theme(strip.background =element_rect(fill="white"),
        axis.title.x = element_blank(), 
        strip.text = element_text(size = 10,
                                  margin = margin(1,1,1,1, "mm")),
        legend.position = "bottom",
        legend.margin=margin(t = -0.15, unit='cm')) ; Tdap_Pert
dev.off()
