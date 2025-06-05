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
current_titer = "FH_Titer"
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

  current_pair  = c("Breastmilk_IgG", "Breastmilk_IgA")
  for (current_pair in list( c("Breastmilk_IgG", "Breastmilk_IgA"),
                             c("Breastmilk_IgG", "Serum_IgG")  )) {
    Pair_data = All_wide#[All_wide$ID != "CPA.06945",]
    
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
  
    
    if (current_titer == "PT_Titer") {
      current_title = "Pertussis Toxin"
    } else if (current_titer == "PRN_Titer") {
      current_title = "Pertactin"
    } else if (current_titer == "FH_Titer") {
      current_title = "Filamentous Hemagglutinin"
    } else {
      current_title = "uh oh"
    }
    
    
    
    correlation_plot = ggplot(Pair_data, aes(x = log10(Pair1), 
                                             y = log10(Pair2),
                                             color = Planned_Time, 
                                             shape = Planned_Time)) + 
      
      ggtitle(current_title) + 
      
      geom_point(size = 1) + 
      geom_smooth(data = Delivery,  aes(x= log10(Pair1), y = log10(Pair2)),
                  color = "#4B0092", fill = "#4B0092", 
                  method = "lm", formula = "y~x") + 
      geom_smooth(data = Month6, aes(x= log10(Pair1), y = log10(Pair2)),
                  color = "#40B0A6", fill = "#40B0A6",
                  method = "lm", formula = "y~x") +
      
      annotate(geom = "label", 
                 label = paste("Del.: rs=", round(Delivery_corr_rho, 2), 
                               ", p",Delivery_corr_pvalue,
                               "\n6m: rs=", round(Month6_corr_rho, 2), 
                               ", p",Month6_corr_pvalue, sep = ""), 
               x = -Inf, y = Inf, hjust = -0.01, vjust = 1.1, 
                 color = "black",size = 3.5) + 
      
      theme_bw() + 
      
      scale_color_manual(name = "Time", 
                         breaks = c("Delivery/Birth", "Birth+6m"), 
                         values = c("#4B0092", "#2a6439")) + 
      scale_shape_manual(name = "Time", 
                         breaks = c("Delivery/Birth", "Birth+6m"), 
                         values = c(16, 17)) +
      
      xlab(paste(gsub("_", " ", current_pair[1]), "[log10]")) + 
      ylab(paste(gsub("_", " ", current_pair[2]), "[log10]")) +   
      
      theme(legend.text = element_text(size=11),
            legend.position = "bottom") + 
      guides(colour = guide_legend(override.aes = list(size=3)))#; plot(correlation_plot)
    
    current_pair_title = paste(current_pair, collapse = "_")
    
    assign(paste(current_titer, current_pair_title,"_correlation", sep = ""), correlation_plot)
    print(paste(current_titer,current_pair_title,  "_correlation", sep = ""))
    
  }
}

BM.IgA_BM.IgG = ggarrange(FH_TiterBreastmilk_IgG_Serum_IgG_correlation, 
                          PRN_TiterBreastmilk_IgG_Serum_IgG_correlation, 
                          PT_TiterBreastmilk_IgG_Serum_IgG_correlation, 
                          FH_TiterBreastmilk_IgG_Breastmilk_IgA_correlation, 
                          PRN_TiterBreastmilk_IgG_Breastmilk_IgA_correlation,
                          PT_TiterBreastmilk_IgG_Breastmilk_IgA_correlation,
                          
                          ncol = 3,nrow = 2, align = "hv",
                          legend = "bottom", common.legend = T)

jpeg(paste("./Figures/Correlation_All_Pertussis_OnlyTdap_V2.jpeg", sep = ""),
     res = 400, height = 2500, width = 4000)
plot(BM.IgA_BM.IgG)
dev.off()

