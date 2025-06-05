#---- Script Metadata #----
# Title: Maternal_HMB_Confounders_3.R
# Author: Trish Milletich, PhD
# Date: 2025-01-13
#--------------------------

library(ggplot2)
library(ggpubr)
library(pheatmap)
library(gridExtra)
library(EnvStats)
library(data.table)
library(patchwork)
library(randomcoloR)

setwd("/Users/pmilletich/OneDrive - University of Maryland School of Medicine/Desktop/PATTI/HuMo/")
#setwd("C:/Users/PMilletich/OneDrive - University of Maryland School of Medicine/Desktop/PATTI")
Metadata = read.csv("../Mother_Infant_All_3.csv")
Metadata$X = NULL 
Metadata = subset(Metadata, Metadata$Treatment_group == "Tdap")
Metadata$Sample_AB = paste(gsub("Breastmilk or colostrum", "BM", Metadata$Sample_Type), 
                           Metadata$Antibody_Type, sep = "_")

Maternal = subset(Metadata, Metadata$AGEU == "YEARS")
Maternal$Difference = Maternal$AAGE2-Maternal$AAGE1
Maternal$Infant_Sex = NA

Infant = subset(Metadata, Metadata$AGEU != "YEARS")
Infant$Infant_Sex = Infant$SEX
Infant$Difference = NA

# Infant Sex 
for (current_id in unique(Maternal$USUBJID)) {
  current_infant = subset(Infant, Infant$RSUBJID == current_id)
  current_infant = unique(current_infant$SEX)
  Maternal[Maternal$USUBJID == current_id, "Infant_Sex"] = current_infant
}



#Default variables and lists 
Antigen_list = c("PT_Titer","FH_Titer","PRN_Titer")
Compare_list = c("AGE","BMI", "IWGTKG", "AAGE1", "Infant_Sex", "Difference", "PARYST", "PREGNN" )
Antigen = "FH_Titer"; Confounder = "BMI"; Ig = "IgA"
current_time = "Birth+6m"; Sample_Type = "Breastmilk or colostrum"
Tdap = Maternal
##########################################################################
#Difference in maternal Serum after Vaccination 
##########################################################################
All_data = data.frame()
Vacc_data = Tdap
summary(Vacc_data[Vacc_data$Planned_Time == "Day 1", "PRN_Titer"])
summary(Vacc_data[Vacc_data$Planned_Time == "Day 1", "PT_Titer"])

current_id = unique(Vacc_data$Subject_ID)[1]
for (current_id in unique(Maternal$Subject_ID)) {
  ID_subset = subset(Vacc_data, Vacc_data$Subject_ID == current_id)
  ID_subset_meta = ID_subset[,c(Compare_list, "Subject_ID")]
  ID_subset_meta = unique(ID_subset_meta)
  
  ID_subset_Titers = ID_subset[,c("Planned_Time", "Sample_AB", "PT_Titer","FH_Titer","PRN_Titer")]
  
  ####################
  #Vaccination uptake 
  Uptake = subset(ID_subset_Titers, ID_subset_Titers$Planned_Time %in% c("Day 1", "Day 31"))
  if (nrow(Uptake) == 2) {
    Uptake_Ratio = Uptake[Uptake$Planned_Time == "Day 31", c("PT_Titer","FH_Titer","PRN_Titer")]/
      Uptake[Uptake$Planned_Time == "Day 1", c("PT_Titer","FH_Titer","PRN_Titer")]
  } else {
    Uptake_Ratio = data.frame("PT_Titer" = NA,  "FH_Titer" = NA,  "PRN_Titer"=NA)
  }
  Uptake_Ratio$Time = "Vaccine Uptake"
  Uptake_Ratio$Sample_AB = "Serum_IgG"
  
  ####################
  #Vaccination Decline
  Decline = subset(ID_subset_Titers, ID_subset_Titers$Planned_Time %in% c("Delivery/Birth", "Birth+6m"))
  Decline_Ratio_Final = data.frame()
  for (ST in unique(Decline$Sample_AB)) {
    Decline_Ratio = subset(Decline, Decline$Sample_AB == ST)
    
    if (nrow(Decline_Ratio) == 2) {
      Decline_Ratio = Decline_Ratio[Decline_Ratio$Planned_Time == "Birth+6m", c("PT_Titer","FH_Titer","PRN_Titer")]/
        Decline_Ratio[Decline_Ratio$Planned_Time == "Delivery/Birth", c("PT_Titer","FH_Titer","PRN_Titer")]
    } else {
      Decline_Ratio = data.frame("PT_Titer" = NA,  "FH_Titer" = NA,  "PRN_Titer"=NA)
    }
    Decline_Ratio$Time = "Vaccine Decline"
    Decline_Ratio$Sample_AB = ST
    Decline_Ratio_Final = rbind(Decline_Ratio_Final, Decline_Ratio)
  }
  
  ####################
  #Cordblood/Maternal 
  CB_Maternal = rbind(ID_subset[ID_subset$Planned_Time == "Delivery/Birth" & ID_subset$Sample_AB == "Serum_IgG",], 
                      Infant[Infant$Planned_Time == "Delivery/Birth" & Infant$RSUBJID == current_id, ])
  if (nrow(CB_Maternal) == 2) {
    Delivery_Ratio = CB_Maternal[CB_Maternal$AGEU != "YEARS", c("PT_Titer","FH_Titer","PRN_Titer")]/
      CB_Maternal[CB_Maternal$AGEU == "YEARS", c("PT_Titer","FH_Titer","PRN_Titer")]
    if (Delivery_Ratio$FH_Titer < 1) {
      print(current_id)
    }
      
  } else {
    Delivery_Ratio = data.frame("PT_Titer" = NA,  "FH_Titer" = NA,  "PRN_Titer"=NA)
  }

  Delivery_Ratio$Time = "Delivery"
  Delivery_Ratio$Sample_AB = "Serum_IgG"
  
  
  ####################
  Final_Ratios = rbind(Uptake_Ratio, Decline_Ratio_Final)
  Final_Ratios = rbind(Final_Ratios, Delivery_Ratio)
  Final_Ratios = merge(Final_Ratios, ID_subset_meta)
  
  All_data = rbind(All_data, Final_Ratios)
  
}

summary(All_data[All_data$Time == "Vaccine Uptake", "FH_Titer"])
summary(All_data[All_data$Time == "Vaccine Uptake", "PT_Titer"])
summary(All_data[All_data$Time == "Vaccine Uptake", "PRN_Titer"])



  All_data$Age_Grouped = ifelse(All_data$AGE <20, "18-19", 
                              ifelse(All_data$AGE <25, "20-24", 
                                     ifelse(All_data$AGE <30, "25-29", 
                                            ifelse(All_data$AGE <35, "30-34", "35-39"))))

All_data$BMI_Grouped = ifelse(All_data$BMI < 21, "17.0-20.9", 
                              ifelse(All_data$BMI <25, "21.0-24.9", 
                                     ifelse(All_data$BMI <30, "25.0-29.9", "30-44")))


summary(All_data$Difference)
All_data$Difference_Grouped = ifelse(All_data$Difference < 16, "11-15", 
                              ifelse(All_data$Difference <21, "16-20", 
                                     ifelse(All_data$Difference <26, "21-25", "26-28")))
table(All_data$Difference_Grouped)

summary(All_data$IWGTKG)
All_data$IWGTKG_Grouped =  ifelse(All_data$IWGTKG < 3, "1.5-3", 
                                  ifelse(All_data$IWGTKG <3.5, "3-3.5", "3.5-4.5"))
                                  

table(All_data$PREGNN)
All_data$PREGNN_Grouped =  ifelse(All_data$PREGNN ==1, "1", 
                                  ifelse(All_data$PREGNN <= 3, "2-3", 
                                         ifelse(All_data$PREGNN <7, "4-6", "7-10")))
table(All_data$PREGNN_Grouped)


Count_table = All_data[,c("Subject_ID", "Age_Grouped", "BMI_Grouped" , "IWGTKG_Grouped", "Difference_Grouped", 
                          "Infant_Sex", "PARYST", "PREGNN_Grouped")]
Count_table = unique(Count_table)
table(Count_table$Age_Grouped)
table(Count_table$BMI_Grouped)
table(Count_table$IWGTKG_Grouped)
table(Count_table$Difference_Grouped)
table(Count_table$PARYST)


All_data_2 = All_data[,c("PT_Titer","FH_Titer","PRN_Titer","Time","Sample_AB",
                         "Age_Grouped", "BMI_Grouped" , "IWGTKG_Grouped", "Difference_Grouped", 
                         "Infant_Sex", "PARYST", "PREGNN_Grouped")]


long <- melt(setDT(All_data_2), id.vars = c("Time","Sample_AB",
                                          "Age_Grouped", "BMI_Grouped" , "IWGTKG_Grouped", "Difference_Grouped", 
                                          "Infant_Sex", "PARYST", "PREGNN_Grouped"), 
             variable.name = "Antigen")

long$Time = factor(long$Time, levels = c("Vaccine Uptake", "Delivery", "Vaccine Decline"))

long$Antigen= ifelse(long$Antigen == "PT_Titer", "Pertussis Toxin", 
                     ifelse(long$Antigen == "FH_Titer", "Filamentous Hemagglutinin", "Pertactin"))
long = subset(long, is.na(long$value) == F)
current_test = 'PARYST'
set.seed(2)
for (current_test in c("Age_Grouped", "BMI_Grouped" , "IWGTKG_Grouped", "Difference_Grouped", 
                        "Infant_Sex", "PARYST", "PREGNN_Grouped"
                       )){
  print(current_test)
  long$test_value = long[,..current_test]
  long$Sample_AB = factor(long$Sample_AB, levels = c("Serum_IgG", "BM_IgG", "BM_IgA"))
  long$Time_SAB = paste(long$Time, long$Sample_AB, sep = "\n")
  long$Time_SAB = factor(long$Time_SAB, levels = c("Vaccine Uptake\nSerum_IgG", "Delivery\nSerum_IgG",
                                                   "Vaccine Decline\nSerum_IgG",
                                                    "Vaccine Decline\nBM_IgG", "Vaccine Decline\nBM_IgA"))

  long$X1 = paste(long$test_value, long$Time_SAB)

  plot_1 = ggplot(long, aes(x = Time_SAB , y = log2(value),fill = test_value)) + 
    geom_hline(yintercept = 0, color = "grey56") + 
    geom_boxplot(color= "black") + 
    ggh4x::facet_grid2(~Antigen, scales = "free", independent = "all") + 
    stat_compare_means(aes(label = paste0("p=", after_stat(p.format)))) + 
    theme_bw() + ylab("Antibody Response [log2]") +
    scale_fill_manual(values = randomColor(count = length(unique(long$test_value)),
                                           luminosity = " "))+
    xlab(gsub("_", " ", current_test)) + labs(fill = gsub("_", "\n", current_test)) + 
    theme(legend.position = "bottom", 
          axis.title.x = element_blank(), 
          strip.background = element_rect(fill="white"),
          legend.box.spacing = margin(0)); suppressWarnings(plot(plot_1))
  
  jpeg(paste("./Figures/", current_test, "_AB.Response_2.jpeg", sep = ""),
       res = 400, height = 2000, width = 6500)
  plot(plot_1)
  dev.off()
  
}
