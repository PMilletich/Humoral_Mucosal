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

##############################################################################
#
Metadata = subset(Metadata, Metadata$Sample_Type != "Serum")

Metadata = Metadata[,c("Treatment_group", "Tetanus_Titer", "Diphtheria_Titer", 
                       "Antibody_Type", "PT_Titer", "FH_Titer", "PRN_Titer",
                       "Subject_ID", "Planned_Time")]


Metadata_long = melt(setDT(Metadata), id.vars = c("Treatment_group", "Subject_ID", "Planned_Time", "Antibody_Type"), 
                     variable.name = "Antigen")

Metadata_long$TG_Ab = paste(Metadata_long$Treatment_group, Metadata_long$Antibody_Type, sep = ": ")
# Metadata_long$Treatment_group = gsub("Tdap", "Vaccinated", Metadata_long$Treatment_group)
# Metadata_long$Treatment_group = gsub("Td", "Unvaccinated", Metadata_long$Treatment_group)

# c("Vaccinated: PT_Titer", "Vaccinated: FH_Titer", "Vaccinated: PRN_Titer",
# "Unvaccinated: PT_Titer", "Unvaccinated: FH_Titer", "Unvaccinated: PRN_Titer")

Metadata_long$Antigen = gsub("PT_Titer", "Pertussis\nToxin", Metadata_long$Antigen)
Metadata_long$Antigen = gsub("FH_Titer", "Filamentous\nHemagglutinin", Metadata_long$Antigen)
Metadata_long$Antigen = gsub("PRN_Titer", "Pertactin", Metadata_long$Antigen)

Metadata_long$Antigen = gsub("Tetanus_Titer", "Tetanus", Metadata_long$Antigen)
Metadata_long$Antigen = gsub("Diphtheria_Titer", "Diphtheria", Metadata_long$Antigen)

table(Metadata_long$Antigen)
Metadata_long$Antigen = factor(Metadata_long$Antigen, 
                               levels =c("Filamentous\nHemagglutinin","Pertactin","Pertussis\nToxin",
                                         "Diphtheria", "Tetanus") )

Metadata_long$Planned_Time = factor(Metadata_long$Planned_Time, 
                                    levels = c("Delivery/Birth", "Birth+42d", "Birth+70d", "Birth+130d", "Birth+6m"))


Days_Labels = data.frame("Days" = c("Day 1","Day 31",
                                    "Delivery/Birth", "Birth+42d","Birth+70d",
                                    "Birth+130d", "Birth+6m"), 
                         "Planned_Time" = c("Day 1","Day 31",
                                            "Delivery/Birth", "Birth+42d","Birth+70d",
                                            "Birth+130d", "Birth+6m"), 
                         "Days.Since.Birth" = c(-90, -60, 0, 42, 70, 130, 180))

Metadata_long = merge(Metadata_long, Days_Labels)

Metadata_long$Days.Since.Birth = factor(Metadata_long$Days.Since.Birth)

BM_plot = ggplot(Metadata_long, aes(x = Days.Since.Birth, y = log2(value), fill = TG_Ab)) + 
  geom_boxplot(color = "black") + 
  facet_grid(Antibody_Type~Antigen) + 
  stat_compare_means(aes(label = after_stat(p.signif)), size = 5, hide.ns = T, 
                     label.y = 6.7,
                     method = "wilcox.test") + 
  theme_bw() + 
  ylab("Antibody Titer [log2]")+
  xlab("Days after Delivery(0)") + 
  scale_fill_manual(breaks = c("Td: IgA", "Td: IgG",
                               "Tdap: IgA", "Tdap: IgG" ), 
                    values = c("dodgerblue", "coral2",
                               "darkblue", "red3"), 
                    name = "Treatment") + 
  theme(legend.position = "bottom",
        legend.margin=margin(c(0,0,0,0)), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 12,
                                  margin = margin(1,1,1,1, "mm"))); BM_plot
#dev.off()



IgG = Metadata_long[Metadata_long$Antibody_Type== "IgG",]

jpeg("./Figures/Breastmilk_Td.Tdap.jpeg", res = 400, height = 2250, width = 5000)
BM_plot
dev.off()
