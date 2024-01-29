library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)


#Adjust to import appropriate file
Df1 <- read_csv("C:/Users/tfoll/OneDrive - Johns Hopkins/Defensins/Flow/Brain_Flow/Brain_Master_Flow.csv")
#Def_D1 <- read_csv("D:/OneDrive - Johns Hopkins/Defensins/Flow\7_9_23\7_9_23_RSTUDIO.csv")
# "C:\Users\tfoll\OneDrive - Johns Hopkins\Defensins\Flow\7_9_23\7_9_23_RSTUDIO.csv"
#C:\Users\tfoll\OneDrive - Johns Hopkins\Defensins\Meningies_Expression\3_15_23\3_15_23_qPCR_WT.csv"
#Make sure file imported properly
View(Df1)

#Calculate Percent of Live (P4), making new columns 

Df1 <-  mutate(Df1, 
          Neutrophils_Perc_Live =  (Df1$P7_Events / Df1$P4_Events)*100)
Df1 <-  mutate(Df1, 
               Monocytes_Perc_Live =  (Df1$P10_Events / Df1$P4_Events)*100)
Df1 <-  mutate(Df1, 
               Mast_Cells_Perc_Live =  (Df1$P13_Events / Df1$P4_Events)*100)
Df1 <-  mutate(Df1, 
               Microglia_Perc_Live =  (Df1$P11_Events / Df1$P4_Events)*100)
Df1 <-  mutate(Df1, 
               CD45_Perc_Live =  (Df1$P5_Events / Df1$P4_Events)*100)
Df1 <-  mutate(Df1, 
               Neutrophils_Perc_CD45 =  (Df1$P7_Events / Df1$P5_Events)*100)
Df1 <-  mutate(Df1, 
               Monocytes_Perc_CD45 =  (Df1$P10_Events / Df1$P5_Events)*100)
Df1 <-  mutate(Df1, 
               Mast_Cells_Perc_CD45 =  (Df1$P13_Events / Df1$P5_Events)*100)
Df1 <-  mutate(Df1, 
               Microglia_Perc_CD45 =  (Df1$P11_Events / Df1$P5_Events)*100)

write.csv(Df1, "C:\\Users\\tfoll\\Desktop\\FlowdataDura.csv", row.names=FALSE)

#Plot by Target the summary data

ggplot(data = Df1, mapping = aes(x = Injection, y = P7_Events)) +
  geom_boxplot()+
  facet_wrap(~ Day, nrow = 3, scales = "free")

ggplot(data = Df1, mapping = aes(x = Day, fill = Genotype)) +
  geom_bar(position = "dodge")

################################################################
# Perform the two-way ANOVA
model <- aov(P7_Events ~ Day * Genotype, data = Df1)
anova_table <- summary(model)
print(anova_table)

# This is good code to generate a plot of the cell-type (no iteration)
summary_data_Neu <- Df1 %>% 
  group_by(Day, Genotype) %>% 
  summarize(mean_output = mean(P7_Events),
            std_error = sd(P7_Events) / sqrt(n()))
view(summary_data_Neu)

# Create a plot with means and error bars using ggplot2
ggplot(summary_data_Neu, aes(x = Day, y = mean_output, color = Genotype)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_output - std_error, ymax = mean_output + std_error),
                width = 0.2) +
  labs(x = "Day", y = "Neutrophils", color = "Genotype") +
  theme_minimal()

######################################################################



#Two-Way ANOVA

mean_Neutrophils <- Df1 %>% 
  group_by(Day, Genotype) %>% 
  summarize(means_neu = mean(P7_Events))

ggplot(mean_Neutrophils, aes(x = Day, y = means_neu, color = Genotype)) +
  geom_point() +
  geom_line() +
  labs(x = "Day", y = "Mean Neutrophils", color = "Genotype") +
  theme_minimal()
###
mean_Monocytes <- Df1 %>% 
  group_by(Day, Genotype) %>% 
  summarize(means_mono = mean(P10_Events))

ggplot(mean_Monocytes, aes(x = Day, y = means_mono, color = Genotype)) +
  geom_point() +
  geom_line() +
  labs(x = "Day", y = "Mean Monocytes", color = "Genotype") +
  theme_minimal()
###
mean_Micro <- Df1 %>% 
  group_by(Day, Genotype) %>% 
  summarize(means_micro = mean(P11_Events))

ggplot(mean_Micro, aes(x = Day, y = means_micro, color = Genotype)) +
  geom_point() +
  geom_line() +
  labs(x = "Day", y = "Mean Microglia", color = "Genotype") +
  theme_minimal()
###
mean_Mast <- Df1 %>% 
  group_by(Day, Genotype) %>% 
  summarize(means_mast = mean(P13_Events))

ggplot(mean_Mast, aes(x = Day, y = means_mast, color = Genotype)) +
  geom_point() +
  geom_line() +
  labs(x = "Day", y = "Mean Mast Cells", color = "Genotype") +
  theme_minimal()

#########################################################################
#It works. Time to iterate

# Group the data by "Day" and "Genotype" and calculate means for the output variables

means_group <- Df1 %>%
  group_by(Day, Genotype) %>%
  summarize(mean_P7_Events = mean(P7_Events),
            mean_P10_Events = mean(P10_Events),
            mean_P11_Events = mean(P11_Events),
            mean_P13_Events = mean(P13_Events))

# View the resulting means for each combination of "Day" and "Genotype"
print(means_group)



#T-test analysis across the targets
"C:\Users\tfoll\Desktop"


# helper function
t_tester2 <- function(TARGET){
  tmp_data <- filter(Def_D1,Target_Name == TARGET)
  t_obj <- t.test(Relative_Expression~Injection,
                  data = tmp_data)
  finalFrame <- data.frame(Target = TARGET,
                           pValue = t_obj$p.val,
                           testStat = t_obj$p.val,
                           LowerCI = t_obj$conf.int[1],
                           UpperCI = t_obj$conf.int[2])
  return(finalFrame)
}

# test function
t_tester2("IL6") # all looks good

# now we iterate
myResults <- lapply(unique(Def_D1$Target_Name),t_tester2) %>% bind_rows(.)

#Post Hoc Test

myResults$pValue %>% p.adjust(.,method="BH")


