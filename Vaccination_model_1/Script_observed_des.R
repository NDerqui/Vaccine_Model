###SARS-CoV-2 Rt Model for Vaccination and Variants###

      ###Model update exploratory script###
    ##Proposal 2: Multiple doses, one variant##

#MRes Biomedical Research - EECID stream - Project 1#



# SET UP ------------------------------------------------------------------

# Packages

library(dplyr)
library(rstudioapi)
library(StanHeaders)
library(ggplot2)
library(rstan)
library(matrixStats)
library(ggpubr)
library(here)
library(tidyr)
library(sjmisc)
library(xlsx)



# LOAD DATA FOR DESCRIPTION -----------------------------------------------

data_model <- readRDS("data_model_for_plots.Rds")

dir.create(here("Figures/Description"),recursive = TRUE)
dir.create(here("Results/Description"),recursive = TRUE) 



# PLOTS --------------------------------------------------------------


Rt_data <- exp(data_model$Rt)

sum_data <- data.frame(Rt_data, 
                       LTLA = data_model$LTLAs,
                       Dose_1 = data_model$First_Prop,
                       Dose_2 = data_model$Second_Prop,
                       Dose_3 = data_model$Third_Prop,
                       date = data_model$date,
                       row.names = paste0("Rt", 1:12726))


#### Vax Description ####

# Vaccination over time

png("Figures/Description/Vax_over_time_good.png", width = 12, height = 8, units = 'in', res = 600)

Vax_Date <- ggplot(data = sum_data) +
  geom_boxplot (mapping = aes(x = date, y = Dose_1, group = date, color = "Dose_1"),
                size = rel(0.5)) +
  geom_boxplot (mapping = aes(x = date, y = Dose_2, group = date, color = "Dose_2"),
                size = rel(0.5)) +
  geom_boxplot (mapping = aes(x = date, y = Dose_3, group = date, color = "Dose_3"),
                size = rel(0.5)) +
  scale_color_manual(name="Number of doses",
                     breaks = c("Dose_1", "Dose_2", "Dose_3"),
                     values = c("Dose_1"="navy", 
                                "Dose_2"="cyan3",
                                "Dose_3"="lightgreen"),
                     labels=c("1 dose", "2 doses", "3 doses")) +
  theme_classic() +
  labs(title = "Distribution of vaccinated individuals proportion across all LTLAs over time",
       x = "Date",
       y = "Total proportion of vaccinated population") +
  theme(
    plot.title = element_text(size = rel(1.7), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(1.5), face="bold"),
    axis.title.y = element_text(size = rel(1.5), face="bold"),
    axis.text = element_text(size=rel(1.2)),
    legend.title = element_text(size = rel(1.5), face="bold"),
    legend.text = element_text(size=rel(1.2)))
Vax_Date

dev.off()

png("Figures/Description/Vax_over_time_line.png", width = 12, height = 8, units = 'in', res = 600)

Vax_Date_line <- ggplot(data = sum_data) +
  geom_line (mapping = aes(x = date, y = Dose_1, group = LTLA,
                           color = "Dose_1")) +
  geom_line (mapping = aes(x = date, y = Dose_2, group = LTLA,
                           color = "Dose_2")) +
  geom_line (mapping = aes(x = date, y = Dose_3, group = LTLA,
                           color = "Dose_3")) +
  scale_color_manual(name="Number of doses",
                     breaks = c("Dose_1", "Dose_2", "Dose_3"),
                     values = c("Dose_1"="navy", 
                                "Dose_2"="cyan3",
                                "Dose_3"="lightgreen"),
                     labels=c("1 dose", "2 doses", "3 doses")) +
  theme_classic() +
  labs(title = "Proportion of vaccinated individuals in all LTLAs over time",
       x = "Date",
       y = "Total proportion of vaccinated population") +
  theme(
    plot.title = element_text(size = rel(1.7), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(1.5), face="bold"),
    axis.title.y = element_text(size = rel(1.5), face="bold"),
    axis.text = element_text(size=rel(1.2)),
    legend.title = element_text(size = rel(1.5), face="bold"),
    legend.text = element_text(size=rel(1.2)))
Vax_Date_line

dev.off()


#### Rt over time ####

# General

png("Figures/Description/Observed_Rt.png", width = 12, height = 8, units = 'in', res = 600)

Rt_Obs_Date <- ggplot(data = sum_data) +
  geom_boxplot (mapping = aes(x = date, y = Rt_data, group = date), size = rel(0.5)) +
  theme_classic() +
  labs(title = "Observed Rt across all LTLAs over time",
       x = "Date",
       y = "Rt Observed") +
  theme(
    plot.title = element_text(size = rel(1.7), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(1.5), face="bold"),
    axis.title.y = element_text(size = rel(1.5), face="bold"),
    axis.text = element_text(size=rel(1.2)),
    legend.title = element_text(size = rel(1.5), face="bold"),
    legend.text = element_text(size=rel(1.2)))
Rt_Obs_Date

dev.off()


#### Vaccination in age groups ####

# Get the data - apply same cleaning criteria

# data_vax <- read.csv("rtm_incoming_vaccination_20211116-173218-f36a1245_vacc_coverage_ltla.csv")

#We want age groups
data_age <- data_vax
data_age <- filter(data_age, !is.na(age_band_min))

#Same date cut-offs
data_age <- rename(data_age, date = vaccination_date)
data_age <- filter(data_age, date >= Date_Start)
data_age <- filter(data_age, date <= Date_End)

# Select the vars of interest
data_age <- select(data_age, "ltla_name", "date", "age_band_min", "age_band_max",
                   "First_Prop", "Second_Prop", "Third_Prop")

# Weekly dates: var for no. of weeks & only one obs per week (per age group!)
data_age$week <- round(as.numeric(floor((data_age$date - Date_Start)/7)), digits = 0)
#Combi for unique combination of LTLA*week*age_group
data_age <- mutate(data_age, combi = paste0(data_age$week, data_age$ltla_name, data_age$age_band_min))
#Remove the duplicates
data_age <- filter(data_age, !duplicated(data_age$combi))

#Var for age group
data_age$age_group <- paste0(data_age$age_band_min, "-", data_age$age_band_max)

#No NA
data_age <- data_age[complete.cases(data_age),]

#Select cols
data_age <- select(data_age, "ltla_name", "date", "week", "age_group",
                   "First_Prop", "Second_Prop", "Third_Prop")

# Plot with age groups

png("Figures/Description/Vax_over_time_with_age.png", width = 15, height = 10, units = 'in', res = 600)

library(forcats)
forcats::fct_relevel

Vax_Date_age <-  ggplot(data = data_age) +
  geom_boxplot (mapping = aes(x = date, y = First_Prop, group = date,
                              color = "First_Prop"), size = rel(0.2), outlier.size = 0.01) +
  geom_boxplot (mapping = aes(x = date, y = Second_Prop, group = date,
                              color = "Second_Prop"), size = rel(0.2), outlier.size = 0.01) +
  geom_boxplot (mapping = aes(x = date, y = Third_Prop, group = date,
                              color = "Third_Prop"), size = rel(0.2), outlier.size = 0.01) +
  scale_color_manual(name="Number of doses",
                     breaks = c("First_Prop", "Second_Prop", "Third_Prop"),
                     values = c("First_Prop"="navy", 
                                "Second_Prop"="cyan3",
                                "Third_Prop"="lightgreen"),
                     labels=c("1 dose", "2 doses", "3 doses")) +
  theme_classic() +
  labs(title = "Distribution of vaccinated individuals proportion across all LTLAs over time",
       x = "Date",
       y = "Age-group proportion of vaccinated population") +
  theme(
    plot.title = element_text(size = rel(1.7), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(1.5), face="bold"),
    axis.title.y = element_text(size = rel(1.5), face="bold"),
    axis.text = element_text(size=rel(1.2)),
    legend.title = element_text(size = rel(1.5), face="bold"),
    legend.text = element_text(size=rel(1.2))) +
  facet_wrap( ~ fct_relevel(age_group,"0-4", "5-9", "10-14", "15-19","20-24",
                            "25-29", "30-34", "35-39", "40-44",
                            "45-49", "50-54","55-59", "60-64", "65-69",
                            "70-74", "75-79", "80-84", "85-89", "90-120")) +
  theme(strip.text.x = element_text(size = rel(1.2)))
Vax_Date_age

dev.off()



# TABLES -------------------------------------------------------------


names(sum_data)
library(tidyr)
library(sjmisc)
library(xlsx)
library(matrixStats)


##### Vax #####

# Dose 1

vax_long_1_mean <- sum_data %>%
  select("LTLA", "date", "Dose_1") %>%
  pivot_wider(names_from = date, values_from = Dose_1) %>%
  colMeans() %>%
  as.data.frame() %>%
  rename(mean = ".") %>%
  slice(-1)
vax_long_1_max <- sum_data %>%
  select("LTLA", "date", "Dose_1") %>%
  pivot_wider(names_from = date, values_from = Dose_1) %>%
  as.matrix() %>%
  colMaxs () %>%
  as.data.frame() %>%
  rename(max = ".") %>%
  slice(-1)
vax_long_1_min <- sum_data %>%
  select("LTLA", "date", "Dose_1") %>%
  pivot_wider(names_from = date, values_from = Dose_1) %>%
  as.matrix() %>%
  colMins () %>%
  as.data.frame() %>%
  rename(min = ".") %>%
  slice(-1)
vax_long_1_sd <- sum_data %>%
  select("LTLA", "date", "Dose_1") %>%
  pivot_wider(names_from = date, values_from = Dose_1) %>%
  as.matrix() %>%
  colSds () %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  slice(-1)
vax_long_1_iqr <- sum_data %>%
  select("LTLA", "date", "Dose_1") %>%
  pivot_wider(names_from = date, values_from = Dose_1) %>%
  as.matrix() %>%
  colQuantiles(probs = seq(from = 0, to = 1, by = 0.25)) %>%
  as.data.frame() %>%
  slice(-1)

vax_dose1 <- round(data.frame(vax_long_1_mean, vax_long_1_sd,
                              vax_long_1_min, vax_long_1_iqr, vax_long_1_max),
                   digits = 4)
write.xlsx(vax_dose1, row.names = TRUE,
           "Results/Description/Dose1_time.xlsx")

# Dose 2

vax_long_2_mean <- sum_data %>%
  select("LTLA", "date", "Dose_2") %>%
  pivot_wider(names_from = date, values_from = Dose_2) %>%
  colMeans() %>%
  as.data.frame() %>%
  rename(mean = ".") %>%
  slice(-1)
vax_long_2_max <- sum_data %>%
  select("LTLA", "date", "Dose_2") %>%
  pivot_wider(names_from = date, values_from = Dose_2) %>%
  as.matrix() %>%
  colMaxs () %>%
  as.data.frame() %>%
  rename(max = ".") %>%
  slice(-1)
vax_long_2_min <- sum_data %>%
  select("LTLA", "date", "Dose_2") %>%
  pivot_wider(names_from = date, values_from = Dose_2) %>%
  as.matrix() %>%
  colMins () %>%
  as.data.frame() %>%
  rename(min = ".") %>%
  slice(-1)
vax_long_2_sd <- sum_data %>%
  select("LTLA", "date", "Dose_2") %>%
  pivot_wider(names_from = date, values_from = Dose_2) %>%
  as.matrix() %>%
  colSds () %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  slice(-1)
vax_long_2_iqr <- sum_data %>%
  select("LTLA", "date", "Dose_2") %>%
  pivot_wider(names_from = date, values_from = Dose_2) %>%
  as.matrix() %>%
  colQuantiles(probs = seq(from = 0, to = 1, by = 0.25)) %>%
  as.data.frame() %>%
  slice(-1)

vax_dose2 <- round(data.frame(vax_long_2_mean, vax_long_2_sd,
                              vax_long_2_min, vax_long_2_iqr, vax_long_2_max),
                   digits = 4)
write.xlsx(vax_dose2, row.names = TRUE,
           "Results/Description/Dose2_time.xlsx")

# Dose 3

vax_long_3_mean <- sum_data %>%
  select("LTLA", "date", "Dose_3") %>%
  pivot_wider(names_from = date, values_from = Dose_3) %>%
  colMeans() %>%
  as.data.frame() %>%
  rename(mean = ".") %>%
  slice(-1)
vax_long_3_max <- sum_data %>%
  select("LTLA", "date", "Dose_3") %>%
  pivot_wider(names_from = date, values_from = Dose_3) %>%
  as.matrix() %>%
  colMaxs () %>%
  as.data.frame() %>%
  rename(max = ".") %>%
  slice(-1)
vax_long_3_min <- sum_data %>%
  select("LTLA", "date", "Dose_3") %>%
  pivot_wider(names_from = date, values_from = Dose_3) %>%
  as.matrix() %>%
  colMins () %>%
  as.data.frame() %>%
  rename(min = ".") %>%
  slice(-1)
vax_long_3_sd <- sum_data %>%
  select("LTLA", "date", "Dose_3") %>%
  pivot_wider(names_from = date, values_from = Dose_3) %>%
  as.matrix() %>%
  colSds () %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  slice(-1)
vax_long_3_iqr <- sum_data %>%
  select("LTLA", "date", "Dose_3") %>%
  pivot_wider(names_from = date, values_from = Dose_3) %>%
  as.matrix() %>%
  colQuantiles(probs = seq(from = 0, to = 1, by = 0.25)) %>%
  as.data.frame() %>%
  slice(-1)

vax_dose3 <- round(data.frame(vax_long_3_mean, vax_long_3_sd,
                              vax_long_3_min, vax_long_3_iqr, vax_long_3_max),
                   digits = 4)
write.xlsx(vax_dose3, row.names = TRUE,
           "Results/Description/Dose3_time.xlsx")


##### Rt #####

rt_mean <- sum_data %>%
  select("LTLA", "date", "Rt_data") %>%
  pivot_wider(names_from = date, values_from = Rt_data) %>%
  colMeans() %>%
  as.data.frame() %>%
  rename(mean = ".") %>%
  slice(-1)
rt_max <- sum_data %>%
  select("LTLA", "date", "Rt_data") %>%
  pivot_wider(names_from = date, values_from = Rt_data) %>%
  as.matrix() %>%
  colMaxs () %>%
  as.data.frame() %>%
  rename(max = ".") %>%
  slice(-1)
rt_min <- sum_data %>%
  select("LTLA", "date", "Rt_data") %>%
  pivot_wider(names_from = date, values_from = Rt_data) %>%
  as.matrix() %>%
  colMins() %>%
  as.data.frame() %>%
  rename(min = ".") %>%
  slice(-1)
rt_sd <- sum_data %>%
  select("LTLA", "date", "Rt_data") %>%
  pivot_wider(names_from = date, values_from = Rt_data) %>%
  as.matrix() %>%
  colSds() %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  slice(-1)
rt_iqr <- sum_data %>%
  select("LTLA", "date", "Rt_data") %>%
  pivot_wider(names_from = date, values_from = Rt_data) %>%
  as.matrix() %>%
  colQuantiles(probs = seq(from = 0, to = 1, by = 0.25)) %>%
  as.data.frame() %>%
  slice(-1)  

rt_table <- round(data.frame(rt_mean, rt_sd,
                             rt_min, rt_iqr, rt_max), digits = 4)
write.xlsx(rt_table, row.names = TRUE,
           "Results/Description/Rt_time.xlsx")


