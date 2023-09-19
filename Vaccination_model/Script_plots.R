###SARS-CoV-2 Rt Model for Vaccination and Variants###

        ###Model update exploratory script###
     ##Proposal 2: Multiple doses, one variant##

#MRes Biomedical Research - EECID stream - Project 1#


### Beter to not save the RData when closing the session ###



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
library(openxlsx)



# INPUT DATA  ------------------------------------------------------------


#### Basic data ####

# data_vax <- read.csv("Data/rtm_incoming_vaccination_20211116-173218-f36a1245_vacc_coverage_ltla.csv")
# data_rt <- read.csv("Data/UK_hotspot_Rt_estimates.csv")
# data_var <- read.csv("Data/vam_by_ltla.csv")
# 
# data_vax$vaccination_date <- as.Date(data_vax$vaccination_date)
# data_rt$date <- as.Date(data_rt$date)

# If needed, run the get_data function with the same specs as model

# data_stan <- get_data(data_vax = data_vax, data_rt = data_rt, data_var = data_var,
#                       covar_var = covar_var, covar_vax = covar_vax,
#                       Date_Start = Date_Start, Date_End = Date_End,
#                       lockdown_steps = lockdown_steps,
#                       DoVariants = 0, DoAge = 0,
#                       DoKnots = 0, Quadratic = 0,
#                       IncludeIntercept = 1, IncludeScaling = 1)


#### Options ####

DoKnots <- 0
Quadratic <- 0

DoAge <- 0
DoVariants <- 0


#### Extract parameters ####

if (DoAge == 0) {
  
  Rt_data <- data_stan[[22]]
  
  LTLA <- data_stan[[20]]
  
  Dose_1 <- data_stan[[23]][,1]
  Dose_2 <- data_stan[[23]][,2]
  Dose_3 <- data_stan[[23]][,3]
  
  date <- max(c(min(data_rt$date), min(data_var$date))) + data_stan[[19]]*7
  
} else {
  
  Rt_data <- data_stan_age[[22]]
  
  LTLA <- data_stan_age[[20]]
  
  Dose_1 <- data_stan_age[[23]][,1]
  Dose_2 <- data_stan_age[[23]][,2]
  Dose_3 <- data_stan_age[[23]][,3]
  
  date <- max(c(min(data_rt$date), min(data_var$date))) + data_stan_age[[19]]*7
  
}

Steps <- c(max(c(min(data_rt$date), min(data_var$date))),
           lockdown_steps[2:length(lockdown_steps)],
           min(c(max(data_rt$date), max(data_var$date))))



# MODEL DATA ------------------------------------------------------------


#### Load data ####

model_name <- "CorrectVar_2000_8_3D"

# fit <- readRDS(paste0(model_name, ".Rds"))
  
# loo_run <- readRDS(paste0("loo_", model_name, ".Rds"))

fit <- readRDS(paste0("C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Analyses/Models_BackUp/", model_name, ".Rds"))


#### Substract parameters ####

model_matrix <- as.matrix(fit)
# rm(fit)
# gc()

Rt_LogP <- colMeans(model_matrix[, grep("LogPredictions", colnames(model_matrix))])
# Rt_LogP <- sum[grep("LogPredictions", rownames(sum)), 1]

RegionalTrends <- colMeans(model_matrix[, grep("RegionalTrends", colnames(model_matrix))])
# RegionalTrends <- sum[grep("RegionalTrends", rownames(sum)), 1]

NationalTrend <- colMeans(model_matrix[, grep('^NationalTrend\\[', colnames(model_matrix))])
# NationalTrend <- sum[grep('^NationalTrend\\[', rownames(sum)), 1]

Lambda <- colMeans(model_matrix[, grep('^lambda\\[', colnames(model_matrix))])
# Lambda <- sum[grep('^lambda\\[', rownames(sum)), 1]

if (DoVariants == 1) {
  VarAdvantage <- colMeans(model_matrix[, grep('^VarAdvantage\\[', colnames(model_matrix))])
  # VarAdvantage <- sum[grep('^VarAdvantage\\[', rownames(sum)), 1]
}

sum_rt <- data.frame(LTLA, date, 
                     Rt_data, Rt_LogP, RegionalTrends, NationalTrend,
                     Dose_1, Dose_2, Dose_3)

if (DoKnots == 1) {
  sum_line <- data.frame(date = rep(Steps, times = length(unique(LTLAs))),
                         Lambda = Lambda)
}



# DIRECTORY ---------------------------------------------------------------

 dir.create(here(paste0("Figures/", model_name)),recursive = TRUE)
 dir.create(here(paste0("Results/", model_name)),recursive = TRUE) 



# TABLES ---------------------------------------------------------------
 
 
#### VaxEffect ####
 
VERedMean <- colMeans(model_matrix[, grep("VaxEffect", colnames(model_matrix))])
# VERedMean <- sum[grep('^VaxEffect\\[', rownames(sum)), 1]

VERedQuan <- colQuantiles(model_matrix[, grep("VaxEffect", colnames(model_matrix))], probs=c(0.025,0.975))
# VERedQuan <- sum[grep('^VaxEffect\\[', rownames(sum)), c(4,8)]
 
sum_ve <- round(data.frame(VERedMean, VERedQuan), digits = 4) %>%
  select("VERedMean", "X2.5.", "X97.5.")
colnames(sum_ve) <- c("VE (1-Red)", " VE 2.5% Q", " VE97.5% Q")

if (nrow(sum_ve) == 3) {
  row.names(sum_ve) <- c("Dose 1", "Dose 2", "Dose 3")

} else {
  
  if(nrow(sum_ve) == 6) {
    row.names(sum_ve) <- c("Alpha1", "Alpha2", "Alpha3",
                           "Delta1", "Delta2", "Delta3")
    } else {
      
      if (DoAge == 1) {
        row.names(sum_ve) <- c("15-49_D1", "15-49_D2", "15-49_D3",
                               "50-69_D1", "50-69_D2", "50-69_D3",
                               "70+_D1", "70+_D2", "70+_D3")
      } else {
        row.names(sum_ve) <- c("PreAl_1", "PreAl_2", "PreAl_3",
                               "Alpha1", "Alpha2", "Alpha3",
                               "Delta1", "Delta2", "Delta3")
      }
  }
}
                
write.xlsx(sum_ve, rowNames = TRUE,
            paste0("Results/", model_name, "/Vax_VE_table.xlsx"))


#### Var Advantage ####

if (DoVariants == 1) {
  
  write.xlsx(as.data.frame(VarAdvantage), rowNames = TRUE,
             paste0("Results/", model_name, "/Var_Ad_table.xlsx"))
}

 
 
# PLOTS -------------------------------------------------------------------
 
 
#### LogPred vs Obs ####
 
# General
 
png(paste0("Figures/", model_name, "/Obs_Pre_Rt.png"),
    width = 9, height = 8, units = 'in', res = 600)
 
ggplot(data = sum_rt) +
   geom_point(mapping = aes(x = Rt_data, y = Rt_LogP), size = rel(0.8)) +
   geom_abline(x = 0, y = 1, col = "red") +
   theme_classic() +
   labs(title = "Observed versus Predicted Rt in all LTLAs and timepoints",
        x = "Rt Observed",
        y = "Rt Predicted") +
   theme(
     plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(0.9), face="bold"),
     axis.title.y = element_text(size = rel(0.9), face="bold"),
     axis.text = element_text(size=rel(0.7)))

dev.off()
 
# As a scatter plot
 
png(paste0("Figures/", model_name, "/Pre_Obs_Rt_time.png"),
    width = 9, height = 5, units = 'in', res = 600)
 
ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = Rt_data, group = date,
                               color = "Rt_data"), size = rel(0.5)) +
  geom_boxplot (mapping = aes(x = date, y = Rt_LogP, group = date,
                               color = "Rt_LogP"), size = rel(0.5)) + 
  scale_color_manual(name="Reproduction number",
                      breaks = c("Rt_data", "Rt_LogP"),
                      values = c("Rt_data"="firebrick", 
                                 "Rt_LogP"="forestgreen"),
                      labels=c("Observed", "Predicted")) +
  theme_classic() +
  labs(title = "Observed and Predicted Rt in all LTLAs over time",
        x = "Date",
        y = "Reproduction number") +
  theme(
     plot.title = element_text(face="bold", hjust = 0.5),
     axis.title.x = element_text(face="bold"),
     axis.title.y = element_text(face="bold"),
     legend.position = "bottom",
     legend.title = element_text(face="bold")) +
  if (DoKnots == 1) {
  geom_vline(xintercept = as.Date(Steps, format = "%d/%m/%Y"),
             color = "darkmagenta")}

dev.off()

 
#### Reg Trend ####
 
png(paste0("Figures/", model_name, "/RegionalTrend.png"),
    width = 10, height = 6, units = 'in', res = 600)
 
ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = RegionalTrends, group = date), 
                 size = rel(0.5)) +
  theme_classic() +
   labs(title = "Regional Trends in all LTLAs over time",
        x = "Date",
        y = "Secular Trend") +
  theme(
     plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(0.9), face="bold"),
     axis.title.y = element_text(size = rel(0.9), face="bold"),
     axis.text = element_text(size=rel(0.7)),
     legend.title = element_text(size = rel(0.9), face="bold"),
     legend.text = element_text(size=rel(0.7)))+
  if (DoKnots == 1) {
    geom_vline(xintercept = as.Date(Steps, format = "%d/%m/%Y"),
               color = "darkmagenta")}

dev.off()
 
 
#### National Trend ####
 
png(paste0("Figures/", model_name, "/NationalTrend.png"),
    width = 10, height = 6, units = 'in', res = 600)
 
ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = NationalTrend, group = date), 
                 size = rel(0.5)) +
  theme_classic() +
  labs(title = "National Trend over time",
        x = "Date",
        y = "National Trend") +
  theme(
     plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(0.9), face="bold"),
     axis.title.y = element_text(size = rel(0.9), face="bold"),
     axis.text = element_text(size=rel(0.7)),
     legend.title = element_text(size = rel(0.9), face="bold"),
     legend.text = element_text(size=rel(0.7)))+
  if (DoKnots == 1) {
    geom_vline(xintercept = as.Date(Steps, format = "%d/%m/%Y"),
               color = "darkmagenta")}

dev.off()
 
# Spline: National Trend plotted with knots 

if(DoKnots == 1) {
  
  png(paste0("Figures/", model_name, "/NationalTrendSpline.png"),
      width = 10, height = 6, units = 'in', res = 300)
  
  ggplot() +
    geom_point (data = sum_line,
                mapping = aes(x = date, y = Lambda, group = date), 
                size = rel(1.2), color = "red") +
    geom_line (data = sum_line,
               mapping = aes(x = date, y = Lambda), 
               color = "red") +
    geom_point (data = sum_rt,
                mapping = aes(x = date, y = NationalTrend, group = date), 
                size = rel(0.8), color = "black") +
    theme_classic() +
    labs(title = "National Trend over time and knots inputted",
         x = "Week",
         y = "National Trend") +
    theme(
      plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
      axis.title.x = element_text(size = rel(0.9), face="bold"),
      axis.title.y = element_text(size = rel(0.9), face="bold"),
      axis.text = element_text(size=rel(0.7)),
      legend.title = element_text(size = rel(0.9), face="bold"),
      legend.text = element_text(size=rel(0.7)))
  
  dev.off()
}

 
#### In one LTLA ####

some_LTLA <- sum_rt %>%
  filter(LTLA == 1 | LTLA == 47 | LTLA == 93 |
         LTLA == 104 | LTLA == 123 | LTLA == 147 |
         LTLA == 208 | LTLA == 179 | LTLA == 221)
 
png(paste0("Figures/", model_name, "/LTLA_sum_plot_all.png"),
    width = 12, height = 8, units = 'in', res = 300)

ggplot(data = some_LTLA) +
  geom_line (mapping = aes(x = date, y = NationalTrend, 
                           color = "Lambda"), linewidth = 1.3) +
  geom_line (mapping = aes(x = date, y = RegionalTrends,
                           color = "RanEffect"), linewidth = 1.3) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), linewidth = 1.3) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), linewidth = 1.3) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "Regional Trend",
                                "Predicted Rt", "Observed Rt")) +
   theme_classic() +
   labs(title = "Parameters in various LTLAs",
        x = "Date",
        y = "Parameter value") +
   theme(
     plot.title = element_text(size = rel(1.5), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(1.3), face="bold"),
     axis.title.y = element_text(size = rel(1.3), face="bold"),
     axis.text = element_text(size=rel(1.2)),
     legend.title = element_text(size = rel(1.3), face="bold"),
     legend.position = "bottom",
     legend.text = element_text(size=rel(1.2))) +
  facet_wrap(LTLA ~ .)

dev.off()

