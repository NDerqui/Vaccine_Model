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
DoVarriants <- 0


#### Extract parameters ####

Rt_data <- exp(data_stan[[20]])

LTLA <- data_stan[[18]]

Dose_1 <- data_stan[[21]][,1]
Dose_2 <- data_stan[[21]][,2]
Dose_3 <- data_stan[[21]][,3]

date <- max(c(min(data_rt$date), min(data_var$date))) + data_stan[[17]]*7

Steps <- c(max(c(min(data_rt$date), min(data_var$date))),
           lockdown_steps[2:length(lockdown_steps)],
           min(c(max(data_rt$date), max(data_var$date))))



# MODEL DATA ------------------------------------------------------------


#### Load data ####

model_name <- "Base_1A"

# fit <- readRDS(paste0(model_name, ".Rds"))
  
# loo_run <- readRDS(paste0("loo_", model_name, ".Rds"))

# [Or directly open from directory]


#### Substract parameters ####

model_matrix <- as.matrix(fit)
rm(fit)
gc()

Rt_LogP <- exp(colMeans(model_matrix[, grep(
  "LogPredictions", colnames(model_matrix))]))

RegionalTrends <- exp(colMeans(model_matrix[, grep(
  "RegionalTrends", colnames(model_matrix))]))

NationalTrend <- exp(colMeans(model_matrix[, grep(
  '^NationalTrend\\[', colnames(model_matrix))]))

Lambda <- exp(colMeans(model_matrix[, grep(
  '^lambda\\[', colnames(model_matrix))]))

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



# VAX TABLE ---------------------------------------------------------------
 
 
#### VaxEffect ####
 
VEMean <- colMeans(model_matrix[, grep(
   "VaxEffect", colnames(model_matrix))])
VEQuan <- colQuantiles(model_matrix[, grep(
   "VaxEffect", colnames(model_matrix))], probs=c(0.025,0.975))

VEMean_exp <-(exp(-VEMean))
VEQuan_exp <-(exp(-VEQuan))
 
sum_ve <- round(data.frame(VEMean, VEQuan, VEMean_exp, VEQuan_exp), digits = 2) %>%
  select("VEMean", "X2.5.", "X97.5.",
         "VEMean_exp", "X97.5..1", "X2.5..1")
colnames(sum_ve) <- c("Mean Effect", "2.5% Q", "97.5% Q",
                      "Exp(Mean Effect)", "Exp(2.5% Q)", "Exp(97.5% Q)")

if (nrow(sum_ve) == 3) {
  row.names(sum_ve) <- c("Dose 1", "Dose 2", "Dose 3")
} else {
  row.names(sum_ve) <- c("Alpha1", "Alpha2", "Alpha3", "Delta1", "Delta2", "Delta3")
}
                
write.xlsx(sum_ve, rowNames = TRUE,
            paste0("Results/", model_name, "/Vax_VE_table.xlsx"))

 
 
# PARAMETERs TABLE --------------------------------------------------------
 
 
#### Log Predictions ####
 
rt_mean <- sum_rt %>%
 select("LTLA", "date", "Rt_LogP") %>%
 pivot_wider(names_from = date, values_from = Rt_LogP) %>%
 colMeans() %>%
 as.data.frame() %>%
 rename(mean = ".")%>%
 slice(-1)
rt_max <- sum_rt %>%
 select("LTLA", "date", "Rt_LogP") %>%
 pivot_wider(names_from = date, values_from = Rt_LogP) %>%
 as.matrix() %>%
 colMaxs () %>%
 as.data.frame() %>%
 rename(max = ".") %>%
 slice(-1)
rt_min <- sum_rt %>%
 select("LTLA", "date", "Rt_LogP") %>%
 pivot_wider(names_from = date, values_from = Rt_LogP) %>%
 as.matrix() %>%
 colMins() %>%
 as.data.frame() %>%
 rename(min = ".") %>%
 slice(-1)
rt_sd <- sum_rt %>%
 select("LTLA", "date", "Rt_LogP") %>%
 pivot_wider(names_from = date, values_from = Rt_LogP) %>%
 as.matrix() %>%
 colSds() %>%
 as.data.frame() %>%
 rename(sd = ".") %>%
 slice(-1)
rt_iqr <- sum_rt %>%
 select("LTLA", "date", "Rt_LogP") %>%
 pivot_wider(names_from = date, values_from = Rt_LogP) %>%
 as.matrix() %>%
 colQuantiles(probs = seq(from = 0, to = 1, by = 0.25)) %>%
 as.data.frame() %>%
 slice(-1)  

rt_table <- round(data.frame(rt_mean, rt_sd,
                          rt_min, rt_iqr, rt_max), digits = 4)
write.xlsx(rt_table, rowNames = TRUE,
          paste0("Results/", model_name, "/Rt_LogP.xlsx"))

 
#### Reg Trends ####
 
ref_mean <- sum_rt %>%
 select("LTLA", "date", "RegionalTrends") %>%
 pivot_wider(names_from = date, values_from = RegionalTrends) %>%
 colMeans() %>%
 as.data.frame() %>%
 rename(mean = ".")%>%
 slice(-1)
ref_max <- sum_rt %>%
 select("LTLA", "date", "RegionalTrends") %>%
 pivot_wider(names_from = date, values_from = RegionalTrends) %>%
 as.matrix() %>%
 colMaxs () %>%
 as.data.frame() %>%
 rename(max = ".") %>%
 slice(-1)
ref_min <- sum_rt %>%
 select("LTLA", "date", "RegionalTrends") %>%
 pivot_wider(names_from = date, values_from = RegionalTrends) %>%
 as.matrix() %>%
 colMins() %>%
 as.data.frame() %>%
 rename(min = ".") %>%
 slice(-1)
ref_sd <- sum_rt %>%
 select("LTLA", "date", "RegionalTrends") %>%
 pivot_wider(names_from = date, values_from = RegionalTrends) %>%
 as.matrix() %>%
 colSds() %>%
 as.data.frame() %>%
 rename(sd = ".") %>%
 slice(-1)
ref_iqr <- sum_rt %>%
 select("LTLA", "date", "RegionalTrends") %>%
 pivot_wider(names_from = date, values_from = RegionalTrends) %>%
 as.matrix() %>%
 colQuantiles(probs = seq(from = 0, to = 1, by = 0.25)) %>%
 as.data.frame() %>%
 slice(-1)  
 
ref_table <- round(data.frame(ref_mean, ref_sd,
                             ref_min, ref_iqr, ref_max), digits = 4)
write.xlsx(ref_table, rowNames = TRUE,
          paste0("Results/", model_name, "/RegionalTrends.xlsx"))
 
 
#### National Trends ####
 
lam_table <- sum_rt %>%
 select("LTLA", "date", "NationalTrend") %>%
 pivot_wider(names_from = date, values_from = NationalTrend) %>%
 colMeans() %>%
 as.data.frame() %>%
 rename(lambda = ".")%>%
 slice(-1)
 
write.xlsx(lam_table, rowNames = TRUE,
          paste0("Results/", model_name, "/NationalTrend.xlsx"))
 
 

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
    width = 12, height = 6, units = 'in', res = 600)
 
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
     plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(0.9), face="bold"),
     axis.title.y = element_text(size = rel(0.9), face="bold"),
     axis.text = element_text(size=rel(0.7)),
     legend.position = "bottom",
     legend.title = element_text(size = rel(0.9), face="bold"),
     legend.text = element_text(size=rel(0.7))) +
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
 
png(paste0("Figures/", model_name, "/LTLA_sum_plot_1.png"), 
    width = 10, height = 6, units = 'in', res = 300)
 
Plot_sum <- ggplot(data = sum_rt[sum_rt$LTLA == 1,]) +
   geom_line (mapping = aes(x = date, y = NationalTrend, 
                            color = "Lambda"), size =rel(1)) +
   geom_line (mapping = aes(x = date, y = RegionalTrends,
                            color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
   geom_line (mapping = aes(x = date, y = Rt_data,
                            color = "RtData"), size =rel(1)) +
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
   labs(title = "Parameters in LTLA 1",
        x = "Date",
        y = "Parameter value") +
   theme(
     plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(0.9), face="bold"),
     axis.title.y = element_text(size = rel(0.9), face="bold"),
     axis.text = element_text(size=rel(0.7)),
     legend.title = element_text(size = rel(0.9), face="bold"),
     legend.text = element_text(size=rel(0.7)))

dev.off()
 
png(paste0("Figures/", model_name, "/LTLA_sum_plot_2.png"),
    width = 10, height = 6, units = 'in', res = 300)
 
Plot_sum_2 <- ggplot(data = sum_rt[sum_rt$LTLA == 147,]) +
  geom_line (mapping = aes(x = date, y = NationalTrend, 
                           color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = RegionalTrends,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
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
   labs(title = "Parameters in LTLA 147",
        x = "Date",
        y = "Parameter value") +
   theme(
     plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(0.9), face="bold"),
     axis.title.y = element_text(size = rel(0.9), face="bold"),
     axis.text = element_text(size=rel(0.7)),
     legend.title = element_text(size = rel(0.9), face="bold"),
     legend.text = element_text(size=rel(0.7)))
 
dev.off()
 
png(paste0("Figures/", model_name, "/LTLA_sum_plot_3.png"),
    width = 10, height = 6, units = 'in', res = 300)
 
Plot_sum_3 <- ggplot(data = sum_rt[sum_rt$LTLA == 208,]) +
  geom_line (mapping = aes(x = date, y = NationalTrend, 
                           color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = RegionalTrends,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
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
   labs(title = "Parameters in LTLA 208",
        x = "Date",
        y = "Parameter value") +
   theme(
     plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(0.9), face="bold"),
     axis.title.y = element_text(size = rel(0.9), face="bold"),
     axis.text = element_text(size=rel(0.7)),
     legend.title = element_text(size = rel(0.9), face="bold"),
     legend.text = element_text(size=rel(0.7)))

dev.off()
 
png(paste0("Figures/", model_name, "/LTLA_sum_plot_4.png"),
    width = 10, height = 6, units = 'in', res = 300)
 
Plot_sum_4 <- ggplot(data = sum_rt[sum_rt$LTLA == 47,]) +
  geom_line (mapping = aes(x = date, y = NationalTrend, 
                           color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = RegionalTrends,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
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
   labs(title = "Parameters in LTLA 47",
        x = "Date",
        y = "Parameter value") +
   theme(
     plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
     axis.title.x = element_text(size = rel(0.9), face="bold"),
     axis.title.y = element_text(size = rel(0.9), face="bold"),
     axis.text = element_text(size=rel(0.7)),
     legend.title = element_text(size = rel(0.9), face="bold"),
     legend.text = element_text(size=rel(0.7)))

dev.off()
 
library(ggpubr)
 
png(paste0("Figures/", model_name, "/LTLA_sum_plot_all.png"),
    width = 10, height = 6, units = 'in', res = 300)
 
ggarrange(Plot_sum,Plot_sum_2, Plot_sum_3, Plot_sum_4,
           common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 2)
dev.off()

