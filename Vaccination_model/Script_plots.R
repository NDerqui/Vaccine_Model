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


# Let's do results' retrieval in a loop
# Extract data and plot in a loop

# Imp: need to adjust for options like model with/without variants
# Loop options account for the variant but not Age or Spline options


# Insert models names

names <- c("Checks_2k_8cha_2D", "Checks_5k_10cha_2D",
           "Checks_10k_10cha_2D", "Checks_10k_20cha_2D", "Checks_20k_10cha_2D")

# Insert Variant status

variants <- c(rep(1, times = 5))


#### Loop

for (i in 1:length(names)) {
  
  
  #### Load data ####
  
  DoVariants <- variants[i]
  
  model_name <- names[i]
  
  list_result <- readRDS(paste0("C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Analyses/Models_BackUp/", model_name, ".Rds"))
  
  
  #### Substract parameters ####
  
  # Parameters
  
  VaxEffect_data <- list_result[[2]]
  VaxEffect <- list_result[[2]][1]
  
  if (DoVariants == 1) {
    VarAdvantage_data <- list_result[[3]]
    VarAdvantage <- list_result[[3]][1]
  } else {
    VarAdvantage_data <- data.frame(Null = NA)
  }
  
  Rt_Predictions_data <- list_result[[4]]
  Rt_Predictions <- list_result[[4]][1]
  
  RegionalTrends_data <- list_result[[5]]
  RegionalTrends <- list_result[[5]][1]
  
  NationalTrend_data <- list_result[[6]]
  NationalTrend <- list_result[[6]][1]
  
  Lambda_data <- list_result[[7]]
  Lambda <- list_result[[7]][1]
  
  sum_rt <- data.frame(LTLA, date, Rt_data,
                       Rt_Predictions_data, RegionalTrends_data, NationalTrend_data,
                       Dose_1, Dose_2, Dose_3)
  
  
  
  # DIRECTORY ---------------------------------------------------------------
  
  dir.create(here(paste0("Figures/", model_name)),recursive = TRUE)
  
  
  
  # TABLES ---------------------------------------------------------------
  
  
  #### VaxEffect ####
  
  # Change rownames for easier interpretation
  
  if (nrow(VaxEffect_data) == 3) {
    row.names(VaxEffect_data) <- c("Dose 1", "Dose 2", "Dose 3")
    
  } else {
    
    if(nrow(VaxEffect_data) == 6) {
      row.names(VaxEffect_data) <- c("Alpha1", "Alpha2", "Alpha3",
                                     "Delta1", "Delta2", "Delta3")
    } else {
      
      if (DoAge == 1) {
        row.names(VaxEffect_data) <- c("15-49_D1", "15-49_D2", "15-49_D3",
                                       "50-69_D1", "50-69_D2", "50-69_D3",
                                       "70+_D1", "70+_D2", "70+_D3")
      } else {
        row.names(VaxEffect_data) <- c("PreAl_1", "PreAl_2", "PreAl_3",
                                       "Alpha1", "Alpha2", "Alpha3",
                                       "Delta1", "Delta2", "Delta3")
      }
    }
  }
  
  
  #### Save ####
  
  table_list <- list("VaxEffect" = VaxEffect_data,
                     "VarAdvantage" = VarAdvantage_data,
                     "RtPredictions" = Rt_Predictions_data,
                     "RegionalTrends" = RegionalTrends_data,
                     "NationalTrends" = NationalTrend_data,
                     "Lambda" = Lambda_data)
  
  write.xlsx(table_list, rowNames = TRUE,
             file = paste0("Results/", model_name, "_results_table.xlsx"))
  
  
  
  
  # PLOTS -------------------------------------------------------------------
  
  
  #### LogPred vs Obs ####
  
  # General
  
  p <- ggplot(data = sum_rt) +
    geom_point(mapping = aes(x = Rt_data, y = Rt), size = rel(0.8)) +
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
  
  
  png(paste0("Figures/", model_name, "/Obs_Pre_Rt.png"),
      width = 9, height = 8, units = 'in', res = 600)
  print(p)
  dev.off()
  
  # As a scatter plot
  
  p <- ggplot(data = sum_rt) +
    geom_boxplot (mapping = aes(x = date, y = Rt_data, group = date,
                                color = "Rt_data"), size = rel(0.5)) +
    geom_boxplot (mapping = aes(x = date, y = Rt, group = date,
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
  
  png(paste0("Figures/", model_name, "/Pre_Obs_Rt_time.png"),
      width = 9, height = 5, units = 'in', res = 600)
  print(p)
  dev.off()
  
  
  #### National Trend ####
  
  # Spline: National Trend plotted with knots 
  
  if(DoKnots == 1) {
    
    p <- ggplot() +
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
    
    png(paste0("Figures/", model_name, "/NationalTrendSpline.png"),
        width = 10, height = 6, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  
  
  #### In one LTLA ####
  
  some_LTLA <- sum_rt %>%
    filter(LTLA == 1 | LTLA == 47 | LTLA == 93 |
             LTLA == 104 | LTLA == 123 | LTLA == 147 |
             LTLA == 208 | LTLA == 179 | LTLA == 221)
  
  
  p <- ggplot(data = some_LTLA) +
    geom_line (mapping = aes(x = date, y = NationalTrend, 
                             color = "Lambda"), linewidth = 1.3) +
    geom_line (mapping = aes(x = date, y = RegionalTrends,
                             color = "RanEffect"), linewidth = 1.3) +
    geom_line (mapping = aes(x = date, y = Rt,
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
  
  png(paste0("Figures/", model_name, "/LTLA_sum_plot_all.png"),
      width = 12, height = 8, units = 'in', res = 300)
  print(p)
  dev.off()
  
  p <- ggplot(data = some_LTLA) +
    geom_line (mapping = aes(x = date, y = NationalTrend, 
                             color = "Lambda")) +
    geom_ribbon(mapping = aes(x = date, ymin = X2.5..2, ymax = X97.5..2, 
                              fill = "Lambda"), alpha = 0.2) +
    geom_line (mapping = aes(x = date, y = RegionalTrends,
                             color = "RanEffect")) +
    geom_ribbon(mapping = aes(x = date, ymin = X2.5..1, ymax = X97.5..1,
                              fill = "RanEffect"), alpha = 0.2) +
    geom_line (mapping = aes(x = date, y = Rt,
                             color = "RtLogP")) +
    geom_ribbon(mapping = aes(x = date, ymin = X2.5., ymax = X97.5.,
                              fill = "RtLogP"), alpha = 0.2) +
    geom_line (mapping = aes(x = date, y = Rt_data,
                             color = "RtData")) +
    scale_color_manual(name = "Parameter",
                       breaks = c("Lambda", "RanEffect",
                                  "RtLogP", "RtData"),
                       values = c("Lambda" = "navy",
                                  "RanEffect" = "lightgreen",
                                  "RtLogP" = "forestgreen",
                                  "RtData" = "firebrick"),
                       labels = c("National Trend", "Regional Trend",
                                  "Predicted Rt", "Observed Rt")) +
    scale_fill_manual(name = "Parameter",
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
  
  png(paste0("Figures/", model_name, "/LTLA_plot.png"),
      width = 12, height = 8, units = 'in', res = 300)
  print(p)
  dev.off()

}


#### Summary VE effects ####

## Now that results have been extracted, we can loop through some models
## to summarise VE estimates together

library(rcartocolor)

all_ve <- data.frame()

models <- c("1A", "1B", "1C", "1D",
            "2A", "2B", "2C", "2D",
            "3A", "3B", "3C", "3D")

for (variation in 1:length(models)) {
  
  add <- read.xlsx(xlsxFile = paste0("Results/CorrectVar_10k_20_sum_", models[variation], "_results_table.xlsx"),
                   sheet = "VaxEffect")
  
  add$model <- paste0(models[variation])
  
  if(models[variation] == "2A" | models[variation] == "2B" |
     models[variation] == "2C" | models[variation] == "2D") {
    
    add[,1] <- rep(c("Dose 1", "Dose 2", "Dose 3"), times = 3)}
  
  colnames(add) <- c("variant", "ve", "low", "upp", "model")
  
  all_ve <- rbind(all_ve, add)
}

all_ve[all_ve == "PreAl_1"] <- "WT1"
all_ve[all_ve == "PreAl_2"] <- "WT2"
all_ve[all_ve == "PreAl_3"] <- "WT3"

### Model 1

png(paste0("Figures/VE_sum_model1.png"),
    width = 8, height = 5, units = 'in', res = 1200)

ggplot(data = filter(all_ve, model == "1A" |  model == "1B" |
                       model == "1C" |  model == "1D")) +
  geom_point(mapping = aes(x = variant, y = ve, color = model),
             position = position_dodge(0.5)) +
  geom_errorbar(mapping = aes(x = variant, y = ve, color = model,
                              ymin = low, ymax = upp),
                position = position_dodge(0.5), width = 0.2) +
  scale_color_manual(values = carto_pal(name = "Safe")) +
  labs(title = "Model 1", y = "VE estimates", color = "Model") +
  theme_bw() +
  theme(title = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom")

dev.off()

### Model 2

png(paste0("Figures/VE_sum_model2.png"),
    width = 8, height = 5, units = 'in', res = 1200)

ggplot(data = filter(all_ve, model == "2A" |  model == "2B" |
                       model == "2C" |  model == "2D")) +
  geom_point(mapping = aes(x = variant, y = ve, color = model),
             position = position_dodge(0.5)) +
  geom_errorbar(mapping = aes(x = variant, y = ve, color = model,
                              ymin = low, ymax = upp),
                position = position_dodge(0.5), width = 0.2) +
  scale_color_manual(values = carto_pal(name = "Safe")) +
  labs(title = "Model 2", y = "VE estimates", color = "Model") +
  theme_bw() +
  theme(title = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom")

dev.off()

### Model 3

png(paste0("Figures/VE_sum_model3.png"),
    width = 10, height = 6, units = 'in', res = 1200)

ggplot(data = filter(all_ve, model == "3A" |  model == "3B" |
                       model == "3C" |  model == "3D")) +
  geom_point(mapping = aes(x = variant, y = ve, color = model),
             position = position_dodge(0.5)) +
  geom_errorbar(mapping = aes(x = variant, y = ve, color = model,
                              ymin = low, ymax = upp),
                position = position_dodge(0.5), width = 0.2) +
  scale_color_manual(values = carto_pal(name = "Safe")) +
  labs(title = "Model 3", y = "VE estimates", color = "Model") +
  theme_bw() +
  theme(title = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom")

dev.off()