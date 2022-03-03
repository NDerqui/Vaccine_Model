###SARS-CoV-2 Rt Model for Vaccination and Variants###

        ###Model update exploratory script###
     ##Proposal 2: Multiple doses, one variant##

#MRes Biomedical Research - EECID stream - Project 1#



# SET UP ------------------------------------------------------------------

rm(list = ls())

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



# DATA: OBS  ------------------------------------------------------------

# data_model <- readRDS("data_model_for_plots.Rds")


#### Overview ####

names(data_model)
View(data_model)
dim(data_model)



# DATA: MODEL ------------------------------------------------------------


#### Load data ####

model_name <- "fit_2000_8_base"

# fit <- readRDS(paste0(model_name, ".Rds"))
  
# loo_run <- readRDS(paste0("loo_", model_name, ".Rds"))


#### Substract parameters ####

model_matrix <- as.matrix(fit)

#dim(model_matrix)
#View(model_matrix)

Rt_data <- exp(data_model$Rt)

Rt_LogP <- exp(colMeans(model_matrix[, grep(
  "LogPredictions", colnames(model_matrix))]))

Ran_Eff <- exp(colMeans(model_matrix[, grep(
  "random_effects", colnames(model_matrix))]))

Lambda <- exp(colMeans(model_matrix[, grep(
  "lambda", colnames(model_matrix))]))

Gamma <- exp(colMeans(model_matrix[, grep(
  "gamma", colnames(model_matrix))]))

Intercept <- exp(colMeans(model_matrix[, grep(
  "intercept", colnames(model_matrix))]))

sum_rt <- data.frame(Rt_data, Rt_LogP, Ran_Eff,
                     Lambda, LTLA = data_model$LTLAs,
                     Dose_1 = data_model$First_Prop,
                     Dose_2 = data_model$Second_Prop,
                     Dose_3 = data_model$Third_Prop,
                     date = data_model$date,
                     row.names = paste0("Rt", 1:12726))



# DIRECTORY ---------------------------------------------------------------

dir.create(here(paste0("Figures/", model_name)),recursive = TRUE)
dir.create(here(paste0("Results/", model_name)),recursive = TRUE) 



# VAX TABLE ---------------------------------------------------------------


#### VaxEffect ####

VEMean <- colMeans(model_matrix[, grep(
  "VaxEffect", colnames(model_matrix))])
VEQuan <- colQuantiles(model_matrix[, grep(
  "VaxEffect", colnames(model_matrix))], probs=c(0.025,0.975))

sum_ve <- round(data.frame(VEMean, VEQuan,
                           row.names = c("Dose 1", "Dose 2", "Dose 3")),
                digits = 4)
colnames(sum_ve) <- c("Mean Effect", "2.5% Q", "97.5% Q")
View(sum_ve)

write.xlsx(sum_ve, row.names = TRUE,
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
write.xlsx(rt_table, row.names = TRUE,
           paste0("Results/", model_name, "/Rt_LogP.xlsx"))


#### RanEffects ####

ref_mean <- sum_rt %>%
  select("LTLA", "date", "Ran_Eff") %>%
  pivot_wider(names_from = date, values_from = Ran_Eff) %>%
  colMeans() %>%
  as.data.frame() %>%
  rename(mean = ".")%>%
  slice(-1)
ref_max <- sum_rt %>%
  select("LTLA", "date", "Ran_Eff") %>%
  pivot_wider(names_from = date, values_from = Ran_Eff) %>%
  as.matrix() %>%
  colMaxs () %>%
  as.data.frame() %>%
  rename(max = ".") %>%
  slice(-1)
ref_min <- sum_rt %>%
  select("LTLA", "date", "Ran_Eff") %>%
  pivot_wider(names_from = date, values_from = Ran_Eff) %>%
  as.matrix() %>%
  colMins() %>%
  as.data.frame() %>%
  rename(min = ".") %>%
  slice(-1)
ref_sd <- sum_rt %>%
  select("LTLA", "date", "Ran_Eff") %>%
  pivot_wider(names_from = date, values_from = Ran_Eff) %>%
  as.matrix() %>%
  colSds() %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  slice(-1)
ref_iqr <- sum_rt %>%
  select("LTLA", "date", "Ran_Eff") %>%
  pivot_wider(names_from = date, values_from = Ran_Eff) %>%
  as.matrix() %>%
  colQuantiles(probs = seq(from = 0, to = 1, by = 0.25)) %>%
  as.data.frame() %>%
  slice(-1)  

ref_table <- round(data.frame(ref_mean, ref_sd,
                              ref_min, ref_iqr, ref_max), digits = 4)
write.xlsx(ref_table, row.names = TRUE,
           paste0("Results/", model_name, "/Ran_Eff.xlsx"))


#### Lambda ####

lam_table <- sum_rt %>%
  select("LTLA", "date", "Lambda") %>%
  pivot_wider(names_from = date, values_from = Lambda) %>%
  colMeans() %>%
  as.data.frame() %>%
  rename(lambda = ".")%>%
  slice(-1)

write.xlsx(lam_table, row.names = TRUE,
           paste0("Results/", model_name, "/Lambda.xlsx"))


#### Gamma ####

gamma_table <- Gamma %>%
  as.data.frame(row.names = paste0("LTLA", 1:303)) %>%
  rename(gamma = ".")
  
write.xlsx(gamma_table, row.names = TRUE,
           paste0("Results/", model_name, "/Gamma.xlsx"))


#### Intercept ####

intercept_table <- Intercept %>%
  as.data.frame(row.names = paste0("LTLA", 1:303)) %>%
  rename(intercept = ".")

write.xlsx(intercept_table, row.names = TRUE,
           paste0("Results/", model_name, "/Intercept.xlsx"))



# PLOTS -------------------------------------------------------------------


#### LogPred vs Obs ####

# General

png(paste0("Figures/", model_name, "/Obs_Pre_Rt.png"), width = 9, height = 8, units = 'in', res = 600)

Rt_Obs_Pre <- ggplot(data = sum_rt) +
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
Rt_Obs_Pre

dev.off()

# As a scatter plot

png(paste0("Figures/", model_name, "/Pre_Obs_Rt_time.png"), width = 10, height = 6, units = 'in', res = 600)

Obs_Pre_Date <- ggplot(data = sum_rt) +
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
  labs(title = "Observed and Predicted Rt across LTLAs over time",
       x = "Date",
       y = "Reproduction number") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Obs_Pre_Date

dev.off()

# Compared to Vaccination proportion

p1 <- ggplot(data = sum_rt, mapping = aes(x = Rt_data, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_1), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "navy", limits = c(0,1)) +
  xlab("Rt Observed") + ylab("Rt Predicted") +
  xlim(0.5, 2)+ ylim(0.5, 2) + theme_classic()
p1
p2 <- ggplot(data = sum_rt, mapping = aes(x = Rt_data, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_2), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "cyan3", limits = c(0,1)) +
  xlab("Rt Observed") + ylab("Rt Predicted") +
  xlim(0.5, 2)+ ylim(0.5, 2) + theme_classic()
p2
p3 <- ggplot(data = sum_rt, mapping = aes(x = Rt_data, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_3), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "lightgreen", limits = c(0,1)) +
  xlab("Rt Observed") + ylab("Rt Predicted") +
  xlim(0.5, 2)+ ylim(0.5, 2) + theme_classic()
p3

#png(paste0("Figures/", model_name, "/Obs_Pre_Rt_Vax.png"), width = 6, height = 4, units = 'in', res = 300)

Rts_VaxProp <- ggarrange(p1, p2, p3,
                         ncol = 3, nrow = 1)
Rts_VaxProp

dev.off()


#### LogPred vs Time ####

# General

png(paste0("Figures/", model_name, "/Pre_Rt_time.png"), width = 10, height = 6, units = 'in', res = 600)

Rt_Pre_Date <- ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = Rt_LogP, group = date), size = rel(0.5)) +
  theme_classic() +
  labs(title = "Predicted Rt across LTLAs over time",
       x = "Date",
       y = "Rt Predicted") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Rt_Pre_Date

dev.off()

# With vaccination proportion

p7 <- ggplot(data = sum_rt, mapping = aes(x = date, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_1), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "navy", limits = c(0,1)) +
  xlab("Date") + ylab("Rt Predicted") +
  ylim(0.5, 2) + 
  theme_classic() +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
p7
p8 <- ggplot(data = sum_rt, mapping = aes(x = date, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_2), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "cyan3", limits = c(0,1)) +
  xlab("Date") + ylab("Rt Predicted") +
  ylim(0.5, 2) + 
  theme_classic() +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
p8
p9 <- ggplot(data = sum_rt, mapping = aes(x = date, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_3), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "palegreen", limits = c(0,1)) +
  xlab("Date") + ylab("Rt Predicted") +
  ylim(0.5, 2) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
p9

#png(paste0("Figures/", model_name, "/Pre_Rt_time_Vax.png"), width = 6, height = 6, units = 'in', res = 300)

PreRts_Date <- ggarrange(p7, p8, p9,
                         ncol = 1, nrow = 3)
PreRts_Date

dev.off()


####Trend####

RE_uniq <- Ran_Eff[data_model$ltla_name == NamesLTLAs[1]]
for(i in 2:length(NamesLTLAs))	{
  RE_uniq = RE_uniq + Ran_Eff[data_model$ltla_name == NamesLTLAs[i]]}
RE_uniq = RE_uniq/length(NamesLTLAs)
#Random Effects normalised to the number of LTLAs

SES <- ggplot()+
  geom_point(data = data.frame(x = data_model$date[data_model$ltla_name==NamesLTLAs[1]],
                               y = Ran_Eff[data_model$ltla_name==NamesLTLAs[1]]),
             aes(x=x,y=y), alpha=0.5, size = rel(0.8)) + 
  theme_classic() +
  labs(x = "Date",
       y = "Secular Effect Size") +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
#Note we only plot unique values for each LTLA

for(i in 2:length(NamesLTLAs)){
  SES <- SES + geom_point(data = data.frame(x = data_model$date[data_model$ltla_name==NamesLTLAs[i]],
                                            y = Ran_Eff[data_model$ltla_name==NamesLTLAs[i]]),
                          aes(x=x,y=y), alpha=0.5, size = rel(0.8))}
#This is the loop to do so

png(paste0("Figures/", model_name, "/SES.png"), width = 10, height = 6, units = 'in', res = 600)

SES = SES + geom_point(
  data = data.frame(x = data_model$date[data_model$ltla_name==NamesLTLAs[i]],
                    y = RE_uniq), aes(x=x,y=y), alpha=1, col = "red", size = rel(0.8)) +
  theme_classic() +
  labs(title = "Secular Trend across all LTLAs over time",
       x = "Date",
       y = "Secular Trend") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
SES

dev.off()


#### Ran Effect ####

png(paste0("Figures/", model_name, "/RanEffect.png"), width = 10, height = 6, units = 'in', res = 600)

RanEff_time <- ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = Ran_Eff, group = date), 
                size = rel(0.5)) +
  theme_classic() +
  labs(title = "Secular Trend across all LTLAs over time",
       x = "Date",
       y = "Secular Trend") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
RanEff_time

dev.off()


#### Lambda ####

png(paste0("Figures/", model_name, "/Lambda.png"), width = 10, height = 6, units = 'in', res = 600)

Lambda_time <- ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = Lambda, group = date), 
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
    legend.text = element_text(size=rel(0.7)))
Lambda_time

dev.off()


#### In one LTLA ####

png(paste0("Figures/", model_name, "/LTLA_sum_plot_1.png"), width = 10, height = 6, units = 'in', res = 300)

Plot_sum <- ggplot(data = sum_rt[sum_rt$LTLA == 1,]) +
  geom_line (mapping = aes(x = date, y = Lambda, 
                           color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[1]),
                           color = "Lambda*Gamma"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
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
Plot_sum

dev.off()

png(paste0("Figures/", model_name, "/LTLA_sum_plot_2.png"), width = 10, height = 6, units = 'in', res = 300)

Plot_sum_2 <- ggplot(data = sum_rt[sum_rt$LTLA == 147,]) +
  geom_line (mapping = aes(x = date, y = Lambda, 
                            color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[147]),
                           color = "Lambda*Gamma"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
  theme_classic() +
  labs(title = "Parameters in LTLA 2",
       x = "Date",
       y = "Parameter value") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Plot_sum_2

dev.off()

png(paste0("Figures/", model_name, "/LTLA_sum_plot_3.png"), width = 10, height = 6, units = 'in', res = 300)

Plot_sum_3 <- ggplot(data = sum_rt[sum_rt$LTLA == 208,]) +
  geom_line (mapping = aes(x = date, y = Lambda, 
                            color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[208]),
                           color = "Lambda*Gamma"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
  theme_classic() +
  labs(title = "Parameters in LTLA 3",
       x = "Date",
       y = "Parameter value") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Plot_sum_3

dev.off()

png(paste0("Figures/", model_name, "/LTLA_sum_plot_4.png"), width = 10, height = 6, units = 'in', res = 300)

Plot_sum_4 <- ggplot(data = sum_rt[sum_rt$LTLA == 299,]) +
  geom_line (mapping = aes(x = date, y = Lambda, 
                            color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[299]),
                           color = "Lambda*Gamma"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
  theme_classic() +
  labs(title = "Parameters in LTLA 4",
       x = "Date",
       y = "Parameter value") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Plot_sum_4

dev.off()

library(ggpubr)

png(paste0("Figures/", model_name, "/LTLA_sum_plot_all.png"), width = 10, height = 6, units = 'in', res = 300)

ggarrange(Plot_sum,Plot_sum_2, Plot_sum_3, Plot_sum_4,
             common.legend = TRUE, legend = "right", ncol = 2, nrow = 2)
dev.off()



# FOR THESIS --------------------------------------------------------------


#### Renaming 4LTLA plot ####

Psum <- ggplot(data = sum_rt[sum_rt$LTLA == 1,]) +
  geom_line (mapping = aes(x = date, y = Lambda, 
                            color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[1]),
                           color = "Lambda*Gamma"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
  theme_classic() +
  labs(title = "i",
       x = "Date",
       y = "Parameter value") +
  theme(
    plot.title = element_text(size = rel(1), face="bold"),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Psum_2 <- ggplot(data = sum_rt[sum_rt$LTLA == 147,]) +
  geom_line (mapping = aes(x = date, y = Lambda, 
                            color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[147]),
                           color = "Lambda*Gamma"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
  theme_classic() +
  labs(title = "ii",
       x = "Date",
       y = "Parameter value") +
  theme(
    plot.title = element_text(size = rel(1), face="bold"),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Psum_3 <- ggplot(data = sum_rt[sum_rt$LTLA == 208,]) +
  geom_line (mapping = aes(x = date, y = Lambda, 
                            color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[208]),
                           color = "Lambda*Gamma"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
  theme_classic() +
  labs(title = "iii",
       x = "Date",
       y = "Parameter value") +
  theme(
    plot.title = element_text(size = rel(1), face="bold"),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Psum_4 <- ggplot(data = sum_rt[sum_rt$LTLA == 299,]) +
  geom_line (mapping = aes(x = date, y = Lambda, 
                            color = "Lambda"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[299]),
                           color = "Lambda*Gamma"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
  theme_classic() +
  labs(title = "iv",
       x = "Date",
       y = "Parameter value") +
  theme(
    plot.title = element_text(size = rel(1), face="bold"),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))

png(paste0("Figures/", model_name, "/LTLA_sum_plot_all_abc.png"), width = 10, height = 6, units = 'in', res = 300)

ggarrange(Psum,Psum_2, Psum_3, Psum_4,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 2)

dev.off()


#### Combine w/ other ####

# Obs vs Pre

Keep_1a <- ggplot(data = sum_rt) +
  geom_point(mapping = aes(x = Rt_data, y = Rt_LogP), size = rel(0.8)) +
  geom_abline(x = 0, y = 1, col = "red") +
  theme_classic() +
  labs(title = "a",
       x = "Rt Observed",
       y = "Rt Predicted") +
  theme(
    plot.title = element_text(size = rel(1), face="bold"),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)))
Keep_1a

saveRDS(Keep_1a, paste0("Figures/Combined_figures/Data/Keep_1", model_name, ".Rds"))

png("Figures/Combined_figures/No_iter_check_Pre_Obs.png", width = 15, height = 6, units = 'in', res = 300)

ggarrange(Keep_1a ,
          Keep_1b + rremove("ylab"),
          Keep_1c + rremove("ylab"),
          #Keep_1d,
          ncol = 3, nrow = 1)
dev.off()

# Scatter

Keep_2a <- ggplot(data = sum_rt) +
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
  labs(title = "a",
       x = "Date",
       y = "Reproduction number") +
  theme(
    plot.title = element_text(size = rel(1), face="bold"),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Keep_2a

saveRDS(Keep_2a, paste0("Figures/Combined_figures/Data/Keep_2", model_name, ".Rds"))

png("Figures/Combined_figures/No_iter_check_Pre_Obs_time.png", width = 18, height = 6, units = 'in', res = 300)

ggarrange(Keep_2a ,
          Keep_2b + rremove("ylab"),
          Keep_2c + rremove("ylab"),
          #Keep_2d,
          ncol = 3, nrow = 1,
          common.legend = TRUE, legend = "bottom")
dev.off()

# Secular Trend

Keep_3a = SES + geom_point(
  data = data.frame(x = data_model$date[data_model$ltla_name==NamesLTLAs[i]],
                    y = RE_uniq), aes(x=x,y=y), alpha=1, col = "red", size = rel(0.8)) +
  theme_classic() +
  labs(title = "a",
       x = "Date",
       y = "Secular Trend") +
  theme(
    plot.title = element_text(size = rel(1), face="bold"),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Keep_3a

saveRDS(Keep_3a, paste0("Figures/Combined_figures/Data/Keep_3", model_name, ".Rds"))

png("Figures/Combined_figures/No_iter_check_SES.png", width = 18, height = 6, units = 'in', res = 300)

ggarrange(Keep_3a ,
          Keep_3b + rremove("ylab"),
          Keep_3c + rremove("ylab"),
          #Keep_3d,
          ncol = 3, nrow = 1)
dev.off()

# Single LTLAs

Keep_4c <- ggarrange(Psum + rremove("xlab"),
                     Psum_2 + rremove("xlab") + rremove("ylab"),
                     Psum_3,
                     Psum_4 + rremove("ylab"),
                     ncol = 2, nrow = 2,
                     legend = "none")
#                     common.legend = TRUE, legend = "right")
Keep_4c

saveRDS(Keep_4c, paste0("Figures/Combined_figures/Data/Keep_4", model_name, ".Rds"))

legend <- ggplot(data = sum_rt[sum_rt$LTLA == 299,]) +
  geom_line (mapping = aes(x = date, y = Lambda, color = "Lambda")) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[299]), color = "Lambda*Gamma")) +
  geom_line (mapping = aes(x = date, y = Ran_Eff, color = "RanEffect")) +
  geom_line (mapping = aes(x = date, y = Rt_LogP, color = "RtLogP")) +
  geom_line (mapping = aes(x = date, y = Rt_data, color = "RtData")) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtLogP", "RtData"),
                     values = c("Lambda" = "navy",
                                "Lambda*Gamma" = "cyan3",
                                "RanEffect" = "lightgreen",
                                "RtLogP" = "forestgreen",
                                "RtData" = "firebrick"),
                     labels = c("National Trend", "National Trend x LTLA factor",
                                "Secular Trend", "Predicted Rt", "Observed Rt")) +
  lims(x = c(0,0), y = c(0,0)) +
  theme_void() +
  theme(
    legend.position = c(0.5,0.5),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))

png("Figures/Combined_figures/No_iter_check_LTLA.png", width = 10, height = 6, units = 'in', res = 300)

ggarrange(annotate_figure(Keep_4a, fig.lab = "a", fig.lab.pos = "top.left",
                          fig.lab.size = 12, fig.lab.face = "bold"),
          annotate_figure(Keep_4b, fig.lab = "b", fig.lab.pos = "top.left",
                          fig.lab.size = 12, fig.lab.face = "bold"),
          annotate_figure(Keep_4c, fig.lab = "c", fig.lab.pos = "top.left",
                          fig.lab.size = 12, fig.lab.face = "bold"),
          legend,
          ncol = 2, nrow = 2)

dev.off()






# ARCHIVE PLOTS -----------------------------------------------------------


#### Descriptive ####

png("Figures/Vax_over_time.png", width = 6, height = 4, units = 'in', res = 300)

Vax_Date <- ggplot(data = sum_data) +
  geom_boxplot (mapping = aes(x = date, y = Dose_1, group = date),
                color = "lightskyblue", size = rel(0.5)) +
  geom_boxplot (mapping = aes(x = date, y = Dose_2, group = date),
                color = "palegreen", size = rel(0.5)) +
  geom_boxplot (mapping = aes(x = date, y = Dose_3, group = date),
                color = "salmon", size = rel(0.5)) +
  theme_classic() +
  labs(title = "Proportion of vaccinated individuals across LTLAs",
       x = "Date",
       y = "Total proportion of vaccinated population") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text = element_text(size=rel(0.9), face="bold"))
Vax_Date

dev.off()


#### Preliminary ####

library(ggnewscale)

Rt_Obs_Pre <- ggplot(data = sum_rt, mapping = aes(x = Rt_data, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_1)) +
  scale_colour_gradient(low = "lightskyblue", high = "darkblue") +
  new_scale_color() +
  geom_point(mapping = aes(colour = Dose_2)) +
  scale_colour_gradient(low = "palegreen", high = "olivedrab") +
  new_scale_color() +
  geom_point(mapping = aes(colour = Dose_3)) +
  scale_colour_gradient(low = "salmon", high = "darkred") +
  xlab("Rt Observed") + ylab("Rt Predicted") + theme_classic()
Rt_Obs_Pre

Rt_Observed <- ggplot(data = sum_rt) +
  geom_point(mapping = aes(x = Rt_data, y = Dose_1, col = "lightskyblue")) +
  geom_point(mapping = aes(x = Rt_data, y = Dose_2, col = "palegreen")) +
  geom_point(mapping = aes(x = Rt_data, y = Dose_3, col = "salmon")) +
  theme_classic() +
  labs(x = "Observed Rt",
       y = "Proportion of vaccinated population") +
  scale_colour_discrete(name  = "Vaccine Doses",
                        labels = c("Dose 1", "Dose 2", "Dose 3")) +
  theme(
    legend.title = element_text(size = rel(1), face="bold"),
    legend.text = element_text(size = rel(1)),
    axis.title.x = element_text(size = rel(1.1), face="bold"),
    axis.title.y = element_text(size = rel(1.1), face="bold"),
    axis.text = element_text(size=rel(0.9), face="bold"))
Rt_Observed

Rt_Predicted <- ggplot(data = sum_rt) +
  geom_point(mapping = aes(x = Rt_LogP, y = Dose_1, col = "lightskyblue")) +
  geom_point(mapping = aes(x = Rt_LogP, y = Dose_2, col = "palegreen")) +
  geom_point(mapping = aes(x = Rt_LogP, y = Dose_3, col = "salmon")) +
  theme_classic() +
  labs(x = "Predicted Rt",
       y = "Proportion of vaccinated population") +
  scale_colour_discrete(name  = "Vaccine Doses",
                        labels = c("Dose 1", "Dose 2", "Dose 3")) +
  theme(
    legend.title = element_text(size = rel(1), face="bold"),
    legend.text = element_text(size = rel(1)),
    axis.title.x = element_text(size = rel(1.1), face="bold"),
    axis.title.y = element_text(size = rel(1.1), face="bold"),
    axis.text = element_text(size=rel(0.9), face="bold"))
Rt_Predicted


#### 20 LTLAs ####

# Rt

png("Figures/Obs_Pre_Rt_SUB.png", width = 6, height = 4, units = 'in', res = 300)

Rt_Obs_Pre_SUB <- ggplot(data = sum_rt_20) +
  geom_point(mapping = aes(x = Rt_data, y = Rt_LogP), size = rel(0.8)) +
  geom_abline(x = 0, y = 1, col = "red") +
  theme_classic() +
  labs(title = "Observed vs Predicted Rt in 20 LTLAs",
       x = "Rt Observed",
       y = "Rt Predicted") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Rt_Obs_Pre_SUB

dev.off()

# SES

data_model_20 <- data_model[1:820,]

Names20LTLAs <- unique(data_model_20$ltla_name)

Ran_Eff_20 <- Ran_Eff[1:820]

RE_uniq_20 <- Ran_Eff_20[data_model_20$ltla_name == Names20LTLAs[1]]
for(i in 2:length(Names20LTLAs))	{
  RE_uniq_20 = RE_uniq_20 + Ran_Eff_20[data_model_20$ltla_name == Names20LTLAs[i]]}
RE_uniq_20 = RE_uniq_20/length(Names20LTLAs)
#Random Effects normalised to the number of LTLAs: in his we case, we know 20!

SES_SUB <- ggplot()+
  geom_point(data = data.frame(x = data_model_20$date[data_model_20$ltla_name==Names20LTLAs[1]],
                               y = Ran_Eff_20[data_model_20$ltla_name==Names20LTLAs[1]]),
             aes(x=x,y=y), alpha=0.5, size = rel(0.8)) + 
  theme_classic() +
  labs(x = "Date",
       y = "Secular Effect Size") +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
#Note we only plot unique values for each LTLA

for(i in 2:length(Names20LTLAs)){
  SES_SUB <- SES_SUB + geom_point(data = data.frame(x = data_model_20$date[data_model_20$ltla_name==Names20LTLAs[i]],
                                                    y = Ran_Eff_20[data_model_20$ltla_name==Names20LTLAs[i]]),
                                  aes(x=x,y=y), alpha=0.5, size = rel(0.8))}
#This is the loop to do so

png("Figures/SES_SUB.png", width = 6, height = 4, units = 'in', res = 300)

SES_SUB = SES_SUB + geom_point(
  data = data.frame(x = data_model_20$date[data_model_20$ltla_name==Names20LTLAs[i]],
                    y = RE_uniq_20), aes(x=x,y=y), alpha=1, col = "red", size = rel(0.8)) +
  theme_classic() +
  labs(title = "SES across 20 LTLAs",
       x = "Date",
       y = "Secular Effect Size") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
SES_SUB

dev.off()
