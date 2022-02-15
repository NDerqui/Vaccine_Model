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



# DATA: OBS  ------------------------------------------------------------

# data_model <- readRDS("data_model_for_plots.Rds")


#### Overview ####

names(data_model)
View(data_model)
dim(data_model)



# DATA: MODEL ------------------------------------------------------------


#### Load data ####

model_name <- ""

# fit <- readRDS(paste0(model_name, ".Rds"))
  
# loo_run <- readRDS(paste0("loo_", model_name, ".Rds"))


#### Substract parameters ####

model_matrix <- as.matrix(fit)

dim(model_matrix)
View(model_matrix)

Rt_data <- exp(data_model$Rt)

Rt_LogP <- exp(colMeans(model_matrix[, grep(
  "LogPredictions", colnames(model_matrix))]))

Ran_Eff <- exp(colMeans(model_matrix[, grep(
  "random_effects", colnames(model_matrix))]))

Lambda <- exp(colMeans(model_matrix[, grep(
  '^lambda\\[', colnames(model_matrix))]))

LambdaParameter <- exp(colMeans(model_matrix[, grep(
  '^lambda_parameters\\[', colnames(model_matrix))]))

Gamma <- exp(colMeans(model_matrix[, grep(
  '^gamma\\[', colnames(model_matrix))]))

Intercept <- exp(colMeans(model_matrix[, grep(
  '^intercept\\[', colnames(model_matrix))]))

Origin <- (colMeans(model_matrix[, grep(
  "origin", colnames(model_matrix))]))

Slope <- (colMeans(model_matrix[, grep(
  "slope", colnames(model_matrix))]))

sum_rt <- data.frame(Rt_data, Rt_LogP, Ran_Eff, LambdaParameter,
                     LTLA = data_model$LTLAs,
                     Dose_1 = data_model$First_Prop,
                     Dose_2 = data_model$Second_Prop,
                     Dose_3 = data_model$Third_Prop,
                     date = data_model$date,
                     week = data_model$week,
                     row.names = paste0("Rt", 1:12726))
sum_line <- data.frame(date = rep(Steps, times = 303),
                       Lambda = Lambda)

sum_rt_1 <- sum_rt[1:42,]



# DIRECTORY ---------------------------------------------------------------

# dir.create(here(paste0("Figures/", model_name)),recursive = TRUE)
# dir.create(here(paste0("Results/", model_name)),recursive = TRUE) 



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

library(formattable)

formattable(sum_ve)



# PLOTS -------------------------------------------------------------------


#### LogPred vs Obs ####

# General

#png("Figures/Obs_Pre_Rt.png", width = 6, height = 4, units = 'in', res = 300)

Rt_Obs_Pre <- ggplot(data = sum_rt) +
  geom_point(mapping = aes(x = Rt_data, y = Rt_LogP), size = rel(0.8)) +
  geom_abline(x = 0, y = 1, col = "red") +
  theme_classic() +
  labs(title = "Observed vs Predicted Rt across all LTLAs",
       x = "Rt Observed",
       y = "Rt Predicted") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Rt_Obs_Pre

dev.off()

# Compared to Vaccination proportion

p1 <- ggplot(data = sum_rt, mapping = aes(x = Rt_data, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_1), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
  xlab("Rt Observed") + ylab("Rt Predicted") +
  xlim(0.5, 2)+ ylim(0.5, 2) + theme_classic()
p1
p2 <- ggplot(data = sum_rt, mapping = aes(x = Rt_data, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_2), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "olivedrab", limits = c(0,1)) +
  xlab("Rt Observed") + ylab("Rt Predicted") +
  xlim(0.5, 2)+ ylim(0.5, 2) + theme_classic()
p2
p3 <- ggplot(data = sum_rt, mapping = aes(x = Rt_data, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_3), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "darkred", limits = c(0,1)) +
  xlab("Rt Observed") + ylab("Rt Predicted") +
  xlim(0.5, 2)+ ylim(0.5, 2) + theme_classic()
p3

#png("Figures/Obs_Pre_Rt_Vax.png", width = 6, height = 4, units = 'in', res = 300)

Rts_VaxProp <- ggarrange(p1, p2, p3,
                         ncol = 3, nrow = 1)
Rts_VaxProp

dev.off()


#### LogPred vs Time ####

# General

#png("Figures/Pre_Rt_time.png", width = 6, height = 4, units = 'in', res = 300)

Rt_Pre_Date <- ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = Rt_LogP, group = date), size = rel(0.5)) +
  theme_classic() +
  labs(title = "Predicted Rt over time",
       x = "Date",
       y = "Rt Predicted") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Rt_Pre_Date

dev.off()

# With vaccination proportion

p7 <- ggplot(data = sum_rt, mapping = aes(x = date, y = Rt_LogP)) +
  geom_point(mapping = aes(colour = Dose_1), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
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
  scale_colour_gradient(low = "white", high = "olivedrab", limits = c(0,1)) +
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
  scale_colour_gradient(low = "white", high = "darkred", limits = c(0,1)) +
  xlab("Date") + ylab("Rt Predicted") +
  ylim(0.5, 2) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
p9

#png("Figures/Pre_Rt_time_Vax.png", width = 6, height = 6, units = 'in', res = 300)

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

#png("Figures/SES.png", width = 6, height = 4, units = 'in', res = 300)

SES = SES + geom_point(
  data = data.frame(x = data_model$date[data_model$ltla_name==NamesLTLAs[i]],
                    y = RE_uniq), aes(x=x,y=y), alpha=1, col = "red", size = rel(0.8)) +
  theme_classic() +
  labs(title = "SES across all LTLAs",
       x = "Date",
       y = "Secular Effect Size") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
SES

dev.off()


#### Ran Effect ####

#png("Figures/RanEffect.png", width = 6, height = 4, units = 'in', res = 300)

RanEff_time <- ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = Ran_Eff, group = date), 
                size = rel(0.5)) +
  theme_classic() +
  labs(title = "Random Effects: Lambda*Gamma*Intercept",
       x = "Date",
       y = "Random Effects") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
RanEff_time

dev.off()


#### Lambda ####

#png("Figures/Lambda.png", width = 6, height = 4, units = 'in', res = 300)

Lambda_time <- ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = Lambda, group = date), 
                size = rel(0.5)) +
  ylim(0, 2) +
  theme_classic() +
  labs(title = "Lambda over time in all LTLAs",
       x = "Date",
       y = "Lambda") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Lambda_time

dev.off()


#### In one LTLA ####

#png("Figures/LTLA_sum_plot.png", width = 6, height = 4, units = 'in', res = 300)

Plot_sum <- ggplot(data = sum_rt_1) +
  geom_point (mapping = aes(x = date, y = Lambda, 
                           color = "Lambda"), size =rel(1.5)) +
  geom_line (mapping = aes(x = date, y = Lambda*(Gamma[1]),
                           color = "Lambda*Gamma"), size =rel(1.5)) +
  geom_line (mapping = aes(x = date, y = Ran_Eff,
                           color = "RanEffect"), size =rel(1.5)) +
  geom_line (mapping = aes(x = date, y = Rt_data,
                           color = "RtData"), size =rel(1.5)) +
  geom_line (mapping = aes(x = date, y = Rt_LogP,
                           color = "RtLogP"), size =rel(1.5)) +
  scale_color_manual(name = "Parameter",
                     breaks = c("Lambda", "Lambda*Gamma", "RanEffect",
                                "RtData", "RtLogP"),
                     values = c("Lambda" = "olivedrab",
                                "Lambda*Gamma" = "darkblue",
                                "RanEffect" = "salmon",
                                "RtData" = "grey",
                                "RtLogP" = "black")) +
  theme_classic() +
  labs(title = "Parameters in example LTLA",
       x = "Date",
       y = "Parameter value") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Plot_sum

dev.off()



#### Playing with linear intercept ####

Lamda_knots <- as.numeric(Lambda[1:5])

Lambda_dummy <- ggplot() +
  geom_point (data = sum_line,
              mapping = aes(x = date, y = Lambda, group = date), 
                size = rel(1.2), color = "red") +
  geom_line (data = sum_line,
              mapping = aes(x = date, y = Lambda), 
              color = "red") +
  theme_classic() +
  labs(title = "Lambda over time in all LTLAs",
       x = "Date",
       y = "Lambda") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Lambda_dummy

Lambda_dummy_2 <- ggplot() +
  geom_point (data = sum_line,
              mapping = aes(x = date, y = Lambda, group = date), 
              size = rel(1.2), color = "red") +
  geom_line (data = sum_line,
             mapping = aes(x = date, y = Lambda), 
             color = "red") +
  geom_point (data = sum_rt,
              mapping = aes(x = date, y = LambdaParameter, group = date), 
              size = rel(0.8), color = "black") +
  theme_classic() +
  labs(title = "Lambda over time in all LTLAs",
       x = "Week",
       y = "Lambda") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Lambda_dummy_2

Lambda_dummy <- ggplot() +
  geom_point (data = sum_line,
              mapping = aes(x = date, y = Lambda, group = date), 
              size = rel(1.2), color = "red") +
  geom_line (data = sum_line,
             mapping = aes(x = date, y = Lambda), 
             color = "red") +
  geom_abline(slope = Slope[1], intercept = Origin[1], color = "blue") +
  geom_abline(slope = Slope[2], intercept = Origin[2], color = "green") +
  geom_abline(slope = Slope[3], intercept = Origin[3], color = "gray") +
  geom_abline(slope = Slope[4], intercept = Origin[4], color = "black") +
  theme_classic() +
  labs(title = "Lambda over time in all LTLAs",
       x = "Week",
       y = "Lambda") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Lambda_dummy



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
