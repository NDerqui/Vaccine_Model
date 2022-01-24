###SARS-CoV-2 Rt Model for Vaccination and Variants###

        ###Model update exploratory script###
     ##Proposal 2: Multiple doses, one variant##

#MRes Biomedical Research - EECID stream - Project 1#



# SET UP ------------------------------------------------------------------

rm(list = ls())

# Data

data_vax <- read.csv("rtm_incoming_vaccination_20211116-173218-f36a1245_vacc_coverage_ltla.csv")
data_rt <- read.csv("UK_hotspot_Rt_estimates.csv")

# Packages

library(dplyr)
library(rstudioapi)
library(rstan)



# OVERVIEW DATA ------------------------------------------------------------


#### Vaccination Data ####

names(data_vax)
dim(data_vax)

# Data on 306 LTLA vaccination, per 20 age groups, per 342 dates

l_code <- length(unique(data_vax$ltla_code))

l_agemin <- length(unique(data_vax$age_band_min))

l_date <- length(unique(data_vax$vaccination_date))

dim(data_vax)
l_code*l_date*l_agemin

# As described before, when age is NA, this is the sum of all age_groups

sum(data_vax$First[is.na(data_vax$age_band_min)],
    data_vax$Second[is.na(data_vax$age_band_min)],
    data_vax$Third[is.na(data_vax$age_band_min)])
sum(data_vax$First[!is.na(data_vax$age_band_min)],
    data_vax$Second[!is.na(data_vax$age_band_min)],
    data_vax$Third[!is.na(data_vax$age_band_min)])

# Set as date

data_vax$vaccination_date <- as.Date(data_vax$vaccination_date)


#### Rt Data ####

names(data_rt)
dim(data_rt)

# Data on 391 LTLS Rt, per 338 dates

l_name <- length(unique(data_rt$area))

l_date_2 <- length(unique(data_rt$date))

dim(data_rt)
l_name*l_date_2

# Set as date

data_rt$date <- as.Date(data_rt$date)



# SET PARAMETERS ----------------------------------------------------------

# Set date parameters

Date_Start <- as.Date("01/01/2021", format = "%d/%m/%Y")
Date_End <- as.Date("31/12/2021", format = "%d/%m/%Y")

# Covariates

covar_var <- c("Var_Alpha", "Var_Delta")
covar_vax <- c("First_Prop", "Second_Prop", "Third_Prop")

# Nested and BackDate (For later Stan Model)

BackDate_Char <- "_BD"

Nested_Char <- ""



# DATA CLEANING & MERGE --------------------------------------------------


# Dummy data for merging

dvax <- data_vax
drt <- data_rt


#### Making Vax/Rt Data comparable ####

# No age groups

dvax <- filter(dvax, is.na(age_band_min))

# Dates of interest

dvax <- rename(dvax, date = vaccination_date)
dvax <- filter(dvax, date >= Date_Start)
dvax <- filter(dvax, date <= Date_End)

drt <- filter(drt, date >= Date_Start)
drt <- filter(drt, date <= Date_End)

# LTLA names

drt <- rename(drt, ltla_name = area)

# Select the vars of interest

dvax <- select(dvax, "ltla_name", "date",
                "First_Prop", "Second_Prop", "Third_Prop")
drt <- select(drt, "ltla_name", "date", "Rt")

# Data from Rt has now more obs

names(dvax)
dim(dvax)
names(drt)
dim(drt)

min(dvax$date)
min(drt$date)
max(dvax$date)
max(drt$date)
#Vax data earlier (January in), but Rt data later (December in)

length(unique(dvax$ltla_name))
length(unique(drt$ltla_name))
#Rt has data on more LTLAs


#### Merge ####

data_merge <- merge(dvax, drt, by = c("ltla_name", "date"), all = TRUE)

# Checks

dim(data_merge)
names(data_merge)


#### Final Cleaning ####

# Weekly dates: var for no. of weeks & only one obs per week

#Var for the number of weeks
data_merge$week <- round(as.numeric(floor((data_merge$date - Date_Start)/7)), digits = 0)

#Combi for unique combination of LTLA*week
data_merge <- mutate(data_merge, combi = paste0(data_merge$week, data_merge$ltla_name))

#Remove the duplicates
data_merge <- filter(data_merge, !duplicated(data_merge$combi))

# No NA

data_merge <- data_merge[complete.cases(data_merge),]

# Select cols

data_model <- select(data_merge, "ltla_name", "date", "week",
                     "First_Prop", "Second_Prop", "Third_Prop", "Rt")



# DATA FOR STAN MODEL ---------------------------------------------------------


#### No. observations ####

dim(data_model)

# No of LTLA

NumLTLAs <- length(unique(data_model$ltla_name))
NumLTLAs

# No of weeks

NumTimepoints <- length(unique(data_model$date))
NumTimepoints
length(unique(data_model$week))

# No of total obs 

NumDatapoints <- NumLTLAs * NumTimepoints
NumDatapoints
nrow(data_model)

# No of LTLAs in the data

NamesLTLAs <- unique(data_model$ltla_name)

data_model$LTLAs <- NA
for(i in 1:NumLTLAs){
  data_model$LTLAs[data_model$ltla_name == NamesLTLAs[i]] = i
}

# No of weeks per LTLA

NumWeeksByLTLA <- rep(NA, NumLTLAs)
for(i in 1: NumLTLAs){   
  
  d_sub = data_model[data_model$ltla_name == NamesLTLAs[i], ]
  
  NumWeeksByLTLA[i] = length(unique(d_sub$week))
}
NumWeeksByLTLA


#### Rt log ####

data_model$Rt <- log(data_model$Rt)


#### Stan Data ####

data_model <- select(data_model,
                     "ltla_name", "LTLAs", "date", "week", covar_vax, "Rt")

data_stan <- list(
          RtVals = data_model$Rt,
          VaxProp = data_model[,covar_vax],
          NumLTLAs = NumLTLAs,
          NumDoses = length(covar_vax),
          NumDatapoints = NumDatapoints,
          LTLAs = data_model$LTLAs,
          NumTimepoints = NumTimepoints
   )



# STAN MODEL --------------------------------------------------------------


#### Model ####

library(here)

ModelChar <- "2stage_LFM_non_centered_Vax"
StanModel <- stan_model(here(paste0("Stan_models/",ModelChar, ".stan")))
cat(paste0("Model compilation done\n"))

# Create and write meta data

ModelMetaData 				= c()
ModelMetaData$iter 			= 6000 #Increase
ModelMetaData$warmup 		= 600 #Increase
ModelMetaData$thin 			= 1
ModelMetaData$chains 		= 8 #Increase
ModelMetaData$adapt_delta 	= 0.9
ModelMetaData$max_treedepth = 15
ModelMetaData$ModelChar 	= ModelChar
ModelMetaData$BackDate_Char = BackDate_Char
ModelMetaData$Nested_Char 	= Nested_Char

ModelMetaData_dummy = as.data.frame(unlist(ModelMetaData))
colnames(ModelMetaData_dummy) = NULL


#### Run ####

memory.limit(size = 1000000)

fit = sampling(StanModel, data = data_stan, 
               iter 	= ModelMetaData$iter, 
               warmup 	= ModelMetaData$warmup, 
               thin 	= ModelMetaData$thin, 
               chains 	= ModelMetaData$chains, 
               pars 	= c("VaxEffect", "LogPredictions", "random_effects", "lambda", "gamma", "intercept"), 
               control = list(adapt_delta = ModelMetaData$adapt_delta,
                              max_treedepth = ModelMetaData$max_treedepth))



# MODEL RESULTS -----------------------------------------------------------

#saveRDS(fit, file = "fit_6000.Rds")
#saveRDS(fit, file = "C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Analyses/Models_BackUp/fit_6000.Rds")

dir.create(here("Figures"),recursive = TRUE)
dir.create(here("Results"),recursive = TRUE) 


#### Data ####

library(matrixStats)

#model_6000C8_iter <- as.matrix(fit)
dim(model_6000C8_iter)
View(model_6000C8_iter)

Rt_data <- exp(data_model$Rt)
Rt_LogP <- exp(colMeans(model_6000C8_iter[, grep(
  "LogPredictions", colnames(model_6000C8_iter))]))

sum_rt <- data.frame(Rt_data, Rt_LogP,
                     LTLA = data_model$LTLAs,
                     Dose_1 = data_model$First_Prop,
                     Dose_2 = data_model$Second_Prop,
                     Dose_3 = data_model$Third_Prop,
                     date = data_model$date,
                     row.names = paste0("Rt", 1:12423))
sum_rt_20 <- sum_rt[1:820,]


#### Vax Description ####

library(ggpubr)

# Vaccination over time

png("Figures/Vax_over_time.png", width = 6, height = 4, units = 'in', res = 300)

Vax_Date <- ggplot(data = sum_rt) +
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
    axis.text=element_text(size=rel(0.9), face="bold"))
Vax_Date

dev.off()


#### LogPredictions ####

# Observed vs Predicted

png("Figures/Obs_Pre_Rt.png", width = 6, height = 4, units = 'in', res = 300)

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

png("Figures/Obs_Pre_Rt_Vax.png", width = 6, height = 4, units = 'in', res = 300)

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

Rts_VaxProp <- ggarrange(p1, p2, p3,
                         ncol = 3, nrow = 1)
Rts_VaxProp

dev.off()


#### Rt Over Time ####

# General

png("Figures/Obs_Rt_time.png", width = 6, height = 4, units = 'in', res = 300)

Rt_Obs_Date <- ggplot(data = sum_rt) +
  geom_boxplot (mapping = aes(x = date, y = Rt_data, group = date), size = rel(0.5)) +
  theme_classic() +
  labs(title = "Observed Rt over time",
       x = "Date",
       y = "Rt Observed") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold"),
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
Rt_Obs_Date

dev.off()

png("Figures/Pre_Rt_time.png", width = 6, height = 4, units = 'in', res = 300)

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

p4 <- ggplot(data = sum_rt, mapping = aes(x = date, y = Rt_data)) +
  geom_point(mapping = aes(colour = Dose_1), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
  xlab("Date") + ylab("Rt Observed") +
  ylim(0.5, 2) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
p4
p5 <- ggplot(data = sum_rt, mapping = aes(x = date, y = Rt_data)) +
  geom_point(mapping = aes(colour = Dose_2), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "olivedrab", limits = c(0,1)) +
  xlab("Date") + ylab("Rt Observed") +
  ylim(0.5, 2) + 
  theme_classic() +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
p5
p6 <- ggplot(data = sum_rt, mapping = aes(x = date, y = Rt_data)) +
  geom_point(mapping = aes(colour = Dose_3), size = rel(0.8)) +
  scale_colour_gradient(low = "white", high = "darkred", limits = c(0,1)) +
  xlab("Date") + ylab("Rt Observed") +
  ylim(0.5, 2) + 
  theme_classic() +
  theme(
    axis.title.x = element_text(size = rel(1), face="bold"),
    axis.title.y = element_text(size = rel(1), face="bold"),
    axis.text=element_text(size=rel(0.9), face="bold"))
p6

png("Figures/Obs_Rt_time_Vax.png", width = 6, height = 6, units = 'in', res = 300)

ObsRts_Date <- ggarrange(p4, p5, p6,
                         ncol = 1, nrow = 3)
ObsRts_Date

dev.off()

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

png("Figures/Pre_Rt_time_Vax.png", width = 6, height = 6, units = 'in', res = 300)

PreRts_Date <- ggarrange(p7, p8, p9,
                         ncol = 1, nrow = 3)
PreRts_Date

dev.off()


#### VaxEffect ####

VEMean <- colMeans(model_6000C8_iter[, grep(
         "VaxEffect", colnames(model_6000C8_iter))])
VEQuan <- colQuantiles(model_6000C8_iter[, grep(
        "VaxEffect", colnames(model_6000C8_iter))], probs=c(0.025,0.975))

sum_ve <- round(data.frame(VEMean, VEQuan,
                     row.names = c("Dose 1", "Dose 2", "Dose 3")),
                digits = 4)
colnames(sum_ve) <- c("Mean Effect", "2.5% Q", "97.5% Q")
View(sum_ve)

library(formattable)

formattable(sum_ve)


####Trend####

Ran_Eff <- exp(colMeans(model_6000C8_iter[, grep(
  "random_effects", colnames(model_6000C8_iter))]))

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

png("Figures/SES.png", width = 6, height = 4, units = 'in', res = 300)

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


#### In a subset of LTLAs ####

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



# ARCHIVE PLOTS -----------------------------------------------------------

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
