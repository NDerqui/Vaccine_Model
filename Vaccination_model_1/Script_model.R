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
library(StanHeaders)
library(ggplot2)
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

# No NA

data_merge <- data_merge[complete.cases(data_merge),]


# Weekly dates: var for no. of weeks & only one obs per week

#Var for the number of weeks
data_merge$week <- round(as.numeric(floor((data_merge$date - Date_Start)/7)), digits = 0)

#Combi for unique combination of LTLA*week
data_merge <- mutate(data_merge, combi = paste0(data_merge$week, data_merge$ltla_name))

#Remove the duplicates
data_merge <- filter(data_merge, !duplicated(data_merge$combi))

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


#### Switches for parameters ####

IncludeIntercept <- 1
IncludeScaling <- 1


#### Stan Data ####

data_model <- select(data_model,
                     "ltla_name", "LTLAs", "date", "week", covar_vax, "Rt")
saveRDS(data_model, "data_model_for_plots.Rds")

data_stan <- list(
          RtVals = data_model$Rt,
          VaxProp = data_model[,covar_vax],
          NumLTLAs = NumLTLAs,
          NumDoses = length(covar_vax),
          NumDatapoints = NumDatapoints,
          LTLAs = data_model$LTLAs,
          NumTimepoints = NumTimepoints,
          IncludeIntercept = IncludeIntercept,
          IncludeScaling = IncludeScaling
   )



# DESCRIPTION -------------------------------------------------------------


## A more detailed description  of data is on the description script ##

Rt_data <- exp(data_model$Rt)

sum_data <- data.frame(Rt_data, 
                     Dose_1 = data_model$First_Prop,
                     Dose_2 = data_model$Second_Prop,
                     Dose_3 = data_model$Third_Prop,
                     date = data_model$date,
                     row.names = paste0("Rt", 1:12423))


#### Vax Description ####

# Vaccination over time

Vax_Date <- ggplot(data = sum_data) +
  geom_boxplot (mapping = aes(x = date, y = Dose_1, group = date,
                              color = "Dose_1"), size = rel(0.5)) +
  geom_boxplot (mapping = aes(x = date, y = Dose_2, group = date,
                              color = "Dose_2"), size = rel(0.5)) +
  geom_boxplot (mapping = aes(x = date, y = Dose_3, group = date,
                              color = "Dose_3"), size = rel(0.5)) +
  scale_color_manual(name="Number of doses",
                     breaks = c("Dose_1", "Dose_2", "Dose_3"),
                     values = c("Dose_1"="lightskyblue", 
                                "Dose_2"="palegreen",
                                "Dose_3"="salmon"),
                     labels=c("1 dose", "2 doses", "3 doses")) +
  theme_classic() +
  labs(title = "Distribution of vaccinated individuals proportion across all LTLAs over time",
       x = "Date",
       y = "Total proportion of vaccinated population") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Vax_Date


#### Rt over time ####

# General

Rt_Obs_Date <- ggplot(data = sum_data) +
  geom_boxplot (mapping = aes(x = date, y = Rt_data, group = date), size = rel(0.5)) +
  theme_classic() +
  labs(title = "Observed Rt over time",
       x = "Date",
       y = "Rt Observed") +
  theme(
    plot.title = element_text(size = rel(1), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(0.9), face="bold"),
    axis.title.y = element_text(size = rel(0.9), face="bold"),
    axis.text = element_text(size=rel(0.7)),
    legend.title = element_text(size = rel(0.9), face="bold"),
    legend.text = element_text(size=rel(0.7)))
Rt_Obs_Date



# STAN MODEL --------------------------------------------------------------


#### Model ####

library(here)

#ModelChar <- "2stage_LFM_non_centered_Vax"
ModelChar <- "2stage_LFM_non_centered_Vax_log_lik"
StanModel <- stan_model(here(paste0("Stan_models/",ModelChar, ".stan")))
cat(paste0("Model compilation done\n"))

# Create and write meta data

ModelMetaData 				= c()
ModelMetaData$iter 			= 10000 #Increase
ModelMetaData$warmup 		= 2000 #Increase
ModelMetaData$thin 			= 1
ModelMetaData$chains 		= 10 #Increase
ModelMetaData$adapt_delta 	= 0.9
ModelMetaData$max_treedepth = 15
ModelMetaData$ModelChar 	= ModelChar
ModelMetaData$BackDate_Char = BackDate_Char
ModelMetaData$Nested_Char 	= Nested_Char

ModelMetaData_dummy = as.data.frame(unlist(ModelMetaData))
colnames(ModelMetaData_dummy) = NULL


#### Run ####

memory.limit(size = 100000000)

fit = sampling(StanModel, data = data_stan, 
               iter 	= ModelMetaData$iter, 
               warmup 	= ModelMetaData$warmup, 
               thin 	= ModelMetaData$thin, 
               chains 	= ModelMetaData$chains, 
               pars 	= c("VaxEffect", "LogPredictions", "random_effects",
                          "lambda", "gamma", "intercept", "log_lik"), 
               control = list(adapt_delta = ModelMetaData$adapt_delta,
                              max_treedepth = ModelMetaData$max_treedepth))



# MODEL RESULTS -----------------------------------------------------------

# Comment with 6000iter-2chains
# 
# Warning messages:
#   1: The largest R-hat is NA, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#r-hat 
# 2: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#bulk-ess 
# 3: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#tail-ess 


#### Model name ####

model_note <- "" # Ie. "_nointercept"

model_name <- paste0("fit_", ModelMetaData$iter, "_", ModelMetaData$chains, model_note)


#### Save data ####

# saveRDS(fit, paste0(model_name, ".Rds"))
# saveRDS(fit, paste0("C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Analyses/Models_BackUp/", model_name, ".Rds"))


#### Loo_cv ####

loo_run = loo(fit)
loo_run

# saveRDS(loo_run, paste0("loo_", model_name, ".Rds")
# saveRDS(loo_run, paste0("C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Analyses/Models_BackUp/", "loo_", model_name,".Rds")
