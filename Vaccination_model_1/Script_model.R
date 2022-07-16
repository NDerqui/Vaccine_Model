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

# Sets of easing lockdown

lockdown_steps <- as.Date(c("05/01/2021", "08/03/2021", "19/04/2021",
                            "17/05/2021", "19/07/2021"), format = "%d/%m/%Y")



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


#### Merge ####

data_merge <- merge(dvax, drt, by = c("ltla_name", "date"), all = TRUE)

# Checks

dim(data_merge)
names(data_merge)


#### Final Cleaning ####

# No NA

data_merge <- data_merge[complete.cases(data_merge),]

# Steps

# First lockdown step is not in the data: take the initial date

lockdown_steps %in% data_merge$date
min(data_merge$date)

Steps <- c(min(data_merge$date), lockdown_steps[2:length(lockdown_steps)], max(data_merge$date))
Steps %in% data_merge$date

# Define knots with a day and week time scale

Knots <- round(as.numeric(floor((Steps - Steps[1]))), digits = 0)
Knots_weeks <- round(as.numeric(floor((Steps - Steps[1])/7)), digits = 0)
#Match the weeks of the Knots to the ones in the data (start on week 4)

# Weekly dates: var for no. of weeks & only one obs per week

#Var for the number of weeks
data_merge$week <- round(as.numeric(floor((data_merge$date - min(data_merge$date))/7)), digits = 0)

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

Timepoints <- data_model$week
Timepoints

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

# No of Knots

Knots
Knots_weeks

NumKnots <- length(Knots)
NumKnots

# No of LTLAs x No of Knots

NumPointsLine <- NumLTLAs*NumKnots
NumPointsLine


#### Rt log ####

data_model$Rt <- log(data_model$Rt)


#### Switches for parameters ####

IncludeIntercept <- 1
IncludeScaling <- 1

# Analyses by steps of lockdown

DoKnots <- 1
Quadratic <- 1

if (DoKnots == 1) {
  NumTrendPar <- NumKnots
} else {
  NumTrendPar <- NumTimepoints
}


#### Stan Data ####

data_model <- select(data_model,
                     "ltla_name", "LTLAs", "date", "week", covar_vax, "Rt")
# saveRDS(data_model, "data_model_for_plots.Rds")

data_stan <- list(
          IncludeIntercept = IncludeIntercept,
          IncludeScaling = IncludeScaling,
          DoKnots = DoKnots,
          Quadratic = Quadratic,
          NumDatapoints = NumDatapoints,
          NumLTLAs = NumLTLAs,
          NumDoses = length(covar_vax),
          NumTimepoints = NumTimepoints,
          NumKnots = NumKnots,
          NumPointsLine = NumPointsLine,
          NumTrendPar = NumTrendPar,
          Knots = Knots_weeks,
          Timepoints = Timepoints,
          LTLAs = data_model$LTLAs,
          RtVals = data_model$Rt,
          VaxProp = data_model[,covar_vax]
)



# STAN MODEL --------------------------------------------------------------


#### Model ####

library(here)

#ModelChar <- "2stage_LFM_non_centered_Vax"
ModelChar <- "2stage_LFM_non_centered_Vax_consolidated"
StanModel <- stan_model(here(paste0("Stan_models/",ModelChar, ".stan")))
#StanModel <- stan_model(here(paste0("Vaccination_model_1/Stan_models/", ModelChar, ".stan"))) ## added this as my here() package playing up and loads to root directory. Switch back to previous line as required.

cat(paste0("Model compilation done\n"))

# Create and write meta data

ModelMetaData 				= c()
ModelMetaData$iter 			= 200 #Increase
ModelMetaData$warmup 		= 50 #Increase
ModelMetaData$thin 			= 1
ModelMetaData$chains 		= 1 #Increase
ModelMetaData$adapt_delta 	= 0.9
ModelMetaData$max_treedepth = 15
ModelMetaData$ModelChar 	= ModelChar
ModelMetaData$BackDate_Char = BackDate_Char
ModelMetaData$Nested_Char 	= Nested_Char

ModelMetaData_dummy = as.data.frame(unlist(ModelMetaData))
colnames(ModelMetaData_dummy) = NULL


#### Run ####

fit = sampling(StanModel, data = data_stan, 
                 iter 	= ModelMetaData$iter, 
                 warmup 	= ModelMetaData$warmup, 
                 thin 	= ModelMetaData$thin, 
                 chains 	= ModelMetaData$chains, 
                 pars 	= c("VaxEffect", "VacEffects_Regional",
                           "LogPredictions", "RegionalTrends",
                           "NationalTrend", "gamma", "intercept", "lambda",
                           "log_lik"), 
                 control = list(adapt_delta = ModelMetaData$adapt_delta,
                                max_treedepth = ModelMetaData$max_treedepth))


# MODEL RESULTS -----------------------------------------------------------


#### Model name ####

model_note <- "_linear" # Ie. "_nointercept"

model_name <- paste0("fit_", ModelMetaData$iter, "_", ModelMetaData$chains, model_note)


#### Save data ####

# saveRDS(fit, paste0(model_name, ".Rds"))
# saveRDS(fit, paste0("C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Analyses/Models_BackUp/", model_name, ".Rds"))


#### Loo_cv ####

loo_run = loo(fit)
loo_run

# saveRDS(loo_run, paste0("loo_", model_name, ".Rds")
# saveRDS(loo_run, paste0("C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Analyses/Models_BackUp/", "loo_", model_name,".Rds"))
