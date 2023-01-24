###SARS-CoV-2 Rt Model for Vaccination and Variants###

        ###Model update exploratory script###
     ##Proposal: Multiple doses, Two variants##

#MRes Biomedical Research - EECID stream - Project 1#



# SET UP ------------------------------------------------------------------

rm(list = ls())

# Data

data_vax <- read.csv("Data/rtm_incoming_vaccination_20211116-173218-f36a1245_vacc_coverage_ltla.csv")
data_rt <- read.csv("Data/UK_hotspot_Rt_estimates.csv")
data_var <- read.csv("Data/vam_by_ltla.csv")

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


#### Var Data ####

names(data_var)
dim(data_var)

# Data on 314 LTLA and 871 dates

l_code2 <- length(unique(data_var$ltlacode))

l_date3 <- length(unique(data_var$date))

dim(data_var)
l_code2*l_date3


#### Clean ####

rm(l_agemin, l_code, l_code2, l_name, l_date, l_date_2, l_date3)



# SET PARAMETERS ----------------------------------------------------------

# Set date parameters

Date_Start <- as.Date("01/01/2021", format = "%d/%m/%Y")
Date_End <- as.Date("31/12/2021", format = "%d/%m/%Y")

# Covariates

# covar_var <- c("Var_PreAlpha", "Var_Alpha", "Var_Delta")
covar_var <- c("Var_Alpha", "Var_Delta")
covar_vax <- c("First_Prop", "Second_Prop", "Third_Prop")

# Nested and BackDate (For later Stan Model)

BackDate_Char <- "_BD"

Nested_Char <- ""

# Sets of easing lockdown

lockdown_steps <- as.Date(c("05/01/2021", "08/03/2021", "19/04/2021",
                            "17/05/2021", "19/07/2021"), format = "%d/%m/%Y")



# DATA GET FUNCTION --------------------------------------------------


get_data <- function(data_vax, data_rt, data_var,
                     covar_var, covar_vax,
                     Date_Start, Date_End,
                     lockdown_steps,
                     DoVariants, DoVaxVariants, DoAge,
                     DoKnots, Quadratic,
                     IncludeIntercept, IncludeScaling) {
  
  # Substract data
  
  dvax <- data_vax
  drt <- data_rt
  dvar <- data_var
  dage <- data_vax
  
  
  ## Making Vax/Var/Rt Data comparable ##
  
  # No/Yes age groups
  
  dvax <- filter(dvax, is.na(age_band_min))
  
  dage <- filter(dage, !is.na(age_band_min))
  
  # Dates of interest
  
  dvax <- rename(dvax, date = vaccination_date)
  dvax <- filter(dvax, date >= Date_Start)
  dvax <- filter(dvax, date <= Date_End)
  
  drt <- filter(drt, date >= Date_Start)
  drt <- filter(drt, date <= Date_End)
  
  dvar <- filter(dvar, date >= Date_Start)
  dvar <- filter(dvar, date <= Date_End)
  
  dage <- rename(dage, date = vaccination_date)
  dage <- filter(dage, date >= Date_Start)
  dage <- filter(dage, date <= Date_End)
  
  # LTLA names
  
  drt <- rename(drt, ltla_name = area)
  dvar <- rename(dvar, ltla_code = ltlacode)
  
  # Three age groups
  
  dage <- dage %>%
    mutate(total = First/First_Prop) %>%
    mutate(group = case_when(age_band_min >= 15 & age_band_min <= 45 ~ "15-49",
                             age_band_min >= 50 & age_band_min <= 65 ~ "50-69",
                             age_band_min >= 70 & age_band_min <= 90 ~ "70plus"))
  
  #To calculate new group proportions, only one weekly entry per age band
  dage$week <- round(as.numeric(floor((dage$date - Date_Start)/7)), digits = 0)
  dage <- dage %>%
    mutate(combi = paste0(week, ltla_name, age_band_min))
  dage <- filter(dage, !duplicated(dage$combi))
  
  
  ## Merge ##
  
  # Select the vars of interest
  
  dvax <- select(dvax, "ltla_code", "ltla_name", "date",
                 "First_Prop", "Second_Prop", "Third_Prop")
  drt <- select(drt, "ltla_name", "date", "Rt")
  dvar <- select(dvar, "ltla_code", "date", "n_all_wildtype_variant",
                 "n_all_alpha_variant", "n_all_delta_variant")
  dage <- select(dage, "ltla_code", "ltla_name", "date", "group",
                 "First", "Second", "Third", "total")
  
  # Merge
  
  data_merge <- merge(dvax, drt, by = c("ltla_name", "date"), all = TRUE)
  data_merge <- merge(data_merge, dvar, by = c("ltla_code", "date"), all = TRUE)
  
  data_merge_age <- merge(dage, drt, by = c("ltla_name", "date"), all = TRUE)
  
  
  ## Final Cleaning ##
  
  # No NA
  
  data_merge <- data_merge[complete.cases(data_merge),]
  data_merge_age <- data_merge_age[complete.cases(data_merge_age),]
  
  # Steps
  
  # First lockdown step is not in the data: take the initial date
  
  lockdown_steps %in% data_merge$date
  
  Steps <- c(min(data_merge$date), lockdown_steps[2:length(lockdown_steps)], max(data_merge$date))
  
  # Define knots with a day and week time scale
  
  Knots <- round(as.numeric(floor((Steps - Steps[1]))), digits = 0)
  Knots_weeks <- round(as.numeric(floor((Steps - Steps[1])/7)), digits = 0)
  
  # Weekly dates: var for no. of weeks & only one obs per week
  
  #Var for the number of weeks
  data_merge$week <- round(as.numeric(floor((data_merge$date - min(data_merge$date))/7)), digits = 0)
  data_merge_age$week <- round(as.numeric(floor((data_merge_age$date - min(data_merge_age$date))/7)), digits = 0)
  
  #Combi for unique combination of LTLA*week
  data_merge <- data_merge %>%
    mutate(combi = paste0(week, "/", ltla_name))
  data_merge_age <- data_merge_age %>%
    mutate(combi = paste0(week, "/", group, "/", ltla_name))
  
  # Age group proportion: prop of vaccinated and age prop in each region
  
  data_merge_age <- data_merge_age %>%
    group_by(combi) %>%
    mutate(FirstDose = sum(First),
           SecondDose = sum(Second),
           ThirdDose = sum(Third)) %>%
    mutate(First_Prop = FirstDose / sum(total),
           Second_Prop = SecondDose / sum(total),
           Third_Prop = ThirdDose / sum(total)) %>%
    ungroup() %>%
    arrange(ltla_name, week, group)
  
  # Calculate age proportion in each LTLA
  
  age_prop <- data_merge_age %>%
    group_by(ltla_name, week) %>%
    mutate(total_ltla = sum(total)) %>%
    ungroup() %>%
    group_by(ltla_name, week, group) %>%
    mutate(age_prop = sum(total)/total_ltla) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(age_prop)
  
  #Remove the duplicates
  
  data_merge <- filter(data_merge, !duplicated(data_merge$combi))
  data_merge_age <- filter(data_merge_age, !duplicated(data_merge_age$combi))
  
  # Create vars for the proportion of alpha vs delta
  
  data_merge <- data_merge %>%
    # mutate(total = (n_all_wildtype_variant + n_all_alpha_variant + n_all_delta_variant)) %>%
    # mutate(Var_PreAlpha = case_when(total == 0 ~ 0,
    #                              total != 0 ~ n_all_wildtype_variant/total)) %>%
    mutate(total = (n_all_alpha_variant + n_all_delta_variant)) %>%
    mutate(Var_Alpha = case_when(total == 0 ~ 1,
                                 total != 0 ~ n_all_alpha_variant/total)) %>%
    mutate(Var_Delta = case_when(total == 0 ~ 0,
                                 total != 0 ~ n_all_delta_variant/total))
  
  # Select only the ltla with max no of weeks
  
  data_merge <- data_merge %>%
    group_by(ltla_name) %>%
    mutate(max_week = length(unique(week))) %>%
    ungroup() %>%
    filter(max_week == length(unique(week)))
  
  # Select cols
  
  data_model <- select(data_merge, "ltla_name", "date", "week", "Rt",
                       "First_Prop", "Second_Prop", "Third_Prop",
                       "Var_Alpha", "Var_Delta")
                       # "Var_PreAlpha", "Var_Alpha", "Var_Delta")
  data_model_age <- select(data_merge_age, "ltla_name", "date", "week", "Rt","group",
                           "First_Prop", "Second_Prop", "Third_Prop")
  
  
  ## Stan list data ##
  
  if (DoAge == 0) {
    
    #### No age: var or not ####
    
    dim(data_model)
    
    # No of LTLA
    
    NumLTLAs <- length(unique(data_model$ltla_name))
    NumLTLAs
    
    # No of LTLAs in the data
    
    NamesLTLAs <- unique(data_model$ltla_name)
    
    data_model$LTLAs <- NA
    for(i in 1:NumLTLAs){
      data_model$LTLAs[data_model$ltla_name == NamesLTLAs[i]] = i
    }
    
    LTLAs <- data_model$LTLAs
    LTLAs
    
    # No of weeks
    
    NumTimepoints <- length(unique(data_model$week))
    NumTimepoints
    
    Timepoints <- data_model$week
    Timepoints
    
    # No of weeks per LTLA
    
    NumWeeksByLTLA <- rep(NA, NumLTLAs)
    for(i in 1: NumLTLAs){   
      
      d_sub = data_model[data_model$ltla_name == NamesLTLAs[i], ]
      
      NumWeeksByLTLA[i] = length(unique(d_sub$week))
    }
    NumWeeksByLTLA
    
    # No groups
    
    NumGroup <- 1
    NumGroup
    
    Groups <- rep(1, times = nrow(data_model))
    Groups
    
    # No of total obs 
    
    NumDatapoints <- NumLTLAs * NumTimepoints * NumGroup
    NumDatapoints
    
    nrow(data_model)
    
    # No of Knots
    
    Knots
    Knots_weeks
    
    NumKnots <- length(Knots)
    NumKnots
    
    # No of LTLAs x No of Knots
    
    NumPointsLine <- NumLTLAs*NumKnots
    NumPointsLine
    
    
    # Rt log #
    
    data_model$Rt <- log(data_model$Rt)
    
    
    # Spline option
    
    if (DoKnots == 1) {
      NumTrendPar <- NumKnots
    } else {
      NumTrendPar <- NumTimepoints
    }
    
    # Variants
    
    if (DoVariants == 1) {
      NumVar <- length(covar_var)
      VarProp <- data_model[,covar_var]
    } else {
      NumVar <- 1
      VarProp <- as.data.frame(matrix(1, nrow = NumDatapoints))
    }
    
    if (DoVaxVariants == 1) {
      NumVaxVar <- NumVar
    } else {
      NumVaxVar <- 1
    }
    
    # Stan list
    
    data_stan <- list(
      IncludeIntercept = IncludeIntercept,
      IncludeScaling = IncludeScaling,
      DoKnots = DoKnots,
      Quadratic = Quadratic,
      DoVariants = DoVariants,
      DoVaxVariants = DoVaxVariants,
      DoAge = DoAge,
      
      NumDatapoints = NumDatapoints,
      NumLTLAs = NumLTLAs,
      NumDoses = length(covar_vax),
      NumVar = NumVar,
      NumVaxVar = NumVaxVar,
      NumGroup = NumGroup,
      NumTimepoints = NumTimepoints,
      NumKnots = NumKnots,
      NumPointsLine = NumPointsLine,
      NumTrendPar = NumTrendPar,
      
      Knots = Knots_weeks,
      Timepoints = Timepoints,
      LTLAs = LTLAs,
      Groups = Groups,
      
      RtVals = data_model$Rt,
      VaxProp = data_model[,covar_vax],
      VarProp = VarProp,
      AgeProp = as.data.frame(matrix(1, nrow = NumDatapoints))
    )  
    
  } else {
    
    #### Age model ####
    
    dim(data_model_age)
    
    # No of LTLA
    
    NumLTLAs <- length(unique(data_model_age$ltla_name))
    NumLTLAs
    
    # LTLAs in the data
    
    NamesLTLAs <- unique(data_model_age$ltla_name)
    
    data_model_age$LTLAs <- NA
    for(i in 1:NumLTLAs){
      data_model_age$LTLAs[data_model_age$ltla_name == NamesLTLAs[i]] = i
    }
    
    LTLAs <- data_model_age$LTLAs
    LTLAs
    
    # No of weeks
    
    NumTimepoints <- length(unique(data_model_age$week))
    NumTimepoints
    
    Timepoints <- data_model_age$week
    Timepoints
    
    # No of weeks per LTLA
    
    NumWeeksByLTLA <- rep(NA, NumLTLAs)
    for(i in 1: NumLTLAs){   
      
      d_sub = data_model_age[data_model_age$ltla_name == NamesLTLAs[i], ]
      
      NumWeeksByLTLA[i] = length(unique(d_sub$week))
    }
    NumWeeksByLTLA
    
    # No groups
    
    NumGroup <- length(unique(data_model_age$group))
    NumGroup
    
    NamesGroups <- unique(data_model_age$group)
    
    data_model_age$Groups <- NA
    for(i in 1:NumGroup){
      data_model_age$Groups[data_model_age$group == NamesGroups[i]] = i
    }
    
    Groups <- data_model_age$Groups
    Groups
    
    # No of total obs 
    
    NumDatapoints <- NumLTLAs * NumTimepoints * NumGroup
    NumDatapoints
    
    nrow(data_model_age)
    
    # No of Knots
    
    Knots
    Knots_weeks
    
    NumKnots <- length(Knots)
    NumKnots
    
    # No of LTLAs x No of Knots
    
    NumPointsLine <- NumLTLAs*NumKnots
    NumPointsLine
    
    
    # Rt log #
    
    data_model_age$Rt <- log(data_model_age$Rt)
    
    
    # Spline option
    
    if (DoKnots == 1) {
      NumTrendPar <- NumKnots
    } else {
      NumTrendPar <- NumTimepoints
    }
    
    # Variants / Age groups option
    
    NumVar <- 1
    VarProp <- as.data.frame(matrix(1, nrow = NumDatapoints))
    NumVaxVar <- 1
    
    # Stan list
    
    data_stan <- list(
      IncludeIntercept = IncludeIntercept,
      IncludeScaling = IncludeScaling,
      DoKnots = DoKnots,
      Quadratic = Quadratic,
      DoVariants = DoVariants,
      DoVaxVariants = DoVaxVariants,
      DoAge = DoAge,
      
      NumDatapoints = NumDatapoints,
      NumLTLAs = NumLTLAs,
      NumDoses = length(covar_vax),
      NumVar = NumVar,
      NumVaxVar = NumVaxVar,
      NumGroup = NumGroup,
      NumTimepoints = NumTimepoints,
      NumKnots = NumKnots,
      NumPointsLine = NumPointsLine,
      NumTrendPar = NumTrendPar,
      
      Knots = Knots_weeks,
      Timepoints = Timepoints,
      LTLAs = LTLAs,
      Groups = Groups,
      
      RtVals = data_model_age$Rt,
      VaxProp = data_model_age[,covar_vax],
      VarProp = VarProp,
      AgeProp = age_prop
    )
    
  }
    
  data_stan
}

# Run function to get data

data_stan <- get_data(data_vax = data_vax, data_rt = data_rt, data_var = data_var,
                      covar_var = covar_var, covar_vax = covar_vax,
                      Date_Start = Date_Start, Date_End = Date_End,
                      lockdown_steps = lockdown_steps,
                      DoVariants = 1, DoVaxVariants = 1, DoAge = 0,
                      DoKnots = 0, Quadratic = 0,
                      IncludeIntercept = 0, IncludeScaling = 0)

data_stan_age <- get_data(data_vax = data_vax, data_rt = data_rt, data_var = data_var,
                      covar_var = covar_var, covar_vax = covar_vax,
                      Date_Start = Date_Start, Date_End = Date_End,
                      lockdown_steps = lockdown_steps,
                      DoVariants = 0, DoAge = 1,
                      DoKnots = 0, Quadratic = 0,
                      IncludeIntercept = 1, IncludeScaling = 1)



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
