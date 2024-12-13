###SARS-CoV-2 Rt Model for Vaccination and Variants###

#MRes Biomedical Research - EECID stream - Project 1#


      ## Plots for paper ##



# SET UP ------------------------------------------------------------------

# Packages

library(tidyverse)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(rcartocolor)
library(here)



# DATA --------------------------------------------------------------------


#### 'Observed' data ####

# Data inputed into the model
# Normally saved in session, use Script_model.R if not

Rt_data <- data_stan[[24]]

LTLA <- data_stan[[22]]

VaxProp <- as.data.frame(data_stan[[25]])

VarProp <- data_stan[[26]]

NamesLTLAs <- data_stan[[28]]

date <- max(c(min(data_rt$date), min(data_var$date))) + data_stan[[21]]*7


#### Model data ####

# Let's do a model we are happy with

## Model with Variants, but not Variant's VE

## Script from before DJL3 branch, and
## before trying to add age

## (before messing up)

list_result <- readRDS("C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Analyses/Models_BackUp/Before DJL3 corrections - Summer 2024/Hipercow_5k_10cha_1B_Old_NoAgeModel.Rds")


#### Substract parameters ####

DoVariants <- 1

DoAge <- 0

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
RegionalTrends_low <- list_result[[5]][2]
RegionalTrends_upp <- list_result[[5]][3]

NationalTrend_data <- list_result[[6]]
NationalTrend <- list_result[[6]][1]

Scaling_data <- list_result[[7]]
Scaling <- list_result[[7]][1]

Rt_NoVax <- (as.matrix(VarProp) %*% as.matrix(VarAdvantage))*RegionalTrends
Rt_NoVax_low <- (as.matrix(VarProp) %*% as.matrix(VarAdvantage))*RegionalTrends_low
Rt_NoVax_upp <- (as.matrix(VarProp) %*% as.matrix(VarAdvantage))*RegionalTrends_upp
Rt_NoVax_data <- data.frame(Rt_NoVax, Rt_NoVax_low, Rt_NoVax_upp)
colnames(Rt_NoVax_data) <- c("Rt_NoVax", "Rt_NoVax_low", "Rt_NoVax_upp")

estimates <- data.frame(LTLA, NamesLTLAs, date, Rt_data, Rt_Predictions_data,
                        Rt_NoVax_data, RegionalTrends_data, NationalTrend_data)


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

list_table <- list("VaxEffect" = VaxEffect_data,
                   "VarAdvantage" = VarAdvantage_data,
                   "Scaling" = Scaling_data,
                   "NationalTrends" = NationalTrend_data,
                   "RegionalTrends" = RegionalTrends_data,
                   "Rt_NoVax_[VarAd_x_RegTrend]" = Rt_NoVax,
                   "RtPredictions" = Rt_Predictions_data,
                   "RtData" = Rt_data)



# MINOR -------------------------------------------------------------------


## How much does vaccination reduce Rt overall

mean(estimates$Rt)

mean(estimates$Rt_NoVax)

mean(estimates$Rt)/mean(estimates$Rt_NoVax)



# PLOTS -------------------------------------------------------------------


#### Figure 1 ####

# Observed Rt

fig_1a <- ggplot(data = estimates,
                 mapping = aes(x = date, y = Rt_data, group = date)) +
          geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
          geom_boxplot(color = carto_pal(name = "Safe")[10],
                       fill = carto_pal(name = "Safe")[10], alpha = 0.2) +
          scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
          theme_bw() +
          labs(y = "Observed Reproduction Number (Rt)") +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
                axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.05)),
                legend.text = element_text(size = rel(1.1)))

# Vaccination

fig_1b <- ggplot(data = estimates,
                 mapping = aes(x = date, group = date)) +
          geom_boxplot(data = VaxProp, mapping = aes(y = V1, fill = "V1", color = "V1"), alpha = 0.2) +
          geom_boxplot(data = VaxProp, mapping = aes(y = V2, fill = "V2", color = "V2"), alpha = 0.2) +
          geom_boxplot(data = VaxProp, mapping = aes(y = V3, fill = "V3", color = "V3"), alpha = 0.2) +
          scale_fill_manual(breaks = c("V1", "V2", "V3"),
                            values = c(carto_pal(name = "Safe")[4],
                                       carto_pal(name = "Safe")[8],
                                       carto_pal(name = "Safe")[3]),
                            labels = c("Dose 1", "Dose 2", "Dose 3"),
                            name = "Dose") +
          scale_color_manual(breaks = c("V1", "V2", "V3"),
                             values = c(carto_pal(name = "Safe")[4],
                                        carto_pal(name = "Safe")[8],
                                        carto_pal(name = "Safe")[3]),
                            labels = c("Dose 1", "Dose 2", "Dose 3"),
                            name = "Dose") +
          scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
          scale_y_continuous(labels = scales::percent) +
          theme_bw() +
          labs(y = "Proportion of vaccinated population") +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
                axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.05)),
                legend.text = element_text(size = rel(1.1)), 
                legend.position = "bottom",
                legend.title = element_blank())

# Variants

VarProp_plot <- cbind(date = estimates$date, VarProp) %>%
  pivot_longer(-date, names_to = "variant", values_to = "prop")

fig_1c <- ggplot(data = VarProp_plot,
                 mapping = aes(x = date, y = prop, fill = variant)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(breaks = c("Var_PreAlpha", "Var_Alpha", "Var_Delta"),
                    values = c(carto_pal(name = "Safe")[4],
                               carto_pal(name = "Safe")[8],
                               carto_pal(name = "Safe")[3]),
                    labels = c("Wild-type", "Alpha", "Delta")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(y = "Proportion of circulating variants") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.05)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "bottom",
        legend.title = element_blank())

### All

fig_1 <- ggarrange(fig_1a, fig_1b, fig_1c, ncol = 1, labels = c("A", "B", "C"))

png(filename = "C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Paper/Plots/Figure1.png",
    height = 12, width = 8, res = 1200, units = "in")
fig_1
dev.off()


#### Figure 2 ####

# Regional Trend

fig_2a <- ggplot(data = estimates,
                  mapping = aes(x = date, group = LTLA)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
  geom_line(mapping = aes(y = RegionalTrends, color = "RanEffect"), alpha = 0.2) +
  geom_line(mapping = aes(y = NationalTrend, color = "Lambda"), linewidth = 1.2) +
  scale_color_manual(breaks = c("Lambda", "RanEffect"),
                     values = c(carto_pal(name = "Safe")[5], carto_pal(name = "Safe")[1]),
                     labels = c("National Trend", "National Trend scaled to Regional level")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  scale_y_continuous(limits = c(0.2, 2.5)) +
  theme_bw() +
  labs(y = "Reproduction Number (Rt)") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "bottom",
        legend.title = element_blank())

# Regional Trend with VarAd

# Calculating National Trend with VarAdvantage as:
# Average prop of each variant at each timepoint x VarAd

VarProp <- cbind(LTLA = estimates$LTLA,
                 date = estimates$date,
                 VarProp) %>%
  group_by(date) %>%
  mutate(ave_pre = mean(Var_PreAlpha, na.rm = TRUE)) %>%
  mutate(ave_alp = mean(Var_Alpha, na.rm = TRUE)) %>%
  mutate(ave_del = mean(Var_Delta, na.rm = TRUE)) %>% 
  ungroup()

estimates <- estimates %>%
  mutate(NationalTrend_VarAd = (as.matrix(VarProp[,c("ave_pre", "ave_alp", "ave_del")]) %*% as.matrix(VarAdvantage))*NationalTrend)
colnames(estimates)[17] <- "NationalTrend_VarAd"

fig_2b <- ggplot(data = estimates,
                 mapping = aes(x = date, group = LTLA)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
  geom_line(mapping = aes(y = Rt_NoVax, color = "WithVarAd", linetype = "WithVarAd"), alpha = 0.2) +
  geom_line(mapping = aes(y = NationalTrend_VarAd, color = "Lambda", linetype = "Lambda"), linewidth = 1.2) +
  scale_color_manual(breaks = c("Lambda", "WithVarAd"),
                     values = c(carto_pal(name = "Safe")[5], carto_pal(name = "Safe")[7]),
                     labels = c("National Rt in the absence of vaccination", "LTLA Rt in the absence of vaccination"),
                     name = "Type") +
  scale_linetype_manual(breaks = c("Lambda", "WithVarAd"),
                        values = c("dotdash", "solid"),
                        labels = c("National Rt in the absence of vaccination", "LTLA Rt in the absence of vaccination"),
                        name = "Type") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  scale_y_continuous(limits = c(0.2, 2.5)) +
  theme_bw() +
  labs(y = "Reproduction Number (Rt)") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "bottom",
        legend.title = element_blank())

### All

fig_2 <- ggarrange(fig_2a, fig_2b, ncol = 1, labels = c("A", "B"))

png(filename = "C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Paper/Plots/Figure2.png",
    height = 5, width = 8, res = 1200, units = "in")
fig_2b
dev.off()


#### Figure 3 ####

## Parameter in a couple of LTLAs

estimates_sum <- estimates %>%
  filter(LTLA %in% c(1,47,93,104,122,147,208,179,221))

fig_3 <- ggplot(data = estimates_sum) +
  geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
  geom_line (mapping = aes(x = date, y = Rt_NoVax, color = "WithVarAd"), linewidth = 1.05) +
  geom_ribbon(mapping = aes(x = date, ymin = Rt_NoVax_low, ymax = Rt_NoVax_upp, fill = "WithVarAd"), alpha = 0.2) +
  geom_line (mapping = aes(x = date, y = Rt, color = "RtLogP"), linewidth = 1.05) +
  geom_ribbon(mapping = aes(x = date, ymin = X2.5., ymax = X97.5., fill = "RtLogP"), alpha = 0.2) +
  geom_line (mapping = aes(x = date, y = Rt_data, color = "RtData"), linewidth = 1.05) +
  scale_fill_manual(breaks = c("WithVarAd", "RtLogP", "RtData"),
                    values = c(carto_pal(name = "Safe")[7], carto_pal(name = "Safe")[4],
                               carto_pal(name = "Safe")[10]),
                    labels = c("Rt in the absence of vaccination",
                               "Model-Predicted Rt", "Observed Rt"),
                    name = "Type") +
  scale_color_manual(breaks = c("WithVarAd", "RtLogP", "RtData"),
                     values = c(carto_pal(name = "Safe")[7], carto_pal(name = "Safe")[4],
                              carto_pal(name = "Safe")[10]),
                     labels = c("Rt in the absence of vaccination",
                                "Model-Predicted Rt", "Observed Rt"),
                     name = "Type") +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y") +
  labs(y = "Reproduction Number (Rt)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.05), angle = 10), axis.text.y = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)), strip.text = element_text(size = rel(1.2)),
        legend.position = "bottom", legend.title = element_blank()) +
  facet_wrap(NamesLTLAs ~ .) + guides(fill = "none")

png(filename = "C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Paper/Plots/Figure3.png",
    height = 8, width = 12, res = 1200, units = "in")
fig_3
dev.off()


#### Figure 4 ####

# VE

VaxEffect_plot <- cbind(dose = c("Dose 1", "Dose 2", "Dose 3"),
                        VaxEffect_data[1:3,])

fig_4a <- ggplot(data = VaxEffect_plot,
                 mapping = aes(x = dose, y = VE, color = dose)) +
  geom_point(shape = 15, size = 3) +
  geom_errorbar(mapping = aes(ymin = X2.5., ymax = X97.5.), width = 0.1) +
  scale_color_manual(values = c(carto_pal(name = "Safe")[4],
                                carto_pal(name = "Safe")[8],
                                carto_pal(name = "Safe")[3]),
                     labels = c("Dose 1", "Dose 2", "Dose 3")) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(y = "Vaccine Effect") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

# VarAd

VarAdvantage_plot <- cbind(variant = c(" Pre-Alpha", "Alpha", "Delta"),
                           VarAdvantage_data)

fig_4b <- ggplot(data = VarAdvantage_plot,
                 mapping = aes(x = variant, y = VarAdvantage, color = variant)) +
  geom_point(shape = 15, size = 3) +
  geom_errorbar(mapping = aes(ymin = X2.5., ymax = X97.5.), width = 0.1) +
  scale_x_discrete(breaks = c(" Pre-Alpha", "Alpha", "Delta"),
                   labels = c("Wild-type", "Alpha", "Delta")) +
  scale_color_manual(breaks = c(" Pre-Alpha", "Alpha", "Delta"),
                     values = c(carto_pal(name = "Safe")[4],
                                carto_pal(name = "Safe")[8],
                                carto_pal(name = "Safe")[3]),
                     labels = c("Wild-type", "Alpha", "Delta")) +
  theme_bw() +
  labs(y = "Variants Advantage") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

## All

fig_4 <- ggarrange(fig_4a, fig_4b, ncol = 2, labels = c("A", "B"))

png(filename = "C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Paper/Plots/Figure4.png",
    height = 4, width = 6, res = 1200, units = "in")
fig_4
dev.off()


#### Supplementary ####

# Model vs observed Rt

fig_s1a <- ggplot(data = estimates) +
  geom_point(mapping = aes(x = Rt_data, y = Rt)) +
  geom_abline(x = 0, y = 1, color = carto_pal(name = "Safe")[9]) +
  theme_bw() +
  labs(x = "Observed Reproduction Number (Rt)",
       y = "Model-Predicted Reproduction Number (Rt)") 

fig_s1b <- ggplot(data = estimates,
                 mapping = aes(x = date, group = date)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
  geom_boxplot(mapping = aes(y = Rt_data, color = "RtData", fill = "RtData"), alpha = 0.2) +
  geom_boxplot(mapping = aes(y = Rt, color = "RtLogP", fill = "RtLogP"), alpha = 0.2) +
  scale_color_manual(breaks = c("RtLogP", "RtData"),
                     values = c(carto_pal(name = "Safe")[4], carto_pal(name = "Safe")[10]),
                     labels = c("Model-Predicted Rt", "Observed Rt"), name = "Rt") +
  scale_fill_manual(breaks = c("RtLogP", "RtData"),
                    values = c(carto_pal(name = "Safe")[4], carto_pal(name = "Safe")[10]),
                    labels = c("Model-Predicted Rt", "Observed Rt"), name = "Rt") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  theme_bw() +
  labs(y = "Reproduction Number (Rt)") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank())

fig_s1 <- ggarrange(fig_s1a, fig_s1b, ncol = 1, labels = c("A", "B"))

png(filename = "C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Paper/Plots/FigureS1.png",
    height = 10, width = 8, res = 1200, units = "in")
fig_s1
dev.off()

# Fig 3 for the rest of LTLAs

estimates_rest <- estimates %>%
  filter(!(LTLA %in% c(1,47,93,104,122,147,208,179,221)))

cut_offs <- list(c(1:20), c(21:40), c(41:60), c(61:80), c(81:100),
                 c(101:120), c(121:140), c(141:160), c(161:180), c(181:200),
                 c(201:212))

cut_titles <- c("FigureS2a", "FigureS2b", "FigureS2c", "FigureS2d", "FigureS2e",
                "FigureS2f", "FigureS2g", "FigureS2h", "FigureS2i", "FigureS2j",
                "FigureS2k")

for (cut in 1:length(cut_offs)) {
  
  p <- ggplot(data = filter(estimates_rest, LTLA %in% cut_offs[[cut]])) +
    geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
    geom_line (mapping = aes(x = date, y = Rt_NoVax, color = "WithVarAd"), linewidth = 1.05) +
    geom_ribbon(mapping = aes(x = date, ymin = Rt_NoVax_low, ymax = Rt_NoVax_upp, fill = "WithVarAd"), alpha = 0.2) +
    geom_line (mapping = aes(x = date, y = Rt, color = "RtLogP"), linewidth = 1.05) +
    geom_ribbon(mapping = aes(x = date, ymin = X2.5., ymax = X97.5., fill = "RtLogP"), alpha = 0.2) +
    geom_line (mapping = aes(x = date, y = Rt_data, color = "RtData"), linewidth = 1.05) +
    scale_fill_manual(breaks = c("WithVarAd", "RtLogP", "RtData"),
                      values = c(carto_pal(name = "Safe")[7], carto_pal(name = "Safe")[4],
                                 carto_pal(name = "Safe")[10]),
                      labels = c("Rt in the absence of vaccination",
                                 "Model-Predicted Rt", "Observed Rt"),
                      name = "Type") +
    scale_color_manual(breaks = c("WithVarAd", "RtLogP", "RtData"),
                       values = c(carto_pal(name = "Safe")[7], carto_pal(name = "Safe")[4],
                                  carto_pal(name = "Safe")[10]),
                       labels = c("Rt in the absence of vaccination",
                                  "Model-Predicted Rt", "Observed Rt"),
                       name = "Type") +
    scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y") +
    labs(y = "Reproduction Number (Rt)") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
          axis.text.x = element_text(size = rel(0.85), angle = 10), axis.text.y = element_text(size = rel(1.1)),
          legend.text = element_text(size = rel(1.1)), strip.text = element_text(size = rel(1.2)),
          legend.position = "bottom", legend.title = element_blank()) +
    facet_wrap(NamesLTLAs ~ .) + guides(fill = "none")
  
  png(filename = paste0("C:/Users/nd1316/OneDrive - Imperial College London/MRes/PROJECT 1/Paper/Plots/", cut_titles[cut], ".png"),
      height = 12, width = 16, res = 1200, units = "in")
  print(p)
  dev.off()
  
}
rm(p, cut)

