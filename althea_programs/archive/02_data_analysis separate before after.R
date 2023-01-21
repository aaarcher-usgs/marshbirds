#' # Data processing
#' 
#' ## Project: Marshbird response to herbicide treatments
#' 
#' ### Original programmer: Nina Hill
#' 
#' ### Secondary programmer: Althea Archer
#' 
#' ## Header
#' 
library(ezknitr)


# Analysis libraries
library(lme4)

#remotes::install_github("palday/coefplot2", subdir = "pkg")
library(coefplot2)

library(interactions)
library(glmmTMB)
library(modEvA)
library(ggplot2)
library(tidyr)
library(dplyr)
library(emmeans)


remove(list = ls())
set.seed(22059747)

#' ## Load data
#' 
#' Point counts in both wide and long format
load(file = "data/processed/NWPOINTR2.Rdata")

#' Separate data into before treatment and after treatment datasets and
#' add in new terms (binary, detection/non-detection)
#' 
NWprePOINTR2 <- NWPOINTR2 %>%
  filter(Year == "2015") %>%
  mutate(PBGRbin = case_when(PBGR == 0 ~ 0, 
                             PBGR != 0 ~1),
         AMBIbin = case_when(AMBI == 0 ~ 0, 
                             AMBI != 0 ~1),
         LEBIbin = case_when(LEBI == 0 ~ 0, 
                             LEBI != 0 ~1),
         SORAbin = case_when(SORA == 0 ~ 0, 
                             SORA != 0 ~1),
         VIRAbin = case_when(VIRA == 0 ~ 0, 
                             VIRA != 0 ~1))

NWpostPOINTR2 <- NWPOINTR2 %>%
  filter(Year != "2015") %>%
  mutate(PBGRbin = case_when(PBGR == 0 ~ 0, 
                             PBGR != 0 ~1),
         AMBIbin = case_when(AMBI == 0 ~ 0, 
                             AMBI != 0 ~1),
         LEBIbin = case_when(LEBI == 0 ~ 0, 
                             LEBI != 0 ~1),
         SORAbin = case_when(SORA == 0 ~ 0, 
                             SORA != 0 ~1),
         VIRAbin = case_when(VIRA == 0 ~ 0, 
                             VIRA != 0 ~1))

NWprePOINTR2.long <- NWPOINTR2.long %>%
  filter(Year == "2015") %>%
  mutate(CountBin = case_when(Count == 0 ~ 0, 
                              Count != 0 ~1))

NWpostPOINTR2.long <- NWPOINTR2.long %>%
  filter(Year != "2015") %>%
  mutate(CountBin = case_when(Count == 0 ~ 0, 
                             Count != 0 ~1))


#' Summarize counts for table in manuscript
#' 
NWPOINTR2.long %>%
  group_by(Species, Year, Treatment) %>%
  summarise(sum = sum(Count)) %>%
  print.data.frame()

NWPOINTR2.long %>%
  group_by(Species) %>%
  summarise(sum = sum(Count)) %>%
  print.data.frame()


NWPOINTR2.long %>%
  group_by(Year, Treatment) %>%
  summarise(sum = sum(Count)) %>%
  print.data.frame()

NWPOINTR2.long %>%
  summarise(sum = sum(Count)) %>%
  print.data.frame()

#' Plot histograms of the counts
#' 
table(NWPOINTR2.long$Species, NWPOINTR2.long$Count)
ggplot(data = NWPOINTR2.long, aes(x = Count))+
  geom_histogram() +
  facet_wrap(~ Species)

#' Set up a data.frame to hold the results of analysis for plotting, printing
#' 
results <- data.frame(
  Year = rep(c(rep(2015:2016, each = 2), rep(2017:2018, 2)), 6),
  Species = rep(c("AMBI", "LEBI", "PBGR", "SORA", "VIRA", "All"), each = 8),
  Treatment = rep(c(rep(c("Control", "Herbicide"), 2), rep(c("Control", "Herbicide"), each = 2)),6),
  effectSize = NA,
  effectSizeSE = NA,
  effectSizeP = NA,
  logEstimate = NA,
  logEstLL = NA,
  logEstUL = NA
)


#' ____________________________________________________________________________
#' ### AMBI Analysis
#' 
#' #### Pre-treatment analysis
#' 
#' 
table(NWprePOINTR2$AMBI)
var(NWprePOINTR2$AMBI)/mean(NWprePOINTR2$AMBI)
AMBI_pre_glmerBin <- glmmTMB(AMBI ~ Treatment + (1|SiteName),
                             data = NWprePOINTR2,
                             family = "poisson")
(AMBIpre <- summary(AMBI_pre_glmerBin))
(AMBIpreEmm <- emmeans(AMBI_pre_glmerBin, ~ Treatment, level = 0.9, type = "response"))
#emmip(AMBI_pre_glmerBin, ~ Treatment, type = "response")

#' Results for printing
#' 
results$effectSize[results$Year == 2015 & 
                     results$Species == "AMBI" ] <- AMBIpre$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year == 2015 & 
                     results$Species == "AMBI"] <- AMBIpre$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year == 2015 & 
                       results$Species == "AMBI"] <- AMBIpre$coefficients$cond[,"Pr(>|z|)"]

results$logEstimate[results$Year == 2015 & 
                      results$Species == "AMBI" ] <- summary(AMBIpreEmm)[,"rate"]

results$logEstLL[results$Year == 2015 & 
                      results$Species == "AMBI" ] <- summary(AMBIpreEmm)[,"lower.CL"]

results$logEstUL[results$Year == 2015 & 
                   results$Species == "AMBI" ] <- summary(AMBIpreEmm)[,"upper.CL"]


#' #### Post-treatment analysis
#' 
#' First, test interaction between treatment and year
#' 
table(NWpostPOINTR2$AMBI)
var(NWpostPOINTR2$AMBI)/mean(NWpostPOINTR2$AMBI)
AMBI_post_glmerInt <- glmmTMB(AMBI ~ Treatment + Year + Treatment*Year + (1|SiteName),
                                 data = NWpostPOINTR2,
                                 family = "poisson")
(AMBI_post <- summary(AMBI_post_glmerInt))
(AMBI_postEmm <- emmeans(AMBI_post_glmerInt, Treatment*Year, level = 0.9, type = "response"))
#cat_plot(AMBI_post_glmerInt, pred = Year, modx = Treatment, interval = T, int.width = 0.9)
#emmip(AMBI_post_glmerInt, Treatment ~ Year, type = "response")

#' Results for plotting
#' 
results$effectSize[results$Year != 2015 & 
                     results$Species == "AMBI" ] <- AMBI_post$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year != 2015 & 
                       results$Species == "AMBI"] <- AMBI_post$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year != 2015 & 
                      results$Species == "AMBI"] <- AMBI_post$coefficients$cond[,"Pr(>|z|)"]

AMBI_postEmmSort <- summary(AMBI_postEmm) %>% 
  slice(1,2,3,5,4,6) 

results$logEstimate[results$Year != 2015 & 
                      results$Species == "AMBI" ] <- AMBI_postEmmSort[,"rate"]


results$logEstLL[results$Year != 2015 & 
                   results$Species == "AMBI" ] <- AMBI_postEmmSort[,"lower.CL"]

results$logEstUL[results$Year != 2015 & 
                   results$Species == "AMBI" ] <- AMBI_postEmmSort[,"upper.CL"]

#' ____________________________________________________________________________
#' ### LEBI Analysis
#' 
#' #### Pre-treatment analysis
#' 
#' 
table(NWprePOINTR2$LEBI)
var(NWprePOINTR2$LEBI)/mean(NWprePOINTR2$LEBI)
LEBI_pre_glmerBin <- glmmTMB(LEBI ~ Treatment + (1|SiteName),
                             data = NWprePOINTR2,
                             family = "poisson")
(LEBIpre <- summary(LEBI_pre_glmerBin))
(LEBIpreEmm <- emmeans(LEBI_pre_glmerBin, ~ Treatment, level = 0.9, type = "response"))
#emmip(LEBI_pre_glmerBin, ~ Treatment, type = "response")

#' Results for printing
#' 
results$effectSize[results$Year == 2015 & 
                     results$Species == "LEBI" ] <- LEBIpre$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year == 2015 & 
                       results$Species == "LEBI"] <- LEBIpre$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year == 2015 & 
                      results$Species == "LEBI"] <- LEBIpre$coefficients$cond[,"Pr(>|z|)"]

results$logEstimate[results$Year == 2015 & 
                      results$Species == "LEBI" ] <- summary(LEBIpreEmm)[,"rate"]

results$logEstLL[results$Year == 2015 & 
                   results$Species == "LEBI" ] <- summary(LEBIpreEmm)[,"lower.CL"]

results$logEstUL[results$Year == 2015 & 
                   results$Species == "LEBI" ] <- summary(LEBIpreEmm)[,"upper.CL"]

#' #### Post-treatment analysis
#' 
#' First, test interaction between treatment and year
#' 
table(NWpostPOINTR2$LEBI)
LEBI_post_glmerBinInt <- glmmTMB(LEBI ~ Treatment + Year + Treatment*Year + (1|SiteName),
                                 data = NWpostPOINTR2,
                                 family = "poisson")
(LEBI_post <- summary(LEBI_post_glmerBinInt))
(LEBI_postEmm <- emmeans(LEBI_post_glmerBinInt, ~ Treatment:Year, level = 0.9, type = "response"))
# cat_plot(LEBI_post_glmerBinInt, pred = Year, modx = Treatment, interval = T, int.width = 0.9)
# emmip(LEBI_post_glmerBinInt, Treatment ~ Year, type = "response") +
#   ylim(c(0,1))

#' Results for plotting
#' 
results$effectSize[results$Year != 2015 & 
                     results$Species == "LEBI" ] <- LEBI_post$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year != 2015 & 
                       results$Species == "LEBI"] <- LEBI_post$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year != 2015 & 
                      results$Species == "LEBI"] <- LEBI_post$coefficients$cond[,"Pr(>|z|)"]

LEBI_postEmmSort <- summary(LEBI_postEmm) %>% 
  slice(1,2,3,5,4,6) 

results$logEstimate[results$Year != 2015 & 
                      results$Species == "LEBI" ] <- LEBI_postEmmSort[,"rate"]


results$logEstLL[results$Year != 2015 & 
                   results$Species == "LEBI" ] <- LEBI_postEmmSort[,"lower.CL"]

results$logEstUL[results$Year != 2015 & 
                   results$Species == "LEBI" ] <- LEBI_postEmmSort[,"upper.CL"]



#' ____________________________________________________________________________
#' ### PBGR Analysis
#' 
#' #### Pre-treatment analysis
#' 
table(NWprePOINTR2$PBGR)
var(NWprePOINTR2$PBGR)/mean(NWprePOINTR2$PBGR)
PBGR_pre_glmer <- glmmTMB(PBGR ~ Treatment + (1|SiteName),
                             data = NWprePOINTR2,
                             family = "poisson")
(PBGRpre <- summary(PBGR_pre_glmer))
(PBGRpreEmm <- emmeans(PBGR_pre_glmer, ~ Treatment, level = 0.9, type = "response"))
#emmip(PBGR_pre_glmer, ~ Treatment, type = "response")

#' Results for printing
#' 
results$effectSize[results$Year == 2015 & 
                     results$Species == "PBGR" ] <- PBGRpre$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year == 2015 & 
                       results$Species == "PBGR"] <- PBGRpre$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year == 2015 & 
                      results$Species == "PBGR"] <- PBGRpre$coefficients$cond[,"Pr(>|z|)"]

results$logEstimate[results$Year == 2015 & 
                      results$Species == "PBGR" ] <- summary(PBGRpreEmm)[,"rate"]

results$logEstLL[results$Year == 2015 & 
                   results$Species == "PBGR" ] <- summary(PBGRpreEmm)[,"lower.CL"]

results$logEstUL[results$Year == 2015 & 
                   results$Species == "PBGR" ] <- summary(PBGRpreEmm)[,"upper.CL"]


#' #### Post-treatment analysis
#' 
#' First, test interaction between treatment and year
#' 
table(NWpostPOINTR2$PBGR)
var(NWpostPOINTR2$PBGR)/mean(NWpostPOINTR2$PBGR)
PBGR_post_glmerInt <- glmmTMB(PBGR ~ Treatment + Year + Treatment*Year + (1|SiteName),
                             data = NWpostPOINTR2,
                             family = "poisson")
(PBGR_post <- summary(PBGR_post_glmerInt))
(PBGR_postEmm <- emmeans(PBGR_post_glmerInt, ~ Treatment:Year, level = 0.9, type = "response"))
#cat_plot(PBGR_post_glmerInt, pred = Year, modx = Treatment, interval = T, int.width = 0.9)


#' Results for plotting
#' 
results$effectSize[results$Year != 2015 & 
                     results$Species == "PBGR" ] <- PBGR_post$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year != 2015 & 
                       results$Species == "PBGR"] <- PBGR_post$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year != 2015 & 
                      results$Species == "PBGR"] <- PBGR_post$coefficients$cond[,"Pr(>|z|)"]

PBGR_postEmmSort <- summary(PBGR_postEmm) %>% 
  slice(1,2,3,5,4,6) 

results$logEstimate[results$Year != 2015 & 
                      results$Species == "PBGR" ] <- PBGR_postEmmSort[,"rate"]


results$logEstLL[results$Year != 2015 & 
                   results$Species == "PBGR" ] <- PBGR_postEmmSort[,"lower.CL"]

results$logEstUL[results$Year != 2015 & 
                   results$Species == "PBGR" ] <- PBGR_postEmmSort[,"upper.CL"]



#' ____________________________________________________________________________
#' ### SORA Analysis
#' 
#' #### Pre-treatment analysis
#' 
#' 
table(NWprePOINTR2$SORA)
var(NWprePOINTR2$SORA)/mean(NWprePOINTR2$SORA)
SORA_pre_glmerBin <- glmmTMB(SORA ~ Treatment + (1|SiteName),
                             data = NWprePOINTR2,
                             family = "poisson")
(SORApre <- summary(SORA_pre_glmerBin))
(SORApreEmm <- emmeans(SORA_pre_glmerBin, ~ Treatment, level = 0.9, type = "response"))
#emmip(SORA_pre_glmerBin, ~ Treatment, type = "response")

#' Results for printing
#' 
results$effectSize[results$Year == 2015 & 
                     results$Species == "SORA" ] <- SORApre$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year == 2015 & 
                       results$Species == "SORA"] <- SORApre$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year == 2015 & 
                      results$Species == "SORA"] <- SORApre$coefficients$cond[,"Pr(>|z|)"]

results$logEstimate[results$Year == 2015 & 
                      results$Species == "SORA" ] <- summary(SORApreEmm)[,"rate"]

results$logEstLL[results$Year == 2015 & 
                   results$Species == "SORA" ] <- summary(SORApreEmm)[,"lower.CL"]

results$logEstUL[results$Year == 2015 & 
                   results$Species == "SORA" ] <- summary(SORApreEmm)[,"upper.CL"]

#' #### Post-treatment analysis
#' 
#' First, test interaction between treatment and year
#' 
table(NWpostPOINTR2$SORA)
SORA_post_glmerBinInt <- glmmTMB(SORA ~ Treatment + Year + Treatment*Year + (1|SiteName),
                                 data = NWpostPOINTR2,
                                 family = "poisson")
(SORA_post <- summary(SORA_post_glmerBinInt))
(SORA_postEmm <- emmeans(SORA_post_glmerBinInt, ~ Treatment:Year, level = 0.9, type = "response"))
#cat_plot(SORA_post_glmerBinInt, pred = Year, modx = Treatment, interval = T, int.width = 0.9)
#emmip(SORA_post_glmerBinInt, Treatment ~ Year, type = "response")


#' Results for plotting
#' 
results$effectSize[results$Year != 2015 & 
                     results$Species == "SORA" ] <- SORA_post$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year != 2015 & 
                       results$Species == "SORA"] <- SORA_post$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year != 2015 & 
                      results$Species == "SORA"] <- SORA_post$coefficients$cond[,"Pr(>|z|)"]

SORA_postEmmSort <- summary(SORA_postEmm) %>% 
  slice(1,2,3,5,4,6) 

results$logEstimate[results$Year != 2015 & 
                      results$Species == "SORA" ] <- SORA_postEmmSort[,"rate"]


results$logEstLL[results$Year != 2015 & 
                   results$Species == "SORA" ] <- SORA_postEmmSort[,"lower.CL"]

results$logEstUL[results$Year != 2015 & 
                   results$Species == "SORA" ] <- SORA_postEmmSort[,"upper.CL"]

#' ____________________________________________________________________________
#' ### VIRA Analysis
#' 
#' #### Pre-treatment analysis
#' 
table(NWprePOINTR2$VIRA)
VIRA_pre_glmer <- glmmTMB(VIRA ~ Treatment + (1|SiteName),
                          data = NWprePOINTR2,
                          family = "poisson")
(VIRApre <- summary(VIRA_pre_glmer))
(VIRApreEmm <- emmeans(VIRA_pre_glmer, ~ Treatment, level = 0.9, type = "response"))


#' Results for printing
#' 
results$effectSize[results$Year == 2015 & 
                     results$Species == "VIRA" ] <- VIRApre$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year == 2015 & 
                       results$Species == "VIRA"] <- VIRApre$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year == 2015 & 
                      results$Species == "VIRA"] <- VIRApre$coefficients$cond[,"Pr(>|z|)"]

results$logEstimate[results$Year == 2015 & 
                      results$Species == "VIRA" ] <- summary(VIRApreEmm)[,"rate"]

results$logEstLL[results$Year == 2015 & 
                   results$Species == "VIRA" ] <- summary(VIRApreEmm)[,"lower.CL"]

results$logEstUL[results$Year == 2015 & 
                   results$Species == "VIRA" ] <- summary(VIRApreEmm)[,"upper.CL"]


#' #### Post-treatment analysis
#' 
#' First, test interaction between treatment and year
#' 
table(NWpostPOINTR2$VIRA)
VIRA_post_glmerInt <- glmmTMB(VIRA ~ Treatment + Year + Treatment*Year + (1|SiteName),
                              data = NWpostPOINTR2,
                              family = "poisson")
(VIRA_post <- summary(VIRA_post_glmerInt))
(VIRA_postEmm <- emmeans(VIRA_post_glmerInt, ~ Treatment:Year, level = 0.9, type = "response"))
#cat_plot(VIRA_post_glmerInt, pred = Year, modx = Treatment, interval = T, int.width = 0.9)
#emmip(VIRA_post_glmerInt, Treatment ~ Year, type = "response")



#' Results for plotting
#' 
results$effectSize[results$Year != 2015 & 
                     results$Species == "VIRA" ] <- VIRA_post$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year != 2015 & 
                       results$Species == "VIRA"] <- VIRA_post$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year != 2015 & 
                      results$Species == "VIRA"] <- VIRA_post$coefficients$cond[,"Pr(>|z|)"]

VIRA_postEmmSort <- summary(VIRA_postEmm) %>% 
  slice(1,2,3,5,4,6) 

results$logEstimate[results$Year != 2015 & 
                      results$Species == "VIRA" ] <- VIRA_postEmmSort[,"rate"]
  

results$logEstLL[results$Year != 2015 & 
                   results$Species == "VIRA" ] <- VIRA_postEmmSort[,"lower.CL"]

results$logEstUL[results$Year != 2015 & 
                   results$Species == "VIRA" ] <- VIRA_postEmmSort[,"upper.CL"]

#' ____________________________________________________________________________
#' ### All species together analysis
#' 
#' #### Pre-treatment analysis
#' 
table(NWprePOINTR2.long$Count)
ALL_pre_glmer <- glmmTMB(Count ~ Treatment + (1|SiteName),
                          data = NWpostPOINTR2.long,
                          family = "poisson")
(ALLpre <- summary(ALL_pre_glmer))
(ALLpreEmm <- emmeans(ALL_pre_glmer, ~ Treatment, level = 0.9, type = "response"))


#' Results for printing
#' 
results$effectSize[results$Year == 2015 & 
                     results$Species == "All" ] <- ALLpre$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year == 2015 & 
                       results$Species == "All"] <- ALLpre$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year == 2015 & 
                      results$Species == "All"] <- ALLpre$coefficients$cond[,"Pr(>|z|)"]

results$logEstimate[results$Year == 2015 & 
                      results$Species == "All" ] <- summary(ALLpreEmm)[,"rate"]

results$logEstLL[results$Year == 2015 & 
                   results$Species == "All" ] <- summary(ALLpreEmm)[,"lower.CL"]

results$logEstUL[results$Year == 2015 & 
                   results$Species == "All" ] <- summary(ALLpreEmm)[,"upper.CL"]

#' #### Post-treatment analysis
#' 
#' First, test interaction between treatment and year
#' 
table(NWpostPOINTR2.long$Count)
ALL_post_glmerInt <- glmmTMB(Count ~ Treatment + Year + Treatment*Year + (1|SiteName) + (1|Species),
                              data = NWpostPOINTR2.long,
                              family = "poisson")
(ALL_post <- summary(ALL_post_glmerInt))
(ALL_postEmm <- emmeans(ALL_post_glmerInt,pairwise  ~ Treatment:Year, level = 0.9, type = "response"))
#cat_plot(ALL_post_glmerInt, pred = Year, modx = Treatment, interval = T, int.width = 0.9)
#emmip(ALL_post_glmerInt, Treatment ~ Year, type = "response") +
#  ylim(c(0,1))



#' Results for plotting
#' 
results$effectSize[results$Year != 2015 & 
                     results$Species == "All" ] <- ALL_post$coefficients$cond[,"Estimate"]

results$effectSizeSE[results$Year != 2015 & 
                       results$Species == "All"] <- ALL_post$coefficients$cond[,"Std. Error"]

results$effectSizeP[results$Year != 2015 & 
                      results$Species == "All"] <- ALL_post$coefficients$cond[,"Pr(>|z|)"]

ALL_postEmmSort <- summary(ALL_postEmm) %>% 
  slice(1,2,3,5,4,6) 

results$logEstimate[results$Year != 2015 & 
                      results$Species == "All" ] <- ALL_postEmmSort[,"rate"]


results$logEstLL[results$Year != 2015 & 
                   results$Species == "All" ] <- ALL_postEmmSort[,"lower.CL"]

results$logEstUL[results$Year != 2015 & 
                   results$Species == "All" ] <- ALL_postEmmSort[,"upper.CL"]


#' ### Save files for printing
#' 
#' Significance test
#' 
results$significance <- ifelse(results$effectSizeP > 0.1,
                               "",
                               ifelse(results$effectSizeP >0.05,
                                      ".",
                                      ifelse(results$effectSizeP > 0.01,
                                             "*",
                                             "**")))
#' 
#' 
save(results, file = "data/processed/results.Rdata")


#' ### Footer
#' 
#' ezspin(file = "althea_programs/02_data_analysis.R", out_dir = "output", fig_dir = "figures",keep_md = F)


