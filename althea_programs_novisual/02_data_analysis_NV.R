#' # Data analysis
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
load(file = "data/processed/NWPOINTR2_NV.Rdata")

#' Convert Treatment to factor
#' 
NWPOINTR2$Treatment <- factor(NWPOINTR2$TreatmentN, 
                               levels = c("Spray", "Control"), 
                               labels = c("Herbicide", "Control"))

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
  group_by(Species, Treatment) %>%
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
  group_by(Treatment) %>%
  summarise(sum = sum(Count)) %>%
  print.data.frame()



#' Plot histograms of the counts
#' 
#+ histogram
table(NWPOINTR2.long$Species, NWPOINTR2.long$Count)
ggplot(data = NWPOINTR2.long, aes(x = Count))+
  geom_histogram() +
  facet_wrap(~ Species)

#' Set up a data.frame to hold the results of analysis for plotting, printing
#' 
results.graph <- data.frame(
  Year = rep(rep(2015:2018, each = 2), 6),
  Species = rep(c("AMBI", "LEBI", "PBGR", "SORA", "VIRA", "All"), each = 8),
  Treatment = rep(c("Herbicide", "Control"), 24),
  logEstimate = NA,
  logEstLL = NA,
  logEstUL = NA
)

#' 
results.table <- data.frame(
  Species = unique(results.graph$Species),
  t15ratioEst = NA,
  t15ratioSE = NA,
  t15t.ratio = NA,
  t15p.value = NA,
  t16ratioEst = NA,
  t16ratioSE = NA,
  t16t.ratio = NA,
  t16p.value = NA,
  t17ratioEst = NA,
  t17ratioSE = NA,
  t17t.ratio = NA,
  t17p.value = NA,
  t18ratioEst = NA,
  t18ratioSE = NA,
  t18t.ratio = NA,
  t18p.value = NA,
  df = NA
)

diff.results.table <- data.frame(
  Species = unique(results.graph$Species),
  diff15Est = NA,
  diff15LL = NA,
  diff15UL = NA,
  diff16Est = NA,
  diff16LL = NA,
  diff16UL = NA,
  diff17Est = NA,
  diff17LL = NA,
  diff17UL = NA,
  diff18Est = NA,
  diff18LL = NA,
  diff18UL = NA
)
diff.results.graph <- data.frame(
  Year = rep(rep(2015:2018), 6),
  Species = rep(c("AMBI", "LEBI", "PBGR", "SORA", "VIRA", "All"), each = 4),
  #Treatment = rep(c("Control", "Herbicide"), 12),
  diffEstimate = NA,
  diffEstLL = NA,
  diffEstUL = NA
)
#' ____________________________________________________________________________
#' ### AMBI Analysis
#' 
#' #### GLMM analysis
#' 
#' 
table(NWPOINTR2$AMBI)
var(NWPOINTR2$AMBI)/mean(NWPOINTR2$AMBI)
AMBI_glmerBin <- glmmTMB(AMBI ~ Treatment + Year + Treatment*Year + (1|SiteName),
                             data = NWPOINTR2,
                             family = "poisson")
(AMBI <- summary(AMBI_glmerBin))
(AMBIEmm <- emmeans(AMBI_glmerBin, pairwise~Treatment*Year, level = 0.9, type = "response"))
(AMBIdiff <- eff_size(object = AMBIEmm, sigma = sigma(AMBI_glmerBin), 
                      edf = 193, level = 0.9, type = "response"))
#emmip(AMBI_glmerBin,Treatment ~ Year, type = "response")

#' Results for plotting
#' 
AMBI_emmeans <- summary(AMBIEmm$emmeans)
results.graph[results.graph$Species == "AMBI",
              c("logEstimate", "logEstLL", "logEstUL")] <-
  AMBI_emmeans[,c("rate", "lower.CL", "upper.CL")]

AMBI_contrasts <- summary(AMBIEmm$contrasts)
AMBI_diff <- summary(AMBIdiff)
#2015
results.table[results.table$Species == "AMBI", 
              c("t15ratioEst", "t15ratioSE", "t15t.ratio", "t15p.value")] <- 
  AMBI_contrasts[AMBI_contrasts$contrast == "Herbicide Year2015 / Control Year2015",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "AMBI" & diff.results.graph$Year == 2015, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  AMBI_diff[AMBI_diff$contrast == "(Herbicide Year2015 - Control Year2015)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2016
results.table[results.table$Species == "AMBI", 
              c("t16ratioEst", "t16ratioSE", "t16t.ratio", "t16p.value")] <- 
  AMBI_contrasts[AMBI_contrasts$contrast == "Herbicide Year2016 / Control Year2016",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "AMBI" & diff.results.graph$Year == 2016, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  AMBI_diff[AMBI_diff$contrast == "(Herbicide Year2016 - Control Year2016)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2017
results.table[results.table$Species == "AMBI", 
              c("t15ratioEst", "t17ratioSE", "t17t.ratio", "t17p.value")] <- 
  AMBI_contrasts[AMBI_contrasts$contrast == "Herbicide Year2017 / Control Year2017",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "AMBI" & diff.results.graph$Year == 2017, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  AMBI_diff[AMBI_diff$contrast == "(Herbicide Year2017 - Control Year2017)", 
            c("effect.size", "lower.CL", "upper.CL")]

#2018
results.table[results.table$Species == "AMBI", 
              c("t18ratioEst", "t18ratioSE", "t18t.ratio", "t18p.value")] <- 
  AMBI_contrasts[AMBI_contrasts$contrast == "Herbicide Year2018 / Control Year2018",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "AMBI" & diff.results.graph$Year == 2018, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  AMBI_diff[AMBI_diff$contrast == "(Herbicide Year2018 - Control Year2018)", 
            c("effect.size", "lower.CL", "upper.CL")]


#' ____________________________________________________________________________
#' ### LEBI Analysis
#' 
#' #### GLMM analysis
#' 
#' 
table(NWPOINTR2$LEBI)
var(NWPOINTR2$LEBI)/mean(NWPOINTR2$LEBI)
LEBI_glmerBin <- glmmTMB(LEBI ~ Treatment + Year + Treatment*Year + (1|SiteName),
                         data = NWPOINTR2,
                         family = "poisson")
(LEBI <- summary(LEBI_glmerBin))
(LEBIEmm <- emmeans(LEBI_glmerBin, pairwise ~ Treatment*Year, level = 0.9, type = "response"))
(LEBIdiff <- eff_size(object = LEBIEmm, sigma = sigma(LEBI_glmerBin), 
                      edf = 193, level = 0.9, type = "reponse"))
#emmip(LEBI_glmerBin,Treatment ~ Year, type = "response")

#' Results for plotting
#' 
LEBI_emmeans <- summary(LEBIEmm$emmeans)
results.graph[results.graph$Species == "LEBI",
              c("logEstimate", "logEstLL", "logEstUL")] <-
  LEBI_emmeans[,c("rate", "lower.CL", "upper.CL")]

LEBI_contrasts <- summary(LEBIEmm$contrasts)
LEBI_diff <- summary(LEBIdiff)
#2015
results.table[results.table$Species == "LEBI", 
              c("t15ratioEst", "t15ratioSE", "t15t.ratio", "t15p.value")] <- 
  LEBI_contrasts[LEBI_contrasts$contrast == "Herbicide Year2015 / Control Year2015",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "LEBI" & diff.results.graph$Year == 2015, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  LEBI_diff[LEBI_diff$contrast == "(Herbicide Year2015 - Control Year2015)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2016
results.table[results.table$Species == "LEBI", 
              c("t16ratioEst", "t16ratioSE", "t16t.ratio", "t16p.value")] <- 
  LEBI_contrasts[LEBI_contrasts$contrast == "Herbicide Year2016 / Control Year2016",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "LEBI" & diff.results.graph$Year == 2016, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  LEBI_diff[LEBI_diff$contrast == "(Herbicide Year2016 - Control Year2016)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2017
results.table[results.table$Species == "LEBI", 
              c("t15ratioEst", "t17ratioSE", "t17t.ratio", "t17p.value")] <- 
  LEBI_contrasts[LEBI_contrasts$contrast == "Herbicide Year2017 / Control Year2017",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "LEBI" & diff.results.graph$Year == 2017, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  LEBI_diff[LEBI_diff$contrast == "(Herbicide Year2017 - Control Year2017)", 
            c("effect.size", "lower.CL", "upper.CL")]

#2018
results.table[results.table$Species == "LEBI", 
              c("t18ratioEst", "t18ratioSE", "t18t.ratio", "t18p.value")] <- 
  LEBI_contrasts[LEBI_contrasts$contrast == "Herbicide Year2018 / Control Year2018",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "LEBI" & diff.results.graph$Year == 2018, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  LEBI_diff[LEBI_diff$contrast == "(Herbicide Year2018 - Control Year2018)", 
            c("effect.size", "lower.CL", "upper.CL")]



#' ____________________________________________________________________________
#' ### PBGR Analysis
#' 
#' #### GLMM analysis
#' 
#' 
table(NWPOINTR2$PBGR)
var(NWPOINTR2$PBGR)/mean(NWPOINTR2$PBGR)
PBGR_glmerBin <- glmmTMB(PBGR ~ Treatment + Year + Treatment*Year + (1|SiteName),
                         data = NWPOINTR2,
                         family = "poisson")
(PBGR <- summary(PBGR_glmerBin))
(PBGREmm <- emmeans(PBGR_glmerBin, pairwise ~ Treatment*Year, level = 0.9, type = "response"))
(PBGRdiff <- eff_size(object = PBGREmm, sigma = sigma(PBGR_glmerBin), 
                      edf = 193, level = 0.9, type = "reponse"))
#emmip(PBGR_glmerBin,Treatment ~ Year, type = "response")

#' Results for plotting
#' 
PBGR_emmeans <- summary(PBGREmm$emmeans)
results.graph[results.graph$Species == "PBGR",
              c("logEstimate", "logEstLL", "logEstUL")] <-
  PBGR_emmeans[,c("rate", "lower.CL", "upper.CL")]

PBGR_contrasts <- summary(PBGREmm$contrasts)
PBGR_diff <- summary(PBGRdiff)
#2015
results.table[results.table$Species == "PBGR", 
              c("t15ratioEst", "t15ratioSE", "t15t.ratio", "t15p.value")] <- 
  PBGR_contrasts[PBGR_contrasts$contrast == "Herbicide Year2015 / Control Year2015",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "PBGR" & diff.results.graph$Year == 2015, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  PBGR_diff[PBGR_diff$contrast == "(Herbicide Year2015 - Control Year2015)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2016
results.table[results.table$Species == "PBGR", 
              c("t16ratioEst", "t16ratioSE", "t16t.ratio", "t16p.value")] <- 
  PBGR_contrasts[PBGR_contrasts$contrast == "Herbicide Year2016 / Control Year2016",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "PBGR" & diff.results.graph$Year == 2016, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  PBGR_diff[PBGR_diff$contrast == "(Herbicide Year2016 - Control Year2016)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2017
results.table[results.table$Species == "PBGR", 
              c("t15ratioEst", "t17ratioSE", "t17t.ratio", "t17p.value")] <- 
  PBGR_contrasts[PBGR_contrasts$contrast == "Herbicide Year2017 / Control Year2017",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "PBGR" & diff.results.graph$Year == 2017, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  PBGR_diff[PBGR_diff$contrast == "(Herbicide Year2017 - Control Year2017)", 
            c("effect.size", "lower.CL", "upper.CL")]

#2018
results.table[results.table$Species == "PBGR", 
              c("t18ratioEst", "t18ratioSE", "t18t.ratio", "t18p.value")] <- 
  PBGR_contrasts[PBGR_contrasts$contrast == "Herbicide Year2018 / Control Year2018",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "PBGR" & diff.results.graph$Year == 2018, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  PBGR_diff[PBGR_diff$contrast == "(Herbicide Year2018 - Control Year2018)", 
            c("effect.size", "lower.CL", "upper.CL")]




#' ____________________________________________________________________________
#' ### SORA Analysis
#' 
#' #### GLMM analysis
#' 
#' 
table(NWPOINTR2$SORA)
var(NWPOINTR2$SORA)/mean(NWPOINTR2$SORA)
SORA_glmerBin <- glmmTMB(SORA ~ Treatment + Year + Treatment*Year + (1|SiteName),
                         data = NWPOINTR2,
                         family = "poisson")
(SORA <- summary(SORA_glmerBin))
(SORAEmm <- emmeans(SORA_glmerBin, pairwise ~ Treatment*Year, level = 0.9, type = "response"))
(SORAdiff <- eff_size(object = SORAEmm, sigma = sigma(SORA_glmerBin), 
                      edf = 198, level = 0.90, type = "response"))
#emmip(SORA_glmerBin,Treatment ~ Year, type = "response")

#' Results for plotting
#' 
SORA_emmeans <- summary(SORAEmm$emmeans)
results.graph[results.graph$Species == "SORA",
              c("logEstimate", "logEstLL", "logEstUL")] <-
  SORA_emmeans[,c("rate", "lower.CL", "upper.CL")]

SORA_contrasts <- summary(SORAEmm$contrasts)
SORA_diff <- summary(SORAdiff)
#2015
results.table[results.table$Species == "SORA", 
              c("t15ratioEst", "t15ratioSE", "t15t.ratio", "t15p.value")] <- 
  SORA_contrasts[SORA_contrasts$contrast == "Herbicide Year2015 / Control Year2015",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "SORA" & diff.results.graph$Year == 2015, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  SORA_diff[SORA_diff$contrast == "(Herbicide Year2015 - Control Year2015)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2016
results.table[results.table$Species == "SORA", 
              c("t16ratioEst", "t16ratioSE", "t16t.ratio", "t16p.value")] <- 
  SORA_contrasts[SORA_contrasts$contrast == "Herbicide Year2016 / Control Year2016",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "SORA" & diff.results.graph$Year == 2016, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  SORA_diff[SORA_diff$contrast == "(Herbicide Year2016 - Control Year2016)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2017
results.table[results.table$Species == "SORA", 
              c("t15ratioEst", "t17ratioSE", "t17t.ratio", "t17p.value")] <- 
  SORA_contrasts[SORA_contrasts$contrast == "Herbicide Year2017 / Control Year2017",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "SORA" & diff.results.graph$Year == 2017, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  SORA_diff[SORA_diff$contrast == "(Herbicide Year2017 - Control Year2017)", 
            c("effect.size", "lower.CL", "upper.CL")]

#2018
results.table[results.table$Species == "SORA", 
              c("t18ratioEst", "t18ratioSE", "t18t.ratio", "t18p.value")] <- 
  SORA_contrasts[SORA_contrasts$contrast == "Herbicide Year2018 / Control Year2018",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "SORA" & diff.results.graph$Year == 2018, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  SORA_diff[SORA_diff$contrast == "(Herbicide Year2018 - Control Year2018)", 
            c("effect.size", "lower.CL", "upper.CL")]





#' ____________________________________________________________________________
#' ### VIRA Analysis
#' 
#' #### GLMM analysis
#' 
#' 
table(NWPOINTR2$VIRA)
var(NWPOINTR2$VIRA)/mean(NWPOINTR2$VIRA)
VIRA_glmerBin <- glmmTMB(VIRA ~ Treatment + Year + Treatment*Year + (1|SiteName),
                         data = NWPOINTR2,
                         family = "poisson")
(VIRA <- summary(VIRA_glmerBin))
(VIRAEmm <- emmeans(VIRA_glmerBin, pairwise ~ Treatment*Year, level = 0.9, type = "response"))
(VIRAdiff <- eff_size(object = VIRAEmm, sigma = sigma(VIRA_glmerBin), 
                      edf = 193, level = 0.9, type = "reponse"))
#emmip(VIRA_glmerBin,Treatment ~ Year, type = "response")

#' Results for plotting
#' 
VIRA_emmeans <- summary(VIRAEmm$emmeans)
results.graph[results.graph$Species == "VIRA",
              c("logEstimate", "logEstLL", "logEstUL")] <-
  VIRA_emmeans[,c("rate", "lower.CL", "upper.CL")]

VIRA_contrasts <- summary(VIRAEmm$contrasts)
VIRA_diff <- summary(VIRAdiff)
#2015
results.table[results.table$Species == "VIRA", 
              c("t15ratioEst", "t15ratioSE", "t15t.ratio", "t15p.value")] <- 
  VIRA_contrasts[VIRA_contrasts$contrast == "Herbicide Year2015 / Control Year2015",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "VIRA" & diff.results.graph$Year == 2015, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  VIRA_diff[VIRA_diff$contrast == "(Herbicide Year2015 - Control Year2015)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2016
results.table[results.table$Species == "VIRA", 
              c("t16ratioEst", "t16ratioSE", "t16t.ratio", "t16p.value")] <- 
  VIRA_contrasts[VIRA_contrasts$contrast == "Herbicide Year2016 / Control Year2016",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "VIRA" & diff.results.graph$Year == 2016, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  VIRA_diff[VIRA_diff$contrast == "(Herbicide Year2016 - Control Year2016)", 
            c("effect.size", "lower.CL", "upper.CL")]
#2017
results.table[results.table$Species == "VIRA", 
              c("t15ratioEst", "t17ratioSE", "t17t.ratio", "t17p.value")] <- 
  VIRA_contrasts[VIRA_contrasts$contrast == "Herbicide Year2017 / Control Year2017",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "VIRA" & diff.results.graph$Year == 2017, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  VIRA_diff[VIRA_diff$contrast == "(Herbicide Year2017 - Control Year2017)", 
            c("effect.size", "lower.CL", "upper.CL")]

#2018
results.table[results.table$Species == "VIRA", 
              c("t18ratioEst", "t18ratioSE", "t18t.ratio", "t18p.value")] <- 
  VIRA_contrasts[VIRA_contrasts$contrast == "Herbicide Year2018 / Control Year2018",
                 c("ratio", "SE", "t.ratio", "p.value")]

diff.results.graph[diff.results.graph$Species == "VIRA" & diff.results.graph$Year == 2018, 
                   c("diffEstimate","diffEstLL","diffEstUL")] <-
  VIRA_diff[VIRA_diff$contrast == "(Herbicide Year2018 - Control Year2018)", 
            c("effect.size", "lower.CL", "upper.CL")]



#' ____________________________________________________________________________
#' ### All Species Analysis
#' 
#' #### GLMM analysis
#' 
#' 
table(NWPOINTR2.long$Count)
var(NWPOINTR2.long$Count)/mean(NWPOINTR2.long$Count)
All_glmerBin <- glmmTMB(Count ~ Treatment + Year + Treatment*Year + (1|SiteName),
                         data = NWPOINTR2.long,
                         family = "poisson")
(All <- summary(All_glmerBin))
(AllEmm <- emmeans(All_glmerBin, pairwise ~ Treatment*Year, level = 0.9, type = "response"))
#emmip(All_glmerBin,Treatment ~ Year, type = "response")


#' Test with species effects too
#' 
All_glmerBin2 <- glmmTMB(Count ~ Treatment*Year*Species + (1|SiteName),
                        data = NWPOINTR2.long,
                        family = "poisson", na.action = "na.omit")
(All2 <- summary(All_glmerBin2))
(AllEmm2 <- emmeans(All_glmerBin2, pairwise ~ Treatment*Year*Species, level = 0.9, type = "response"))

#' Results for plotting
#' 
All_emmeans <- summary(AllEmm$emmeans)
results.graph[results.graph$Species == "All",
              c("logEstimate", "logEstLL", "logEstUL")] <-
  All_emmeans[,c("rate", "lower.CL", "upper.CL")]

All_contrasts <- summary(AllEmm$contrasts)
#2015
results.table[results.table$Species == "All", 
              c("t15ratioEst", "t15ratioSE", "t15t.ratio", "t15p.value")] <- 
  All_contrasts[All_contrasts$contrast == "Control Year2015 / Spray Year2015",
                 c("ratio", "SE", "t.ratio", "p.value")]
#2016
results.table[results.table$Species == "All", 
              c("t16ratioEst", "t16ratioSE", "t16t.ratio", "t16p.value")] <- 
  All_contrasts[All_contrasts$contrast == "Control Year2016 / Spray Year2016",
                 c("ratio", "SE", "t.ratio", "p.value")]
#2017
results.table[results.table$Species == "All", 
              c("t17ratioEst", "t17ratioSE", "t17t.ratio", "t17p.value")] <- 
  All_contrasts[All_contrasts$contrast == "Control Year2017 / Spray Year2017",
                 c("ratio", "SE", "t.ratio", "p.value")]
#2018
results.table[results.table$Species == "All", 
              c("t18ratioEst", "t18ratioSE", "t18t.ratio", "t18p.value")] <- 
  All_contrasts[All_contrasts$contrast == "Control Year2018 / Spray Year2018",
                 c("ratio", "SE", "t.ratio", "p.value")]
results.table$df[results.table$Species == "All"] <- All_contrasts$df[1]





#' ### Save files for printing
#' 
#' 
#' 
save(results.graph, results.table, diff.results.graph, file = "data/processed/results_NV.Rdata")


#' ### Footer
#' 
#' ezspin(file = "althea_programs_novisual/02_data_analysis_NV.R", out_dir = "output", fig_dir = "figures",keep_md = F)


