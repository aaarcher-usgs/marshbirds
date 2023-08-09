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
library(tidyverse)


remove(list = ls())
set.seed(22059747)


#' ## Load data
#' 
#' 
NWBIRD.ALL <- read.csv(file = "data/NWBIRD.csv") 

NWSURV.orig <- read.csv("data/NWSURV.csv") 

#' ## Set global parameters
#' 
bird.broadcast <- c("AMBI","LEBI","PBGR","SORA","VIRA")
years <- c("2015","2016", "2017", "2018")
PTdropped <- c("MNS2", "MNS4", "MNC5", "MNC6", "MNS2", "MNS4", "MNC5", "MNC6","TLS4")

#' remove visual records
#' 
column_names <- colnames(NWBIRD.ALL[,grepl("MIN", colnames(NWBIRD.ALL))])

auditory_filter_function<- function(x){
  case_when(x == "0" ~ "not auditory")
}


NWBIRD.ALL_TO_REMOVE <- NWBIRD.ALL |>
  mutate(across(matches("MIN"), ~ case_when(. == "v" ~ "not auditory",
                                            . == "0" ~ "not auditory",
                                            . == "" ~ "not auditory",
                                             TRUE ~ .))) |>
  filter(if_all(matches("MIN"), ~ . ==  "not auditory"))
recordIDS_to_remove <- NWBIRD.ALL_TO_REMOVE$DetectionRecID


#' ## Clean up data
#' 
#' Point count data
NWBIRD.SURV <- NWBIRD.ALL %>%
  dplyr::filter(OutTimeONLY != "Y" & OutTarget != "Y" & Previous !="Y") %>%
  dplyr::filter(!PointID %in% PTdropped) %>% 
  dplyr::filter(Species %in% bird.broadcast) %>%
  dplyr::filter(is.na(Distance) | Distance <= 400) %>%
  filter(! DetectionRecID %in% recordIDS_to_remove) %>%
  ### droplevels() %>% # drops unused levels. like AMCO, BLTE from Species list
  dplyr::select(SurveyID, DetectionRecID, Year, Round, PointID, 
                Species, Distance, OutTimeONLY, OutTarget, Previous, Time)
# Subset for only birds detected DURING the survey time, at the TARGET wetland, NOT previously detected
# 1252 obs, now 1495 w 2017 data added; 1635 w 2018 data
# 1635 ALL; 1051 SURV

#' 
#' Merge NWBIRD.SURV and NWSURV to create NWDATA.SURV
NWDATA.SURV <- NWSURV.orig %>%
  filter(!SiteName == "Manston") %>%
  left_join(NWBIRD.SURV, by=c("SurveyID", "Year", "PointID", "Round")) %>%
  droplevels() %>%   
  mutate("Year" = factor(Year)) %>% # Recode factor(SurveyID) and factor(year)
  mutate("SurveyID" = factor(SurveyID)) %>%
  dplyr::select(SurveyID, Year, Round, PointID, Treatment, Species, SiteName, Time)
# OutTimeONLY, OutTarget, Previous)  
# is.na(Species) were surveys when no broadcast species detected

#'
#'  Tally bird detections by species for each visit
NWPOINT <- NWDATA.SURV %>%
  pivot_wider(#data = NWDATA.SURV,
              #id_cols = c(PointID, Treatment, SiteName, Year, Round, Species), 
              names_from = Species, 
              names_sort = TRUE,
              values_from = Species,
              values_fn = list(Species = length),
              values_fill = 0) 



    

#' Determine when each plot ID was visited for surveys
visits <- NWSURV.orig %>%
  dplyr::select(PointID, Year, Round, SurveyID) %>%
  pivot_wider(id_cols = PointID, names_from=c(Year, Round), 
              values_from = SurveyID, 
              values_fn = list(SurveyID = length))

# Remove the MN points 
visits <- visits[-grep("MN", visits$PointID),]

#' How many times were each surveyed by year?
colnames(visits) <- c("PointID", "R20151", "R20152", "R20161", "R20162", "R20172", "R20182")
colSums(visits[,c("R20151", "R20152", "R20161", "R20162", "R20172", "R20182")], 
        na.rm = T)
# 321 surveys, 313 without Manston

#' How many times were each surveyed by time of day?
#' 



#' Summarize the counts by survey ID and species
#' 
counts <- NWDATA.SURV %>%
  dplyr::select(Species, PointID, Year, Round, SurveyID) %>%
  pivot_wider(id_cols = PointID, names_from=c(Year, Round, Species),
              values_from = SurveyID,
              values_fn = list(SurveyID = length))

#' Compile counts 
#' 
aa <- colSums(counts[,2:ncol(counts)], na.rm = T)
aa
sum(aa, na.rm = T)


#' Compile data for GLMM format
#' 
#' First, make the treatment, pointID, and Year into factors
NWPOINT <- NWPOINT %>% 
  dplyr::mutate(TreatmentN = as.factor(Treatment)) %>% 
  dplyr::mutate(PointIDN = as.factor(PointID)) %>% 
  dplyr::mutate(YearN = as.integer(Year)) 

#' Subset Round 2 data
#' 
NWPOINTR2 <- NWPOINT %>% 
  dplyr::filter(Round == 2)
    

#' Also, convert NWPOINTR2 to long version so "species" can be a variable in GLMM model
#' 
NWPOINTR2.long <- NWPOINTR2 %>%
  pivot_longer(cols = c(AMBI, LEBI, PBGR, SORA, VIRA),
               names_to = "Species",
               values_to = "Count")


#' View Surveys by each point ID
#' 
(table(NWPOINTR2$PointID, NWPOINTR2$Year)) # Year and Round visits to each PointID

#' Save data
#' 
#' 
save(NWPOINTR2, NWPOINTR2.long, file = "data/processed/NWPOINTR2_NV.Rdata")


    

#' ### Footer
#' 
#' ezspin(file = "althea_programs_novisual/01_data_processing_NV.R", out_dir = "output", fig_dir = "figures",keep_md = F)


