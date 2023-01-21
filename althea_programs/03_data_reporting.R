#' # Results reporting
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
#' Model output results
load(file = "data/processed/results.Rdata")

#' Point counts in both wide and long format
#' 
load(file = "data/processed/NWPOINTR2.Rdata")
NWPOINTR2.long$Treatment <- factor(NWPOINTR2$TreatmentN, 
                                   levels = c("Spray", "Control"), 
                                   labels = c("Herbicide", "Control"))

#' Table 2
#' 
#' Summarize counts 
#' 
NWPOINTR2.long %>%
  group_by(Species, Year, Treatment) %>%
  summarise(sum = sum(Count)) %>%
  print.data.frame()

#' Table 2 statistical results
#' 
results.table %>%
  mutate_if(is.numeric, ~round(.,3))

#' Figure 3
#' 
#+ modelResults, fig.height = 4, fig.width = 8
(p1 <- ggplot(data = results.graph[results.graph$Species != "All",], aes(y = logEstimate, x = Year, group = Treatment)) + 
  #geom_point(data = NWPOINTR2.long, 
  #           aes(y = Count, x = YearN+2014, shape = Treatment),
   #          position = position_jitter(width = 0.1, height = 0.1),
   #          size = 2,
   #          color = "#cccccc", alpha = 0.4)+
  geom_linerange(aes(ymin = logEstLL, ymax = logEstUL, color = Treatment),
                position = position_dodge(width = 0.15)) +
  geom_point(aes(y = logEstimate, color = Treatment, shape = Treatment),
            position = position_dodge(width = 0.15),
            size = 3)+
  facet_wrap(~Species, ncol = 6)+
  ylab("Expected Mean Counts")+
  scale_color_manual(values = c("#969696", "#525252"))+
  theme_bw()+
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        text = element_text(size = 12),
        axis.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12)))

#' Figure 4
#' 
#+ modelResultsDIFF, fig.height = 4, fig.width = 8
(p2 <- ggplot(data = diff.results.graph[diff.results.graph$Species != "All",], aes(y = diffEstimate, x = Year)) + 
  #geom_point(data = NWPOINTR2.long, 
  #           aes(y = Count, x = YearN+2014, shape = Treatment),
  #          position = position_jitter(width = 0.1, height = 0.1),
  #          size = 2,
  #          color = "#cccccc", alpha = 0.4)+
  geom_linerange(aes(ymin = diffEstLL, ymax = diffEstUL),
                 position = position_dodge(width = 0.15)) +
  geom_point(aes(y = diffEstimate),
             position = position_dodge(width = 0.15),
             size = 2.1)+
  facet_wrap(~Species, ncol = 5)+
  geom_hline(yintercept = 0)+
  ylim(c(-2.7, 5.2))+
  ylab("Herbicide-Control")+
  theme_bw()+
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        text = element_text(size = 12),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12)))


#' Figure 4
#' 
#+ modelResultsCombo, fig.height = 5, fig.width = 8
cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(0.6,0.28)) +
  cowplot::draw_text("*", x = 0.791, y = 0.8) + 
  cowplot::draw_text("*", x = 0.979, y = 0.8) #+
  #cowplot::draw_text("*", x = 0.791, y = 0.405) + 
  #cowplot::draw_text("*", x = 0.979, y = 0.405)
#' ### Footer
#' 
#' ezspin(file = "althea_programs/03_data_reporting.R", out_dir = "output", fig_dir = "figures",keep_md = F)


