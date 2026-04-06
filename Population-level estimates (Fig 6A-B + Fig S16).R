# This script was used to assess how different titer estimation methods
# influence RSV serological measures at the population level, 
# including age-stratified seroprevalence, geometric mean titers (GMTs), fold-rises and seroconversion.

# Clear workspace
rm(list = ls())

# Load required libraries
library(readxl)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Import source data ----
GMT <- read_excel("Source_data.xlsx", sheet = "Figure 6A") %>% data.frame()
Seroprevalence <- read_excel("Source_data.xlsx", sheet = "Figure 6B") %>% data.frame()
Fold_rise <- read_excel("Source_data.xlsx", sheet = "Figure S16") %>% data.frame()

GMT$Method <- factor(GMT$Method, 
                     levels = c("Karber Method", 
                                "4PL Model",
                                "BHM-Unadjusted", 
                                "BHM-Adjusted"))
Seroprevalence$Method <- factor(Seroprevalence$Method, 
                                levels = c("Karber Method", 
                                           "4PL Model",
                                           "BHM-Unadjusted",
                                           "BHM-Adjusted"))
Fold_rise$Method <- factor(Fold_rise$Method, 
                           levels = c("Karber Method", 
                                      "4PL Model",
                                      "BHM-Unadjusted",
                                      "BHM-Adjusted")) 
# Figure 6A-B ---- 
## GMTs ----
GMT %>% 
  ggplot(aes(x = Age_grp, color = Method)) +
  geom_point(aes(y = GMT), position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI),
                position = position_dodge(0.4), width = 0) + 
  geom_hline(yintercept = log2(40), linetype = "dashed", color = "grey") +
  scale_color_manual("",
                     breaks = c("Karber Method", "4PL Model", 
                                "BHM-Unadjusted", "BHM-Adjusted"),
                     labels = c("Kärber Method", "4PL Model", 
                                "BHM-Unadjusted", "BHM-Adjusted"),
                     values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age groups") +
  scale_y_continuous("\nGeometric Mean Titers", limits = c(3.5,9.2), breaks = seq(4, 9, 1),
                     labels = c(expression(2^4), expression(2^5), expression(2^6), 
                                expression(2^7), expression(2^8), expression(2^9))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 7)) -> p1

## Seroprevalence ----
ggplot(Seroprevalence, aes(x = Age_grp, y = P_POS, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.2) +
  geom_point(aes(color = Method), position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI, color = Method),
                position = position_dodge(width = 0.9), width = 0.2) +
  scale_color_manual("",
                     breaks = c("Karber Method", "4PL Model", 
                                "BHM-Unadjusted", "BHM-Adjusted"),
                     labels = c("Kärber Method", "4PL Model", 
                                "BHM-Unadjusted", "BHM-Adjusted"),
                     values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_fill_manual("",
                    breaks = c("Karber Method", "4PL Model", 
                               "BHM-Unadjusted", "BHM-Adjusted"),
                    labels = c("Kärber Method", "4PL Model", 
                               "BHM-Unadjusted", "BHM-Adjusted"),
                    values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age groups") +
  scale_y_continuous("\nSeroprevalence (%)", 
                     limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 7)) -> p2

## Output ----
tiff("Fig_6A_B.tiff", width = 6, height = 5, units = "in", res = 600)
ggarrange(p1, p2, nrow = 2, align = "hv",
          labels = c(" ", " "),
          font.label = list(size = 8), 
          common.legend = T, legend = "top")
dev.off()

# Figure S16 ----
#### Seroconvension ----
Fold_rise %>% ggplot() +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -0.08, ymax = -0.02,
           fill = "grey", alpha = 0.2) +
  geom_text(aes(x = Visit, y = -0.05, label = N),
            color = "black", size = 2.5) +
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  geom_bar(aes(x = Visit, y = P_Conv, fill = Method),
           stat = "identity", position = position_dodge(width = 0.9), alpha = 0.2) +
  geom_point(aes(x = Visit, y = P_Conv, fill = Method, colour = Method), 
             position = position_dodge(0.9)) +
  geom_errorbar(aes(x = Visit, colour = Method,
                    ymin = P_Conv_Lower, ymax = P_Conv_Upper), 
                position = position_dodge(0.9), width = 0) +
  scale_color_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_fill_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age at follow-up visits", labels = c("2m", "4m", "6m", "1y", "2y", "3y", "Overall")) +
  scale_y_continuous("\nSeroconversion (%)", 
                     breaks = seq(0, 0.6, 0.1),
                     labels = seq(0, 60, 10)) +
  coord_cartesian(ylim = c(-0.05, 0.6), clip = "off") +
  annotate("segment", x = 0, xend = 8, y = 0, yend = 0, color = "grey", size = 0.3) +
  annotate("segment", x = 0, xend = 8, y = 0.6, yend = 0.6, color = "grey", size = 0.3) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 0.6, color = "grey", size = 0.3) +
  annotate("segment", x = 8, xend = 8, y = 0, yend = 0.6, color = "grey", size = 0.3) +
  annotate("text", x = 0.5, y = -0.05, label = "N =", hjust = 1, size = 2.5) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7)) -> p3

#### 2-Fold titer rise ----
Fold_rise %>% ggplot() +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -0.08, ymax = -0.02,
           fill = "grey", alpha = 0.2) +
  geom_text(aes(x = Visit, y = -0.05, label = N),
            color = "black", size = 2.5) +
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  geom_bar(aes(x = Visit, y = P_2fold, fill = Method),
           stat = "identity", position = position_dodge(width = 0.9), alpha = 0.2) +
  geom_point(aes(x = Visit, y = P_2fold, fill = Method, colour = Method), 
             position = position_dodge(0.9)) +
  geom_errorbar(aes(x = Visit, colour = Method,
                    ymin = P_2fold_Lower, ymax = P_2fold_Upper), 
                position = position_dodge(0.9), width = 0) +
  scale_color_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_fill_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age at follow-up visits", labels = c("2m", "4m", "6m", "1y", "2y", "3y", "Overall")) +
  scale_y_continuous("\nTwo-fold titer rise (%)", 
                     breaks = seq(0, 0.6, 0.1),
                     labels = seq(0, 60, 10)) +
  coord_cartesian(ylim = c(-0.05, 0.6), clip = "off") +
  annotate("segment", x = 0, xend = 8, y = 0, yend = 0, color = "grey", size = 0.3) +
  annotate("segment", x = 0, xend = 8, y = 0.6, yend = 0.6, color = "grey", size = 0.3) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 0.6, color = "grey", size = 0.3) +
  annotate("segment", x = 8, xend = 8, y = 0, yend = 0.6, color = "grey", size = 0.3) +
  annotate("text", x = 0.5, y = -0.05, label = "N =", hjust = 1, size = 2.5) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7)) -> p4

#### 4-Fold titer rise ----
Fold_rise %>% ggplot() +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -0.08, ymax = -0.02,
           fill = "grey", alpha = 0.2) +
  geom_text(aes(x = Visit, y = -0.05, label = N),
            color = "black", size = 2.5) +
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  geom_bar(aes(x = Visit, y = P_4fold, fill = Method),
           stat = "identity", position = position_dodge(width = 0.9), alpha = 0.2) +
  geom_point(aes(x = Visit, y = P_4fold, fill = Method, colour = Method), 
             position = position_dodge(0.9)) +
  geom_errorbar(aes(x = Visit, colour = Method,
                    ymin = P_4fold_Lower, ymax = P_4fold_Upper), 
                position = position_dodge(0.9), width = 0) +
  scale_color_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_fill_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age at follow-up visits", labels = c("2m", "4m", "6m", "1y", "2y", "3y", "Overall")) +
  scale_y_continuous("\nFour-fold titer rise (%)", 
                     breaks = seq(0, 0.6, 0.1),
                     labels = seq(0, 60, 10)) +
  coord_cartesian(ylim = c(-0.05, 0.6), clip = "off") +
  annotate("segment", x = 0, xend = 8, y = 0, yend = 0, color = "grey", size = 0.3) +
  annotate("segment", x = 0, xend = 8, y = 0.6, yend = 0.6, color = "grey", size = 0.3) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 0.6, color = "grey", size = 0.3) +
  annotate("segment", x = 8, xend = 8, y = 0, yend = 0.6, color = "grey", size = 0.3) +
  annotate("text", x = 0.5, y = -0.05, label = "N =", hjust = 1, size = 2.5) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7)) -> p5

#### Output ----
tiff("Fig_S16.tiff", width = 10, height = 2.8, units = "in", res = 600)
ggarrange(p3, p4, p5, nrow = 1, ncol = 3, common.legend = T, 
          labels = c("A ", "B ", "C "),
          font.label = list(size = 8))
dev.off()