# In this script, we applied the proposed model to our real RSV FRNT experimental data
# It consists of the following four parts:
# 1) Data cleaning,
# 2) Applying standard methods — including the Kärber method and the 4PL model — to the real data,
# 3) Fitting Bayesian hierarchical model (BHM) to the real data,
# 4) Creating figures to visualize and interpret the results.

# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Source custom functions and model scripts
source("1_Func_Karber_and_4PL.R")
source("2_Model_Stan.R")

# Set seed ----
set.seed(20250411)

# 1. Clean experimental data ----
source("5_Clean_experimental_data.R")
data$stock <- as.numeric(factor(data$vd, levels = c(200, 300, 330, 400), 1:4))
data_VC$stock <- as.numeric(factor(data_VC$vd, levels = c(200, 300, 330, 400), 1:4))
VC_data = data_VC[, c("batch", "stock", "Count")]

# 2. Apply standard methods to the experimental data ----
which_sto_per_bat <- VC_data %>% group_by(batch, stock) %>% summarise(mean_obs_vc = mean(Count))
which_bat_per_sam <- data %>% group_by(Serum_id, batch) %>% slice(1) %>% 
  dplyr::select(Serum_id, batch)
data <- left_join(data, which_sto_per_bat)

# Kärber method
Karber_data <- Karber_function(data = data,
                               sam_id_var_name = "Serum_id",
                               bat_id_var_name = "batch",
                               dilution_var_name = "dilution", 
                               Count_var_name = "Count", 
                               mean_vc_var_name = "mean_obs_vc",
                               max_dilution = 87480,
                               dilution_interval = 3)

# 4PL Model
FRNT50 <- DRM_function(data = data,
                       sam_id_var_name = "Serum_id",
                       bat_id_var_name = "batch",
                       dilution_var_name = "dilution", 
                       Count_var_name = "Count", 
                       mean_vc_var_name = "mean_obs_vc",
                       max_dilution = 87480,
                       min_dilution = 40,
                       dilution_interval = 3)

# Merge Kärber and 4PL results
Karber_and_4PL <- left_join(Karber_data, FRNT50[, c(1,2,7)])
Karber_and_4PL[Karber_and_4PL$Nt.Karber <= log2(40/3), ]$Nt.Karber <-  log2(40/3)

# Check results
# plot(Karber_and_4PL$Nt.Karber, Karber_and_4PL$Nt.4pl)
# abline(0, 1)

rm(Karber_data, FRNT50)
rm(Karber_function, DRM_function)

# 3. Fit BHM Model ----
#### Estimate effect of each working virus stock -----
data1_stan <- list(num_sto_vc = 4,
                   num_bat_vc = length(unique(data_VC$batch)),
                   num_obs_vc = nrow(data_VC),
                   which_sto_vc = data_VC$stock,
                   which_bat_vc = data_VC$batch,
                   count_vc = data_VC$Count,
                   
                   which_sto_per_bat = which_sto_per_bat$stock,
                   # Priors
                   delta_prior = c(48, 44, 40, 25),
                   sample_from_prior = 0)

# model <- stan_model(model_code=model1_code)
# fit1 <- stan(model_code=model1_code, data=data1_stan,
#              chains=4, cores=4, iter=3000,
#              control = list(adapt_delta = 0.95))
# saveRDS(fit1, file = "../3. Output/Real_fit1.rds")
fit1 <- readRDS("../3. Output/Real_fit1.rds")
posterior1 <- extract(fit1)

# Check results
# plot(which_sto_per_bat$mean_obs_vc, 
#      apply(posterior1$mu_vc_per_bat, 2, median))
# abline(0, 1)

## Calculate Foci reduction data (FR) ----
sto_param <- data.frame(stock = 1:4,
                        post_sto_effect = apply(posterior1$delta, 2, median))

data <- left_join(data, sto_param)
data$FR <- 1 - data$Count/data$post_sto_effect

df_PC <- data %>% filter(Serum_id == "PC")
df_IS1 <- data %>% filter(Serum_id == "IS500")
df_IS2 <- data %>% filter(Serum_id == "IS1000")
df_sam <- data %>% filter(!Serum_id %in% c("PC", "IS500", "IS1000"))

rm(sto_param)

## Remove 'clear negative' samples ----
pos_sam <- df_sam %>% filter(dilution == 40) %>% 
  group_by(Serum_id, batch) %>% 
  summarise(max_y = mean(FR)) %>% 
  filter(max_y > 0.3) %>% 
  arrange(batch, Serum_id) 
pos_sam$sam_id <- 1:nrow(pos_sam)

df_sam <- df_sam %>% filter(Serum_id %in% pos_sam$Serum_id)
df_sam <- left_join(df_sam, pos_sam[, c("Serum_id", "sam_id")])

## Estimate batch effects and nAb titers ----
data2_stan <- list(num_bat_total = 28,
                   
                   num_obs_pc = nrow(df_PC),
                   which_bat_pc = df_PC$batch,
                   x_log_pc = log(df_PC$dilution),
                   y_pc = df_PC$FR,
                   
                   num_obs_IS1 = nrow(df_IS1),
                   which_bat_IS1 = df_IS1$batch,
                   x_log_IS1 = log(df_IS1$dilution),
                   y_IS1 = df_IS1$FR,
                   
                   num_obs_IS2 = nrow(df_IS2),
                   which_bat_IS2 = df_IS2$batch,
                   x_log_IS2 = log(df_IS2$dilution),
                   y_IS2 = df_IS2$FR,
                   
                   num_sam = length(unique(df_sam$sam_id)),
                   num_obs = nrow(df_sam),
                   which_bat = df_sam$batch,
                   which_sam = df_sam$sam_id,
                   x_log = log(df_sam$dilution),
                   y = df_sam$FR,
                   
                   f_pc_prior = c(1, 0, 1, 6),
                   f_IS1_prior = c(1, 0, 1, 6),   
                   f_IS2_prior = c(1, 0, 1, 6),   
                   f_pop_prior = c(1, 0, 1, 6),   
                   
                   f_lower = c(0, -1, 0, log(40)),
                   f_upper = c(2,  1, 1, log(87480)),
                   
                   which_bat_per_sam = pos_sam$batch,
                   
                   sample_from_prior = 0)

# model <- stan_model(model_code=model2_code)
# fit2 <- stan(model_code=model2_code, data=data2_stan,
#              chains=4, cores=4, iter=3000,
#              control = list(adapt_delta = 0.95))
# saveRDS(fit2, file = "../3. Output/Real_fit2.rds")
fit2 <- readRDS("../3. Output/Real_fit2.rds")
posterior2 <- extract(fit2)

# rm(model)
rm(data1_stan, data2_stan, model1_code, model2_code)

# ## Dist plot: distribution of the posterior titer distribution ----
# df_f_pop <- data.frame(posterior2$f_pop)[3001:6000, ]
# colnames(df_f_pop) <- c("b.pop", "c.pop", "d.pop", "e.pop")
# df_sd_pop <- data.frame(posterior2$sd_f_sam_bat)[3001:6000, ]
# colnames(df_sd_pop) <- c("sd.b.pop", "sd.c.pop", "sd.d.pop", "sd.e.pop")
# df_f_pop <- cbind(df_f_pop, df_sd_pop)
# rownames(df_f_pop) <- 1:3000
# df_f_bat <- posterior2$f_bat
# 
# for (i in 1:3000){
#   
#   b_bat <- rep(df_f_bat[3000+i, , 1], each = 50)
#   c_bat <- rep(df_f_bat[3000+i, , 2], each = 50)
#   d_bat <- rep(df_f_bat[3000+i, , 3], each = 50)
#   e_bat <- rep(df_f_bat[3000+i, , 4], each = 50)
#   
#   b <- rnorm(1400, df_f_pop[i, "b.pop"], df_f_pop[i,  "sd.b.pop"]) - b_bat
#   c <- rnorm(1400, df_f_pop[i, "c.pop"], df_f_pop[i,  "sd.c.pop"]) - c_bat
#   d <- rnorm(1400, df_f_pop[i, "d.pop"], df_f_pop[i,  "sd.d.pop"]) - d_bat
#   e <- rnorm(1400, df_f_pop[i, "e.pop"], df_f_pop[i,  "sd.e.pop"]) - e_bat
#   
#   tmp <- data.frame(iter = i, b = b, c = c, d = d, e = e)
#   
#   if (i == 1){
#     output <- tmp
#   } else(
#     output <- rbind(output, tmp)
#   )
# }
# 
# fit2_prior <- readRDS("../3. Output/Simu_fit2_prior.rds")
# posterior2_prior <- extract(fit2_prior)
# 
# df_f_pop_prior <- data.frame(posterior2_prior$f_pop)[3001:6000, ]
# colnames(df_f_pop_prior) <- c("b.pop", "c.pop", "d.pop", "e.pop")
# df_sd_pop_prior <- data.frame(posterior2_prior$sd_f_sam_bat)[3001:6000, ]
# colnames(df_sd_pop_prior) <- c("sd.b.pop", "sd.c.pop", "sd.d.pop", "sd.e.pop")
# df_f_pop_prior <- cbind(df_f_pop_prior, df_sd_pop_prior)
# rownames(df_f_pop) <- 1:3000
# df_f_bat_prior <- posterior2_prior$f_bat
# 
# for (i in 1:3000){
#   
#   b_bat <- rep(df_f_bat_prior[3000+i, , 1], each = 50)
#   c_bat <- rep(df_f_bat_prior[3000+i, , 2], each = 50)
#   d_bat <- rep(df_f_bat_prior[3000+i, , 3], each = 50)
#   e_bat <- rep(df_f_bat_prior[3000+i, , 4], each = 50)
#   
#   b <- rnorm(1400, df_f_pop_prior[i, "b.pop"], df_f_pop[i,  "sd.b.pop"]) - b_bat
#   c <- rnorm(1400, df_f_pop_prior[i, "c.pop"], df_f_pop[i,  "sd.c.pop"]) - c_bat
#   d <- rnorm(1400, df_f_pop_prior[i, "d.pop"], df_f_pop[i,  "sd.d.pop"]) - d_bat
#   e <- rnorm(1400, df_f_pop_prior[i, "e.pop"], df_f_pop[i,  "sd.e.pop"]) - e_bat
#   
#   tmp <- data.frame(iter = i, b = b, c = c, d = d, e = e)
#   
#   if (i == 1){
#     output_prior <- tmp
#   } else(
#     output_prior <- rbind(output_prior, tmp)
#   )
# }
# 
# 
# output <- data.table(output)
# output_prior <- data.table(output_prior)
# 
# output$titer <- as.numeric()
# 
# output[c >= 0.5, titer := log(87480)]
# output[d > 0.5 & c < 0.5, titer := log((d - c) / (0.5 - c) - 1) / b + e,]
# 
# output_prior[b <= 0, titer := log(40/3)]
# output_prior[c >= 0.5, titer := log(87480)]
# output_prior[d <= 0.5, titer := log(40/3)]
# output_prior[d > 0.5 & c < 0.5, titer := log((d - c) / (0.5 - c) - 1) / b + e,]
# 
# output_prior[titer < log(40/3), titer := log(40/3)]
# output_prior[titer > log(87480), titer := log(87480)]
# output[titer < log(40/3), titer := log(40/3)]
# output[titer > log(87480), titer := log(87480)]
# 
# df_tmp <- data.frame(titer = log2(exp(apply(posterior2$Nt_sam_adj, 2, median))))
# df_tmp <- df_tmp %>% mutate(titer = case_when(
#   titer < log2(40/3) ~ log2(40/3),TRUE ~ titer))
# 
# ggplot() +
#   geom_density(data=output_prior, aes(x = log2(exp(titer)), group = iter,
#                                       color ="Draws from prior distribution"),
#                alpha= 0.5, linewidth = 0.2) +
#   geom_density(data=output, aes(x = log2(exp(titer)), group = iter, 
#                                 color = "Draws from posterior distribution"),
#                alpha= 0.5, linewidth = 0.2) +
#   geom_density(data = df_tmp, aes(x = titer,
#                    color = "Posterior median estimates")) +
#   geom_vline(xintercept = log2(40/3), color = "grey", linetype = "dashed") +
#   geom_vline(xintercept = log2(87480), color = "grey", linetype = "dashed") +
#   scale_color_manual("",values = c( "Draws from prior distribution" = "#9ECAE1", 
#                                     "Draws from posterior distribution" = "#FDAE61",
#                                     "Posterior median estimates" = "#D53E4F")) +
#   scale_x_continuous("RSV nAb titers", breaks = seq(4, 16, 2),
#                      labels = c(expression(2^4), expression(2^6), expression(2^8), 
#                                 expression(2^10), expression(2^12), expression(2^14),
#                                 expression(2^16))) +
#   scale_y_continuous("Density") +
#   theme(panel.background = element_blank(),
#         panel.border = element_rect(color = "black", fill = "transparent"),
#         legend.position = c(0.7, 0.8)) -> p
# 
# 
# tiff("../3. Output/Dist-pop-titer.tiff", width = 5, height = 4 , units = "in", res = 600)
# p
# dev.off()
# 
# rm(fit2_prior, posterior2_prior, df_f_pop_prior, df_sd_pop_prior,
#    df_f_pop, df_sd_pop, output, output_prior, p)
# rm(b_bat, c_bat, d_bat, e_bat, b, c, d, e, tmp, i)

## Output nAb titer estimates ---- 
pos_sam$Nt.unadj <- log2(exp(apply(posterior2$Nt_sam, 2, median)))
pos_sam$Nt.adj <- log2(exp(apply(posterior2$Nt_sam_adj, 2, median)))
df <- left_join(Karber_and_4PL, pos_sam)
df <- df %>% dplyr::select(Serum_id, sam_id, batch, 
                           Nt.Karber, Nt.4pl, Nt.unadj, Nt.adj)
df <- data.table(df)
df[Nt.unadj <= log2(40/3), Nt.unadj := log2(40/3)]
df[Nt.adj <= log2(40/3), Nt.adj := log2(40/3)]
df[(!Serum_id %in% c("PC", "IS500", "IS1000")) & (is.na(Nt.unadj)), Nt.unadj := log2(40/3)]
df[(!Serum_id %in% c("PC", "IS500", "IS1000")) & (is.na(Nt.adj)), Nt.adj := log2(40/3)]

Nt_all <- df[!Serum_id %in% c("PC", "IS500", "IS1000"), 
             c("Serum_id", "sam_id", "batch", "Nt.Karber", "Nt.4pl", "Nt.unadj", "Nt.adj")]
Nt_all$IU <- 2^Nt_all$Nt.adj/(exp(median(posterior2$Nt_IS2))/2000)
write.csv(Nt_all, "../3. Output/Final_Nt.csv")

# 4. Figures ----
x_log <- log(unique(df_sam$dilution))

#### VC ----
batch_effect_vc <- data.frame(batch = 1:length(unique(VC_data$batch)),
                              mu_vc_per_bat = apply(posterior1$mu_vc_per_bat, 2, median),
                              mu_vc_lower = apply(posterior1$mu_vc_per_bat, 2, quantile, 0.025),
                              mu_vc_upper = apply(posterior1$mu_vc_per_bat, 2, quantile, 0.975))

ggplot() + 
  geom_jitter(data = data_VC, width = 0.1, 
              aes(x = factor(batch), y = Count, color = factor(stock)), 
              alpha = 0.3, fill = "white", shape = 1) +
  geom_point(data = data_VC %>% group_by(stock, batch) %>% summarise(mean= mean(Count)),
             aes(x = batch, y = mean, color = factor(stock)), position = position_nudge(-0.1)) +
  geom_errorbar(data = data_VC %>% group_by(stock, batch) %>% 
                  summarise(mean= mean(Count), n = n(), sd = sd(Count)) %>%
                  mutate(lower = mean - 1.96 * sd /sqrt(n), 
                         upper = mean + 1.96 * sd /sqrt(n)),
                aes(x = batch, ymin = lower, ymax = upper, 
                    color = factor(stock)), 
                position = position_nudge(-0.1), width = 0) +
  geom_point(data = batch_effect_vc, aes(x = batch, y = mu_vc_per_bat),
             color = "#D7191C", position = position_nudge(0.1)) +
  geom_errorbar(data = batch_effect_vc, 
                aes(x = batch, ymin = mu_vc_lower, ymax = mu_vc_upper), 
                width = 0, color = "#D7191C", position = position_nudge(0.1)) +
  scale_color_manual("Virus working dilutions",
                     breaks = 1:4, 
                     labels = c("1:200", "1:300", "1:330", "1:400"),
                     values = c("#FC8D62","#66C2A5","#999999","#8DA0CB")) +
  scale_x_discrete("Batch id", expand = c(0.02, 0.02)) +
  scale_y_continuous("Count of foci\n", limits = c(2.5, 72.5),
                     breaks = seq(10, 70, 10)) +
  annotate("text", x = 2.5, y = 5.5, label = "Virus Control", size = 3) +
  theme_bw() + theme(legend.position = c(0.5, 0.9),
                     legend.direction = "horizontal",
                     # legend.text = element_text(size = 8),
                     legend.title = element_text(size = 8),
                     legend.key.height = unit(8, "point"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 8),
                     plot.margin = unit(c(5, 5, 5, 5), "point"),
                     panel.grid = element_blank()) -> plot_vc

#### PC ----
df_PC$mu_y <- apply(posterior2$mu_y_pc, 2, median)
df_PC$mu_y_lower <- apply(posterior2$mu_y_pc, 2, quantile, 0.025)
df_PC$mu_y_upper <- apply(posterior2$mu_y_pc, 2, quantile, 0.975)

df_PC$legend <- "BHM-Unadjusted"
df_PC$strip <- paste("Batch ", df_PC$batch, sep="")

PC_adj <- data.frame(x_log = x_log, 
                     y = apply(posterior2$y_pc_adj, 2, median),
                     y_l = apply(posterior2$y_pc_adj, 2, quantile, 0.025),
                     y_u = apply(posterior2$y_pc_adj, 2, quantile, 0.975),
                     legend = "BHM-Adjusted")

ggplot() +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "grey", linewidth = 0.25) +
  geom_line(data = df_PC, 
            aes(x = log(dilution), y = mu_y, group = batch,
                linetype = legend), color = "#FDAE61", linewidth = 0.1) +
  annotate("segment", x = 6.1, xend = 7.2,
           y = 0.5, yend = 0.5,  color = "black", linewidth = 0.5) +
  annotate("segment", x = 6.1, xend = 6.1,
           y = 0.475, yend = 0.525,  color = "black", linewidth = 0.5) +
  annotate("segment", x = 7.2, xend = 7.2,
           y = 0.475, yend = 0.525,  color = "black", linewidth = 0.5) +
  geom_line(data = PC_adj, aes(x = x_log, y = y, linetype = legend), 
            color = "#D7191C", linewidth = 0.5) +
  geom_ribbon(data = PC_adj, aes(x = x_log, ymin = y_l, ymax = y_u), 
              fill = "#D7191C", alpha = 0.1) +  
  scale_linetype_manual("", values = c("BHM-Adjusted" = 1,
                                       "BHM-Unadjusted" = 1)) +
  scale_x_continuous("",limits = c(3.5, 11.56),
                     breaks = x_log, labels = unique(df_PC$dilution)) +
  scale_y_continuous("Foci reduction (%)",
                     limits = c(-0.25, 1), 
                     breaks = seq(-0.2, 1, 0.2),
                     labels = 100 * seq(-0.2, 1, 0.2)) +
  annotate("text", x = log(13000), y = 0.89, 
           label = "Positive Control", size = 3) +
  annotate("point", x = median(posterior2$Nt_PC), y = 0.5, color = "black", size = 1.25) +
  theme_bw() + theme(panel.grid = element_blank(),
                     plot.margin = unit(c(5, 0, 5, 5), "point"),
                     legend.position = c(0.25, 0.2),
                     legend.key.height = unit(11, "point"),
                     legend.background = element_blank(),
                     # legend.text = element_text(size = 8),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                     legend.spacing.y = unit(0, "cm"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 7)) -> plot_PC

df_PC <- df_PC[df_PC$batch %in% c(11, 19, 21), ]
ggplot() +
  geom_point(data = df_PC, aes(x = log(dilution), y = FR, color = factor(batch)), 
             alpha = 1, size = 1, shape = 1) +
  geom_line(data = df_PC, aes(x = log(dilution), y = mu_y, group = batch,
                              color = factor(batch), linetype = legend), linewidth = 0.25) +
  geom_ribbon(data = df_PC, aes(x = log(dilution), ymin = mu_y_lower, group = batch,
                                ymax = mu_y_upper, fill = factor(batch)), alpha = 0.1) +    
  geom_line(data = PC_adj, aes(x = x_log, y = y, linetype = legend), color = "#D7191C", linewidth = 0.5) +
  geom_ribbon(data = PC_adj, aes(x = x_log, ymin = y_l, ymax = y_u), fill = "#D7191C", alpha = 0.1) +  
  scale_color_manual("", values = brewer.pal(4, "Set1")[2:4], 
                     labels = c("Batch 11", "Batch 19", "Batch 21")) +
  scale_fill_manual("", values = brewer.pal(4, "Set1")[2:4], 
                    labels = c("Batch 11", "Batch 19", "Batch 21")) +
  scale_linetype_manual("", values = c("BHM-Adjusted" = "solid",
                                       "BHM-Unadjusted" =  "dashed")) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "grey", linewidth = 0.25) +
  scale_x_continuous("Serial dilution",
                     breaks = x_log, labels = unique(df_PC$dilution)) +
  scale_y_continuous("Foci reduction (%)",
                     limits = c(-0.5, 1), 
                     breaks = seq(-0.5, 1, 0.25),
                     labels = 100 * seq(-0.5, 1, 0.25)) +
  annotate("text", x = log(23000), y = 0.89, 
           label = "Positive Control", size = 3) +
  annotate("point", x = median(posterior2$Nt_PC), y = 0.5, color = "#D7191C") +
  theme_bw() + theme(panel.grid = element_blank(),
                     plot.margin = unit(c(15, 0, 5, 5), "point"),
                     legend.position = "bottom",
                     legend.key.height = unit(11, "point"),
                     legend.background = element_blank(),
                     # legend.text = element_text(size = 8),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                     legend.spacing.y = unit(0, "cm"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 7)) -> plot_PC_supp

#### IS500 ----
df_IS1$mu_y <- apply(posterior2$mu_y_IS1, 2, median)
df_IS1$mu_y_lower <- apply(posterior2$mu_y_IS1, 2, quantile, 0.025)
df_IS1$mu_y_upper <- apply(posterior2$mu_y_IS1, 2, quantile, 0.975)
df_IS1$legend <- "BHM-Unadjusted"

IS1_adj <- data.frame(x_log = x_log, 
                      y = apply(posterior2$y_IS1_adj, 2, median),
                      y_l = apply(posterior2$y_IS1_adj, 2, quantile, 0.025),
                      y_u = apply(posterior2$y_IS1_adj, 2, quantile, 0.975),
                      legend = "BHM-Adjusted")

ggplot() +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "grey", linewidth = 0.25) +
  geom_line(data = df_IS1, 
            aes(x = log(dilution), y = mu_y, 
                group = factor(batch), linetype = legend), 
            linewidth = 0.1, color = "#FDAE61") +
  annotate("segment", x = 6.03, xend = 6.85,
           y = 0.5, yend = 0.5,  color = "black", linewidth = 0.5) +
  annotate("segment", x = 6.03, xend = 6.03,
           y = 0.475, yend = 0.525,  color = "black", linewidth = 0.5) +
  annotate("segment", x = 6.85, xend = 6.85,
           y = 0.475, yend = 0.525,  color = "black", linewidth = 0.5) +
  geom_line(data = IS1_adj, aes(x = x_log, y = y, linetype = legend), 
            color = "#D7191C", linewidth = 0.5) +
  geom_ribbon(data = IS1_adj, aes(x = x_log, ymin = y_l, ymax = y_u), 
              fill = "#D7191C", alpha = 0.1) +  
  scale_linetype_manual("", values = c("BHM-Adjusted" = 1,
                                       "BHM-Unadjusted" = 1)) +
  scale_x_continuous("Serial dilution", limits = c(3.5, 11.56),
                     breaks = x_log, labels = unique(df_IS1$dilution)) +
  scale_y_continuous("", limits = c(-0.25, 1), breaks = seq(-0.2, 1, 0.2)) +
  annotate("text", x = log(8000), y = 0.875, 
           label = "International Standard\n(Concentration=500 IU)", size = 3) +
  annotate("point", x = median(posterior2$Nt_IS1), y = 0.5, color = "black", size = 1.25) +
  theme_bw() + theme(panel.grid = element_blank(),
                     plot.margin = unit(c(5, 0, 5, 0), "point"),
                     legend.position = c(0.25, 0.2),
                     legend.key.height = unit(11, "point"),
                     legend.background = element_blank(),
                     # legend.text = element_text(size = 8),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                     legend.spacing.y = unit(0, "cm"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 7)) -> plot_IS1


#### IS1000 ----
df_IS2$mu_y <- apply(posterior2$mu_y_IS2, 2, median)
df_IS2$mu_y_lower <- apply(posterior2$mu_y_IS2, 2, quantile, 0.025)
df_IS2$mu_y_upper <- apply(posterior2$mu_y_IS2, 2, quantile, 0.975)
df_IS2$legend <- "BHM-Unadjusted"

IS2_adj <- data.frame(x_log = x_log, 
                      y = apply(posterior2$y_IS2_adj, 2, median),
                      y_l = apply(posterior2$y_IS2_adj, 2, quantile, 0.025),
                      y_u = apply(posterior2$y_IS2_adj, 2, quantile, 0.975),
                      legend = "BHM-Adjusted")

ggplot() +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "grey", linewidth = 0.25) +
  geom_line(data = df_IS2, 
            aes(x = log(dilution), y = mu_y, 
                group = factor(batch), linetype = legend), 
            linewidth = 0.1, color = "#FDAE61") +
  annotate("segment", x = 6.4, xend = 7.35,
           y = 0.5, yend = 0.5,  color = "black", linewidth = 0.5) +
  annotate("segment", x = 6.4, xend = 6.4,
           y = 0.475, yend = 0.525,  color = "black", linewidth = 0.5) +
  annotate("segment", x = 7.35, xend = 7.35,
           y = 0.475, yend = 0.525,  color = "black", linewidth = 0.5) +
  geom_line(data = IS2_adj, aes(x = x_log, y = y, linetype = legend), color = "#D7191C", linewidth = 0.5) +
  geom_ribbon(data = IS2_adj, aes(x = x_log, ymin = y_l, ymax = y_u), fill = "#D7191C", alpha = 0.1) +  
  scale_linetype_manual("", values = c("BHM-Adjusted" = 1,
                                       "BHM-Unadjusted" = 1)) +
  scale_x_continuous("", limits = c(3.5, 11.56),
                     breaks = x_log, labels = unique(df_IS2$dilution)) +
  scale_y_continuous("", limits = c(-0.25, 1), breaks = seq(-0.2, 1, 0.2)) +
  annotate("text", x = log(8000), y = 0.875, 
           label = "International Standard\n(Concentration=1000 IU)", size = 3) +
  annotate("point", x = median(posterior2$Nt_IS2), y = 0.5, color = "black", size = 1.25) +
  theme_bw() + theme(panel.grid = element_blank(),
                     plot.margin = unit(c(5, 5, 5, 0), "point"),
                     legend.position = c(0.25, 0.2),
                     legend.key.height = unit(11, "point"),
                     legend.background = element_blank(),
                     # legend.text = element_text(size = 8),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                     legend.spacing.y = unit(0, "cm"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 7)) -> plot_IS2

df_IS2 <- df_IS2[df_IS2$batch %in% c(11, 19, 21), ]
ggplot() +
  geom_point(data = df_IS2, aes(x = log(dilution), y = FR, color = factor(batch)), 
             alpha = 1, size = 1, shape = 1) +
  geom_line(data = df_IS2, aes(x = log(dilution), y = mu_y, group = batch,
                                 color = factor(batch), linetype = legend), linewidth = 0.25) +
  geom_ribbon(data = df_IS2, aes(x = log(dilution), ymin = mu_y_lower, group = batch,
                                   ymax = mu_y_upper, fill = factor(batch)), alpha = 0.1) +    
  geom_line(data = IS2_adj, aes(x = x_log, y = y, linetype = legend), color = "#D7191C", linewidth = 0.5) +
  geom_ribbon(data = IS2_adj, aes(x = x_log, ymin = y_l, ymax = y_u), fill = "#D7191C", alpha = 0.1) +  
  scale_color_manual("", values = brewer.pal(4, "Set1")[2:4], 
                     labels = c("Batch 11", "Batch 19", "Batch 21")) +
  scale_fill_manual("", values = brewer.pal(4, "Set1")[2:4], 
                    labels = c("Batch 11", "Batch 19", "Batch 21")) +
  scale_linetype_manual("", values = c("BHM-Adjusted" = "solid",
                                       "BHM-Unadjusted" =  "dashed")) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "grey", linewidth = 0.25) +
  scale_x_continuous("Serial dilution",
                     breaks = x_log, labels = unique(df_IS2$dilution)) +
  scale_y_continuous("", limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.25)) +
  annotate("text", x = log(13000), y = 0.85, 
           label = "International Standard\n(Concentration=1000 IU)", size = 3) +
  annotate("point", x = median(posterior2$Nt_IS2), y = 0.5, color = "#D7191C") +
  theme_bw() + theme(panel.grid = element_blank(),
                     plot.margin = unit(c(5, 5, 5, 0), "point"),
                     legend.position = "bottom",
                     legend.key.height = unit(11, "point"),
                     legend.background = element_blank(),
                     # legend.text = element_text(size = 8),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                     legend.spacing.y = unit(0, "cm"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 7)) +
  guides(linetype = guide_legend(order = 2),
         color = guide_legend(order = 1),
         fill = guide_legend(order = 1)) -> plot_IS2_supp

#### Sample ----
x_log1 <- seq(x_log[1], x_log[8], 0.1)
df_sam$FR1 <- 1 - df_sam$Count/df_sam$mean_obs_vc
df_sam$mu_y_BHM_unadj <- apply(posterior2$mu_y, 2, median)
mu_y_BHM_adj <- apply(posterior2$y_sam_adj, c(2, 3), median)

df1 <- Nt_all %>% filter(sam_id %in% c(378, 372, 387, 
                                       801, 826, 838, 
                                       904, 907, 936))
df1 <- df1 %>% arrange(batch, sam_id)
df1$id <- 1:nrow(df1)
df1$legend.adj <- "BHM-Adjusted"
df1$legend.unadj <- "BHM-Unadjusted"
df1$legend.4pl <- "4PL Model"
df1$legend.karber <- "Kärber Method"
df1$strip <- paste("Sample ", df1$sam_id, ", Batch ", df1$batch, sep="")

plot_data <- df_sam %>% filter((sam_id %in% df1$sam_id))
plot_data <- left_join(plot_data, df1)

for (j in unique(plot_data$sam_id)){

  k = df1[df1$sam_id == j, ]$batch
  
  # 4PL model
  tmp <- df_sam[df_sam$sam_id == j, ] %>% 
    group_by(dilution) %>% 
    summarise(FR1 = mean(FR1))
  model <- drm(data = tmp, FR1 ~ dilution, fct = LL.4())
  FR_4pl.j <- predict(model, newdata = data.frame(dilution = exp(x_log1)))
  FR_4pl.j <- data.frame(sam_id = j, 
                         batch = k,
                         x_log = x_log1,
                         y_m = FR_4pl.j)
  # BHM-Adjusted
  BHM_adj.j <- data.frame(sam_id = j, 
                          batch = k,
                          x_log = x_log,
                          y_m = mu_y_BHM_adj[, j])
  if (j == unique(plot_data$sam_id)[1]){
    FR_4pl <- FR_4pl.j
    BHM_adj <- BHM_adj.j
  } else{
    FR_4pl <- rbind(FR_4pl, FR_4pl.j)
    BHM_adj <- rbind(BHM_adj, BHM_adj.j)
  }
}

FR_4pl <- left_join(FR_4pl, df1)
BHM_adj <- left_join(BHM_adj, df1)

ggplot() + 
  facet_wrap(~strip, ncol = 3)+
  scale_color_manual("", values = c("Kärber Method" = "#2B83BA",
                                    "4PL Model" =  "#41AB5D",
                                    "BHM-Unadjusted" = "#FDAE61",
                                    "BHM-Adjusted" = "#D7191C")) +
  # Kärber Method
  geom_point(data = df1, aes(x = log(2^Nt.Karber), y = -0.4, color = legend.karber),
             shape = 6, size = 1) +
  
  # 4PL Model
  geom_point(data = plot_data,
             aes(x = log(dilution), y = FR1), color = "#41AB5D",
             size = 0.8, alpha = 0.8, shape = 1) +
  geom_line(data = FR_4pl,
            aes(x = x_log, y = y_m, color = legend.4pl), linewidth = 0.5) +
  geom_point(data = df1, aes(x = log(2^Nt.4pl), y = -0.4, color = legend.4pl),
             shape = 4, size = 1) +
  geom_segment(data = df1, linetype = "dashed", linewidth = 0.25,
               aes(x = log(2^Nt.4pl), xend = log(2^Nt.4pl),
                   y = -0.3, yend = 0.5, color = legend.4pl)) +
  
  # BHM-Unadjusted
  geom_point(data = plot_data,
             aes(x = log(dilution), y = FR), color = "#FDAE61",
             size = 0.8, alpha = 0.8, shape = 1) +
  geom_line(data = plot_data,
            aes(x = log(dilution), y = mu_y_BHM_unadj, color = legend.unadj), linewidth = 0.5) +
  geom_point(data = df1, aes(x = log(2^Nt.unadj), y = -0.4, color = legend.unadj),
             shape = 7, size = 1) +
  geom_segment(data = df1, linetype = "dashed", linewidth = 0.25,
               aes(x = log(2^Nt.unadj), xend = log(2^Nt.unadj),
                   y = -0.3, yend = 0.5, color = legend.unadj)) +
  
  # BHM-Adjusted
  geom_line(data = BHM_adj,
            aes(x = x_log, y = y_m, color = legend.adj), linewidth = 0.5) +
  geom_point(data = df1, aes(x = log(2^Nt.adj), y = -0.4, color = legend.adj),
             shape = 3, size = 1) +
  geom_segment(data = df1, linetype = "dashed", linewidth = 0.25,
               aes(x = log(2^Nt.adj), xend = log(2^Nt.adj),
                   y = -0.3, yend = 0.5, color = legend.adj)) +
  
  annotate("segment", x = log(40/3), xend = log(87480), y = 0.5, yend = 0.5, 
           linetype = "dashed", linewidth = 0.25, color = "grey") +
  scale_x_continuous("Serial dilution",
                     breaks = x_log, labels = unique(plot_data$dilution)) +
  scale_y_continuous("Foci reduction (%)", limits = c(-0.5, 1),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) + 
  
  theme_bw() + theme(strip.background = element_blank(),
                     panel.grid = element_blank(),
                     plot.margin = unit(c(5, 5, 15, 5), "point"),
                     legend.position = "top",
                     axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) -> curve_sam

#### Titers ----
res <- cor.test(Nt_all$Nt.adj, Nt_all$Nt.Karber, method = "spearman", exact = FALSE)
rho <- round(res$estimate, 2)
rmse <- round(sqrt(mean((Nt_all$Nt.adj - Nt_all$Nt.Karber)^2)), 2)

Nt_all %>% ggplot(aes(x = Nt.adj)) + 
  geom_hline(yintercept = log2(40), linetype = 2, color = "grey") +
  geom_vline(xintercept = log2(40), linetype = 2, color = "grey") +
  geom_point(aes(y = Nt.Karber), color = "#2B83BA", alpha = 0.3) +
  scale_x_continuous("", limits = c(2.5,14), breaks = seq(4, 14, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14))) +
  scale_y_continuous("NAb titer estimates", limits = c(2.5,14), breaks = seq(4, 14, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  annotate("text",  x = 5.7, y = 13.8, hjust = 0, vjust = 1, label = "Kärber Method", 
           size = 2.5, fontface = "bold") +
  annotate("text", x = 9.4, y = 6.5, hjust = 0, vjust = 1, size = 2.5,
           label = paste0("Spearman~italic(rho)==", rho), parse = TRUE) +
  annotate("text", x = 9.4, y = 4.7, hjust = 0, vjust = 1, size = 2.5,
           label = paste0("RMSE==", rmse), parse = TRUE) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(5, 0, 5, 7.5), "point"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7)) -> p1


res <- cor.test(Nt_all$Nt.adj, Nt_all$Nt.4pl, method = "spearman", exact = FALSE)
rho <- round(res$estimate, 2)
rmse <- round(sqrt(mean((Nt_all$Nt.adj - Nt_all$Nt.4pl)^2)), 2)


Nt_all %>% ggplot(aes(x = Nt.adj)) +
  geom_hline(yintercept = log2(40), linetype = 2, color = "grey") +
  geom_vline(xintercept = log2(40), linetype = 2, color = "grey") +
  geom_point(aes(y = Nt.4pl), color = "#41AB5D", alpha = 0.3) +
  scale_x_continuous("BHM adjusted titers", limits = c(2.5,14), breaks = seq(4, 14, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14))) +
  scale_y_continuous("", limits = c(2.5,14), breaks = seq(4, 14, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  annotate("text",  x = 5.7, y = 13.8, hjust = 0, vjust = 1, label = "4PL Model", 
           size = 2.5, fontface = "bold") +
  annotate("text", x = 9.4, y = 6.5, hjust = 0, vjust = 1, size = 2.5,
           label = paste0("Spearman~italic(rho)==", rho), parse = TRUE) +
  annotate("text", x = 9.4, y = 4.7, hjust = 0, vjust = 1, size = 2.5,
           label = paste0("RMSE==", rmse), parse = TRUE) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(5, 0, 5, 0), "point"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7)) -> p2

res <- cor.test(Nt_all$Nt.adj, Nt_all$Nt.unadj, method = "spearman", exact = FALSE)
rho <- round(res$estimate, 2)
rmse <- round(sqrt(mean((Nt_all$Nt.adj - Nt_all$Nt.unadj)^2)), 2)

Nt_all %>% ggplot(aes(x = Nt.adj)) + 
  geom_hline(yintercept = log2(40), linetype = 2, color = "grey") +
  geom_vline(xintercept = log2(40), linetype = 2, color = "grey") +
  geom_point(aes(y = Nt.unadj), color = "#FDAE61", alpha = 0.3) +
  scale_x_continuous("", limits = c(2.5,14), breaks = seq(4, 14, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14))) +
  scale_y_continuous("", limits = c(2.5,14), breaks = seq(4, 14, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  annotate("text", x = 5.7, y = 13.8, hjust = 0, vjust = 1, label = "BHM-Unadjusted", 
           size = 2.5, fontface = "bold") +
  annotate("text", x = 9.4, y = 6.5, hjust = 0, vjust = 1, size = 2.5,
           label = paste0("Spearman~italic(rho)==", rho), parse = TRUE) +
  annotate("text", x = 9.4, y = 4.8, hjust = 0, vjust = 1, size = 2.5,
           label = paste0("RMSE==", rmse), parse = TRUE) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(5, 5, 5, 0), "point"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7)) -> p3

#### Output ----
tiff("../3. Output/Plot_real_curve.tiff", width = 7, height = 8, units = "in", res = 300)
l1 <- ggarrange(plot_PC_supp, plot_IS2_supp, ncol=2, 
                widths = c(1, 1), common.legend = T, legend = "top")
l1 <- annotate_figure(l1,
                      top = text_grob(" "),  
                      fig.lab.face = "plain")
ggarrange(l1, curve_sam, nrow = 2, heights = c(2.7, 5.3),
          labels = c("A. ","B. "),
          font.label = list(size = 10))
dev.off()

tiff("../3. Output/Plot_real_agreement.tiff", width = 8.25, height = 6.4, units = "in", res = 300)
l2 <- ggarrange(plot_PC, plot_IS1, plot_IS2, p1, p2, p3, ncol = 3, nrow = 2,  
                labels = c("B. "," "," ", "C. "," "," "),
                font.label = list(size = 8),
                common.legend = T, legend = "top")
ggarrange(plot_vc, l2, nrow = 2, heights = c(0.9, 1.6),
          labels = c("A. ", " ", ""),
          font.label = list(size = 8))
dev.off()

# 5. Pred. v.s. Obs -----
#### 5_1. PC ----
obs <- data[data$Serum_id == "PC", c("dilution", "FR")]
obs$pred <- apply(posterior2$mu_y_pc, 2, median)
res <- cor.test(obs$pred, obs$FR, method = "spearman", exact = FALSE)
rho <- round(res$estimate, 2)
rmse <- round(sqrt(mean((obs$pred - obs$FR)^2, na.rm = TRUE)), 2)

ggplot(data = obs) +
  geom_point(aes(x = pred, y = FR, color = dilution),
             alpha = 0.5) +
  scale_color_gradientn(
    name = "Dilution",
    trans = "log",
    breaks = unique(obs$dilution),
    labels = unique(obs$dilution),
    colors = c("#D7301F","#FE9929","#FEE391")) +
  scale_x_continuous("Model-fitted foci reduction (%)",
                     limits = c(-0.5, 1.05),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) +
  scale_y_continuous("Observed foci reduction (%)",
                     limits = c(-0.5, 1.05),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) +
  geom_abline(slope = 1, color = "grey", linetype = "dashed") +
  annotate("text", x = -0.5, y = 1, hjust = 0, vjust = 1,
    label = paste0("Spearman~italic(rho)==", rho), parse = TRUE) +
  annotate("text", x = -0.5, y = 0.85, hjust = 0, vjust = 1,
           label = paste0("RMSE==", rmse), parse = TRUE) +
  annotate("text", x = 0.5, y = -0.45, hjust = 0, vjust = 0,
           label = "Positive Control") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.key.height = unit(30, "pt")) -> p1

#### 5_2. IS500 ----
obs <- data[data$Serum_id == "IS500", c("dilution", "FR")]
obs$pred <- apply(posterior2$mu_y_IS1, 2, median)
res <- cor.test(obs$pred, obs$FR, method = "spearman", exact = FALSE)
rho <- round(res$estimate, 2)
rmse <- round(sqrt(mean((obs$pred - obs$FR)^2, na.rm = TRUE)), 2)

ggplot(data = obs) +
  geom_point(aes(x = pred, y = FR, color = dilution),
             alpha = 0.5) +
  scale_color_gradientn(
    name = "Dilution",
    trans = "log",
    breaks = unique(obs$dilution),
    labels = unique(obs$dilution),
    colors = c("#D7301F","#FE9929","#FEE391")) +
  scale_x_continuous("Model-fitted foci reduction (%)",
                     limits = c(-0.5, 1.05),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) +
  scale_y_continuous("Observed foci reduction (%)",
                     limits = c(-0.5, 1.05),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) +
  geom_abline(slope = 1, color = "grey", linetype = "dashed") +
  annotate("text", x = -0.5, y = 1, hjust = 0, vjust = 1,
           label = paste0("Spearman~italic(rho)==", rho), parse = TRUE) +
  annotate("text", x = -0.5, y = 0.85, hjust = 0, vjust = 1,
           label = paste0("RMSE==", rmse), parse = TRUE) +
  annotate("text", x = 0.3, y = -0.45, hjust = 0, vjust = 0,
           label = "International Standard\n(Concentration=500)",) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.key.height = unit(30, "pt")) -> p2

#### 5_3. IS1000 ----
obs <- data[data$Serum_id == "IS1000", c("dilution", "FR")]
obs$pred <- apply(posterior2$mu_y_IS2, 2, median)
res <- cor.test(obs$pred, obs$FR, method = "spearman", exact = FALSE)
rho <- round(res$estimate, 2)
rmse <- round(sqrt(mean((obs$pred - obs$FR)^2, na.rm = TRUE)), 2)

ggplot(data = obs) +
  geom_point(aes(x = pred, y = FR, color = dilution),
             alpha = 0.5) +
  scale_color_gradientn(
    name = "Dilution",
    trans = "log",
    breaks = unique(obs$dilution),
    labels = unique(obs$dilution),
    colors = c("#D7301F","#FE9929","#FEE391")) +
  scale_x_continuous("Model-fitted foci reduction (%)",
                     limits = c(-0.5, 1.05),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) +
  scale_y_continuous("Observed foci reduction (%)",
                     limits = c(-0.5, 1.05),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) +
  geom_abline(slope = 1, color = "grey", linetype = "dashed") +
  annotate("text", x = -0.5, y = 1, hjust = 0, vjust = 1,
           label = paste0("Spearman~italic(rho)==", rho), parse = TRUE) +
  annotate("text", x = -0.5, y = 0.85, hjust = 0, vjust = 1,
           label = paste0("RMSE==", rmse), parse = TRUE) +
  annotate("text", x = 0.3, y = -0.45, hjust = 0, vjust = 0,
           label = "International Standard\n(Concentration=1000)",) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.key.height = unit(30, "pt")) -> p3

#### 5_4. Serum sample ----
obs <- df_sam[, c("dilution", "FR")]
obs$pred <- apply(posterior2$mu_y, 2, median)
res <- cor.test(obs$pred, obs$FR, method = "spearman", exact = FALSE)
rho <- round(res$estimate, 2)
rmse <- round(sqrt(mean((obs$pred - obs$FR)^2, na.rm = TRUE)), 2)

ggplot(data = obs) +
  geom_point(aes(x = pred, y = FR, color = dilution),
             alpha = 0.05) +
  scale_color_gradientn(
    name = "Dilution",
    trans = "log",
    breaks = unique(obs$dilution),
    labels = unique(obs$dilution),
    colors = c("#D7301F","#FE9929","#FEE391")) +
  scale_x_continuous("Model-fitted foci reduction (%)",
                     limits = c(-0.5, 1.05),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) +
  scale_y_continuous("Observed foci reduction (%)",
                     limits = c(-0.5, 1.05),
                     breaks = seq(-0.5, 1, 0.5),
                     labels = 100 * seq(-0.5, 1, 0.5)) +
  geom_abline(slope = 1, color = "grey", linetype = "dashed") +
  annotate("text", x = -0.5, y = 1, hjust = 0, vjust = 1,
           label = paste0("Spearman~italic(rho)==", rho), parse = TRUE) +
  annotate("text", x = -0.5, y = 0.85, hjust = 0, vjust = 1,
           label = paste0("RMSE==", rmse), parse = TRUE) +
  annotate("text", x = 0.5, y = -0.45, hjust = 0, vjust = 0,
           label = "Serum samples",) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.key.height = unit(30, "pt")) -> p4
  
tiff("../3. Output/Real_pred vs obs FR.tiff", width = 8.7, height = 7, units = "in", res = 300)
ggarrange(p4, p1, p2, p3, ncol=2, nrow = 2,
          labels = c("A.", "B.", "C.", "D."),
          common.legend = T, legend = "right")
dev.off()

# 6. Model check ----
#### 6.1 Number of divergence transitions ----
sampler_params <- get_sampler_params(fit2, inc_warmup = FALSE)
divergent_per_chain <- sapply(sampler_params, function(x) sum(x[,"divergent__"]))
divergent_per_chain
#### 6.2 Rhat ----
fit_summary <- summary(fit2)$summary
summary(fit_summary[, "Rhat"])

#### 6.3 Trace plot ----
# Pop-level parameters 
f_pop_samples <- rstan::extract(fit2, pars = c("f_pop[1]", "f_pop[2]", "f_pop[3]",  "f_pop[4]"), permuted = FALSE)
mcmc_trace(f_pop_samples, facet_args = list(nrow = 1)) +
  theme_bw() + 
  theme(plot.margin = unit(c(60, 15, 20, 5), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank()) -> f_pop_plot
# Batch-level parameters 
pars_to_check0 <- c("f_bat[1,1]", "f_bat[1,2]", "f_bat[1,3]", "f_bat[1,4]")
f_bat_samples <- rstan::extract(fit2, pars = pars_to_check0, permuted = FALSE)
mcmc_trace(f_bat_samples, facet_args = list(nrow = 1)) +
  theme_bw() + 
  theme(plot.margin = unit(c(60, 15, 20, 5), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank()) -> f_batch_plot
# Sample-level parameters
pos_sam[pos_sam$sam_id %in% 1:2, ]$batch
pars_to_check1 <- c("f_sam_bat[1,1]", "f_sam_bat[1,2]", "f_sam_bat[1,3]", "f_sam_bat[1,4]")
pars_to_check2 <- c("f_sam_bat[2,1]", "f_sam_bat[2,2]", "f_sam_bat[2,3]", "f_sam_bat[2,4]")
f_sam_bat_samples_1 <- rstan::extract(fit2, pars = pars_to_check1, permuted = FALSE)
f_sam_bat_samples_2 <- rstan::extract(fit2, pars = pars_to_check2, permuted = FALSE)
diff_samples_1 <- f_sam_bat_samples_1 - f_bat_samples 
diff_samples_2 <- f_sam_bat_samples_2 - f_bat_samples 

mcmc_trace(diff_samples_1, facet_args = list(nrow = 1)) +
  theme_bw() + 
  theme(plot.margin = unit(c(60, 15, 20, 5), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank()) -> f_sample_plot1

mcmc_trace(diff_samples_2, facet_args = list(nrow = 1)) +
  theme_bw() + 
  theme(plot.margin = unit(c(20, 15, 10, 5), "pt"),
        strip.background = element_blank(),
        strip.text = element_blank()) -> f_sample_plot2

tiff("../3. Output/Real_trace_plot.tiff", width = 12, height = 10, units = "in", res = 300)
ggarrange(f_pop_plot, f_batch_plot, f_sample_plot1, f_sample_plot2,
          nrow = 4, align = "v", heights = c(1, 1, 1, 0.7),
          labels = c("", "", "", ""),
          common.legend = T, legend = "bottom")
dev.off()
