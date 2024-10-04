############################# MA SternbergSIM RT and ACC barplots
library(ggsignif)
library(ggplot2)
library(beeswarm)
library(ggbeeswarm)
library(effectsize)
library(rstatix)
library(readxl)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(car) 
setwd("~/Documents/R/MA figures")

# Load data
sternberg_file_path <- "/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_ALL.xlsx"
dataStern <- read_excel(sternberg_file_path)
dataSternALL <- read_excel("/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg.xlsx")

# Extract accuracy and reaction time data for each load
SternAcc2 <- dataSternALL$`Correct2 [%]`
SternAcc4 <- dataSternALL$`Correct4 [%]`
SternAcc6 <- dataSternALL$`Correct6 [%]`
SternAcc8 <- dataSternALL$`Correct8 [%]`
SternRT2 <- dataSternALL$`RT2 [ms]`
SternRT4 <- dataSternALL$`RT4 [ms]`
SternRT6 <- dataSternALL$`RT6 [ms]`
SternRT8 <- dataSternALL$`RT8 [ms]`

# Calculate Standard Error of the Mean (SEM) for accuracy and reaction times
sem_SternAcc2 <- sd(SternAcc2) / sqrt(length(SternAcc2))
sem_SternAcc4 <- sd(SternAcc4) / sqrt(length(SternAcc4))
sem_SternAcc6 <- sd(SternAcc6) / sqrt(length(SternAcc6))
sem_SternAcc8 <- sd(SternAcc8) / sqrt(length(SternAcc8))
sem_SternRT2 <- sd(SternRT2) / sqrt(length(SternRT2))
sem_SternRT4 <- sd(SternRT4) / sqrt(length(SternRT4))
sem_SternRT6 <- sd(SternRT6) / sqrt(length(SternRT6))
sem_SternRT8 <- sd(SternRT8) / sqrt(length(SternRT8))

# Reshape the data for accuracy
dataSternALL_long_acc <- dataSternALL %>%
  pivot_longer(cols = starts_with("Correct"),
               names_to = "Condition",
               values_to = "Accuracy") %>%
  mutate(Condition = factor(gsub("Correct", "", Condition)))

# Reshape the data for reaction times
dataSternALL_long_rt <- dataSternALL %>%
  pivot_longer(cols = starts_with("RT"),
               names_to = "Condition",
               values_to = "ReactionTime") %>%
  mutate(Condition = factor(gsub("RT", "", Condition)))

#############################  Create the beeswarm plot for Accuracy
dataSternALL_long_acc <- dataSternALL_long_acc %>%
  mutate(Participant = rep(1:(n()/4), each = 4))

beeswarm(Accuracy ~ Condition, data = dataSternALL_long_acc, pch = 20, 
         col = c("blue", "green", "black", "red"), xlab = "WM load", ylab = "Accuracy [%]", xaxt = "n")
axis(1, at = 1:4, labels = c("2", "4", "6", "8"))

# Add error bars for Accuracy
# Note: Adjust the values in 'at' argument to match the positions of your beeswarm plot
arrows(x0 = c(1, 2, 3, 4), y0 = c(mean(SternAcc2) - sem_SternAcc2, mean(SternAcc4) - sem_SternAcc4, mean(SternAcc6) - sem_SternAcc6, mean(SternAcc8) - sem_SternAcc8),
       x1 = c(1, 2, 3, 4), y1 = c(mean(SternAcc2) + sem_SternAcc2, mean(SternAcc4) + sem_SternAcc4, mean(SternAcc6) + sem_SternAcc6, mean(SternAcc8) + sem_SternAcc8),
       angle = 90, code = 3, length = 0.05, col = "black")

# Add mean value horizontal lines for Accuracy
segments(x0 = 1:4 - 0.1, y0 = c(mean(SternAcc2), mean(SternAcc4), mean(SternAcc6), mean(SternAcc8)),
         x1 = 1:4 + 0.1, y1 = c(mean(SternAcc2), mean(SternAcc4), mean(SternAcc6), mean(SternAcc8)), 
         col = "black")

# Add mean value labels for Accuracy
text(x = 1:4 + 0.1, y = c(mean(SternAcc2), mean(SternAcc4), mean(SternAcc6), mean(SternAcc8)), 
     labels = sprintf("%.2f", c(mean(SternAcc2), mean(SternAcc4), mean(SternAcc6), mean(SternAcc8))), 
     pos = 4)

############################# Create the beeswarm plot for Reaction Times
dataSternALL_long_rt <- dataSternALL_long_rt %>%
  mutate(Participant = rep(1:(n()/4), each = 4))

beeswarm(ReactionTime ~ Condition, data = dataSternALL_long_rt, pch = 20, 
         col = c("blue", "green", "black", "red"), xlab = "WM load", ylab = "Reaction Time [ms]", xaxt = "n")
axis(1, at = 1:4, labels = c("2", "4", "6", "8"))

# Add error bars for Reaction Times
arrows(x0 = c(1, 2, 3, 4), y0 = c(mean(SternRT2) - sem_SternRT2, mean(SternRT4) - sem_SternRT4, mean(SternRT6) - sem_SternRT6, mean(SternRT8) - sem_SternRT8),
       x1 = c(1, 2, 3, 4), y1 = c(mean(SternRT2) + sem_SternRT2, mean(SternRT4) + sem_SternRT4, mean(SternRT6) + sem_SternRT6, mean(SternRT8) + sem_SternRT8),
       angle = 90, code = 3, length = 0.05, col = "black")

mean_SternRT2 <- mean(SternRT2)
mean_SternRT4 <- mean(SternRT4)
mean_SternRT6 <- mean(SternRT6)
mean_SternRT8 <- mean(SternRT8)
# Add mean value horizontal lines for Reaction Times
segments(x0 = 1:4 - 0.1, y0 = c(mean_SternRT2, mean_SternRT4, mean_SternRT6, mean_SternRT8),
         x1 = 1:4 + 0.1, y1 = c(mean_SternRT2, mean_SternRT4, mean_SternRT6, mean_SternRT8), 
         col = "black")

# Add mean value labels for Reaction Times
text(x = 1:4 + 0.1, y = c(mean_SternRT2, mean_SternRT4, mean_SternRT6, mean_SternRT8), 
     labels = sprintf("%.2f", c(mean_SternRT2, mean_SternRT4, mean_SternRT6, mean_SternRT8)), 
     pos = 4)

#################### Perform the pairwise comparisons

# Pairwise comparisons for Accuracy
pairwise_acc_results <- dataSternALL_long_acc %>%
  pairwise_wilcox_test(Accuracy ~ Condition, paired = TRUE, p.adjust.method = "BH")

# Viewing the results
pairwise_acc_results

# Pairwise Wilcoxon tests for Reaction Times
pairwise_rt_results <- dataSternALL_long_rt %>%
  pairwise_wilcox_test(ReactionTime ~ Condition, paired = TRUE, p.adjust.method = "BH")

# Viewing the results for Reaction Times
pairwise_rt_results

############## Plots with MEDIANS

# Custom function to calculate median and MAD
median_mad_fun <- function(x) {
  med <- median(x)
  mad <- mad(x, constant = 1) # constant adjusts the scaling to be consistent with the standard deviation for normal distributions
  return(c(y = med, ymin = med - mad, ymax = med + mad))
}

# Calculate medians for the conditions
medians_acc <- aggregate(Accuracy ~ Condition, dataSternALL_long_acc, median)

# Create a data frame for median values
medians_df_acc <- data.frame(Condition = c("2 [%]", "4 [%]", "6 [%]", "8 [%]"),
                             Median = medians_acc$Accuracy)

# Define the significance levels based on the p-values
signif_levels <- c("***", "**", "*", "n.s.")

# Start the accuracy plot with theme adjustments for axis lines
accuracy_plot <- ggplot(dataSternALL_long_acc, aes(x = Condition, y = Accuracy)) +
  geom_quasirandom(aes(color = Condition), groupOnX = FALSE, size = 2) +
  scale_color_manual(values = c("blue", "green", "black", "red")) +
  labs(x = "WM load", y = "Accuracy [%]") +
  stat_summary(fun.data = median_mad_fun, geom = "errorbar", width = 0.2) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 95, aes(group = Condition)) +
  geom_line(data = medians_df_acc, aes(x = Condition, y = Median), group = 1, colour = "black", size = 0.5) +
  #geom_text(data = medians_df_acc, aes(x = Condition, y = Median + 1, label = sprintf("%.2f", Median)), vjust = 0.85, hjust = -0.35) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14),
    axis.ticks = element_line(color = "black"),
    axis.line.x = element_line(color = "black", size = 0.5), # Add x-axis line
    axis.line.y = element_line(color = "black", size = 0.5),  # Add y-axis line
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  geom_signif(comparisons = list(c("2 [%]", "4 [%]")), annotation = signif_levels[2], y_position = 100, tip_length = 0.02) +
  geom_signif(comparisons = list(c("2 [%]", "6 [%]")), annotation = signif_levels[2], y_position = 102, tip_length = 0.02) +
  geom_signif(comparisons = list(c("2 [%]", "8 [%]")), annotation = signif_levels[2], y_position = 104, tip_length = 0.02) +
  geom_signif(comparisons = list(c("4 [%]", "6 [%]")), annotation = signif_levels[2], y_position = 106, tip_length = 0.02) +
  geom_signif(comparisons = list(c("4 [%]", "8 [%]")), annotation = signif_levels[2], y_position = 108, tip_length = 0.02) +
  #geom_signif(comparisons = list(c("6 [%]", "8 [%]")), annotation = signif_levels[4], y_position = 111, tip_length = 0.02) +
  scale_x_discrete(labels = c("2", "4", "6", "8")) +
  guides(color = guide_legend(title = "WM load"), shape = guide_legend(title = "WM load")) +
  theme(legend.text = element_text(size = 12))

# Print the accuracy plot
print(accuracy_plot)


################# Reaction Time MEDIANS

# Custom function to calculate median and MAD
median_mad_fun <- function(x) {
  med <- median(x)
  mad <- mad(x, constant = 1) # constant adjusts the scaling to be consistent with the standard deviation for normal distributions
  return(c(y = med, ymin = med - mad, ymax = med + mad))
}

# Calculate medians for the conditions
medians_rt <- aggregate(ReactionTime ~ Condition, dataSternALL_long_rt, median)

# Create a data frame for median values
medians_df_rt <- data.frame(Condition = c("2 [ms]", "4 [ms]", "6 [ms]", "8 [ms]"),
                            Median = medians_rt$ReactionTime)

# Define the significance levels based on the p-values
signif_levels_rt <- c("***", "**", "*", "n.s.")

# Start the reaction time plot with theme adjustments for axis lines
rt_plot <- ggplot(dataSternALL_long_rt, aes(x = Condition, y = ReactionTime)) +
  geom_quasirandom(aes(color = Condition), groupOnX = FALSE, size = 2) +
  scale_color_manual(values = c("blue", "green", "black", "red")) +
  labs(x = "WM load", y = "Reaction Time [ms]") +
  stat_summary(fun.data = median_mad_fun, geom = "errorbar", width = 0.2) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 95, aes(group = Condition)) +
  geom_line(data = medians_df_rt, aes(x = Condition, y = Median), group = 1, colour = "black", size = 0.5) +
  #geom_text(data = medians_df_rt, aes(x = Condition, y = Median + 20, label = sprintf("%.2f", Median)), vjust = 4.25, hjust = -0.35) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14),
    axis.ticks = element_line(color = "black"),
    axis.line.x = element_line(color = "black", size = 0.5), # Add x-axis line
    axis.line.y = element_line(color = "black", size = 0.5),  # Add y-axis line
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  geom_signif(comparisons = list(c("2 [ms]", "4 [ms]")), annotation = signif_levels_rt[2], y_position = 540, tip_length = 0.02) +
  geom_signif(comparisons = list(c("2 [ms]", "6 [ms]")), annotation = signif_levels_rt[2], y_position = 560, tip_length = 0.02) +
  geom_signif(comparisons = list(c("2 [ms]", "8 [ms]")), annotation = signif_levels_rt[2], y_position = 580, tip_length = 0.02) +
  geom_signif(comparisons = list(c("4 [ms]", "6 [ms]")), annotation = signif_levels_rt[3], y_position = 600, tip_length = 0.02) +
  geom_signif(comparisons = list(c("4 [ms]", "8 [ms]")), annotation = signif_levels_rt[2], y_position = 620, tip_length = 0.02) +
  geom_signif(comparisons = list(c("6 [ms]", "8 [ms]")), annotation = signif_levels_rt[3], y_position = 650, tip_length = 0.02) +
  scale_x_discrete(labels = c("2", "4", "6", "8")) +
  guides(color = guide_legend(title = "WM load"), shape = guide_legend(title = "WM load")) +
  theme(legend.text = element_text(size = 12)) +
  ylim(200, 700)

# Print the RT plot
print(rt_plot)

############################## OUTPUT STATS RESULTS

# Custom function to calculate median and MAD, and return them as a named vector
median_mad_fun_detailed_acc <- function(x) {
  med <- median(x)
  mad_val <- mad(x, constant = 1)
  return(c(Median = med, MAD = mad_val))
}

# Apply the function to each condition group for Accuracy
medians_mad_acc <- aggregate(Accuracy ~ Condition, dataSternALL_long_acc, median_mad_fun_detailed_acc)

# Print the results for Accuracy
print("Medians and MAD for Accuracy:")
print(medians_mad_acc)


# Custom function to calculate median and MAD, and return them as a named vector
median_mad_fun_detailed_rt <- function(x) {
  med <- median(x)
  mad_val <- mad(x, constant = 1)
  return(c(Median = med, MAD = mad_val))
}

# Apply the function to each condition group for Reaction Time
medians_mad_rt <- aggregate(ReactionTime ~ Condition, dataSternALL_long_rt, median_mad_fun_detailed_rt)

# Print the results for Reaction Time
print("Medians and MAD for Reaction Time:")
print(medians_mad_rt)

# Extract p-values for Accuracy
p_values_acc <- pairwise_acc_results$p.value
formatted_p_values_acc <- format.pval(p_values_acc, method = "holm", digits = 3)
p_values_acc

# Extract p-values for Reaction Times
p_values_rt <- pairwise_rt_results$p.value
formatted_p_values_rt <- format.pval(p_values_rt, method = "holm", digits = 3)
p_values_rt

############# CIs
# Constants
mad_to_sd_factor <- 1.4826
z_value <- 1.96  # For 95% CI

# Sample sizes for each condition (replace these with your actual sample sizes)
n_2 <- 10  # Example sample size for condition 2
n_4 <- 10  # Example sample size for condition 4
n_6 <- 10  # Example sample size for condition 6
n_8 <- 10  # Example sample size for condition 8

# Function to calculate CI
calculate_ci <- function(median, mad, n) {
  sd_approx <- mad * mad_to_sd_factor
  ci_half_width <- z_value * sd_approx / sqrt(n)
  lower_ci <- median - ci_half_width
  upper_ci <- median + ci_half_width
  return(c(lower_ci, upper_ci))
}

# Calculate CIs for Accuracy
ci_acc_2 <- calculate_ci(97.65625, 1.953125, n_2)
ci_acc_4 <- calculate_ci(91.12347, 2.242762, n_4)
ci_acc_6 <- calculate_ci(75.976375, 6.273129, n_6)
ci_acc_8 <- calculate_ci(76.334039, 1.288863, n_8)

# Calculate CIs for Reaction Times
ci_rt_2 <- calculate_ci(336.51492, 38.33802, n_2)
ci_rt_4 <- calculate_ci(376.64246, 62.14245, n_4)
ci_rt_6 <- calculate_ci(382.08133, 49.91853, n_6)
ci_rt_8 <- calculate_ci(403.02869, 50.17007, n_8)

# Print the CIs
print(ci_acc_2)
print(ci_acc_4)
print(ci_acc_6)
print(ci_acc_8)

print(ci_rt_2)
print(ci_rt_4)
print(ci_rt_6)
print(ci_rt_8)

########## pairwise comparison results

pairwise_acc_results
pairwise_rt_results

# Print medians and MAD for Accuracy and Reaction Times
print("Medians and MAD for Accuracy:")
print(medians_mad_acc)

print("Medians and MAD for Reaction Times:")
print(medians_mad_rt)

# Print the correction method used for p-values
print("Correction Method for p-values: Benjamini-Hochberg (BH)")

# Function to calculate Cliff's delta for each pair
cliffs_delta <- function(data, var, group) {
  combn(unique(data[[group]]), 2, function(x) {
    d <- effsize::cliff.delta(data[data[[group]] == x[1], var, drop = TRUE],
                              data[data[[group]] == x[2], var, drop = TRUE])
    return(data.frame(group1 = x[1], group2 = x[2], cliffs_delta = d$estimate))
  }, simplify = FALSE) %>% do.call(rbind, .)
}

# Calculate Cliff's delta for Accuracy
effect_size_acc <- cliffs_delta(dataSternALL_long_acc, "Accuracy", "Condition")

# Calculate Cliff's delta for Reaction Times
effect_size_rt <- cliffs_delta(dataSternALL_long_rt, "ReactionTime", "Condition")

# Viewing the effect sizes for Accuracy
effect_size_acc

# Viewing the effect sizes for Reaction Times
effect_size_rt
