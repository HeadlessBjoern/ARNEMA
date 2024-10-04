############################# Load data
library(ggsignif)
library(ggplot2)
library(beeswarm)
library(ggbeeswarm)
library(rstatix)
library(readxl)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(car) 
library(beeswarm)
library(rstatix)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(car) # Load the 'car' package
setwd("~/Documents/R/MA figures")

# Load data
file_path <- "/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_nback.xlsx" # ALL file means means over all participants
dataNbackALL <- read_excel(file_path) # This file has all participants listed with means over conditions

# Extract accuracy and reaction time data for each load
NbackAcc1 <- dataNbackALL$`corrN1`
NbackAcc2 <- dataNbackALL$`corrN2` 
NbackAcc3 <- dataNbackALL$`corrN3`
NbackRT1 <- dataNbackALL$`RTN1`
NbackRT2 <- dataNbackALL$`RTN2` 
NbackRT3 <- dataNbackALL$`RTN3`

############################# Calculate Standard Error of the Mean (SEM)

# Calculate SEM for each load
sem_NbackAcc1 <- sd(NbackAcc1) / sqrt(length(NbackAcc1))
sem_NbackAcc2 <- sd(NbackAcc2) / sqrt(length(NbackAcc2)) # For load 2
sem_NbackAcc3 <- sd(NbackAcc3) / sqrt(length(NbackAcc3))
sem_NbackRT1 <- sd(NbackRT1) / sqrt(length(NbackRT1))
sem_NbackRT2 <- sd(NbackRT2) / sqrt(length(NbackRT2)) # For load 2
sem_NbackRT3 <- sd(NbackRT3) / sqrt(length(NbackRT3))

############################# Reshape the data

# Reshape the data for accuracy
dataNbackALL_long_acc <- dataNbackALL %>%
  pivot_longer(cols = starts_with("corrN"),
               names_to = "Condition",
               values_to = "Accuracy") %>%
  mutate(Condition = factor(gsub("corrN", "", Condition), levels = c("1", "2", "3")))

# Reshape the data for reaction times
dataNbackALL_long_rt <- dataNbackALL %>%
  pivot_longer(cols = starts_with("RTN"),
               names_to = "Condition",
               values_to = "ReactionTime") %>%
  mutate(Condition = factor(gsub("RTN", "", Condition), levels = c("1", "2", "3")))

############################# Beeswarm plots

# Create the beeswarm plot for Accuracy
beeswarm(Accuracy ~ Condition, data = dataNbackALL_long_acc, pch = 20, 
         col = c("blue", "green", "red"), xlab = "N-back load", ylab = "Accuracy [%]", xaxt = "n")
axis(1, at = 1:3, labels = c("1", "2", "3"))

# Add error bars for Accuracy
arrows(x0 = c(1, 2, 3), y0 = c(mean(NbackAcc1) - sem_NbackAcc1, mean(NbackAcc2) - sem_NbackAcc2, mean(NbackAcc3) - sem_NbackAcc3),
       x1 = c(1, 2, 3), y1 = c(mean(NbackAcc1) + sem_NbackAcc1, mean(NbackAcc2) + sem_NbackAcc2, mean(NbackAcc3) + sem_NbackAcc3),
       angle = 90, code = 3, length = 0.05, col = "black")

# Add mean value horizontal lines for Accuracy
segments(x0 = 1:3 - 0.1, y0 = c(mean(NbackAcc1), mean(NbackAcc2), mean(NbackAcc3)),
         x1 = 1:3 + 0.1, y1 = c(mean(NbackAcc1), mean(NbackAcc2), mean(NbackAcc3)), 
         col = "black")

# Add mean value labels for Accuracy
text(x = 1:3 + 0.2, y = c(mean(NbackAcc1), mean(NbackAcc2), mean(NbackAcc3)), 
     labels = sprintf("%.2f", c(mean(NbackAcc1), mean(NbackAcc2), mean(NbackAcc3))), 
     pos = 4)

# Create the beeswarm plot for Reaction Times
beeswarm(ReactionTime ~ Condition, data = dataNbackALL_long_rt, pch = 20, 
         col = c("blue", "green", "red"), xlab = "N-back load", ylab = "Reaction Time [ms]", xaxt = "n")
axis(1, at = 1:3, labels = c("1", "2", "3"))

# Add error bars for Reaction Times
arrows(x0 = c(1, 2, 3), y0 = c(mean(NbackRT1) - sem_NbackRT1, mean(NbackRT2) - sem_NbackRT2, mean(NbackRT3) - sem_NbackRT3),
       x1 = c(1, 2, 3), y1 = c(mean(NbackRT1) + sem_NbackRT1, mean(NbackRT2) + sem_NbackRT2, mean(NbackRT3) + sem_NbackRT3),
       angle = 90, code = 3, length = 0.05, col = "black")

# Add mean value horizontal lines for Reaction Times
segments(x0 = 1:3 - 0.1, y0 = c(mean(NbackRT1), mean(NbackRT2), mean(NbackRT3)),
         x1 = 1:3 + 0.1, y1 = c(mean(NbackRT1), mean(NbackRT2), mean(NbackRT3)), 
         col = "black")

# Add mean value labels for Reaction Times
text(x = 1:3 + 0.2, y = c(mean(NbackRT1), mean(NbackRT2), mean(NbackRT3)), 
     labels = sprintf("%.2f", c(mean(NbackRT1), mean(NbackRT2), mean(NbackRT3))), 
     pos = 4)

############################# Pairwise comparisons

# Pairwise comparisons for Accuracy
pairwise_acc_results <- dataNbackALL_long_acc %>%
  pairwise_wilcox_test(Accuracy ~ Condition, paired = TRUE, p.adjust.method = "BH")

# Viewing the results
pairwise_acc_results

# Pairwise Wilcoxon tests for Reaction Times
pairwise_rt_results <- dataNbackALL_long_rt %>%
  pairwise_wilcox_test(ReactionTime ~ Condition, paired = TRUE, p.adjust.method = "BH")

# Viewing the results for Reaction Times
pairwise_rt_results

###################### MEDIAN

############ Accuracy plot
# Custom function to calculate median and MAD
median_mad_fun <- function(x) {
  med <- median(x)
  mad <- mad(x, constant = 1) # constant adjusts the scaling to be consistent with the standard deviation for normal distributions
  return(c(y = med, ymin = med - mad, ymax = med + mad))
}

# Calculate medians for the conditions
median_NbackAcc1 <- median(dataNbackALL_long_acc$Accuracy[dataNbackALL_long_acc$Condition == "1"])
median_NbackAcc2 <- median(dataNbackALL_long_acc$Accuracy[dataNbackALL_long_acc$Condition == "2"])
median_NbackAcc3 <- median(dataNbackALL_long_acc$Accuracy[dataNbackALL_long_acc$Condition == "3"])

# Create a data frame for median values
medians_df_acc <- data.frame(Condition = c("1", "2", "3"),
                             Median = c(median_NbackAcc1, median_NbackAcc2, median_NbackAcc3))

# Start the accuracy plot with theme adjustments for axis lines
accuracy_plot <- ggplot(dataNbackALL_long_acc, aes(x = Condition, y = Accuracy)) +
  geom_quasirandom(aes(color = Condition), groupOnX = FALSE, size = 2) +
  scale_color_manual(values = c("blue", "green", "red")) +
  labs(x = "N-back load", y = "Accuracy [%]") +
  stat_summary(fun.data = median_mad_fun, geom = "errorbar", width = 0.2) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 95, aes(group = Condition)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5), # Add x-axis line
    axis.line.y = element_line(color = "black", size = 0.5)  # Add y-axis line
  ) +
  geom_signif(comparisons = list(c("1", "2")), annotation = "*", y_position = 100, tip_length = 0.02) +
  geom_signif(comparisons = list(c("1", "3")), annotation = "**", y_position = 101.25, tip_length = 0.02) +
  geom_signif(comparisons = list(c("2", "3")), annotation = "*", y_position = 102.75, tip_length = 0.02) +
  scale_x_discrete(labels = c("1", "2", "3")) +
  geom_line(data = medians_df_acc, aes(x = Condition, y = Median), group = 1, colour = "black", size = 0.5) +
  #geom_text(data = medians_df_acc, aes(x = Condition, y = Median + 1, label = sprintf("%.2f", Median)), vjust = 3, hjust = -0.65) +
  guides(color = guide_legend(title = "N-back load"), shape = guide_legend(title = "N-back load")) +
  theme(legend.text = element_text(size = 12))

# Print the accuracy plot
print(accuracy_plot)

######################## Reaction Times


# Custom function to calculate median and MAD
median_mad_fun <- function(x) {
  med <- median(x)
  mad <- mad(x, constant = 1) # constant adjusts the scaling to be consistent with the standard deviation for normal distributions
  return(c(y = med, ymin = med - mad, ymax = med + mad))
}

# Calculate medians for the conditions
median_NbackRT1 <- median(dataNbackALL_long_rt$ReactionTime[dataNbackALL_long_rt$Condition == "1"])
median_NbackRT2 <- median(dataNbackALL_long_rt$ReactionTime[dataNbackALL_long_rt$Condition == "2"])
median_NbackRT3 <- median(dataNbackALL_long_rt$ReactionTime[dataNbackALL_long_rt$Condition == "3"])

# Create a data frame for median values
medians_df_rt <- data.frame(Condition = c("1", "2", "3"),
                            Median = c(median_NbackRT1, median_NbackRT2, median_NbackRT3))

# Start the reaction time plot with theme adjustments for axis lines
rt_plot <- ggplot(dataNbackALL_long_rt, aes(x = Condition, y = ReactionTime)) +
  geom_quasirandom(aes(color = Condition), groupOnX = FALSE, size = 2) +
  scale_color_manual(values = c("blue", "green", "red")) +
  labs(x = "N-back load", y = "Reaction Time [ms]") +
  stat_summary(fun.data = median_mad_fun, geom = "errorbar", width = 0.2) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 95, aes(group = Condition)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5), # Add x-axis line
    axis.line.y = element_line(color = "black", size = 0.5)  # Add y-axis line
  ) +
  geom_signif(comparisons = list(c("1", "2")), annotation = "n.s.", y_position = 650, tip_length = 0.02) +
  geom_signif(comparisons = list(c("1", "3")), annotation = "n.s.", y_position = 670, tip_length = 0.02) +
  geom_signif(comparisons = list(c("2", "3")), annotation = "n.s.", y_position = 690, tip_length = 0.02) +
  scale_x_discrete(labels = c("1", "2", "3")) +
  geom_line(data = medians_df_rt, aes(x = Condition, y = Median), group = 1, colour = "black", size = 0.5) +
  #geom_text(data = medians_df_rt, aes(x = Condition, y = Median + 20, label = sprintf("%.2f", Median)), vjust = 2.75, hjust = -0.3) +
  guides(color = guide_legend(title = "N-back load"), shape = guide_legend(title = "N-back load")) +
  ylim(200, 700) +
  theme(legend.text = element_text(size = 12))

# Print the reaction time plot
print(rt_plot)


############################## OUTPUT STATS RESULTS

# Custom function to calculate median and MAD, and return them as a named vector
median_mad_fun_detailed_acc <- function(x) {
  med <- median(x)
  mad_val <- mad(x, constant = 1)
  return(c(Median = med, MAD = mad_val))
}

# Apply the function to each condition group for Accuracy
medians_mad_acc <- aggregate(Accuracy ~ Condition, dataNbackALL_long_acc, median_mad_fun_detailed_acc)

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
medians_mad_rt <- aggregate(ReactionTime ~ Condition, dataNbackALL_long_rt, median_mad_fun_detailed_rt)

# Print the results for Reaction Time
print("Medians and MAD for Reaction Time:")
print(medians_mad_rt)

############# CIs
# Constants
mad_to_sd_factor <- 1.4826
z_value <- 1.96  # For 95% CI

# Sample sizes for each condition (replace these with your actual sample sizes)
n_1 <- 10  # Example sample size for condition 1
n_2 <- 10  # Example sample size for condition 2
n_3 <- 10  # Example sample size for condition 3

# Function to calculate CI
calculate_ci <- function(median, mad, n) {
  sd_approx <- mad * mad_to_sd_factor
  ci_half_width <- z_value * sd_approx / sqrt(n)
  lower_ci <- median - ci_half_width
  upper_ci <- median + ci_half_width
  return(c(lower_ci, upper_ci))
}

medians_mad_acc

# Calculate CIs for Accuracy
ci_acc_1 <- calculate_ci(98.989899, 1.010101, n_1)
ci_acc_2 <- calculate_ci(95.454545, 2.525253, n_2)
ci_acc_3 <- calculate_ci(86.868687, 5.050505, n_3)

medians_mad_rt

# Calculate CIs for Reaction Times
ci_rt_1 <- calculate_ci(417.28825, 59.36354, n_1)
ci_rt_2 <- calculate_ci(397.89271, 61.82604, n_2)
ci_rt_3 <- calculate_ci(369.26190, 56.58442, n_3)

# Print the CIs
print(ci_acc_1)
print(ci_acc_2)
print(ci_acc_3)

print(ci_rt_1)
print(ci_rt_2)
print(ci_rt_3)

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
effect_size_acc <- cliffs_delta(dataNbackALL_long_acc, "Accuracy", "Condition")

# Calculate Cliff's delta for Reaction Times
effect_size_rt <- cliffs_delta(dataNbackALL_long_rt, "ReactionTime", "Condition")

# Viewing the effect sizes for Accuracy
effect_size_acc

# Viewing the effect sizes for Reaction Times
effect_size_rt
