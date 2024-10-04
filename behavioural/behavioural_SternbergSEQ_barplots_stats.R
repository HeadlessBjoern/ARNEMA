############################# MA SternbergSEQ RT and ACC barplots
library(rstatix)
library(beeswarm)
library(ggbeeswarm)
library(tidyr)
library(ggplot2)
library(effsize)
library(ggsignif)
library(dplyr)
library(RColorBrewer)
library(readxl)
library(car) # Load the 'car' package
setwd("~/Documents/R/MA figures")

############################# Load data

sternberg_file_path <- "/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ_ALL.xlsx" #ALL file means means over all participants
dataStern <- read_excel(sternberg_file_path)
dataSternALL <- read_excel("/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ.xlsx") #This file has all participants listed with means over conditions

# Assuming you have the correct columns for conditions 1, 4, and 7 in your Excel files
SternAcc1 <- dataSternALL$`Correct1 [%]`
SternAcc4 <- dataSternALL$`Correct4 [%]`
SternAcc7 <- dataSternALL$`Correct7 [%]`
SternRT1 <- dataSternALL$`RT1 [ms]`
SternRT4 <- dataSternALL$`RT4 [ms]`
SternRT7 <- dataSternALL$`RT7 [ms]`

############################# Calculate Standard Error of the Mean (SEM)

sem_SternAcc1 <- sd(SternAcc1) / sqrt(length(SternAcc1))
sem_SternAcc4 <- sd(SternAcc4) / sqrt(length(SternAcc4))
sem_SternAcc7 <- sd(SternAcc7) / sqrt(length(SternAcc7))
sem_SternRT1 <- sd(SternRT1) / sqrt(length(SternRT1))
sem_SternRT4 <- sd(SternRT4) / sqrt(length(SternRT4))
sem_SternRT7 <- sd(SternRT7) / sqrt(length(SternRT7))

############################# Sternberg Accuracy

# Create the barplot for Accuracy
SternAccFigure <- barplot(height = c(mean(SternAcc1), mean(SternAcc4), mean(SternAcc7)), 
                          names = c(1, 4, 7), col = c("blue", "green", "red"), 
                          xlab = "WM load", ylab = "Accuracy [%]", ylim = c(0, 110))

# Add error bars for Accuracy
arrows(SternAccFigure, c(mean(SternAcc1) - sem_SternAcc1, mean(SternAcc4) - sem_SternAcc4, mean(SternAcc7) - sem_SternAcc7), 
       SternAccFigure, c(mean(SternAcc1) + sem_SternAcc1, mean(SternAcc4) + sem_SternAcc4, mean(SternAcc7) + sem_SternAcc7), 
       angle = 90, code = 3, length = 0.1, col = "black")

############################# Sternberg Reaction Times

# Create the barplot for Reaction Times
SternRTFigure <- barplot(height=c(mean(SternRT1), mean(SternRT4), mean(SternRT7)), 
                         names=c(1, 4, 7), col = c("blue", "green", "red"), 
                         xlab="WM load", ylab="Reaction Time [ms]", ylim = c(0, max(SternRT1, SternRT4, SternRT7)))

# Add error bars for Reaction Times
arrows(SternRTFigure, c(mean(SternRT1) - sem_SternRT1, mean(SternRT4) - sem_SternRT4, mean(SternRT7) - sem_SternRT7), 
       SternRTFigure, c(mean(SternRT1) + sem_SternRT1, mean(SternRT4) + sem_SternRT4, mean(SternRT7) + sem_SternRT7), 
       angle=90, code=3, length=0.1, col="black")

############################# Test

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


###################### Pairwise comparisons 
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

####################### Accuracy WITH MEDIAN
# Custom function to calculate median and MAD
median_mad_fun <- function(x) {
  med <- median(x)
  mad <- mad(x, constant = 1) # The constant adjusts the scaling to make it consistent with the standard deviation for normal distributions
  return(c(y = med, ymin = med - mad, ymax = med + mad))
}

# Calculate medians for the conditions
median_SternAcc1 <- median(dataSternALL_long_acc$Accuracy[dataSternALL_long_acc$Condition == "1 [%]"])
median_SternAcc4 <- median(dataSternALL_long_acc$Accuracy[dataSternALL_long_acc$Condition == "4 [%]"])
median_SternAcc7 <- median(dataSternALL_long_acc$Accuracy[dataSternALL_long_acc$Condition == "7 [%]"])

# Create a data frame for median values
medians_df <- data.frame(Condition = factor(c("1 [%]", "4 [%]", "7 [%]"), levels = c("1 [%]", "4 [%]", "7 [%]")),
                         Median = c(median_SternAcc1, median_SternAcc4, median_SternAcc7))

# Highest y-value for the vertical line to end
highest_y_value <- max(c(median_SternAcc1, median_SternAcc4, median_SternAcc7))

# Start the plot
accuracy_plot <- ggplot(dataSternALL_long_acc, aes(x = Condition, y = Accuracy)) +
  geom_quasirandom(aes(color = Condition), groupOnX = FALSE, size = 2) +
  scale_color_manual(values = c("blue", "green", "red")) +
  labs(x = "WM load", y = "Accuracy [%]") +
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
  )

# Add a line for medians
accuracy_plot <- accuracy_plot +
  geom_line(data = medians_df, aes(x = Condition, y = Median), group = 1, colour = "black", size = 0.5)

# Define the significance levels based on the p-values
signif_levels <- c("***", "**", "*", "n.s.", "??")

# Add the significance annotations for Accuracy
#accuracy_plot <- accuracy_plot +
#  geom_signif(comparisons = list(c("1 [%]", "4 [%]")), annotation = signif_levels[4], y_position = 100, tip_length = 0.02) +
#  geom_signif(comparisons = list(c("1 [%]", "7 [%]")), annotation = signif_levels[4], y_position = 103, tip_length = 0.02) +
#  geom_signif(comparisons = list(c("4 [%]", "7 [%]")), annotation = signif_levels[4], y_position = 106, tip_length = 0.02)

# Add x-axis labels and median values for conditions 1, 4, and 7
accuracy_plot <- accuracy_plot +
  scale_x_discrete(labels = c("1", "4", "7")) 
#  geom_text(data = data.frame(x = c(1, 2, 3), 
 #                             y = c(median_SternAcc1 + 2, median_SternAcc4 + 2, median_SternAcc7 + 2), 
  #                            label = c(median_SternAcc1, median_SternAcc4, median_SternAcc7)),
   #         aes(x, y, label = sprintf("%.2f", label)), vjust = 2.50, hjust = -0.55)

# Add legend for point size with custom title and increased font size
accuracy_plot <- accuracy_plot +
  scale_shape_discrete(name = "WM load", labels = c("Condition 1", "Condition 2", "Condition 3", "Condition 4")) +
  theme(legend.text = element_text(size = 12)) +
  ylim(75, 100)

# Print the accuracy plot
print(accuracy_plot)

####################### RT WITH MEDIAN
# Custom function to calculate median and MAD
median_mad_fun <- function(x) {
  med <- median(x)
  mad <- mad(x, constant = 1) # constant adjusts the scaling to be consistent with the standard deviation for normal distributions
  return(c(y = med, ymin = med - mad, ymax = med + mad))
}

# Calculate medians for the conditions
median_SternRT1 <- median(dataSternALL_long_rt$ReactionTime[dataSternALL_long_rt$Condition == "1 [ms]"])
median_SternRT4 <- median(dataSternALL_long_rt$ReactionTime[dataSternALL_long_rt$Condition == "4 [ms]"])
median_SternRT7 <- median(dataSternALL_long_rt$ReactionTime[dataSternALL_long_rt$Condition == "7 [ms]"])

# Create a data frame for median values
medians_df_rt <- data.frame(Condition = factor(c("1 [ms]", "4 [ms]", "7 [ms]"), levels = c("1 [ms]", "4 [ms]", "7 [ms]")),
                            Median = c(median_SternRT1, median_SternRT4, median_SternRT7))

# Highest y-value for the vertical line to end
highest_y_value_rt <- max(dataSternALL_long_rt$ReactionTime, na.rm = TRUE)

# Start the plot with theme adjustments for axis lines
rt_plot <- ggplot(dataSternALL_long_rt, aes(x = Condition, y = ReactionTime)) +
  geom_quasirandom(aes(color = Condition), groupOnX = FALSE, size = 2) +
  scale_color_manual(values = c("blue", "green", "red")) +
  labs(x = "WM load", y = "Reaction Time [ms]") +
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
    axis.line.y = element_line(color = "black", size = 0.5)
    # y-axis line is added as a segment below to end at the highest y-value
  )

# Add a line for medians
rt_plot <- rt_plot +
  geom_line(data = medians_df_rt, aes(x = Condition, y = Median), group = 1, colour = "black", size = 0.5)

# Define the significance levels based on BH-adjusted p-values
signif_levels_rt <- c("***", "**", "*", "n.s.") 

# Add the significance annotations for Reaction Times
rt_plot <- rt_plot +
  geom_signif(comparisons = list(c("1 [ms]", "4 [ms]")), annotation = signif_levels_rt[2], y_position = 575, tip_length = 0.02) +
  geom_signif(comparisons = list(c("1 [ms]", "7 [ms]")), annotation = signif_levels_rt[2], y_position = 595, tip_length = 0.02) 

# Add x-axis labels and median values for conditions 1, 4, and 7
rt_plot <- rt_plot +
  scale_x_discrete(labels = c("1", "4", "7")) 
  #geom_text(data = data.frame(x = c(1, 2, 3), 
  #                            y = c(median_SternRT1, median_SternRT4, median_SternRT7) + 20, # adjust the y position as needed
  #                           label = c(median_SternRT1, median_SternRT4, median_SternRT7)),
  #         aes(x, y, label = sprintf("%.2f", label)), vjust = 3.5, hjust = -0.75)

# Add legend for point size with custom title and increased font size
rt_plot <- rt_plot +
  scale_shape_discrete(name = "WM load", labels = c("Condition 1", "Condition 2", "Condition 3", "Condition 4")) +
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
p_values_acc <- pairwise_acc_results$p
formatted_p_values_acc <- format.pval(p_values_acc, method = "holm", digits = 3)
p_values_acc

# Extract p-values for Reaction Times
p_values_rt <- pairwise_rt_results$p
formatted_p_values_rt <- format.pval(p_values_rt, method = "holm", digits = 3)
p_values_rt

############# CIs
# Constants
mad_to_sd_factor <- 1.4826
z_value <- 1.96  # For 95% CI

# Sample sizes for each condition (replace these with your actual sample sizes)
n_1 <- 7  # Example sample size for condition 1
n_4 <- 7  # Example sample size for condition 4
n_7 <- 7  # Example sample size for condition 7

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
ci_acc_1 <- calculate_ci(98.4848485, 1.5151515, n_1)
ci_acc_4 <- calculate_ci(97.9166667, 0.6222944, n_4)
ci_acc_7 <- calculate_ci(92.5132275, 5.3373016, n_7)

medians_mad_rt

# Calculate CIs for Reaction Times
ci_rt_1 <- calculate_ci(302.93750, 16.51540, n_1)
ci_rt_4 <- calculate_ci(364.75000, 51.10880, n_4)
ci_rt_7 <- calculate_ci(449.68871, 41.24489, n_7)

# Print the CIs
print(ci_acc_1)
print(ci_acc_4)
print(ci_acc_7)

print(ci_rt_1)
print(ci_rt_4)
print(ci_rt_7)

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
