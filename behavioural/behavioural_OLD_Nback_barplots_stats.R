# MA Nback RT and ACC barplots

library(RColorBrewer)
library(readxl)
library(car) # Load the 'car' package

############################# Load data

# Load data
file_path <- "/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/behavioral_nback_ALL.xlsx" #ALL file means means over all participants
dataNback <- read_excel(file_path)
dataNbackALL <- read_excel("/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/behavioral_nback.xlsx") #This file has all participants listed with means over conditions
NbAccALL <- c(dataNbackALL$`corrN1`, dataNbackALL$`corrN3`)
NbRTALL <- c(dataNbackALL$`RTN1`, dataNbackALL$`RTN3`)
NbAcc <- c(dataNback$`corrN1`, dataNback$`corrN3`)
NbRT <- c(dataNback$`RTN1`, dataNback$`RTN3`)

############################# Separate data into two conditions

# Separate data into two conditions
NbackAcc1 <- dataNbackALL$`corrN1`
NbackAcc3 <- dataNbackALL$`corrN3`
NbackRT1 <- dataNbackALL$`RTN1`
NbackRT3 <- dataNbackALL$`RTN3`

############################# Calculate Standard Error of the Mean (SEM)

# Calculate SEM
sem_NbackAcc1 <- sd(NbackAcc1) / sqrt(length(NbackAcc1))
sem_NbackAcc3 <- sd(NbackAcc3) / sqrt(length(NbackAcc3))
sem_NbackRT1 <- sd(NbackRT1) / sqrt(length(NbackRT1))
sem_NbackRT3 <- sd(NbackRT3) / sqrt(length(NbackRT3))

############################# Histogram and Q-Q
# Create histograms
hist(NbackAcc1, breaks = 10, ylab = "Participants", xlab = "Accuracy [%]")
hist(NbackAcc3, breaks = 10)
hist(NbackRT1, breaks = 10)
hist(NbackRT3, breaks = 10)

############################# Wilcoxon test

# Perform Wilcoxon signed-rank test for Accuracy and RT
ttest_result_acc <- t.test(NbackAcc1, NbackAcc3, paired = TRUE)
wilcox_result_rt <- wilcox.test(NbackRT1, NbackRT3, paired = TRUE)

############################# N-back Accuracy

# Create the barplot
NbAccFigure <- barplot(height = c(mean(NbackAcc1), mean(NbackAcc3)), 
                       names = c(1, 2), col = c("blue", "red"), 
                       xlab = "WM load", ylab = "Accuracy [%]", ylim = c(0, 110))

# Add error bars
arrows(NbAccFigure, c(mean(NbackAcc1) - sem_NbackAcc1, mean(NbackAcc3) - sem_NbackAcc3), 
       NbAccFigure, c(mean(NbackAcc1) + sem_NbackAcc1, mean(NbackAcc3) + sem_NbackAcc3), 
       angle = 90, code = 3, length = 0.1, col = "black")

# Add max values per plot
text_y <- NbAcc * 0.025
text(NbAccFigure, text_y, labels = round(NbAcc, 2), pos = 3, col = "black", cex = 1.1)

############################# N-back Reaction Times

NbRTFigure <- barplot(height=c(mean(NbackRT1), mean(NbackRT3)), 
                      names=c(1, 2), col=c("blue", "red"), 
                      xlab="WM load", ylab="Reaction Time [ms]", ylim = c(0, 500))

# Add error bars
arrows(NbRTFigure, c(mean(NbackRT1) - sem_NbackRT1, mean(NbackRT3) - sem_NbackRT3), 
       NbRTFigure, c(mean(NbackRT1) + sem_NbackRT1, mean(NbackRT3) + sem_NbackRT3), 
       angle=90, code=3, length=0.1, col="black")

# Add max values per plot
text_y <- NbAcc * 0.025
text(NbRTFigure, text_y, labels = round(NbRT, 2), pos = 3, col = "black", cex = 1.1)


# Significance levels
alpha <- 0.05

# Function to determine the number of asterisks based on p-value
get_asterisks <- function(p_value, alpha) {
  if (p_value < alpha) {
    num_asterisks <- sum(p_value < c(0.001, 0.01, 0.05))
    return(paste(rep("*", num_asterisks), collapse = ""))
  } else {
    return("")
  }
}

############################# Check significance levels

# Significance levels
alpha <- 0.05

# Function to determine the number of asterisks based on p-value
get_asterisks <- function(p_value, alpha) {
  if (p_value < alpha) {
    num_asterisks <- sum(p_value < c(0.001, 0.01, 0.05))
    return(paste(rep("*", num_asterisks), collapse = ""))
  } else {
    return("")
  }
}

# Check if the Wilcoxon test for Accuracy is significant and add asterisks
stars_acc <- get_asterisks(ttest_result_acc$p.value, alpha)

# Check if the Wilcoxon test for RT is significant and add asterisks
stars_rt <- get_asterisks(wilcox_result_rt$p.value, alpha)

# Print the results with asterisks
cat("Wilcoxon test results for Accuracy:", stars_acc, "\n")
cat("Wilcoxon test results for RT:", stars_rt, "\n")

