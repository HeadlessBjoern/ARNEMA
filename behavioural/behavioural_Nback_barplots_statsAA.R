# MA Nback RT and ACC barplots

library(RColorBrewer)
library(readxl)
library(car) # Load the 'car' package

# Load data
file_path <- "/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_nback_ALL.xlsx" #ALL file means means over all participants
dataNback <- read_excel(file_path)
dataNbackALL <- read_excel("/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_nback.xlsx") #This file has all participants listed with means over conditions
NbAccALL <- c(dataNbackALL$`corrN1`, dataNbackALL$`corrN3`)
NbRTALL <- c(dataNbackALL$`RTN1`, dataNbackALL$`RTN3`)
NbAcc <- c(dataNback$`corrN1`, dataNback$`corrN3`)
NbRT <- c(dataNback$`RTN1`, dataNback$`RTN3`)

# Calculate SEM for Nback data
sem_NbAcc <- sd(NbAcc) / sqrt(length(NbAcc))
sem_NbRT <- sd(NbRT) / sqrt(length(NbRT))

# Nback Accuracy
NbackAccFigure <- barplot(height=NbAcc, names=c(1, 3), col=c("blue", "red"), 
                          xlab="WM load", ylab="Accuracy [%]", ylim = c(0, 110))
text_y <- NbAcc*0.025
text(NbackAccFigure, text_y, labels = round(NbAcc, 2), pos = 3, col = "black")

# Add error bars for Nback Accuracy
arrows(NbackAccFigure, NbAcc - sem_NbAcc, NbackAccFigure, NbAcc + sem_NbAcc, angle=90, code=3, length=0.1, col="black")

# Nback Reaction Times
NbackRTFigure <- barplot(height=NbRT, names=c(1, 3), col=c("blue", "red"), 
                         xlab="WM load", ylab="Reaction Time [ms]", ylim = c(0, 500))
text_y <- NbRT*0.025
text(NbackRTFigure, text_y, labels = round(NbRT, 2), pos = 3, col = "black")

# Add error bars for Nback Reaction Times
arrows(NbackRTFigure, NbRT - sem_NbRT, NbackRTFigure, NbRT + sem_NbRT, angle=90, code=3, length=0.1, col="black")

############################# Histogram and Q-Q
# Create histograms
hist(NbackAcc1, breaks = 10, ylab = "Participants", xlab = "Accuracy [%]")
hist(NbackAcc3, breaks = 10)
hist(NbackRT1, breaks = 10)
hist(NbackRT3, breaks = 10)

# Create Q-Q plots
qqnorm(NbackAcc1)
qqline(NbackAcc1)
qqnorm(NbackAcc3)
qqline(NbackAcc3)

############################# Shapiro-Wilk Test
# Shapiro-Wilk test for normality
shapiro.test(NbackAcc1)
shapiro.test(NbackAcc3)

############################# Boxplot
# Create boxplots
boxplot(NbackAcc1, NbackAcc3, names = c("1-back", "3-back"))

############################# T-Test

# Separate data into two conditions
NbackAcc1 <- dataNbackALL$`corrN1`
NbackAcc3 <- dataNbackALL$`corrN3`
NbackRT1 <- dataNbackALL$`RTN1`
NbackRT3 <- dataNbackALL$`RTN3`

# Calculate SEM for Sternberg data
sem_NbackAcc1 <- sd(NbackAcc1) / sqrt(length(NbackAcc1))
sem_NbackAcc3 <- sd(NbackAcc3) / sqrt(length(NbackAcc3))
sem_NbackRT1 <- sd(NbackRT1) / sqrt(length(NbackRT1))
sem_NbackRT3 <- sd(NbackRT3) / sqrt(length(NbackRT3))

