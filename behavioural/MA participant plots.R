# MA participant plots


# demographics

# Load the knitr package
library(knitr)

# Create a data frame with your data
data <- data.frame(
  Participants = 10,
  Mean_Age = 25.9,
  SD = 2.13,
  Gender = "50%"
)

# Create a table using kable
table <- kable(data, format = "markdown", col.names = c("Participants", "Mean Age (M)", "Standard Deviation (SD)", "Gender"))

# Print the table
cat(table)
