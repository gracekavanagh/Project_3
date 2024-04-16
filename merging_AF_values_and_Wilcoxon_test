# Script for merging AF values for hybrid population and pure populations and comparing them using a Wilcoxon test
# Required _AF.txt files e.g. OCH_AF.txt for all cohort populations except PIZ

#install package to be able to run normality test
install.packages("nortest")

# Read in all the tetraploids from the cohort

BZD <- read.delim("BZD_AF.txt",sep ="\t")
GYE <- read.delim("GYE_AF.txt",sep ="\t")
KEH <- read.delim("KEH_AF.txt",sep ="\t")
MOD <- read.delim("MOD_AF.txt",sep ="\t")
OCH <- read.delim("OCH_AF.txt",sep ="\t")
PEK <- read.delim("PEK_AF.txt",sep ="\t")

# Merge the AF files for the pure arenosa and lyrata populations
AF_all_pure <- merge(KEH, MOD, by =1, all = TRUE)

# Merge the AF files for the hybrid populations
AF_all_hybrid <- merge(BZD, GYE, by =1, all = TRUE)
AF_all_hybrid <- merge(AF_all_hybrid, OCH, by = 1, all = TRUE)
AF_all_hybrid <- merge(AF_all_hybrid, PEK, by = 1, all = TRUE)

# Select all columns except the first one for each combined file
AF_a <- AF_all_pure[, -1]   
AF_b <- AF_all_hybrid[, -1]  

# Remove NAs
AF_a <- na.omit(AF_a)
AF_b <- na.omit(AF_b)

# Get the column means
means_pure <- colMeans(AF_a)
means_hybrid <- colMeans(AF_b)


#########varience and normality############


# Test for Variance
variance_pure <- apply(AF_all_pure[, -1], 2, var)
variance_hybrid <- apply(AF_all_hybrid[, -1], 2, var)

# Print the variance of the combined data
print(variance_pure)
print(variance_hybrid)

# Test for Normality
# Load the nortest package for Anderson-Darling test


library(nortest)

# Function to perform Anderson-Darling test for each column
anderson_darling_test <- function(data) {
  result <- lapply(data, function(col) ad.test(col))
  return(result)
}

# Test for Normality using Anderson-Darling test for pure populations
anderson_darling_pure <- anderson_darling_test(AF_all_pure[, -1])

# Test for Normality using Anderson-Darling test for hybrid populations
anderson_darling_hybrid <- anderson_darling_test(AF_all_hybrid[, -1])

# Print the results
print(anderson_darling_pure)
print(anderson_darling_hybrid)




############t-test#################


# Perform Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(means_pure, means_hybrid)

# Print the test result
print(wilcox_test_result)


