#STATISTICAL METHODS FOR BIOINFORMATICS - PROJECT: CALORIES BURNING
#Tutorial Group 4
#Alejandro Correa Rojo - r0731524
#Alicia Hernández Giménez - r0734014
#Ana Nogal Macho - r0736547
#Nahdah Sholiha - r0726859
#Felix Anto Ophelia Eusebious Rajan - r0732623
#Patrycja Owczarek - r0737157

#
# PACKAGES
#
library(mice)
library(VIM)
library(Amelia)
library(lattice)
library(ipw)
library(ggplot2)
library(BaylorEdPsych)

#
# DESCRIPTIVE STATISTICS AND EXPLORATORY DATA ANALYSIS        
#

# Dataset Calories Burning - calories, calhour and weight
dataset <- read.table("muscle-incomplete.txt", header = TRUE)
head(dataset)
summary(dataset)
# weight and calhour are expressed as numeric but
# represent categorical variables since they take
# discrete values in the dataset. calories, as
# response variable is continuous.
dataset$weight <- as.factor(dataset$weight)
dataset$calhour <- as.factor(dataset$calhour)
summary(dataset)

# Plots - Density, boxplot
densityplot(dataset$calories, xlab="Heat production [calories]", type = "density", ylab = "Density")
ggplot(dataset, aes(x=factor(weight), y=calories)) + geom_boxplot(colour="blue", fill="white") + xlab("Weight [kgs]") + ylab("Heat production [calories]")
ggplot(dataset, aes(x=factor(calhour), y=calories)) + geom_boxplot(colour="red", fill="white") + xlab("Work levels [calories/hour]") + ylab("Heat production [calories]") + ylim(210,360)
ggplot(dataset, aes(x=calories)) + geom_histogram(aes(y = ..density..),colour = "black", fill = "white", bins= 10) + geom_density(alpha=.2, fill="#FF6666") + xlab("Heat production [calories]") + ylab("Density")

# Exploring the missing values
md.pattern(dataset)
dataset.aggr=aggr(dataset,numbers=TRUE, prop=TRUE, 
                  sortVars = TRUE, ylab=c("Proportion of missing values" ,"Pattern"))
marginplot(dataset[,c("calories", "weight")], col = mdc(1:2), cex.numbers = 1.2, pch = 20, xlab = "Heat production [calories]", ylab = "Weight [kgs]")
marginplot(dataset[,c("calories", "calhour")], cex.numbers = 1.2, pch = 20, xlab = "Heat production [calories]", ylab = "Work levels [calories/hour]")
barMiss(dataset[,c("calhour", "calories")])
barMiss(dataset[,c("weight", "calories")])

# From the exploring data analysis, we see that there are 8 missing values in the dataset.
# there is a pattern in the missing values, only for the calhour equal to 13 and 19
# based on this, we are not dealing with MNAR missing data. We run a hypothesis test to 
# see if it is MAR or MCAR using LittleMCAR by BaylorEdPsych

LittleMCAR(dataset) # Ho: MCAR vs H1: MAR -> P < 0.05 reject Ho

#
# COMPLETE CASE ANALYSIS                        
#
data_cc <- lm(calories ~ weight + calhour, data = dataset)
summary(data_cc)

#
# MULTIPLE IMPUTATION ANALYSIS - mice                  
#

# Imputation methods applied: "pmm", "norm", "sample"
data_pmm <- mice(data = dataset, m = 10 , printFlag = F, method = "pmm")
data_norm <- mice(data = dataset, m = 10, printFlag = F, method = "norm")
data_sample <- mice(data = dataset, m = 10 , printFlag = F, method = "sample")

# Diagnostic plots
densityplot(data_pmm, xlab = "Heat production [calories]", ylab = "Density [Predictive mean matching]")
densityplot(data_norm, xlab = "Heat production [calories]", ylab = "Density [Normal model]")
densityplot(data_sample,  xlab = "Heat production [calories]", ylab = "Density [Random sample]")

stripplot(data_pmm, pch = 20, cex = 1.8, xlab = "Imputation number", ylab = "Heat production [calories] - pmm")
stripplot(data_norm, pch = 20, cex = 1.8, xlab = "Imputation number", ylab = "Heat production [calories] - norm")
stripplot(data_sample, pch = 20, cex = 1.8, xlab = "Imputation number", ylab = "Heat production [calories] - sample")

xyplot(data_pmm, calories ~ calhour + weight , pch = 20, cex = 1.4, main = "Method: 'pmm'")
xyplot(data_norm, calories ~ calhour + weight , pch = 20, cex = 1.4, main = "Method: 'norm'")
xyplot(data_sample, calories ~ calhour + weight , pch = 20, cex = 1.4, main = "Method: 'sample'")

# Pooling
model_pmm <- with(data_pmm, lm(calories ~ calhour + weight))
summary(pool(model_pmm))
pool.r.squared(model_pmm)

model_norm <- with(data_norm, lm(calories ~ calhour + weight))
summary(pool(model_norm))
pool.r.squared(model_norm)

model_sample <- with(data_sample, lm(calories ~ calhour + weight))
summary(pool(model_sample))
pool.r.squared(model_sample)

#
# INVERSE PROBABILITY WEIGHTS                        
#

# Dummy variable r for calories
r <- rep(1,24)
data.ipw <- data.frame(dataset, r)
for (i in 1:nrow(data.ipw)){
  if (is.na(data.ipw$calories[i]) == TRUE) {
    data.ipw$r[i] <- 0
  }
}

# Complete case analysis - IPW
data.ipw.cc <- lm(r ~ weight + calhour, data = data.ipw) 
summary(data.ipw.cc)

# Weights
data.ipw$w <- 1/fitted(data.ipw.cc)

# Complete model
data.ipw.results <- lm(calories ~ weight + calhour, 
                        data = data.ipw,weights = data.ipw$w)
summary(data.ipw.results)
