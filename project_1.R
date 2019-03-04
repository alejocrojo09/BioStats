#Project - Statistic Methods for Bioinformatics
data <- read.table("muscle-incomplete.txt", header = TRUE)
#Packages
library(skimr)
library(ggplot2)
library(plotly)
library(mice)
library(VIM)

###################################################################
#     DESCRIPTIVE STATISTICS AND EXPLORATORY DATA ANALYSIS        #
###################################################################
#Summary - Using Summary function
summary(data)
#Summary - Using skimr function
skim(data)
#Normality of data
weight <- data$weight
calhour <- data$calhour
calories <- data$calories #Calories: There are 8 missing values

shapiro.test(weight) #P-Value: 0.00833 -> Data not normally distributed
shapiro.test(calhour) #P-Value: 0.0339 -> Data not normally distributed
shapiro.test(calories) #P-Value: 0.172 -> Check for the influence of missing values

#Graphics
#Histograms
#Weight
w <- ggplot(data,aes(x=weight))
w + geom_histogram()
#Calhour
calh <- ggplot(data,aes(x=calhour))
calh + geom_histogram()
#Calories
cal <- ggplot(data,aes(x=calories))
cal + geom_histogram()

#Density plots
#Weight
w + geom_density()
#Calhour
calh + geom_density()
#Calories
cal + geom_density()

#Missing values
md.pattern(data) #Missing values only in calories data
aggr_plot <- aggr(data, col=c('navyblue','red'), 
                  numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern")) #Missing values in calories -> 33.3%

#complete case (CC) analysis
#the theory: 
#also called listwise deletion, the NA are deleted
#derive estimates from the complete/filled ones
#basically just analyse the ones that are completed

#fitting regression mode for CC
data <- data.frame(muscle.incomplete)
View(data)
data.omit <- lm(calories~weight+calhour, data=data)
summary(data.omit)
#can be seen 8 observations deleted due to missingness
#R deletes all incomplete observations by default
#intercept is not significant IDK why
data.omit2 <- lm(calories~calhour, data=data)
summary(data.omit2)
#if we omit weight, the intercept became significant 

#making correlation with complete case 
data_cc <- data[complete.cases(data),]
View(data_cc)
weight <- data_cc$weight
calhour <- data_cc$calhour
calories <- data_cc$calories
shapiro.test(weight)
shapiro.test(calhour)
shapiro.test(calories) #only calories normality is accepted
combine<- data.frame(weight, calhour, calories)
#cor with spearman when not normal
cor(combine, method="spearman")

