#Project - Statistic Methods for Bioinformatics
data <- read.table("muscle-incomplete.txt", header = TRUE)
#Packages
library(skimr)
library(ggplot2)
library(plotly)
library(mice)
library(VIM)
library(Amelia)
library(Zelig)
library(ipw)

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
ggplot(data=data,aes(data$weight)) + geom_histogram(aes(y=..density..), col="black",fill="white",binwidth=5)+
  geom_density(col=1, alpha=0.2, fill="#FF6666")+labs(title="Histogram for Weight") +
  labs(x="Weight", y="Count")

#Calhour
ggplot(data=data,aes(data$calhour)) + geom_histogram(aes(y=..density..), col="black",fill="white",binwidth=5)+
  geom_density(col=1, alpha=0.2, fill="#FF6666")+labs(title="Histogram for Calhour") +
  labs(x="Calhour", y="Count")
#Calories
ggplot(data=data,aes(data$calories)) + geom_histogram(aes(y=..density..), col="black",fill="white",binwidth=25)+
  geom_density(col=1, alpha=0.2, fill="#FF6666")+labs(title="Histogram for Calories") +
  labs(x="Calories", y="Count")

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

###################################################################
#             LINEAR REGRESSION - COMPLETE CASE                   #
###################################################################
#Fitting a linear regression model for the complete case
data_cc <- lm(calories ~ calhour + weight, data = data)
summary(data_cc)
anova(data_cc)

###################################################################
#           MULTIPLE IMPUTATION ANALYSIS - mice                   #
###################################################################
data_mice <- mice(data = data, m = 100, method = "pmm")
data_mice
data_mice$imp$calories
#Imputation diagnostics
densityplot(data_mice)
com <- complete(data_mice, "long", inc = T)
stripplot(calories~.imp, data = com, pch = 20, cex = 1.2)
xyplot(data_mice, calories ~ calhour + weight , pch = 20, cex = 1.2)
#Pooling
model_mice <- with(data_mice, lm(calories ~ calhour + weight, data))
summary(pool(model_mice))
pool.r.squared(model_mice)
###################################################################
#           MULTIPLE IMPUTATION ANALYSIS - VIM                    #
###################################################################
data_irmi <- irmi(data, eps = 5, maxit = 50)
vars <- c("weight", "calories", "calories_imp")
pbox(data_irmi[,vars], delimiter = "imp", alpha = 0.5)
model_irmi <- lm(calories ~ calhour + weight, data = data_irmi)
summary(model_irmi)

data_knn <- kNN(data, variable = colnames(data))
marginplot(data_knn[,vars], delimiter = "imp")
model_knn <- lm(calories ~ calhour + weight, data = data_knn)
summary(model_knn)

###################################################################
#         MULTIPLE IMPUTATION ANALYSIS - AMELIA II                #
###################################################################
missmap(data)
complete_data <- amelia(data, m = 50)
write.amelia(obj = complete_data, file.stem = "out_amelia")
model_amelia <- lm(calories ~ calhour + weight, data = complete_data$imputations$imp50)
summary(model_amelia)
anova(model_amelia)
compare.density(complete_data, var = "calories")
overimpute(complete_data, var = "calories")
###################################################################
#              INVERSE PROBABILITY WEIGHTS                        #
###################################################################
level <- rep(1,24)
data_prob <- data.frame(data, level)
for (i in 1:nrow(data_prob)){
  if (is.na(data_prob$calories[i]) == TRUE) {
    data_prob$level[i] <- 0
  }
}
#Complete case
prob_comp <- glm(level ~ calhour + weight , data = data_prob, family = binomial) 
summary(prob_comp)  
#IPW
#Full Model
full <- ipwpoint(exposure = level, family = "binomial", numerator = ~ 1, denominator =
                   ~ calhour + weight, link = "logit", data = data_prob)
summary(temp$ipw.weights)
summary(temp$num.mod)
summary(temp$den.mod)

#Weight
weight <- ipwpoint(exposure = level, family = "binomial", numerator = ~ 1, denominator =
                   ~ weight, link = "logit", data = data_prob)
summary(weight$ipw.weights)
summary(weight$num.mod)
summary(weight$den.mod)

#Calhour
calhour <- ipwpoint(exposure = level, family = "binomial", numerator = ~ 1, denominator =
                   ~ calhour, link = "logit", data = data_prob)
summary(calhour$ipw.weights)
summary(calhour$num.mod)
summary(calhour$den.mod)

