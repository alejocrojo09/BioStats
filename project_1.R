#Project - Statistic Methods for Bioinformatics
data <- read.table("muscle-incomplete.txt", header = TRUE)
#Packages
library(skimr)
library(ggplot2)
library(plotly)
library(mice)
library(VIM)
library(pastecs)
library(Zelig)
library(ipw)
library(ggpubr)
library(BaylorEdPsych)
library(mvnmle)

###################################################################
#     DESCRIPTIVE STATISTICS AND EXPLORATORY DATA ANALYSIS        #
###################################################################
#Summary - Using Summary function
summary(data)
#Summary - Using skimr function
skim(data)
#descriptive stats from ch.1
desc.data <- stat.desc(data[,c("weight", "calhour", "calories")], basic=TRUE, desc=TRUE)
desc.data

#Normality of data
weight <- data$weight
calhour <- data$calhour
calories <- data$calories #Calories: There are 8 missing values

shapiro.test(weight) #P-Value: 0.00833 -> Data not normally distributed
shapiro.test(calhour) #P-Value: 0.0339 -> Data not normally distributed
shapiro.test(calories) #P-Value: 0.172 -> Check for the influence of missing values

combine<- data.frame(weight, calhour, calories)
pairs(combine)

#correlation with spearman when not normal
cor(combine, method="spearman")
cor.test(calories, calhour, method="spearman") #seems high correlation
cor.test(calories, weight, method="spearman") #low correlation
cor.test(weight,calhour, method="spearman") #low correlation

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

#histogram pointing the mean
gghistogram(data, x = "calories", bins = 9, 
            add = "mean")
gghistogram(data, x = "weight", bins = 9, 
            add = "mean")
gghistogram(data, x = "calhour", bins = 9, 
            add = "mean")

#boxplots
par(mfrow = c(1,3))

boxplot(calories, col = "blue", xlab="calories")
boxplot(weight, col = "red", xlab="weight")
boxplot(calhour, col = "pink", xlab="calhour")


#weight seems to be like in categories
#so make boxplots based on categories
weight2 <- as.factor(weight)
summary(data)
data.kat <- data.frame(calories, weight2, calhour)
View(data.kat)
summary(data.kat)
ggplot(data.kat, aes(x=factor(weight2), y=calories)) + geom_boxplot(colour="green", fill="yellow")



#Density plots
par(mfrow = c(1,2))
plot(density(calories), xlab="calories")
plot(density(weight2), xlab="weight")
plot(density(weight), xlab="weight")
plot(density(calhour), xlab="calhour")

#index plots
par(mfrow = c(1,3))
plot(sort(calories), ylab = "", main = "calories", pch =20)
plot(sort(weight), ylab = "", main = "weight", pch =20)
plot(sort(calhour), ylab = "", main = "calhour", pch =20)
plot(sort(weight2), ylab = "", main = "weight", pch =20)


#Missing values
md.pattern(data) #Missing values only in calories data
aggr_plot <- aggr(data, col=c('navyblue','red'), 
                  numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern")) #Missing values in calories -> 33.3%

#combinations
aggr(data, combined=TRUE, numbers = TRUE, prop = TRUE, cex.numbers=0.87, varheight = FALSE)

#amount of missingness of calories in each group?
par(mfrow = c(1,2))
barMiss(data[,c("calhour","calories")])
barMiss(data[,c("weight","calories")])
histMiss(data)

#use tidyverse to filter the NA
themiss <- filter(data, is.na(calories))
View(themiss)
#so turns out those with calhour 13 and 19 miss their calories data


###################################################################
#             LINEAR REGRESSION - COMPLETE CASE                   #
###################################################################
#Fitting a linear regression model for the complete case
data_cc <- lm(calories ~ calhour + weight, data = data)
summary(data_cc)
anova(data_cc)
#second model with interaction
data.omit2 <- lm(calories~weight+calhour+I(weight*calhour))
summary(data.omit2)
#so its better with interaction!!
data.AIC <- step(lm(calories~1), scope=~weight+calhour+weight*calhour, direction = "forward")
summary(data.AIC)
#AIC agrees

#checking underlying assumptions
fit <- fitted(data.omit2)
rs <- rstandard(data.omit2)
plot(rs~fit)
#no pattern = linear

###################################################################
#             CHECK IF MISSING AT RANDOM                          #
###################################################################
LittleMCAR(data)
#Ho = mcar
#H1 = MAR
#p < 0.05 then reject Ho then it is MAR

###################################################################
#           MULTIPLE IMPUTATION ANALYSIS - mice                   #
###################################################################
#patterns on missingness
md.pattern(data)
#16 samples are complete
#8 samples missing in calories
pairs<-md.pairs(data)
marginplot(data[c(3,1)]) #weight
marginplot(data[c(3,2)]) #calhour
pairs
#can be seen it is at random (MAR) for calhour

#impute with mice
data_mice <- mice(data = data, m = 100, method = "pmm")
data_mice
data_mice$imp$calories

#Imputation diagnostics
densityplot(data_mice)
com <- complete(data_mice, "long", inc = T)
col <- rep(c("blue", "red")[1+as.numeric(is.na(data_mice$data$calories))],101)
stripplot(calories~.imp, data = com, jit=TRUE,fac=0.8, col=col,pch=20,cex=1.4, xlab="Imputation number")
stripplot(data_mice, pch=20, cex=1.2)
xyplot(data_mice, calories ~ calhour + weight , pch = 20, cex = 1.2)
#blue=observed, red=NA

#Pooling/analysis of the imputed data
fit<- with(data=data_mice, exp=data.omit2)
est <- pool(fit)
summary(est)
summary(data.omit2)

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

