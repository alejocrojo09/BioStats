# BioStats
Statistical Methods for Bioinformatics Project 2018-2019

## First things first!
Please download the following the following packages into your RStudio:
- mice -- Imputation analysis
- VIM -- Imputation analysis
- Amelia -- Imputation analysis
- ggplot2 -- Data visualization
- Plotly -- Data visualization
- skimr -- Descriptive statistics

When you make an update, **PLEASE COMMENT IT!** so everybody knows the changes

## Project workflow
The idea is to make pairs to answer the questions of the project.

### Imputation analysis
I propose to every pair work with the packages of imputation described before, since they are different in algorithms. First, do all the analysis to classify if the missing data are not at random (MNAR), random (MAR) or completely at random (MCAR). Then, since there is no categorical variables in our dataset do all the imputation using the models that does not take account factor levels and logistic regression, e.g., in mice: pmm,norm, norm.nob, mean and sample. Finally, do all the diagnostic checking and compare the models used, Which is better?

#### Pairs

- Alicia and Ana: mice
- Ophelia and Patrycja: VIM
- Alejandro and Nahdah: Amelia

### Inverse Probability

## Update
Here is the initial code with descriptive statistics of the dataset. Please check the results and feel free to change
it.
