################################################################################
# 00_data_curation.R  —  Psychosocial exposome variable curation (BHRCS)
#
# PURPOSE (documentation of the manual + programmatic curation that produced the
# curated exposome dataset used in this study). This file is the FAITHFUL curation
# record extracted from the original preprocessing pipeline (setup, data loading,
# harmonization/de-duplication, recoding, composite construction, and a-priori
# exclusions). It documents how raw BHRCS wave-0 variables were reduced to the
# curated candidate exposome (~118 variables; data dictionary in Supplementary
# Table S1). The final n=70 analytic exposures result from the programmatic
# filters applied in 01_exwas_enet_cPES.R (CBCL-overlap removal, hard-coded
# exclusion of raw items superseded by recodes, near-zero-variance removal, and
# complete-case restriction on the outcome and covariates).
#
# NOTE: the underlying BHRCS individual-level data are not publicly shareable
# (controlled-access cohort); this script documents the analysis, it is not a
# one-click reproduction. See README.md and the Methods/Supplement.
################################################################################

################################################################################
###############################  CLINICAL EXPOSOME #############################
###############################  BRAZILIAN HIGH-RISK ###########################
################################################################################
# Packages                                                                   ----------------
# importing and working with the dataset
devtools::install_github("strengejacke/strengejacke")
library(strengejacke)
library(openxlsx)
library(car)
library(dplyr)
library(MASS)
library(huge)
library(foreign)
library(sjPlot)
library(foreign)
suppressPackageStartupMessages({
  library("SHAPforxgboost"); library("ggplot2"); library("xgboost")
  library("data.table"); library("here")
})
#library(cvms) 
library(groupdata2) 
library(ebal)
library(lsr)
library(patchwork)
library(glmnet)
library(ggpubr)
library(SHAPforxgboost)
library(boot)
library(ggcorrplot)
#library(dominanceanalysis)
library(magrittr)
library(expss)
library(tableone)
library(tableone)
library(dplyr) 
library(igraph)
library(tcltk)
#library(rgl)
library(ape)
library(groupdata2)
#library(reprtree)
library(MatchIt)
#library(xpectr)
library(hydroGOF) # rmse()
library (e1071)
library(xgboost)
library(adabag)
library(gbm)
library(caTools)
library(class)
library(randomForest)
library(rpart)
library(rpart.plot)
library(gdata)
library(rpart)
library(ranger)
library(cluster)
library(psy)
library(psych)
library(factoextra)
library(dendextend)
library(ggcorrplot)
library(performance)
library(EFAtools)
library(corpcor)
library(EFAtools)
library(GPArotation)
#library(ggradar)
library(gdata)
library(foreign)
library(PMCMRplus)
library(readxl)
library(reshape)
library(reshape2)
library(magrittr)
# filter
library(plyr)
library(dplyr)
library(tidyr)
# gather: wide to ling
# spread: long to wide
# descriptive/exploratiry analysis
library(psych)
library(pastecs)
#install.packages("gdata")
#library(R.basic)
# plot
library(ggplot2)
library(lattice)
library(corrplot)
library(Hmisc)
library(plotly)
# library(gdata)
library(moments)
library(pastecs)
# plot
# library(easyGgplot2)
library(RColorBrewer)
#library(sjPlot)
# required packages for loading multilevel linear growth models
#library(car)
library(nlme)
library(lme4)
library(MASS)
library(QuantPsyc)
library(multilevel)
#library(sjstats)
library(MuMIn)
library(ppcor)
library(multcomp)
library(emmeans)
library(MuMIn)
library(PMCMR)
library(data.table)
library(mltools)
library(car)
library(corrplot) # We'll use corrplot later on in this example too.
library(visreg) # This library will allow us to show multivariate graphs.
library(knitr)
library(scatterplot3d)
library(stringr)
library(reshape2)
# library(ggplot2)
library(ggthemes)
library(pander)
library(ARTool)
library(matrixStats)
library(pROC) # AUC
library(haven)
library(fmsb)
library(Rmisc)
# SEM
library(lavaan)
library(lavaanPlot)
library(semPlot)
# High Density Scatterplot with Binning
library(hexbin) 
# This function draws a (grouped) boxplot-like plot with with kernel density estimators.
library(REdaS)
# plots 
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(RColorBrewer)
# Logistic regression
library(ResourceSelection)
library(rsq)
library(rcompanion)
library(psychometric)
library(ltm)
library(eRm)
library(mokken)
#library(KernSmoothIRT)
library(mirt)
library(gtheory)
library(MPsychoR)
library(gridExtra)
# network 
require(huge)
require(networktools)
require(qgraph)
require(parcor)
require(Matrix)
require(relaimpo)
require(igraph)
require(car)
require(ppcor)
require(bootnet)
require(IsingFit)
require(IsingSampler)
require(semPlot)
require(mlVAR)
require(graphicalVAR)
require(mgm)
require(NetworkComparisonTest)
require(corrplot)
require(digest)
require(devtools)
require(xlsx)
require(bnlearn)
require(Rgraphviz)
require(relaimpo)
require(ape)
require(matrixcalc)
library(Biobase)
library(Hmisc)
# tables
library(knitr)
library(kableExtra)
# BHRCS packages
library(rio)
library(tidyverse)
library(sjmisc)
library(sjlabelled)
library(semPlot)
library(semTools)
library(psych)
#install_formats()
################################################################################
# Directory                                                                  ------ 
setwd("<PROJECT_ROOT>")
# Setting seed                                                               ------
set.seed(241094)
# Cleanning environment                                                      -------------------
rm(list = ls())
graphics.off()
#xpectr::set_test_seed(1)

# Functions                                                                  ------
highleverage <- function(fit) {
  p <- length(coefficients(fit))
  n <- length(fitted(fit))
  ratio <-p/n
  plot(hatvalues(fit), main="Index Plot of Ratio")
  abline(h=c(2,3)*ratio, col="red", lty=2)
  identify(1:n, hatvalues(fit), names(hatvalues(fit)))
}


rows.in.a1.that.are.not.in.a2  <- function(a1,a2)
{
  a1.vec <- apply(a1, 1, paste, collapse = "")
  a2.vec <- apply(a2, 1, paste, collapse = "")
  a1.without.a2.rows <- a1[!a1.vec %in% a2.vec,]
  return(a1.without.a2.rows)
}
#rows.in.a1.that.are.not.in.a2(a1,a2)
spss2date <- function(x) as.Date(x/86400, origin = "1582-10-14")

# install.packages("rstatix")
library(rstatix)

rFromWilcox <- function (wilcoxModel, N){
  z<-qnorm(wilcoxModel$p.value/2)
  r <- z/sqrt(N)
  cat(wilcoxModel$data.name, "Effect Size, r = ", r)
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

plotPCA <- function(x, nGroup) {
  n <- ncol(x) 
  if(!(n %in% c(2,3))) { # check if 2d or 3d
    stop("x must have either 2 or 3 columns")
  }
  
  fit <- hclust(dist(x), method="complete") # cluster
  groups <- cutree(fit, k=nGroup)
  
  if(n == 3) { # 3d plot
    plot3d(x, col=groups, type="s", size=1, axes=F)
    axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=FALSE)
    grid3d("x")
    grid3d("y")
    grid3d("z")
  } else { # 2d plot
    maxes <- apply(abs(x), 2, max)
    rangeX <- c(-maxes[1], maxes[1])
    rangeY <- c(-maxes[2], maxes[2])
    plot(x, col=groups, pch=19, xlab=colnames(x)[1], ylab=colnames(x)[2], xlim=rangeX, ylim=rangeY)
    lines(c(0,0), rangeX*2)
    lines(rangeY*2, c(0,0))
  }
}

residual.stats<-function(matrix){
  residuals<-as.matrix(matrix[upper.tri(matrix)])
  large.resid<-abs(residuals) > 0.05
  numberLargeResids<-sum(large.resid)
  propLargeResid<-numberLargeResids/nrow(residuals)
  rmsr<-sqrt(mean(residuals^2))
  cat("Root means squared residual = ", rmsr, "\n")
  cat("Number of absolute residuals > 0.05 = ", numberLargeResids, "\n")
  cat("Proportion of absolute residuals > 0.05 = ", propLargeResid, "\n")
  hist(residuals)
}

factor.structure <- function(fa, cut = 0.2, decimals = 2){
  structure.matrix <- fa.sort(fa$loadings %*% fa$Phi)
  structure.matrix <- data.frame(ifelse(abs(structure.matrix) < cut, "", round(structure.matrix, decimals)))
  return(structure.matrix)
}

rotate_x <- function(data, column_to_plot, labels_vec, rot_angle) {
  plt <- barplot(data[[column_to_plot]], col='steelblue', xaxt="n")
  text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=0.6) 
}

get_box_stats <- function(y, upper_limit = max(df$mpg) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}
# Special commands from BHRCS -----
# get_label()  # to get variable labels;
# get_labels() # to get value labels;
# as_label()   # to converts a numeric variable into a factor and sets attribute-value-labels as factor levels.
################################################################################
# Data Harmonization                                                         -------
# Loading data                                                               ------
# OFFICE
w0_BHRCS <- read_rds("w0_selecionado_significantes (2).rds")
Data <- w0_BHRCS
BHRC_all <- read_sav("<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav", encoding = "latin1")

w0 <- subset(BHRC_all, BHRC_all$redcap_event_name == "wave0_arm_1")
  summary(w0$reg5)
  myvars <- c("ident", "reg5") # school address
  w0 <- w0[myvars]
  library("writexl")
  write_xlsx(w0,"ClinicalExposome_Data_BHRCS_schools.xlsx")
  rm (w0, myvars)
  
list <- c("ident",
          "reg10"	, #	Childs Skin color (race)
          "reg22"	, #	How does/did the child and the biological mother get along when they are/were together?
          "reg31a"	, #	If the biological mother and father are divorced, did the breakup occur in what terms?
          "reg37"	, #	How do/did (childs name) and the biological father get along when they are/were together?
          "reg38"	, #	In the moment, the biological father is in touch?
          "reg39"	, #	Address - Biological father
          "fat3"	, #	How did the biological mother felt when she found out about the pregnancy?
          "fat4"	, #	Did the parents feel suported during the childs first years?
          "fat14a"	, #	Did anyone often smoked close to the child during his/her childhood?
          "fat15"	, #	Did the biological mother drink any alcohol beveraged during pregnancy?
          "fat16c"	, #	Did the biological mother experienced two or more episodes of urinary tract infections?
          "fat17"	, #	Stress level experienced during pregnancy
          "fat20a_y"	, #	If yes, For how long did the child breastfed? - years
          "fat22"	, #	Has the child ever been bullied?
          "fat22a"	, #	If yes, has the child ever been bullied this year?
          "fat22c"	, #	How much did the bullying bother the child?
          "fat23"	, #	Did the child ever bullied?
          "fat27"	, #	Has the child ever been suspended from school?
          "fat28"	, #	Has the child ever been diagnosed with encephalitis, meningitis and had to be hospitalized for treatment?
          "fes3"	, #	You fight a lot in the family.
          "fes14"	, #	Sometimes your family gets so angry that they throw things
          "fes20r"	, #	Your family members usually loses their heads easily
          "fes26"	, #	Your family members criticize each other frequentily
          "fes32"	, #	Your family members sometimes beat each other
          "fes36"	, #	The group spirit in your family is very low
          "fes44"	, #	Your family members usually try to win in no matter what
          "fat35"	, #	Has the child ever happened to does not have anything to eat in the house and/or have to wear dirty or torn clothes, in a moment of need?
          "fat37"	, #	Has the child ever been cursed as an ass, idiot, stupid, by an adult?
          "cl2"	, #	Drink alcohol without parental permission.
          "cl39"	, #	Walk in bad company
          "cl55"	, #	Is fat, overweight
          "cl56a"	, #	Aches (except for head and stomach)
          "cl56b"	, #	Headaches
          "cl56c"	, #	Nausea, feeling sick
          "cl56d"	, #	Problems with the eyes (which are not corrected with the use of glasses)
          "cl56e"	, #	Skin Problems
          "cl56f"	, #	Stomach or belly pain
          "cl56g"	, #	Vomiting
          "cl76"	, #	Sleep less than most children
          "cl77"	, #	Sleep more than most children during the day and / or at night
          "cl99"	, #	Smokes cigarette or chews tobacco
          "cl100"	, #	Have sleep problems
          "cl105"	, #	Uses drugs (exclude alcohol and tobacco)
          "cl106"	, #	Commits acts of vandalism
          "p_fat39a"	, #	Does your father smile at you, speak in a gentle, caring voice, praise the things you do right, try to understand the problems, and help you solve them?
          "p_fat39b"	, #	Does your father give him freedom, letting him dress as he wants, making his own decisions?
          "p_fat39c"	, #	Does your father treat you like a baby, that is, is very overprotective, not wanting him to grow up?
          "np_cl37"	, #	Gets into many fights
          "np_cl55"	, #	Is fat, overweight
          "p1e1"	, #	PTSD: Exceptionally stressful event (Parent1)
          "p1e2e"	, #	PTSD: Physical abuse (Parent1)
          "p1e2f"	, #	PTSD: Sexual abuse (Parent1)
          "p1e2g"	, #	PTSD: Rape (Parent1)
          "fes_confl"	, #	FES: conflict score
          "pphyab"	, #	Has your child been seriously picked up from an adult (including yourself), to the point of leaving marks in his body?
          "pphyneg"	, #	Has it ever happened to your child not having something to eat and/or having to wear dirty or torn clothes?
          "psexab"	, #	Has somebody ever did sexual things with your child or have treatened your child if he didnt do sexual things with him?
          "pemoab"	, #	Has your child ever been cursed by some adult (including yourself), with words like ass, idiot, stupid, or being yelled that you were no good?
          "cphysab"	, #	Have you been seriously picked up from an adult to the point of leaving marks in your body?
          "cemoab"	, #	Have you ever been cursed by some adult, with words like ass, idiot, stupid, or being yelled that you were no good?
          "traumap_fl"	, #	Parent Rep: Exposure to Child Abuse/Neglect Factor Loading
          "traumac_fl"	, #	Child Rep: Exposure to Child Abuse/Neglect Factor Loading
          "trauma_fl"	, #	Exposure to Child Abuse/Neglect Factor Loading
          "traumap_exp"	, #	Parent Rep: Exposure Child Abuse/Neglect
          "traumac_exp"	, #	Child Rep: Exposure Child Abuse/Neglect
          "trauma_exp"	, #	Exposure to Child Abuse/Neglect
          "trauma_inc"	, #	Inc. Exposure to Child Abuse/Neglect
          "p_bulliedhim_l"	, #	Has the child ever been bullied in his life?
          "p_bulliedhim"	, #	Has the child ever been bullied this year?
          "p_bulliedhim_fr"	, #	If yes, how often has the child been bullied this year?
          "p_bulliedhim_imp"	, #	How much did the bullying bother the child?
          "c_bulliedhim_l"	, #	Have you ever been bullied in your life?
          "c_bulliedhim"	, #	Have you been bullied in the last 12 months?
          "c_bulliedhim_fr"	, #	If yes, how often do you suffered bullying in the last 12 months?
          "c_bulliedhim_imp"	, #	How much did it bothered you?
          "ghelp"	, #	Does anyone who lives with the child receive financial government help?
          "ghelp_income"	, #	What is the monthly value of this benefit? - Reais
          "religion_which"	, #	Which religion is predominant?
          "mrelation"	, #	How does/did the child and the biological mother get along when they are/were together?
          "mmarital"	, #	What is the marital status of the biological mother?
          "divorce"	, #	If the biological mother and father are divorced, did the breakup occur in what terms?
          "frelation"	, #	How do/did (childs name) and the biological father get along when they are/were together?
          "father_sameaddress"	, #	Address - Biological father
          "ffigure" #	Was there any Father figure, ie a WOMAN who cared for the child as a father would do?
)

# Creating Data frame                                                        -------
Data <- Data[,which(colnames(Data) %in% list)]
get_label(Data)  # to get variable labels;
get_labels(Data) # to get value labels;
as_label(Data)   # to converts a numeric variable into a factor and sets attribute-value-labels as factor levels.
rm (w0_BHRCS)
# Harmonizing the variables                                                  ------
# Race 
"reg10"	#	Childs Skin color (race)
hist (Data$reg10)
summary(Data$reg10)
table(Data$reg10)
Data$reg10 <- NULL

# Parents relationship 
"reg22"	#	How does/did the child and the biological mother get along when they are/were together? #DUPLICATE !
"mrelation"	 #	How does/did the child and the biological mother get along when they are/were together? #DUPLICATE !
table(Data$reg22)
table(Data$mrelation)
Data$reg22 <- NULL

summary(Data$mrelation)
table(Data$mrelation)
proporcao <- table(Data$mrelation, useNA = "no")
Data$mrelation[is.na(Data$mrelation)] <- sample(c(1, 2, 3, 4, 5), 
                                                    sum(is.na(Data$mrelation)), 
                                                    replace = TRUE, 
                                                    prob = proporcao / sum(proporcao))
table(Data$mrelation)
summary(Data$mrelation)
Data <- Data %>%
  mutate(mrelation_VeryWell = ifelse(mrelation == 1, 1, 0),
         mrelation_Well = ifelse(mrelation == 2, 1, 0),
         mrelation_Bad = ifelse(mrelation == 4 | mrelation == 5, 1, 0))

table(Data$mrelation_Bad)
Data$mrelation <- NULL

"reg37"	#	How do/did (childs name) and the biological father get along when they are/were together? #DUPLICATE !
"frelation"	#	How do/did (childs name) and the biological father get along when they are/were together? #DUPLICATE !

table(Data$reg37)
table(Data$frelation)

Data$frelation <- NULL

Data <- Data %>%
  mutate(frelation_VeryWell = ifelse(reg37 == 1, 1, 0),
         frelation_Well = ifelse(reg37 == 2, 1, 0),
         frelation_Bad = ifelse(reg37 == 3 |reg37 == 4 | reg37 == 5 | reg37 == 99, 1, 0))

table(Data$frelation_Bad)
table(Data$reg37)
Data$reg37 <- NULL

"mmarital"	 #	What is the marital status of the biological mother?
summary(Data$mmarital)
table(Data$mmarital)
proporcao <- table(Data$mmarital, useNA = "no")
Data$mmarital[is.na(Data$mmarital)] <- sample(c(1, 2, 3, 4, 5), 
                                                      sum(is.na(Data$mmarital)), 
                                                      replace = TRUE, 
                                                      prob = proporcao / sum(proporcao))
summary(Data$mmarital)
Data <- Data %>%
  mutate(mmarital_MarriedWithBioFather = ifelse(mmarital == 1, 1, 0),
         mmarital_MarriedWithOtherPartner = ifelse(mmarital == 2, 1, 0),
         mmarital_Single = ifelse(mmarital == 3, 1, 0),
         mmarital_SeparatedDivorced = ifelse(mmarital == 4, 1, 0),
         mmarital_Widow = ifelse(mmarital == 5, 1, 0))
table(Data$mmarital)
table(Data$mmarital_MarriedWithBioFather)
Data$mmarital <- NULL

"reg31a"	 #	If the biological mother and father are divorced, did the breakup occur in what terms? #DUPLICATE !
"divorce"	 #	If the biological mother and father are divorced, did the breakup occur in what terms? #DUPLICATE !
table(Data$reg31a)
table(Data$divorce)
Data$reg31a <- NULL

summary(Data$divorce)
Data$divorce[is.na(Data$divorce)] <- 0
table (Data$divorce)

Data <- Data %>%
  mutate(divorce_FriendlyNoConflicts = ifelse(divorce == 1, 1, 0),
         divorce_FriendlyFewConflicts = ifelse(divorce == 2, 1, 0),
         divorce_ConflictingNotWitnessed = ifelse(divorce == 3, 1, 0),
         divorce_VeryConflictualWitnessed = ifelse(divorce == 4, 1, 0), 
         divorce_dich = ifelse(divorce == 0, 0, 1))

table (Data$divorce_FriendlyNoConflicts)
table (Data$divorce_dich)
Data$divorce <- NULL

"father_sameaddress"	 #	Address - Biological father #DUPLICATE !
"reg39"	 #	Address - Biological father #DUPLICATE !
Data$father_sameaddress <- NULL
Data$reg39 <- NULL

"ffigure" #	Was there any Father figure, ie a WOMAN who cared for the child as a father would do?
"reg38"	 #	In the moment, the biological father is in touch?

table(Data$ffigure)
Data$ffigure <- NULL

table(Data$reg38)
Data$reg38 <- NULL

# Prenatal                                                                   --------
"fat3"	 #	How did the biological mother felt when she found out about the pregnancy? # UNIFY !!!
"fat17"	 #	Stress level experienced during pregnancy # UNIFY !!!
table (Data$fat3)
Data$fat3 <- ifelse(Data$fat3 == 99, NA, Data$fat3)
table (Data$fat3)
proporcao <- table(Data$fat3, useNA = "no")
Data$fat3[is.na(Data$fat3)] <- sample(c(1, 2, 3, 4, 5), 
                                              sum(is.na(Data$fat3)), 
                                              replace = TRUE, 
                                              prob = proporcao / sum(proporcao))


table (Data$fat3)
summary(Data$fat3)

Data <- Data %>%
  mutate(fat3_PositiveFeelings = ifelse(fat3 == 1, 1, 0),
         fat3_PassiveOrLightStress = ifelse(fat3 == 2, 1, 0),
         fat3_ModerateStress = ifelse(fat3 == 3, 1, 0),
         fat3_SevereStress = ifelse(fat3 == 4, 1, 0),
         fat3_AbortionConsidered = ifelse(fat3 == 5, 1, 0))
table (Data$fat3_AbortionConsidered)
Data$fat3 <- NULL

table (Data$fat17)
  
"fat16c"	 #	Did the biological mother experienced two or more episodes of urinary tract infections?
table (Data$fat16c)

Data <- Data %>%
  mutate(fat16c = ifelse(fat16c == 99, 0, fat16c))
table (Data$fat16c)

# Substance Use (Parents and child)
# Maybe unify the child usage
"fat14a"	 #	Did anyone often smoked close to the child during his/her childhood?
table (Data$fat14a)

Data <- Data %>%
  mutate(fat14a = ifelse(fat14a == 0 | fat14a == 99, 0, 1))
table (Data$fat14a)

"fat15"	 #	Did the biological mother drink any alcohol beveraged during pregnancy?
table (Data$fat15)
Data$fat15 <- ifelse(Data$fat15 == 99, 0, Data$fat15)
table (Data$fat15)
summary(Data$fat15)

Data <- Data %>%
  mutate(fat15_Drink = ifelse(fat15 != 0, 1, 0))
table (Data$fat15_Drink)
Data$fat15 <- NULL

"cl99"	#	Smokes cigarette or chews tobacco
table (Data$cl99)
summary(Data$cl99)

"cl2"	 #	Drink alcohol without parental permission
table (Data$cl2)
summary(Data$cl2)

"cl105"	 #	Uses drugs (exclude alcohol and tobacco)
table (Data$cl105)
summary(Data$cl105)

Data <- Data %>%
  mutate(AnySubstanceUse = ifelse(cl99 != 0 | cl2 != 0 | cl105 != 0, 1, 0))
table (Data$AnySubstanceUse)

Data$cl99 <- NULL
Data$cl2 <- NULL
Data$cl105 <- NULL
 
# Trauma                                                                     --------
# Early trauma = retrospective and better than prospective - we can use both
# NOTE: USE INDIVIDUAL VARIABLES OR FACTOR LOADING (I BELIEVE FACTOR IS BETTER)
"fat22"	 #	Has the child ever been bullied?  #DUPLICATE !
"p_bulliedhim_l"	 #	Has the child ever been bullied in his life? #DUPLICATE !
"c_bulliedhim_l"	#	Have you ever been bullied in your life?
table (Data$fat22)
table (Data$p_bulliedhim_l)
table (Data$c_bulliedhim_l)

Data$fat22 <- NULL

summary(Data$p_bulliedhim_l)
Data$p_bulliedhim_l[is.na(Data$p_bulliedhim_l)] <- 0

"fat22a"	 #	If yes, has the child ever been bullied this year? #DUPLICATE !
"p_bulliedhim"	 #	Has the child ever been bullied this year? #DUPLICATE !
"c_bulliedhim"	#	Have you been bullied in the last 12 months?
summary(Data$fat22a)
table (Data$fat22a)

summary(Data$p_bulliedhim)
table (Data$p_bulliedhim)

Data$fat22a <- NULL
Data$p_bulliedhim[is.na(Data$p_bulliedhim)] <- 0
summary(Data$p_bulliedhim)

summary(Data$c_bulliedhim)
table (Data$c_bulliedhim)
Data$c_bulliedhim[is.na(Data$c_bulliedhim)] <- 0
summary(Data$c_bulliedhim)

"fat22c"	 #	How much did the bullying bother the child? #DUPLICATE !
"p_bulliedhim_imp"	 #	How much did the bullying bother the child? #DUPLICATE !
"c_bulliedhim_imp"	 #	How much did it bothered you?
summary(Data$fat22c)
table (Data$fat22c)

summary(Data$p_bulliedhim_imp)
table (Data$p_bulliedhim_imp)

Data$fat22c <- NULL
Data$p_bulliedhim_imp[is.na(Data$p_bulliedhim_imp)] <- 0
table (Data$p_bulliedhim_imp)

"p_bulliedhim_fr"	#	If yes, how often has the child been bullied this year?
"c_bulliedhim_fr"	 #	If yes, how often do you suffered bullying in the last 12 months?
summary(Data$p_bulliedhim_fr)
table (Data$p_bulliedhim_fr)
Data$p_bulliedhim_fr[is.na(Data$p_bulliedhim_fr)] <- 0

summary(Data$c_bulliedhim_fr)
table (Data$c_bulliedhim_fr)
Data$c_bulliedhim_fr[is.na(Data$c_bulliedhim_fr)] <- 0 # maybe delete!!!

"fat35"	 #	Has the child ever happened to does not have anything to eat in the house and/or have to wear dirty or torn clothes, in a moment of need? 
"pphyneg"	 #	Has it ever happened to your child not having something to eat and/or having to wear dirty or torn clothes? 
summary(Data$fat35)
table (Data$fat35)

summary(Data$pphyneg)
table (Data$pphyneg)
Data$fat35 <- NULL

Data <- Data %>%
  mutate(pphyneg_exposed = ifelse(pphyneg %in% c(1, 2, 3), 1, 0))
table (Data$pphyneg_exposed)

Data$pphyneg <- NULL


"fat37"	 #	Has the child ever been cursed as an ass, idiot, stupid, by an adult? #DUPLICATE !
"pemoab"	 #	Has your child ever been cursed by some adult (including yourself), with words like ass, idiot, stupid, or being yelled that you were no good? #DUPLICATE !
"cemoab"	 #	Have you ever been cursed by some adult, with words like ass, idiot, stupid, or being yelled that you were no good?
table (Data$fat37)

table (Data$pemoab)

Data <- Data %>%
  mutate(pemoab_exposed = ifelse(pemoab %in% c(1, 2, 3), 1, 0))
table (Data$pemoab_exposed)
summary (Data$pemoab_exposed)
Data$pemoab <- NULL

table (Data$cemoab)

Data <- Data %>%
  mutate(cemoab_exposed = ifelse(cemoab %in% c(1, 2, 3), 1, 0))
table (Data$cemoab_exposed)
summary (Data$cemoab_exposed)
Data$cemoab <- NULL

Data$fat37 <- NULL

"pphyab"	#	Has your child been seriously picked up from an adult (including yourself), to the point of leaving marks in his body?
"cphysab" #	Have you been seriously picked up from an adult to the point of leaving marks in your body?
table (Data$pphyab)

Data <- Data %>%
  mutate(pphyab_exposed = ifelse(pphyab %in% c(1, 2, 3), 1, 0))
table (Data$pphyab_exposed)
Data$pphyab <- NULL

table (Data$cphysab)
Data <- Data %>%
  mutate(cphysab_exposed = ifelse(cphysab %in% c(1, 2, 3), 1, 0))
table (Data$cphysab_exposed)
Data$cphysab <- NULL

"psexab"	 #	Has somebody ever did sexual things with your child or have treatened your child if he didnt do sexual things with him?
table (Data$psexab)

Data <- Data %>%
  mutate(psexab_exposed = ifelse(psexab %in% c(1, 2, 3), 1, 0))
table (Data$psexab_exposed)
Data$psexab <- NULL

"traumap_fl"	 #	Parent Rep: Exposure to Child Abuse/Neglect Factor Loading
hist (Data$traumap_fl)
summary(Data$traumap_fl)

"traumac_fl"	 #	Child Rep: Exposure to Child Abuse/Neglect Factor Loading
hist (Data$traumac_fl)
summary(Data$traumac_fl)

"trauma_fl"	 #	Exposure to Child Abuse/Neglect Factor Loading
hist (Data$trauma_fl)
summary(Data$trauma_fl)
Data$trauma_fl <- NULL

"traumap_exp"	 #	Parent Rep: Exposure Child Abuse/Neglect 
hist (Data$traumap_exp)
summary(Data$traumap_exp)
Data$traumap_exp <- NULL

"traumac_exp"	 #	Child Rep: Exposure Child Abuse/Neglect
hist (Data$traumac_exp)
summary(Data$traumac_exp)
Data$traumac_exp <- NULL

"trauma_exp"	 #	Exposure to Child Abuse/Neglect
hist (Data$trauma_exp)
table (Data$trauma_exp)
summary(Data$trauma_exp)
Data$trauma_exp <- NULL

"trauma_inc"	 #	Inc. Exposure to Child Abuse/Neglect
hist (Data$trauma_inc)
table (Data$trauma_inc)
summary(Data$trauma_inc)
Data$trauma_inc <- NULL

# Diseases
"fat28"	#	Has the child ever been diagnosed with encephalitis, meningitis and had to be hospitalized for treatment?
hist (Data$fat28)
table (Data$fat28)
summary(Data$fat28)
Data$fat28 <- NULL

"cl56a"	 #	Aches (except for head and stomach)
Data$cl56a <- NULL

"cl56b"	 #	Headaches
Data$cl56b <- NULL

"cl56c"	 #	Nausea, feeling sick
Data$cl56c <- NULL

"cl56d"	 #	Problems with the eyes (which are not corrected with the use of glasses)
hist (Data$cl56d)
table (Data$cl56d)
get_labels(Data$cl56d) 
Data$cl56d <- NULL

"cl56e"	 #	Skin Problems
hist (Data$cl56e)
table (Data$cl56e)
get_labels(Data$cl56e)
Data$cl56e <- NULL

"cl56f"	#	Stomach or belly pain
hist (Data$cl56f)
table (Data$cl56f)
get_labels(Data$cl56f)
Data$cl56f <- NULL

"cl56g"	 #	Vomiting
hist (Data$cl56g)
table (Data$cl56g)
Data$cl56g <- NULL

"cl55"	 #	Is fat, overweight #DUPLICATE !
"np_cl55"	 #	Is fat, overweight #DUPLICATE !
hist (Data$cl55)
hist (Data$np_cl55)
table (Data$cl55)
table (Data$np_cl55)
get_labels(Data$cl55)
get_labels(Data$np_cl55)
summary(Data$cl55)

Data <- Data %>%
  mutate(Overweight = ifelse(cl55 == 0, 0, 1))
table (Data$Overweight)
Data$cl55 <- NULL
Data$np_cl55 <- NULL

# Sleep - exclude, since it is part of the any-disorder outcome 
"cl76"	 #	Sleep less than most children 
"cl77"	 #	Sleep more than most children during the day and / or at night 
"cl100"	 #	Have sleep problems 
table (Data$cl76)
table (Data$cl77)
table (Data$cl100)

Data$cl76 <- NULL
Data$cl77 <- NULL
Data$cl100 <- NULL

# $$$$
"ghelp"	 #	Does anyone who lives with the child receive financial government help?
"ghelp_income"	 #	What is the monthly value of this benefit? - Reais (EXCLUDE !!!)
hist (Data$ghelp)
table (Data$ghelp)
summary(Data$ghelp)

proporcao <- table(Data$ghelp, useNA = "no")
Data$ghelp[is.na(Data$ghelp)] <- sample(c(0, 1), 
                                                    sum(is.na(Data$ghelp)), 
                                                    replace = TRUE, 
                                                    prob = proporcao / sum(proporcao))
table(Data$ghelp)

hist (Data$ghelp_income)
summary(Data$ghelp_income)
Data$ghelp_income <- NULL

"religion_which"	 #	Which religion is predominant? 
hist (Data$religion_which)
table (Data$religion_which)
summary(Data$religion_which)
get_labels(Data$religion_which)

Data$religion_which <- NULL

# parental trauma history (0 or 1 with the vars below)
"p1e1"	 #	PTSD: Exceptionally stressful event (Parent1)
"p1e2e"	 #	PTSD: Physical abuse (Parent1)
"p1e2f"	 #	PTSD: Sexual abuse (Parent1)
"p1e2g"	 #	PTSD: Rape (Parent1)
hist (Data$p1e1)
table (Data$p1e1)

hist (Data$p1e2e)
table (Data$p1e2e)

hist (Data$p1e2f)
table (Data$p1e2f)

hist (Data$p1e2g)
table (Data$p1e2g)
Data$p1e1 <- NULL
Data$p1e2e <- NULL
Data$p1e2f <- NULL
Data$p1e2g <- NULL

# Externalizing behavior (EXCLUDE = not an exposure)
"cl106"	 #	Commits acts of vandalism (EXCLUDE !!!)
"np_cl37"	 #	Gets into many fights (EXCLUDE !!!)
"fat23"	 #	Did the child ever bullied? (EXCLUDE !!!)
"fat27"	 #	Has the child ever been suspended from school? (EXCLUDE !!!)
table (Data$cl106)
table (Data$np_cl37)
table (Data$fat23)
table (Data$fat27)

Data$cl106 <- NULL
Data$np_cl37 <- NULL
Data$fat23 <- NULL
Data$fat27 <- NULL

# Social Support 
"fat4"	 #	Did the parents feel suported during the childs first years?
"cl39"	 #	Walk in bad company
table (Data$fat4)
summary(Data$fat4)
get_labels(Data$fat4)

Data <- Data %>%
  mutate(fat4_SupportedWheneverNeeded = ifelse(fat4 == 1, 1, 0),
         fat4_SupportedMostTimes = ifelse(fat4 == 2, 1, 0),
         fat4_SupportedSomeMoments = ifelse(fat4 == 3, 1, 0),
         fat4_NoSupport = ifelse(fat4 == 4 | fat4 == 5 | fat4 == 99, 1, 0))  # Unifying the categories
table (Data$fat4_NoSupport)       
table (Data$fat4)
Data$fat4 <- NULL


table (Data$cl39)
get_labels(Data$cl39)

Data <- Data %>%
  mutate(BadCo = ifelse(cl39 == 1 | cl39 == 2, 1, 0))
table (Data$BadCo)
Data$cl39 <- NULL

# Postnatal ---------
# RESEARCH THE WHO-RECOMMENDED MINIMUM TIME (MINIMUM FOR HEALTH, THEN CONVERT TO YES/NO)
"fat20" # Did the child ever breastfed?
"fat20a_y"	#	If yes, For how long did the child breastfed? - years
Df0 <- subset(BHRC_all, BHRC_all$redcap_event_name == "wave0_arm_1")
table (Df0$fat20)
table (Df0$fat20a_m)
table (Data$fat20a_y)
summary(Data$fat20a_y)
get_labels(Data$fat20a_y)

myvars <- c("ident", "fat20", "fat20a_y", "fat20a_m", "fat20b_y", "fat20b_m", 
            "gender", "fat13", "fat10", "fat12", "fat15")
df <- Df0[myvars] # complete W0 dataset

df$fat20a_y[is.na(df$fat20a_y)] <- 0
summary(df$fat20a_y)

df$fat20_y2m <- df$fat20a_y * 12
summary(df$fat20_y2m)

df$fat20_monthsTotal <- rowSums(cbind(df$fat20_y2m, df$fat20a_m), na.rm = TRUE)

summary(df$fat20_monthsTotal)

summary(df$fat20b_y)
table(df$fat20b_y)
summary(df$fat20b_m)
table(df$fat20b_m)

df <- df %>%
  mutate(onlybreastfeeded = ifelse(!is.na(fat20b_y) | !is.na(fat20b_m), 1, 0))
table (df$onlybreastfeeded)

# PRENATAL
table(df$fat13)
summary(df$fat13)

df <- df %>%
  mutate(fat13_Smooking_Gest = ifelse(fat13 == 1 | fat13 == 2 | fat13 == 3 | fat13 == 4, 1, 0),
         fat13_Smooking_Daily_Gest = ifelse(fat13 == 4, 1, 0))  
table (df$fat13_Smooking_Gest) 


""
table(df$fat12)
summary(df$fat12)

df <- df %>%
  mutate(neonatal_UTI = ifelse(fat12 == 1, 1, 0))  
table (df$neonatal_UTI)

"fat15"
table(df$fat15)
summary(df$fat15)
df <- df %>%
  mutate(fat13_Drinking_Gest = ifelse(fat15 == 1 | fat15 == 2 | fat15 == 3 | fat15 == 4, 1, 0))  
table (df$fat13_Drinking_Gest)

# Selecting only the main variables ------

myvars <- c("ident", "fat20", "fat20_monthsTotal", "onlybreastfeeded", "gender", 
            "neonatal_UTI", "fat13_Drinking_Gest", "fat13_Smooking_Gest")
Df0 <- df[myvars] 

# Integrate with the Exposome dataset                                             ----- 
Data_final <- merge (Data, Df0, by = "ident")
rm (df,Data)
Data <- Data_final
rm (Data_final)

Data$fat20a_y <- NULL
# Parenting (KEEP) ------
"p_fat39a"	 #	Does your father smile at you, speak in a gentle, caring voice, praise the things you do right, try to understand the problems, and help you solve them?
"p_fat39b"	 #	Does your father give him freedom, letting him dress as he wants, making his own decisions?
"p_fat39c"	 #	Does your father treat you like a baby, that is, is very overprotective, not wanting him to grow up?
"fes44"	 #	Your family members usually try to win in no matter what (delete!!!)
"fes3"	 #	You (parents) fight a lot in the family.
summary (Data$p_fat39a)
table (Data$p_fat39a)

# Dataset with all paternal variables for score creation 
fathars <- c("ident",
  "reg36", # Is the child's biological father known?
  "reg36b1",	# Did the biological father live with the child at the child's first month of child's life?
  "reg36b2" ,	# Did the biological father live with the child from 2nd to 6th month of child's life?
  "reg36b3",	#	Did the biological father live with the child from 7th to 11th month of child's life?
  "reg36b4",	#	Did the biological father live with the child from 1st to 3rd year of child's life?
  "reg36b5",	#	Did the biological father live with the child from 4th to 6th year of child's life?
  "reg36b6",	#	Did the biological father live with the child from 7th to 9th year of child's life?
  "reg36c1",	#	Did the biological father keet in touch during the first month of child's life?
  "reg36c2",	#	Did the biological father keep in touch between the 2nd to 6th month of child's life?
  "reg36c3",	#	Did the biological father keep in touch between the 7th to 11th month of child's life?
  "reg36c4",	#	Did the biological father keep in touch between the 1st to 3rd year of child's life?
  "reg36c5",	#	Did the biological father keep in touch between the 4th to 6th year of child's life?
  "reg36c6" ,	#	Did the biological father keep in touch between the 7th to 9th year of child's life?
  "reg37" ,	#	How do/did (child's name) and the biological father get along when they are/were together?
  "reg39b",	#	Does the biological father pay alimony?
  "reg47"	,#	What is the marital status of the biological father?
  "fat1"	,#	For how long were the biological mother and father together in a relationship, when she found out to be pregnant?
  "fat38a",	#	Did/Does the biological mother smile at the child, speak in a gentle, friendly voice, praise the things the child did/does (correctly), trying to understand the problems and help solve them?
  "fat38b",	#	Did/Does the biological mother give the child freedom, letting him decides what to wear?
  "fat38c",	#	Did/Does the biological mother treat the child like a baby, very overprotective, not wanting him to grow?
  "fat39a",	#	Did/Does the biological father smile at the child, speak in a gentle, friendly voice, praise the things the child did/does (correctly), trying to understand the problems and help solve them?
  "fat39b",	#	Did/Does the biological father give the child freedom, letting him decides what to wear?
  "fat39c"	#	Did/Does the biological father treat the child like a baby, very overprotective, not wanting him to grow?
  )


Father <- Df0[fathars]
summary(Father$reg36)
table(Father$reg36)
Father <- Father %>%
  mutate(father_Known = ifelse(reg36 == 2 | reg36 == 3, 1, 0),
         father_Unknown = ifelse(reg36 == 1, 1, 0),
         father_Deceased = ifelse(reg36 == 3, 1, 0))
table (Father$father_Known)
Father$reg36 <- NULL


cols_to_replace <- c("reg36b1", "reg36b2", "reg36b3", "reg36b4", "reg36b5", "reg36b6")
Father[cols_to_replace] <- lapply(Father[cols_to_replace], function(x) { x[is.na(x)] <- 1; return(x) })
# Create the score by summing the variables
Father$liveWfather_score <- rowSums(Father[, c("reg36b1", "reg36b2", 
                                               "reg36b3", "reg36b4", 
                                               "reg36b5", "reg36b6")], na.rm = TRUE)

table (Father$liveWfather_score)

# Replace values 99 with 0 in the "in contact" variables
cols_to_replace <- c("reg36c1", "reg36c2", "reg36c3", "reg36c4", "reg36c5", "reg36c6")
Father[cols_to_replace] <- lapply(Father[cols_to_replace], function(x) { x[is.na(x)] <- 1; return(x) })

# Create the "in contact" score by summing the variables
Father$in_touch_father_score <- rowSums(Father[, c("reg36c1", "reg36c2", 
                                                  "reg36c3", "reg36c4", 
                                                  "reg36c5", "reg36c6")], na.rm = TRUE)

table (Father$in_touch_father_score)

table (Father$reg37)
Father <- Father %>%
  mutate(frelation_VeryWell = ifelse(reg37 == 1, 1, 0),
         frelation_Well = ifelse(reg37 == 2, 1, 0),
         frelation_Bad = ifelse(reg37 == 4 | reg37 == 5 | reg37 == 99, 1, 0))

table(Father$frelation_Bad)

table(Father$reg39b)
Father <- Father %>%
  mutate(alimony_problems = ifelse(reg39b == 0 | reg39b == 1, 1, 0))

table(Father$alimony_problems)
Father$alimony_problems[is.na(Father$alimony_problems)] <- 0
table(Father$alimony_problems)

reg47
table(Father$reg47)
Father <- Father %>%
  mutate(fmarital_MarriedWithBioMother = ifelse(reg47 == 1, 1, 0),
         fmarital_MarriedWithOtherPartner = ifelse(reg47 == 2, 1, 0),
         fmarital_Single = ifelse(reg47 == 3, 1, 0),
         fmarital_SeparatedDivorced = ifelse(reg47 == 4, 1, 0),
         fmarital_Widow = ifelse(reg47 == 5, 1, 0))
table(Father$fmarital_MarriedWithBioMother)
summary(Father$fmarital_MarriedWithBioMother)

fat1
summary(Father$fat1)
table(Father$fat1)


summary(Father$fat38a)
table(Father$fat38a)
Father <- Father %>%
  mutate(fat38a_NeverSometimes = ifelse(fat38a == 0 | fat38a == 1 | fat38a == 99, 1, 0),
         fat38a_Often = ifelse(fat38a == 2 | fat38a == 3, 1, 0))
table (Father$fat38a_NeverSometimes)

summary(Father$fat38b)
table(Father$fat38b)
Father <- Father %>%
  mutate(fat38b_NeverSometimes = ifelse(fat38b == 0 | fat38b == 1 | fat38b == 99, 1, 0),
         fat38b_Often = ifelse(fat38b == 2 | fat38b == 3, 1, 0))
table (Father$fat38b_NeverSometimes)

summary(Father$fat38c)
table(Father$fat38c)
Father <- Father %>%
  mutate(fat38c_NeverSometimes = ifelse(fat38c == 0 | fat38c == 1 | fat38c == 99, 1, 0),
         fat38c_Often = ifelse(fat38c == 2 | fat38c == 3, 1, 0)) # Often is the biggest problem here!
table (Father$fat38c_Often)


summary(Father$fat39a)
table(Father$fat39a)
Father <- Father %>%
  mutate(fat39a_NeverSometimes = ifelse(fat39a == 0 | fat39a == 1 | fat39a == 99, 1, 0),
         fat39a_Often = ifelse(fat39a == 2 | fat39a == 3, 1, 0))
table (Father$fat39a_NeverSometimes)

summary(Father$fat39b)
table(Father$fat39b)
Father <- Father %>%
  mutate(fat39b_NeverSometimes = ifelse(fat39b == 0 | fat39b == 1 | fat39b == 99, 1, 0),
         fat39b_Often = ifelse(fat39b == 2 | fat39b == 3, 1, 0))
table (Father$fat39b_NeverSometimes)

summary(Father$fat39c)
table(Father$fat39c)
Father <- Father %>%
  mutate(fat39c_NeverSometimes = ifelse(fat39c == 0 | fat39c == 1 | fat39c == 99, 1, 0),
         fat39c_Often = ifelse(fat39c == 2 | fat39c == 3, 1, 0)) # Often is the biggest problem here!
table (Father$fat39c_NeverSometimes)

# Including the main new variables -----
names(Father)
myvars <- c(
  "ident"     ,                       
  "father_Known"     ,               
  "father_Unknown"      ,             
  "father_Deceased"   ,
  "liveWfather_score"     ,          
  "in_touch_father_score"   ,         
  "frelation_VeryWell"      ,
  "frelation_Well"          ,        
  "frelation_Bad"            ,        
  "alimony_problems"      ,           
  "fmarital_MarriedWithBioMother" ,  
  "fmarital_MarriedWithOtherPartner", 
  "fmarital_Single"                ,  
  "fmarital_SeparatedDivorced"  ,    
  "fmarital_Widow"                   ,
  "fat38a_NeverSometimes"           , 
#  "fat38a_Often"                 ,   
  "fat38b_NeverSometimes"            ,
#  "fat38b_Often"                     ,
#  "fat38c_NeverSometimes"         ,  
  "fat38c_Often"                     ,
  "fat39a_NeverSometimes"            ,
#  "fat39a_Often"                   , 
  "fat39b_NeverSometimes"            ,
#  "fat39b_Often"                     ,
#  "fat39c_NeverSometimes"           ,
  "fat39c_Often"  )

df_father <- Father[myvars]

Data_old <- Data

Data <- merge(Data, df_father, by = "ident")
rm (df_father, Father)
rm (fathars, list, myvars, proporcao)
# Family  
"fes14"	 #	Sometimes your family gets so angry that they throw things (EXCLUDE !!!)
"fes20r"	#	Your family members usually loses their heads easily (EXCLUDE !!!)
"fes26"	 #	Your family members criticize each other frequentily (EXCLUDE !!!)
"fes32"	 #	Your family members sometimes beat each other (EXCLUDE !!!)
"fes36"	#	The group spirit in your family is very low (EXCLUDE !!!)
"fes_confl"	#	FES: conflict score # Latent of the variables above + some others that did not pass the filter (KEEP ONLY THIS ONE)

Data$fes14 <- NULL
Data$fes20r <- NULL
Data$fes26 <- NULL
Data$fes32 <- NULL
Data$fes36 <- NULL

summary (Data$fes_confl)
table (Data$fes_confl)


# Inclusion of the outcome variable (model outcome) -------
Df0 <- subset(BHRC_all, BHRC_all$redcap_event_name == "wave0_arm_1")
myvars <- c("ident", "dcany")
df <- Df0[myvars] # complete W0 dataset

table (df$dcany)
df$dcany <- ifelse(df$dcany == 2, 1, df$dcany)
table (df$dcany)

# SAVE DATASET                                                      -------
Data <- merge (Data, df, by ="ident")

write_sav (Data, "BHRC_exposome_df_2024_preproc.sav")
# ADD SOCIOECONOMIC DATA ------
################################################################################
# Loading Data ------
Data <- read_sav("BHRC_exposome_df_2024_preproc.sav")
# Re-checking variables   I                        ------ 
# Order the variables (columns) by name
df_ordenado <- Data[, sort(names(Data))]
names(df_ordenado)

# ID
Data$ident
Data$gender

# Father problems 
Data$alimony_problems    
  summary(Data$alimony_problems)
  table(Data$alimony_problems)
Data$in_touch_father_score <- NULL
Data$liveWfather_score <- NULL
  {
    BHRC_all <- read_sav("<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav", encoding = "latin1")
    Father <- subset(BHRC_all, BHRC_all$redcap_event_name == "wave0_arm_1")
  Father[, c("reg36b1", "reg36b2", "reg36b3", "reg36b4", "reg36b5", "reg36b6")] <- 
    lapply(Father[, c("reg36b1", "reg36b2", "reg36b3", "reg36b4", "reg36b5", "reg36b6")], 
           function(x) replace(x, x == 99, 1))
  # Replace NA with 1 in the specified variables
  cols_to_replace <- c("reg36b1", "reg36b2", "reg36b3", "reg36b4", "reg36b5", "reg36b6")
  Father[cols_to_replace] <- lapply(Father[cols_to_replace], function(x) { x[is.na(x)] <- 1; return(x) })
  
  # Create the score by summing the variables
  Father$liveWfather_score <- rowSums(Father[, c("reg36b1", "reg36b2", 
                                                 "reg36b3", "reg36b4", 
                                                 "reg36b5", "reg36b6")], na.rm = TRUE)
  
  
  table (Father$liveWfather_score)
  
  # Replace NA with 1 in the specified variables
  cols_to_replace <- c("reg36c1", "reg36c2", "reg36c3", "reg36c4", "reg36c5", "reg36c6")
  Father[cols_to_replace] <- lapply(Father[cols_to_replace], function(x) { x[is.na(x)] <- 1; return(x) })
  
  Father[, c("reg36c1", "reg36c2", "reg36c3", "reg36c4", "reg36c5", "reg36c6")] <- 
    lapply(Father[, c("reg36c1", "reg36c2", "reg36c3", "reg36c4", "reg36c5", "reg36c6")], 
           function(x) replace(x, x == 99, 1))
  # Create the "in contact" score by summing the variables
  Father$in_touch_father_score <- rowSums(Father[, c("reg36c1", "reg36c2", 
                                                     "reg36c3", "reg36c4", 
                                                     "reg36c5", "reg36c6")], na.rm = TRUE)
  
  table (Father$in_touch_father_score)
  myvars <- c("ident", "liveWfather_score", "in_touch_father_score")
  Father <- Father[myvars]
  Data <- merge(Data, Father, by = "ident")
  }
Data$in_touch_father_score
  summary(Data$in_touch_father_score)
  table(Data$in_touch_father_score)
  Data <- Data %>%
    mutate(always_in_touch_father_score = ifelse(in_touch_father_score == 6, 1, 0))
  table (Data$always_in_touch_father_score)
  
Data$liveWfather_score 
  summary(Data$liveWfather_score)
  table(Data$liveWfather_score)
  Data <- Data %>%
    mutate(never_liveWfather_score = ifelse(liveWfather_score == 0 |liveWfather_score == 1, 1, 0))
  table (Data$never_liveWfather_score)
Data$p_fat39a <- NULL
Data$p_fat39b <- NULL
Data$p_fat39c <- NULL
Data$fat39a_NeverSometimes
  summary(Data$fat39a_NeverSometimes)
  table(Data$fat39a_NeverSometimes)
  table(Data$fat39a_NeverSometimes, Data$never_liveWfather_score)
  Data$fat39a_NeverSometimes <- ifelse(is.na(Data$fat39a_NeverSometimes) & Data$never_liveWfather_score == 1, 
                                     1, Data$fat39a_NeverSometimes)
  table(Data$fat39a_NeverSometimes)
  summary(Data$fat39a_NeverSometimes)
  Data$fat39a_NeverSometimes <- ifelse(is.na(Data$fat39a_NeverSometimes) & Data$father_Deceased == 1 | Data$father_Unknown == 1, 
                                       1, Data$fat39a_NeverSometimes)
  table(Data$fat39a_NeverSometimes)
  summary(Data$fat39a_NeverSometimes)

Data$fat39b_NeverSometimes
  summary(Data$fat39b_NeverSometimes)
  table(Data$fat39b_NeverSometimes)
  table(Data$fat39b_NeverSometimes, Data$never_liveWfather_score)
  Data$fat39b_NeverSometimes <- ifelse(is.na(Data$fat39b_NeverSometimes) & Data$never_liveWfather_score == 1, 
                                     1, Data$fat39b_NeverSometimes)
  table(Data$fat39b_NeverSometimes)
  summary(Data$fat39b_NeverSometimes)
  Data$fat39b_NeverSometimes <- ifelse(is.na(Data$fat39b_NeverSometimes) & Data$father_Deceased == 1 | Data$father_Unknown == 1, 
                                       1, Data$fat39b_NeverSometimes)
  table(Data$fat39b_NeverSometimes)
  summary(Data$fat39b_NeverSometimes)

Data$fat39c_Often
  summary(Data$fat39c_Often)
  table(Data$fat39c_Often)
  table(Data$fat39c_Often, Data$never_liveWfather_score)
  Data$fat39c_Often <- ifelse(is.na(Data$fat39c_Often) & Data$never_liveWfather_score == 1, 
                                     0, Data$fat39c_Often)
  table(Data$fat39c_Often)
  summary(Data$fat39c_Often)
  Data$fat39c_Often <- ifelse(is.na(Data$fat39c_Often) & Data$father_Deceased == 1 | Data$father_Unknown == 1, 
                                       0, Data$fat39c_Often)
  table(Data$fat39c_Often)
  summary(Data$fat39c_Often)
  
Data$frelation_Bad.x
  summary(Data$frelation_Bad.x)
  table(Data$frelation_Bad.x)
  Data$frelation_Bad.x <- ifelse(is.na(Data$frelation_Bad.x) & Data$father_Deceased == 1 | Data$father_Unknown == 1, 
                              1, Data$frelation_Bad.x)
  table(Data$frelation_Bad.x)
  summary(Data$frelation_Bad.x)
Data$frelation_Bad.y <- NULL
Data$frelation_VeryWell.x
  summary(Data$frelation_VeryWell.x)
  table(Data$frelation_VeryWell.x)
  Data$frelation_VeryWell.x <- ifelse(is.na(Data$frelation_VeryWell.x) & Data$father_Deceased == 1 | Data$father_Unknown == 1, 
                                 0, Data$frelation_VeryWell.x)
  table(Data$frelation_VeryWell.x)
  summary(Data$frelation_VeryWell.x)
Data$frelation_VeryWell.y <- NULL
Data$frelation_Well.x
  summary(Data$frelation_Well.x)
  table(Data$frelation_Well.x)
  Data$frelation_Well.x <- ifelse(is.na(Data$frelation_Well.x) & Data$father_Deceased == 1 | Data$father_Unknown == 1, 
                                      0, Data$frelation_Well.x)
  table(Data$frelation_Well.x)
  summary(Data$frelation_Well.x)
Data$frelation_Well.y <- NULL

Data$father_Deceased
  table(Data$father_Deceased)
  summary(Data$father_Deceased)
  
Data$father_Known
  table(Data$father_Known)
  summary(Data$father_Known)
  
Data$father_Unknown
  table(Data$father_Unknown)
  summary(Data$father_Unknown)
  
# Mother Problems 
Data$fat38a_NeverSometimes
  table(Data$fat38a_NeverSometimes)
  summary(Data$fat38a_NeverSometimes)

Data$fat38b_NeverSometimes
  table(Data$fat38b_NeverSometimes)
  summary(Data$fat38b_NeverSometimes)

Data$fat38c_Often
  table(Data$fat38c_Often)
  summary(Data$fat38c_Often)

Data$mrelation_Bad
  table(Data$mrelation_Bad)
  summary(Data$mrelation_Bad)

Data$mrelation_VeryWell
  table(Data$mrelation_VeryWell)
  summary(Data$mrelation_VeryWell)

Data$mrelation_Well
  table(Data$mrelation_Well)
  summary(Data$mrelation_Well)

# Parents relationship 
Data$divorce_ConflictingNotWitnessed
  table(Data$divorce_ConflictingNotWitnessed)
  summary(Data$divorce_ConflictingNotWitnessed)
Data$divorce_dich
  table(Data$divorce_dich)
  summary(Data$divorce_dich)
Data$divorce_FriendlyFewConflicts
  table(Data$divorce_FriendlyFewConflicts)
  summary(Data$divorce_FriendlyFewConflicts)
Data$divorce_FriendlyNoConflicts
Data$divorce_VeryConflictualWitnessed

Data$fmarital_MarriedWithBioMother
  table(Data$fmarital_MarriedWithBioMother)
  summary(Data$fmarital_MarriedWithBioMother)
  Data$fmarital_MarriedWithBioMother <- NULL
Data$fmarital_MarriedWithOtherPartner <- NULL
Data$fmarital_SeparatedDivorced <- NULL
Data$fmarital_Single <- NULL
Data$fmarital_Widow <- NULL

Data$mmarital_MarriedWithBioFather
  table(Data$mmarital_MarriedWithBioFather)
  summary(Data$mmarital_MarriedWithBioFather)
Data$mmarital_MarriedWithOtherPartner
  table(Data$mmarital_MarriedWithOtherPartner)
  summary(Data$mmarital_MarriedWithOtherPartner)
Data$mmarital_SeparatedDivorced
  table(Data$mmarital_SeparatedDivorced)
  summary(Data$mmarital_SeparatedDivorced)
Data$mmarital_Single
  table(Data$mmarital_Single)
  summary(Data$mmarital_Single)
Data$mmarital_Widow
  table(Data$mmarital_Widow)
  summary(Data$mmarital_Widow)

# Pregnancy 
Data$fat13_Drinking_Gest
  table(Data$fat13_Drinking_Gest)
  summary(Data$fat13_Drinking_Gest)
Data$fat15_Drink <- NULL
Data$fat13_Smooking_Gest
  table(Data$fat13_Smooking_Gest)
  summary(Data$fat13_Smooking_Gest)
Data$fat16c
  table(Data$fat16c)
  summary(Data$fat16c)
Data$fat17
  table(Data$fat17)
  summary(Data$fat17)
Data$fat3_AbortionConsidered
  table(Data$fat3_AbortionConsidered)
  summary(Data$fat3_AbortionConsidered)
Data$fat3_ModerateStress
Data$fat3_PassiveOrLightStress
Data$fat3_PositiveFeelings
Data$fat3_SevereStress

# Family 
Data$fat14a
Data$fat4_NoSupport
Data$fat4_SupportedMostTimes
Data$fat4_SupportedSomeMoments
Data$fat4_SupportedWheneverNeeded
Data$fes_confl
  table(Data$fes_confl)
  summary(Data$fes_confl)
Data$fes3 <- NULL
Data$fes44 <- NULL
Data$ghelp
  table(Data$ghelp)
  summary(Data$ghelp)

# Children 
Data$AnySubstanceUse
  table(Data$AnySubstanceUse)
  summary(Data$AnySubstanceUse)
Data$BadCo
Data$onlybreastfeeded 
  table(Data$onlybreastfeeded)
  summary(Data$onlybreastfeeded)
Data$fat20
  table(Data$fat20)
  summary(Data$fat20)
Data$fat20_monthsTotal
  table(Data$fat20_monthsTotal)
  summary(Data$fat20_monthsTotal)
  Data$min_breast <- ifelse(Data$fat20_monthsTotal >= 6, 1, 0)
  table(Data$min_breast)
Data$neonatal_UTI
  table(Data$neonatal_UTI)
  summary(Data$neonatal_UTI)
Data$Overweight
  table(Data$Overweight)
  summary(Data$Overweight)

# Trauma 
Data$c_bulliedhim
  table(Data$c_bulliedhim)
  summary(Data$c_bulliedhim)
Data$c_bulliedhim_fr
  table(Data$c_bulliedhim_fr)
  summary(Data$c_bulliedhim_fr)
Data$c_bulliedhim_imp
  table(Data$c_bulliedhim_imp)
  summary(Data$c_bulliedhim_imp)
  Data$c_bulliedhim_imp <- NULL
Data$c_bulliedhim_l
  table(Data$c_bulliedhim_l)
  summary(Data$c_bulliedhim_l)
  Data$c_bulliedhim_l <- NULL

Data$cemoab_exposed
  table(Data$cemoab_exposed)
  summary(Data$cemoab_exposed)
Data$cphysab_exposed
  table(Data$cphysab_exposed)
  summary(Data$cphysab_exposed)
Data$traumac_fl
table(Data$traumac_fl)
summary(Data$traumac_fl)

# Parents report (trauma)
Data$p_bulliedhim
  table(Data$p_bulliedhim)
  summary(Data$p_bulliedhim)
Data$p_bulliedhim_fr
  table(Data$p_bulliedhim_fr)
  summary(Data$p_bulliedhim_fr)
Data$p_bulliedhim_imp
  table(Data$p_bulliedhim_imp)
  summary(Data$p_bulliedhim_imp)
Data$p_bulliedhim_l
  table(Data$p_bulliedhim_l)
  summary(Data$p_bulliedhim_l)
Data$pemoab_exposed
Data$pphyab_exposed
Data$pphyneg_exposed
Data$psexab_exposed
Data$traumap_fl

# Outcome 
Data$dcany

# Re-checking variables   II                        ------ 
table (Df0$reg3) # IF THE MOTHER ANSWERED, THEN IT ENTERS = EXCLUDE 2 !!
myvars <- c("ident", "reg3", "age")
df <- Df0[myvars]

Data <- merge (Data, df, by ="ident")

Data_FINAL <- subset(Data, Data$reg3 != 2)
Data_FINAL$reg3 <- NULL


# Checking for missing values
missing_percentage <- sapply(Data_FINAL, function(x) sum(is.na(x)) / length(x) * 100)
missing_df <- data.frame(Variable = names(missing_percentage), MissingPercentage = missing_percentage)


# Selecting missing values 
missing_vars <- missing_df[missing_df$MissingPercentage > 0, ]
a <- missing_vars$Variable
print(a)

# Adjusting Data frame + Input data -------
Data_FINAL$fat39a_NeverSometimes[is.na(Data_FINAL$fat39a_NeverSometimes)] <- 0
Data_FINAL$fat39b_NeverSometimes[is.na(Data_FINAL$fat39b_NeverSometimes)] <- 0
Data_FINAL$fat39c_Often[is.na(Data_FINAL$fat39c_Often)] <- 0

rm (a, cols_to_replace, missing_percentage, myvars)
rm (missing_df, missing_vars, Father, df_ordenado, df, Data)

# Checking new data frame ------
df <- Data_FINAL
missing_percentage <- sapply(df, function(x) sum(is.na(x)) / length(x) * 100)
missing_df <- data.frame(Variable = names(missing_percentage), MissingPercentage = missing_percentage)
rm (missing_df, missing_percentage)

################################################################################
################################################################################
################################################################################
# HACK - MODIFICATION OF THE M/F RELATION VARIABLE ------
Df0 <- Df0 %>%
  mutate(mrelation_VeryWell = ifelse(mrelation == 1, 1, 0),
         mrelation_Well = ifelse(mrelation == 2, 1, 0),
         mrelation_Bad = ifelse(mrelation == 4 | mrelation == 5, 1, 0))

Df0 <- Df0 %>%
  mutate(frelation_VeryWell = ifelse(reg37 == 1, 1, 0),
         frelation_Well = ifelse(reg37 == 2, 1, 0),
         frelation_Bad = ifelse(reg37 == 4 | reg37 == 5 | reg37 == 99, 1, 0))

table(Df0$frelation_Bad)

Df0 <- Df0 %>%
  mutate(father_Known = ifelse(reg36 == 2 | reg36 == 3, 1, 0),
         father_Unknown = ifelse(reg36 == 1, 1, 0),
         father_Deceased = ifelse(reg36 == 3, 1, 0))
table (Df0$father_Known)
Df0$frelation_Bad <- ifelse(is.na(Df0$frelation_Bad) & Df0$father_Deceased == 1 | Df0$father_Unknown == 1, 
                            1, Df0$frelation_Bad)
table(Df0$frelation_Bad)
summary(Df0$frelation_Bad)

myvars <- c("ident", "frelation_Bad", "frelation_VeryWell", "mrelation_Bad", "mrelation_VeryWell")
Gamb <- Df0[myvars]
df$frelation_Bad.x <- NULL
df$frelation_VeryWell.x <- NULL
df$frelation_Well.x <- NULL
df$mrelation_Bad <- NULL
df$mrelation_VeryWell <- NULL
df$mrelation_Well <- NULL

df <- merge(df, Gamb, by ="ident")
teste <- merge(Data_FINAL, Gamb, by ="ident")

# Checking for missing values
missing_percentage <- sapply(teste, function(x) sum(is.na(x)) / length(x) * 100)
missing_df <- data.frame(Variable = names(missing_percentage), MissingPercentage = missing_percentage)


# Selecting missing values 
missing_vars <- missing_df[missing_df$MissingPercentage > 0, ]
a <- missing_vars$Variable
print(a)

# Adjusting Data frame + Input data -------
teste$frelation_VeryWell[is.na(teste$frelation_VeryWell)] <- 0
teste$mrelation_Bad.y[is.na(teste$mrelation_Bad.y)] <- 0
teste$mrelation_VeryWell.y[is.na(teste$mrelation_VeryWell.y)] <- 0
################################################################################
################################################################################
# NOTE - SUGGESTION: ----------
# LASSO/ELASTIC - REDUCE DATASET 
# STEP TO CHECK IN THE LITERATURE: ACCURACY WILL BE BIASED (DIFFERENT PROPORTIONS IN THE OUTCOME)
# TOOK 70% TRAIN / 30% TEST - SET TO MATCH THE OUTCOME NUMBER (50% / 50%)

# START OF ANALYSES ------
tmp <- teste
tmp <- na.omit(tmp[1:length(tmp)])
tmp <- tmp[, sort(names(tmp))]
names(tmp)

# Double check of variables and exclusions ------
# 75 initial variables (not counting ident)
tmp$age
age_means <- tmp %>%
  group_by(dcany) %>%
  summarise(mean_age = mean(age, na.rm = TRUE))
  print(age_means)
tmp$alimony_problems
  table (tmp$alimony_problems)
tmp$always_in_touch_father_score
  table (tmp$always_in_touch_father_score)
tmp$AnySubstanceUse
  table (tmp$AnySubstanceUse)
tmp$BadCo
tmp$c_bulliedhim
tmp$c_bulliedhim_fr <- NULL
tmp$cemoab_exposed
tmp$cphysab_exposed
tmp$dcany
tmp$divorce_ConflictingNotWitnessed
  table (tmp$divorce_ConflictingNotWitnessed)
tmp$divorce_dich
  table (tmp$divorce_dich)
tmp$divorce_FriendlyFewConflicts
  table (tmp$divorce_FriendlyFewConflicts)
tmp$divorce_FriendlyNoConflicts
  table (tmp$divorce_FriendlyNoConflicts)
tmp$divorce_VeryConflictualWitnessed
tmp$fat13_Drinking_Gest
tmp$fat13_Smooking_Gest
tmp$fat14a
tmp$fat16c
tmp$fat17
tmp$fat20
tmp$fat20_monthsTotal
  summary(tmp$fat20_monthsTotal)
tmp$min_breast
  table(tmp$min_breast)
tmp$fat3_AbortionConsidered
  table(tmp$fat3_AbortionConsidered)
tmp$fat3_ModerateStress
tmp$fat3_PassiveOrLightStress
tmp$fat3_PositiveFeelings
tmp$fat3_SevereStress
tmp$fat38a_NeverSometimes
  table(tmp$fat3_AbortionConsidered)
tmp$fat38b_NeverSometimes
tmp$fat38c_Often
tmp$fat39a_NeverSometimes
tmp$fat39b_NeverSometimes
tmp$fat39c_Often
tmp$fat4_NoSupport
tmp$fat4_SupportedMostTimes
tmp$fat4_SupportedSomeMoments
tmp$fat4_SupportedWheneverNeeded
tmp$father_Deceased
tmp$father_Known
tmp$father_Unknown
tmp$fes_confl
tmp$frelation_Bad
  summary(tmp$frelation_Bad)
  table(tmp$frelation_Bad)
tmp$frelation_Bad.x
  summary(tmp$frelation_Bad.x)
  table(tmp$frelation_Bad.x)
  tmp$frelation_Bad.x <- NULL
tmp$frelation_VeryWell
  summary(tmp$frelation_VeryWell)
  table(tmp$frelation_VeryWell)
tmp$frelation_VeryWell.x
  summary(tmp$frelation_VeryWell.x)
  table(tmp$frelation_VeryWell.x)
  tmp$frelation_VeryWell.x <- NULL
tmp$frelation_Well.x
  summary(tmp$frelation_Well.x)
  table(tmp$frelation_Well.x)
  table(tmp$frelation_Well.x, tmp$frelation_VeryWell)
tmp$gender
  table(tmp$gender)
tmp$ghelp
  table(tmp$ghelp)
tmp$ident
tmp$in_touch_father_score
tmp$liveWfather_score
tmp$mmarital_MarriedWithBioFather
  table(tmp$mmarital_MarriedWithBioFather)
tmp$mmarital_MarriedWithOtherPartner
tmp$mmarital_SeparatedDivorced
tmp$mmarital_Single
tmp$mmarital_Widow
tmp$mrelation_Bad.x
  table(tmp$mrelation_Bad.x)
tmp$mrelation_Bad.y
  table(tmp$mrelation_Bad.y)
  tmp$mrelation_Bad.x <- NULL
tmp$mrelation_VeryWell.x
  table(tmp$mrelation_VeryWell.x)
tmp$mrelation_VeryWell.y
  table(tmp$mrelation_VeryWell.y)
  tmp$mrelation_VeryWell.x <- NULL
tmp$mrelation_Well
  table(tmp$mrelation_Well)
tmp$neonatal_UTI
  table(tmp$neonatal_UTI)
tmp$never_liveWfather_score
tmp$onlybreastfeeded
tmp$Overweight
tmp$p_bulliedhim
tmp$p_bulliedhim_fr <- NULL
tmp$p_bulliedhim_imp <- NULL
tmp$p_bulliedhim_l <- NULL
tmp$pemoab_exposed
tmp$pphyab_exposed
tmp$pphyneg_exposed
tmp$psexab_exposed
tmp$traumac_fl <- NULL
tmp$traumap_fl <- NULL

# 66 selected variables 

write_sav (tmp, "BHRC_exposome_df_2024_INITIAL_ANALYSIS_withID.sav")
tmp$ident <- NULL
write_sav (tmp, "BHRC_exposome_df_2024_INITIAL_ANALYSIS_noID.sav")

# OUTCOME VARIABLE POSITION = BE CAREFUL TO ALWAYS KEEP IT AS THE FIRST VARIABLE IN THE DATASET  ------
tmp <- read_sav("BHRC_exposome_df_2024_INITIAL_ANALYSIS_noID.sav")
tmp <- tmp %>%
  select(dcany, everything())
table(tmp$dcany)
