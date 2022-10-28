library(naniar)
library(missForest)
library(dplyr)
library(lme4)

#################################################################################
# PCA dimension reduction
#################################################################################

df = read.csv("/Users/antonioqin/Desktop/D/2022 FALL/ML for Health/p3_full.csv")
# visualize missingness
vis_miss(df)
gg_miss_var(df, show_pct = TRUE)
plot(df$X.NA, ylab = "missing values", main = "Missing Values for Each Participant")
abline(a = 30,b = 0,col = "red")

df1= df[which(complete.cases(df)),]
df1 = df1[,6:92]
df.mis = prodNA(df1,noNA = 0.2)

iris.imp <- missForest(df.mis, xtrue = df1, verbose = TRUE)

#remove anyone who identified sex as other than male and female
df$sex_other = NULL
df = df[!is.na(df$sex),]
#remove anyone who has over 30% missingness (30 NA)
df$X.NA.1 = NULL
df = df[df$X.NA < 30,]
#remove employee_status_other 
df$employ_status_other = NULL

#CRT group of variables
# 1 = correct, 2 = intuitive error, 3 = error, 4 = Missing
# CRT1 correct (0.05, 5, 25, 0.25)
df$CRT1.t = 3
df$CRT1.t[which(df$CRT1 == 0.05| df$CRT1 == 5| df$CRT1 == 25| df$CRT1 == 0.25)] = 1
# CRT1 intuitive error (10, 0.1, 50, 0.5)
df$CRT1.t[which(df$CRT1 == 10| df$CRT1 == 0.1| df$CRT1 == 50| df$CRT1 == 0.5)] = 2
# Missing = 4
df$CRT1.t[which(is.na(df$CRT1))] = 4
# CRT2 correct(3, 5)
df$CRT2.t = 3
df$CRT2.t[which(df$CRT2 == 3| df$CRT2 == 5)] = 1
# CRT2 intuitive error (100, 300)
df$CRT2.t[which(df$CRT2 == 100| df$CRT2 == 300)] = 2
# Missing = 4
df$CRT2.t[which(is.na(df$CRT2))] = 4
# CRT3 correct(7, 47)
df$CRT3.t = 3
df$CRT3.t[which(df$CRT3 == 7| df$CRT3 == 47)] = 1
# CRT3 intuitive error(4,24)
df$CRT3.t[which(df$CRT3 == 4| df$CRT3 == 24)] = 2
# Missing = 4
df$CRT3.t[which(is.na(df$CRT3))] = 4
# Remove redundant variables
df$CRT1 = NULL
df$CRT2 = NULL
df$CRT3 = NULL
df$X.NA = NULL
gg_miss_var(df, show_pct = TRUE)

md.pattern(df)
# Mice imputation does not wor because of perfect colinear varibles. Still gives NA
df_before = df
# Use random forest to impute without ID and date
df_after = df[,-c(1,2,4,5)]
imp = missForest(df_after, verbose = FALSE)
df_after = cbind(df[,c(1,2,4,5)], imp$ximp)

ps = function(x){
  k = prcomp(x)
  return(as.matrix(x) %*% k$rotation[,1])
}

#df_after = df[which(complete.cases(df)),]

head(df_after[,6:10])
physical_contact = ps(df_after[,6:10])
head(df_after[,11:15])
physical_hygiene = ps(df_after[,11:15])
head(df_after[,16:20])
policy_support = ps(df_after[,16:20])
head(df_after[,21:23])
generosity = ps(df_after[,21:23])
head(df_after[,25:26])
psych_wellbeing = ps(df_after[,25:26])
head(df_after[,27:29])
collective_narcis = ps(df_after[,27:29])
head(df_after[,30:31])
national_identity = ps(df_after[,30:31])
head(df_after[,32:35])
Conspiracy_theories = ps(df_after[,32:35])
head(df_after[,36:41])
open_mindedness = ps(df_after[,36:41])
head(df_after[,42:48])
morality_as_cooperat = ps(df_after[,42:48])
head(df_after[,49:50])
trait_optimism = ps(df_after[,49:50])
head(df_after[,51:54])
social_belonging = ps(df_after[,51:54])
head(df_after[,55:58])
trait_self.control = ps(df_after[,55:58])
head(df_after[,60:65])
Narcissism = ps(df_after[,60:65])
head(df_after[,66:75])
Moral_ID = ps(df_after[,66:75])
head(df_after[,76:77])
risk_perception = ps(df_after[,76:77])
head(df_after[,90:92])
CRT = ps(df_after[,90:92])

final = cbind(df_after[,1:5], 
              physical_contact,
              physical_hygiene,
              policy_support,
              generosity,
              df_after[,24],
              psych_wellbeing,
              collective_narcis,
              national_identity,
              Conspiracy_theories,
              open_mindedness,
              morality_as_cooperat,
              trait_optimism,
              social_belonging,
              trait_self.control,
              df_after[,59],
              Narcissism,
              Moral_ID,
              risk_perception,
              df_after[,78:89],
              CRT)

# remove index
df$X = NULL
# recoding y value
df$tested_positive = as.numeric(df$tested_positive) - 1
# Logistic regression model
logi = glm(tested_positive ~., data = select(df, -Country), family = 'binomial')
k = summary(logi)

#################################################################################
# Mixed effect model
#################################################################################

# Accounting for country

mixedmodel = lmer(tested_positive ~ physical_contact + physical_hygiene + policy_support + generosity + psych_wellbeing + collective_narcis + national_identity 
                  + Conspiracy_theories + open_mindedness + morality_as_cooperat + trait_optimism + social_belonging + trait_self.control + Narcissism + Moral_ID + risk_perception
                  + political_ideology + moral_circle + health_cond + sex + age + marit_status + children + employ_status + Ladder+ urban + tested_positive + know_tested_positive + CRT
                  + (1 | Country), data = df)

l = summary(mixedmodel)
tstat=summary(mixedmodel)$coef[,3]
pvalue=2*(1-pnorm(abs(tstat),0,1))
k2 = cbind(l$coefficients,pvalue)
write.csv(k2, "/Users/antonioqin/Desktop/D/2022 FALL/ML for Health/p3.csv")



