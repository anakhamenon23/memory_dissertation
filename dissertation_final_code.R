# loading the necessary packages

library(dplyr)
library(tidyr)
library(tidyverse)
library(lme4) # for linear mixed effects model
library(pwr) # for calculating power
library(effectsize) # for effect size calculation
# install.packages("rstatix")
library(rstatix)
#install.packages("ggpubr")
library(ggpubr)
library(cowplot)
library(ggplot2)
#install.packages("hrbrthemes")
library(hrbrthemes)
#install.packages('Rmisc', dependencies = TRUE)
library(Rmisc)
#install.packages("nlme")
library(nlme)
#install.packages("effects")
library(effects)
#install.packages("report")
library(report)
library(introdataviz)
#install.packages("coin")
library(coin)
# loading raw data
raw_data <- readr::read_csv("fMRI_Data.csv")

# removing irrelevant sections
mydata <- raw_data|>select(-c(HitorMiss, Ret_OL, DG))

write.csv(mydata, file = "mydata.csv")

# making the tables wide so encoding and retrieval datas are seperate [subject to change]
# mydata_wider <- pivot_wider(mydata, names_from= Enc_or_Ret, values_from = c(CA1, CA3), names_prefix = c("Encoding", "Retrieval"))

mydata_long <- pivot_longer(mydata,
                            cols = CA1:CA3,
                            names_to = "Activity")
# checking for outliers CA1
mydata|>
  group_by(Enc_or_Ret) |>
  identify_outliers(CA1) |>
print(n= 63) 
# results suggest SDT04 has high extreme outliers

# checking for outliers CA3
mydata|>
  group_by(Enc_or_Ret) |>
  identify_outliers(CA3) |>
print(n= 63) 

# checking for normality using shapiro test for CA1
mydata|>
  group_by(CA) |>
  shapiro_test(CA1) |>
  print()
# both enc and ret are <0.05

# checking for normality using shapiro test for CA3
mydata|>
  group_by(Enc_or_Ret) |>
  shapiro_test(CA3) |>
  print()
# both enc and ret are <0.05



# checking for normality using qqplots
ggqqplot(mydata, "CA1", facet.by = "Enc_or_Ret")

ggqqplot(mydata, "CA3", facet.by = "Enc_or_Ret")

ggplot(data = mydata_long,
       aes(Trial, Enc_or_Ret, colour = Activity)) +
        geom_smooth(formula = y ~ x, method="lm")

str(mydata_long)

# turning doubles into factors
mydata_long$Enc_or_Ret <- factor(mydata_long$Enc_or_Ret)
str(mydata_long)

mydata_long$Activity <- factor(mydata_long$Activity)
str(mydata_long)

mydata_long$PID <- factor(mydata_long$PID)

# checking the structure
str(mydata_long)


# summarising data
sum_mydata_long <- summarySE(mydata_long, measurevar = "value", groupvars = c("Activity", "Enc_or_Ret"))
sum_mydata_long

# subseting data for visualisation
data_encoding_long <- subset(mydata_long, Enc_or_Ret == 1)
data_retrieval_long <- subset(mydata_long, Enc_or_Ret == 2)

data_encoding_long$PID <- factor(data_encoding_long$PID)
is.factor(data_encoding_long$PID)

data_encoding_long$Activity <- factor(data_encoding_long$Activity)

#  dividing mydata into encoding and retrieval
mydata_encoding <- subset(mydata, Enc_or_Ret == 1)
mydata_retrieval <- subset(mydata, Enc_or_Ret == 2)

bartlett_retrieval <- bartlett.test(CA3 ~ CA1, data = mydata_retrieval)
print(bartlett_retrieval)

write.csv(mydata_encoding, file = "mydata_encoding.csv")
# making PID factor 
mydata_encoding$PID <- factor(mydata_encoding$PID)
is.factor(mydata_encoding$PID)
View(levels(data_encoding_long$PID))


# extracting encoding data for both regions
encoding_CA1 <- subset(mydata_encoding, select = c(PID, Trial, Enc_or_Ret, CA1))
encoding_CA3 <- subset(mydata_encoding, select = c(PID, Trial, Enc_or_Ret, CA3))

# extracting retrieval data for both regions
retrieval_CA1 <- subset(mydata_retrieval, select = c(PID, Trial, Enc_or_Ret, CA1))
retrieval_CA3 <- subset(mydata_retrieval, select = c(PID, Trial, Enc_or_Ret, CA3))


# t test

wilcox_encoding <- wilcox.test(encoding_CA1$CA1, encoding_CA3$CA3, paired = TRUE)
print(wilcox_encoding)
report(wilcox_encoding)
wilcox_retrieval <- wilcox.test(retrieval_CA1$CA1, retrieval_CA3$CA3, paired = TRUE)
print(wilcox_retrieval)
report(wilcox_retrieval)

mydata_encoding |>
  group_by(Enc_or_Ret) |>
  wilcox_effsize(CA1~CA3)

t_test_result_encoding <- t.test(encoding_CA1$CA1, encoding_CA3$CA3, paired = TRUE, alternative = "two.sided")
print(t_test_result_encoding)
report(t_test_result_encoding)

t_test_result_retrieval <- t.test(retrieval_CA1$CA1, retrieval_CA3$CA3, paired = TRUE, alternative = "two.sided")
print(t_test_result_retrieval)
report(t_test_result_retrieval)

# lmer for CA1
lmm_CA1 <- lmer(CA1 ~ Enc_or_Ret + (1|PID), data = mydata)
summary(lmm_CA1)
confint(lmm_CA1, level = 0.95)
report(lmm_CA1)

#lmer for CA3
lmm_CA3 <- lmer(CA3 ~ Enc_or_Ret + (1|PID), data = mydata)
summary(lmm_CA3)
confint(lmm_CA3, level = 0.95)
report(lmm_CA3)

# Create a coefficient plot
effects_plot <- allEffects(lmm_CA1)
plot(effects_plot)

plot(lmm_CA1, which = 1)

# raincloud plot

raincloudplot_mydata_long <- ggplot(mydata_long, aes(x =Enc_or_Ret, y = value, fill = Activity)) +
  geom_flat_violin(aes(fill = Activity),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(Enc_or_Ret)-.15, y = value, colour = Activity),position = position_jitter(width = .05), size = 1, shape = 20)+
  geom_boxplot(aes(x = Enc_or_Ret, y = value, fill = Activity),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  labs(y = "Bold Activity")+
  labs(x = "Condition")+
  scale_x_discrete(breaks = c(1, 2), labels = c("Encoding", "Retrieval"))+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"))+
  ggtitle("Figure 3: BOLD Activity in CA1 and CA3 Regions During Encoding and Retrieval")+
  coord_flip()


raincloudplot_mydata_long

