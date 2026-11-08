# Using RPACT

# Install RPACT and GGPLOT2

library(rpact)

# If necessary:
#
# install.packages("ggplot2")

# Example 2: Survival data
#
# Two-arm trial, new treatment vs control
#
# Observations follow a proportional hazards model with
# log hazard ratio theta, where
#
#    hazard ratio > 1 and theta > 0
#
# imply the new treatment is superior
#
# To test H_0: theta =< 0 vs theta > 0 with one-sided type I error
# probability alpha=0.025 and power 1-beta=0.8 at
#
#   theta = delta = 0.5 = log(1.65)

alpha=0.025
beta=0.2
delta=0.5
ifix=(qnorm(1-alpha)+qnorm(1-beta))^2/delta^2
ifix

# Fixed sample size design
# We shall analyse data using the logrank statistic, assuming
# a normal approximation, so give the command
#   normalApproximation = TRUE
ss_mean <- getSampleSizeMeans(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  alternative = 0.5,
  stDev = 1,
  normalApproximation = TRUE
)
summary(ss_mean)
ss_mean$maxNumberOfSubjects/4

# "Number of subjects=125.6" in a 2 treatment comparison with
# normal responses of variance 1 implies that the required
# information is 125.6/4 = 31.4
#
# For a logrank test, this corresponds, approximately,
# to 4*31.4 = 125.6 events

# Group sequential design

# Design a group sequential trial with 5 analyses
#
# Error spending design, non-binding futility boundary, rho-family
# error spending function, rho=2 for efficacy and futility boundaries
#
# In RPACT the rho-family is called the Kim and DeMets alpha
# spending function and must be specified for both the
# efficacy and the futility boundaries.
#
# Define a design object with 5 analyses and a parameter
# of 2 for the rho family
des_surviv <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
# See properties of this design, including efficacy
# and futility boundaries on the Z-scale
summary(des_surviv)
# These boundaries on the Z-scale apply to all cases with
# power 1-beta at the treatment effect under which type 2
# error is spent.
# We shall find the information needed to give power 1-beta
# at our specified alternative theta=delta.
# Based on our group sequential design, we first find the
# sample size required in a 2 treatment comparison with
# normal responses of variance 1
ss_gsd <- getSampleSizeMeans(
  design = des_surviv,
  alternative = delta,
  stDev = 1,
  normalApproximation = T
)
ss_gsd$maxNumberOfSubjects

# "Number of subjects=142.3" in a 2 treatment comparison with
# normal responses of variance 1 implies that the required
# information is 142.3/4 = 35.58
#
# For a logrank test, this corresponds, approximately,
# to 4*35.58 = 142.3 events

# Extract design characteristics to get inflation factor
# (and other information) for the group sequential design
dc_des_surviv <- getDesignCharacteristics(des_surviv)

dc_des_surviv$inflationFactor
dc_des_surviv$.design$alphaSpent
dc_des_surviv$.design$betaSpent

# Note the inflation factor is equal to the ratio of numbers
# of events in fixed and sequential designs
#  142.3/125.6 =  1.133

plot(des_surviv)

# store the maximum information for this design
iMax <- ifix * dc_des_surviv$inflationFactor
iMax

# Under theta=0 (hazard ratio=1)
asn <- getPowerAndAverageSampleNumber(
  design = des_surviv,
  theta = 0,
  nMax = iMax
)
# Probabilities of rejecting the null hypothesis per stage
asn$rejectPerStage
# Probabilities of stopping for futility per stage
asn$futilityPerStage
# overall probability of rejecting the null hypothesis
asn$calculatedPower

# Note that, due to the no-binding futility boundary, the
# actual type 2 error rate is 0.0232


# Example data set

# We explore the consequences of observed sample sizes and
# observed information fractions that deviate from the
# originally planned design

# Suppose we have cumulative information levels

iobs=c(5.43,12.58,21.11,30.55,33.28)

# and treatment effect estimates

theta.hat=c(0.10,0.06,0.21,0.31)

# Efficacy boundary

# Analysis 1

# The observed information fraction is iobs[1]/35.58 = 5.43/35.58
#
# We need to state the information fractions for later
# analyses. We can just put in the values that were
# originally planned -- these will not affect the
# first boundary value
designUpdate1 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[1]/35.58,0.4,0.6,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate1)
# We see the first efficacy boundary point is at Z=3.247

# Continuing to the next analysis
# Now we give the observed information at analyses 1 and 2
designUpdate2 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[c(1,2)]/35.58,0.6,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate2)
# We see the second efficacy boundary point is at Z=2.782

# Continuing to analysis 3
# Now we give the observed information at analyses 1, 2 and 3
designUpdate3 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[c(1,2,3)]/35.58,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate3)
# We see the third efficacy boundary point is at Z=2.461

designUpdate4 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[c(1,2,3,4)]/35.58,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate4)
# The fourth efficacy boundary point is at Z=2.201

# The final analysis is handled differently to allow a
# final information level that differs from the target.
#
# In calculating the efficacy boundary, it is important that
# the RATIOS of information levels are correct. These ratios
# are preserved if we sample sizes by the final observed sample
# size to make the last value in informationRates equal to 1.
#
# However, we need to preserve the type I error probabilities
# already seen at earlier analyses and these depended on the
# observed sample sizes. This is achieved by adding the commands
#  typeOfDesign = "asUser",
#  userAlphaSpending = designUpdate4$alphaSpent
# Now, the remainder of alpha will be spent at the final analysis.
designUpdate5 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[c(1,2,3,4,5)]/iobs[5]), 
  typeOfDesign = "asUser",
  userAlphaSpending = designUpdate4$alphaSpent
)
summary(designUpdate5)
# The final efficacy boundary point is at Z=2.088

designUpdate5$criticalValues

# These values agree very nicely with those from CJ's code

# Futility boundary

# RPACT does not yet provide an error spending FUTILITY
# boundary. Given the negative outcome of this trial, we
# would like to know the futility boundary values.
# From CJ's code, the futility boundary values are
#
#   -1.4351 -0.2322 0.7497 1.6382 2.0187

des_surv_1 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1]/35.58,0.4,0.6,0.8,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
des_surv_1$futilityBounds

des_surv_2 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1:2]/35.58,0.6,0.8,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
des_surv_2$futilityBounds

des_surv_3 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1:3]/35.58,0.8,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
des_surv_3$futilityBounds

des_surv_4 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1:4]/35.58,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
des_surv_4$futilityBounds



# Analysis of Cox model output

# The initial design is the same as for the logrank statistic.
# It is the observed information levels and Z-statistics that
# change.

# So, we repeat the above calculations with

iobs=c(4.11,10.89,19.23,28.10,30.96)

# Efficacy boundary

designUpdate1 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[1]/35.58,0.4,0.6,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate1)
# We see the efficacy boundary is at Z=3.403

designUpdate2 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[c(1,2)]/35.58,0.6,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate2)
# We see the second efficacy boundary point is at Z=2.865

designUpdate3 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[c(1,2,3)]/35.58,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate3)
# We see the third efficacy boundary point is at Z=2.521

designUpdate4 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[c(1,2,3,4)]/35.58,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate4)
# The fourth efficacy boundary point is at Z=2.267

designUpdate5 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.2,
  informationRates = c(iobs[c(1,2,3,4,5)]/iobs[5]), 
  typeOfDesign = "asUser",
  userAlphaSpending = designUpdate4$alphaSpent
)
summary(designUpdate5)
# The final efficacy boundary point is at Z=2.049

designUpdate5$criticalValues

# Futility boundary

# RPACT does not yet provide an error spending FUTILITY
# boundary. Given the negative outcome of this trial, we
# would like to know the futility boundary values.
# From CJ's code, the futility boundary values are
#
#   -1.7722 -0.4675 0.5546 1.4137 1.6138

# For an approximate answer, we can follow the approach
# illustrated in Example 1.

des_surv_1 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1]/35.58,0.4,0.6,0.8,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
des_surv_1$futilityBounds

des_surv_2 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1:2]/35.58,0.6,0.8,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
des_surv_2$futilityBounds

des_surv_3 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1:3]/35.58,0.8,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
des_surv_3$futilityBounds

des_surv_4 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1:4]/35.58,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
des_surv_4$futilityBounds
