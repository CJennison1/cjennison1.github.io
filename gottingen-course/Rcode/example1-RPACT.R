# Using RPACT

# Install RPACT and GGPLOT2

library(rpact)

# If necessary:
#
# install.packages("ggplot2")

# Example 1: Normal responses
#
# Two-arm trial, new treatment vs control
#
# Observations distributed as
#
#   X~N(mu_c,sigma^2)  on the control arm
#
#   X~N(mu_t,sigma^2)  on the new treatment arm
#
# Treatment effect theta = mu_t - mu_c
#
# To test H_0: theta =< 0 vs theta > 0 with one-sided type I error
# probability alpha=0.025 and power 1-beta=0.9 at
#
#   theta = delta = 0.4
#
# when sigma^2 = 0.64

alpha=0.025
beta=0.1
sigma=0.8
delta=0.4
ifix=(qnorm(1-alpha)+qnorm(1-beta))^2/delta^2

# Fixed sample size design
# By default, RPACT uses the t-distribution in sample size
# calculations
# We shall assume a known variance, so give the command
#   normalApproximation = TRUE
ss_mean <- getSampleSizeMeans(
  alpha = 0.025,
  beta = 0.1,
  sided = 1,
  alternative = 0.4,
  stDev = 0.8,
  normalApproximation = TRUE
)
summary(ss_mean)

# Group sequential design

# Design a group sequential trial with 5 analyses
#
# Error spending design, non-binding futility boundary,
# rho-family for both spending functions, rho=2 for efficacy
# and futility boundaries
#
# In RPACT the rho-family is called the Kim and DeMets alpha
# spending function and must be specified for both the
# efficacy and the futility boundaries.
#
# Define a design object with 5 analyses and a parameter
# of 2 for the rho family
des_rho5 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.1,
  sided = 1,
  kMax = 5,
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
# See properties of this design, including efficacy
# and futility boundaries on the Z-scale
summary(des_rho5)
# These boundaries on the Z-scale apply to cases with
# power 1-beta at the treatment effect under which type 2
# error is spent.

# The sample size and hence, the absolute information 
# levels, need to be chosen to give power 1-beta
# at our specified alternative theta=delta.
# Based on this group sequential design, we find the
# sample size required (assuming known variance)
ss_gsd <- getSampleSizeMeans(
  design = des_rho5,
  alternative = delta,
  stDev = sigma,
  normalApproximation = T
)

ss_gsd$maxNumberOfSubjects
ss_gsd$maxNumberOfSubjects1
ss_gsd$maxNumberOfSubjects2

# Compare with fixed sample size
ss_gsd$nFixed

# extract design characteristics to get inflation factor
# (and other information)
dc_des_rho5 <- getDesignCharacteristics(des_rho5)

dc_des_rho5$inflationFactor
dc_des_rho5$.design$alphaSpent
dc_des_rho5$.design$betaSpent

# To see everything:
dc_des_rho5

# Information is calculated for power 1-beta at theta=1
# For information needed to give power at theta=delta, type
dc_des_rho5$information/delta^2

# We san also see the inflation factor as
ss_gsd$maxNumberOfSubjects/ss_gsd$nFixed

# We can plot the boundary

plot(des_rho5)

# Probability of stopping early under different values of theta

# store the information max from the design
iMax <- ifix * dc_des_rho5$inflationFactor
iMax

# get the prob. of stopping under theta = 0
asn <- getPowerAndAverageSampleNumber(
  design = des_rho5,
  theta = 0,
  nMax = iMax
)

# Focusing on sub parts, we have probabilities of:

# rejecting the null hypothesis per stage
asn$rejectPerStage

# stopping for futility per stage
asn$futilityPerStage

# overall probability of rejecting the null hypothesis
asn$calculatedPower

# probability of not rejecting the null hypothesis
1 - asn$calculatedPower

# expected number of patients per treatment group
asn$averageSampleNumber

# We can repeat this analysis for theta=delta and
# theta=delta/2
asn <- getPowerAndAverageSampleNumber(
  design = des_rho5,
  theta = c(delta, delta/2),
  nMax = iMax
)

asn$rejectPerStage
asn$futilityPerStage
asn$calculatedPower
asn$averageSampleNumber
# Results headed [,1] are for theta=delta
# Results headed [,2] are for theta=delta/2

# Expected sample size as a function of theta

# While there is an inflation of the maximum number of
# patients compared to a fixed sample design, the
# average sample size is reduced due to the chance
# of early stopping.

# treatment effects
theta.vec <- seq(-0.1, 0.5, by = 0.02)

# store the information max from the design
iMax <- ifix * dc_des_rho5$inflationFactor

# expected information for each theta
asn <- getPowerAndAverageSampleNumber(
  design = des_rho5,
  theta = theta.vec,
  nMax = iMax
)

# expected information as percent of fixed design
einf.percent.vec <- 100 *
  asn$averageSampleNumber / ifix

# plot
plot(
  theta.vec,
  einf.percent.vec,
  ylim = c(0, 100),
  pch = 20,
  main = "E(N) as percent of fixed sample size",
  ylab = "E(N)",
  xlab = "theta"
)

# Example data set

# We explore the consequences of observed sample sizes and
# observed information fractions that deviate from the
# originally planned design

# Suppose we have cumulative sample sizes on treatment
# and control

nt=c(20,40,60,82,95)
nc=c(20,40,60,82,95)

iobs=1/(sigma^2/nt+sigma^2/nc)

# and treatment effect estimates

theta.hat=c(0.10,0.06,0.21,0.31)

# RPACT can calculate the error spending EFFICACY boundary
# stage by stage

# Since the futility boundary is non-binding, we do not
# need too know what that boundary is when calculating
# the efficacy boundary.

# Analysis 1

# We have observed 40 subjects out of a target of
# ss_gsd$maxNumberOfSubjects = 190.434
# So, the observed information fraction is 40/190.434
#
# We need to provide information fractions for later
# analyses. We can just put in the values that were
# originally planned -- these will not affect the
# first boundary value
designUpdate1 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.1,
  informationRates = c(40/ss_gsd$maxNumberOfSubjects,0.4,0.6,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate1)
# We see the efficacy boundary is at Z=3.061

# Continuing to the next analysis
# Now we give the observed sample sizes at analyses 1 and 2
designUpdate2 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.1,
  informationRates = c(c(40,80)/ss_gsd$maxNumberOfSubjects,0.6,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate2)
# We see the second efficacy boundary point is at Z=2.681

# Continuing to the analysis 3
# Now we give the observed sample sizes at analyses 1, 2 and 3
designUpdate3 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.1,
  informationRates = c(c(40,80,120)/ss_gsd$maxNumberOfSubjects,0.8,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate3)
# We see the third efficacy boundary point is at Z=2.436

designUpdate4 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.1,
  informationRates = c(c(40,80,120,164)/ss_gsd$maxNumberOfSubjects,1), 
  typeOfDesign = "asKD",
  gammaA = 2
)
summary(designUpdate4)
# The fourth efficacy boundary point is at Z=2.213

# The final analysis is handled differently to allow a
# final information level that differs from the target, iMax.
#
# In calculating the efficacy boundary, what is important is 
# that the RATIOS of information levels are correct. These
# ratios are preserved if we divide sample sizes by the final
# observed sample size to make the last value in
# informationRates equal to 1.
#
# However, we need to preserve the type I error probabilities
# already seen at earlier analyses and these depended on the
# observed sample sizes. This is achieved by adding the commands
#   typeOfDesign = "asUser",
#   userAlphaSpending = designUpdate4$alphaSpent
# Now, the remainder of alpha will be spent at the final analysis.
designUpdate5 <- getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.1,
  informationRates = c(c(40,80,120,164,190)/190), 
  typeOfDesign = "asUser",
  userAlphaSpending = designUpdate4$alphaSpent
)
summary(designUpdate5)
# The final efficacy boundary point is at Z=2.135

# For more decimal places, see
designUpdate5$criticalValues

# These values agree very nicely with those from CJ's code

# RPACT does not yet provide an error spending FUTILITY
# boundary.

# From CJ's code, the futility boundary values at
# analyses 1 to 4 should be
#
#   -1.0377 0.07157 0.8866 1.6531

# Suppose we type the following:

des_rho5_1 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.1,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[1]/iMax,c(0.4,0.6,0.8,1)),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
summary(des_rho5_1)
des_rho5_1$futilityBounds
# This gives almost what we want for the first futility
# boundary value. However, RPACT has spent type 2 error
# at a value of theta such that the efficacy and
# futility boundaries will meet up at analysis 5,
# spending type 2 error exactly equal to beta=0.1.

# Continuing to analysis 2, etc

des_rho5_2 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.1,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[c(1,2)]/iMax,0.6,0.8,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
summary(des_rho5_2)
des_rho5_2$futilityBounds
# Note the new value for the futility boundary at analysis 1.
# This occurs because the value of theta at which type 2
# error is spent has changed again.

des_rho5_3 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.1,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[c(1,2,3)]/iMax,0.8,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
summary(des_rho5_3)
des_rho5_3$futilityBounds

des_rho5_4 <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.1,
  sided = 1,
  kMax = 5,
  informationRates = c(iobs[c(1,2,3,4)]/iMax,1),
  typeOfDesign = "asKD",
  typeBetaSpending = "bsKD",
  gammaA = 2,
  gammaB = 2
)
summary(des_rho5_4)
des_rho5_4$futilityBounds

# There is no need to calculate a futility boundary
# at analysis 5 as the efficacy boundary values is all
# that matters there.


# We could run the whole sequence of analyses to design
# the trial and process the observed data as follows:

design <- getDesignGroupSequential(
  alpha = 0.025,
  beta = 0.1,
  kMax = 5,
  typeOfDesign = "asKD",
  gammaA = 2,
  typeBetaSpending = "bsKD",
  gammaB = 2,
  bindingFutility = FALSE
)
sampleSizes <- getSampleSizeMeans(
  design = design,
  alternative = 0.4,
  stDev = 0.8,
  normalApproximation = TRUE
)
power <- getPowerMeans(
  design = design,
  alternative = c(-0.1, 0.5),
  stDev = 0.8,
  maxNumberOfSubjects = sampleSizes$maxNumberOfSubjects,
  normalApproximation = TRUE
)
summary(sampleSizes)
summary(power)
plot(sampleSizes, type = 1)
plot(power, type = 9) 
dataObserved <- getDataset(
  cumMeans1 = c(0.1, 0.06, 0.21, 0.31),
  cumMeans2= c(0, 0, 0, 0),
  cumN1 = c(20,40,60,82),
  cumN2 = c(20,40,60,82),
  cumStDev1 = 0.8*c(1, 1, 1, 1),
  cumStDev2 = 0.8*c(1, 1, 1, 1)
)
designDummy <- getDesignGroupSequential(
  alpha = 0.025,
  typeOfDesign = "asKD",
  gammaA = 2
) 
getAnalysisResults(
  # informationEpsilon = 0,  #  can be entered in case of underrunning
  maxInformation = 190,
  design = designDummy,
  dataInput = dataObserved,
  normalApproximation = TRUE
) |> summary()
