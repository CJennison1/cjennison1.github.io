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

# Load functions to be called later

source("~/My files/deming_dec_2025_demos/standards.R")
source("~/My files/deming_dec_2025_demos/ersptest.R")
source("~/My files/deming_dec_2025_demos/findroot.R")

alpha=0.025
beta=0.1
sigma=0.8
delta=0.4

# Design a fixed sample size trial

ifix=(qnorm(1-alpha)+qnorm(1-beta))^2/delta^2
ifix

# nfix= number of observations per treatment

nfix=ifix*(2*sigma^2)
nfix

# ifix=65.67
# nfix=84.06, which rounds up to 85 observations per treatment

# Design a group sequential trial with 5 analyses
#
# Error spending design, non-binding futility boundary, rho-family
# error spending function, rho=2 for efficacy and futility boundaries
#
# Call function find_rfac_ersp1b to find the inflation factor
  
r=16         #  Gives high accuracy for numerical computations
na=5         #  Number of analyses
alpha=0.025  #  Type I error rate
beta=0.1     #  Type II error rate
rho=2        #  Determines rate of error spending

rfac=find_rfac_ersp1b(r,na,alpha,beta,rho)
rfac

# Required inflation factor is rfac=1.133

imax=rfac*ifix
imax

nmax=imax*(2*sigma^2)
nmax

# imax=95.21712
# Maximum sample size is just over 95 per treatment

# Inspect this design with cumulative information levels as planned

iobs=imax*c(1:5)/5
iobs

# Define cumulative error probabilities cume
#
# cume[2;1:na]  cumulative type 1 error probabilities under theta=0
# cume[1;1:na]  cumulative type 2 error probabilities under theta=delta
#
# If we observe data with iobs[na] > imax, we cap type I error spent at
# alpha and the type II error spent at beta
#
# If we observe data with iobs[na] < imax, we set the final type I error
# to be alpha, so cume[2,na=alpha] 

rho=2
cume=array(0,c(2,na))
cume[2,1:na]=pmin(alpha*(iobs/imax)^rho,alpha)
cume[1,1:na]=pmin(beta*(iobs/imax)^rho,beta)
cume

out=ersp1b(r,na,iobs,delta,cume)
out  # Boundary, type I error probability, type II error probability

zbdy=out[[1]]
zbdy

plot(c(1:na),zbdy[1,],col="blue",ylim=c(-2,4),xlim=c(0.5,5.2),pch=20,
     main="Error spending boundary",ylab="Z",xlab="Analysis")
points(c(1:na),zbdy[2,],col="red",pch=20)

# Under theta=0

theta=0
out=gst1(r,na,iobs,zbdy,theta)
out

# Output: pstop(1:3,1:na)  pstop(1,k)=P(Trial stops by crossing the lower
#                                       boundary at analysis k)
#                          pstop(2,k)=P(Trial stops by crossing the upper
#                                       boundary at analysis k)
#                          pstop(3,k)=pstop(1,k)+pstop(2,k)
#         pu               P(Exit by upper boundary)
#         pl               P(Exit by lower boundary)
#         einf             E(information on termination)

einf=out[[4]]
100*einf/ifix # Expected sample size as a percentage of the fixed sample size

# Under theta=delta

theta=delta
out=gst1(r,na,iobs,zbdy,theta)
out
einf=out[[4]]
100*einf/ifix

# Under theta=delta/2

theta=delta/2
out=gst1(r,na,iobs,zbdy,theta)
out
einf=out[[4]]
100*einf/ifix

theta.vec=seq(-0.1,0.5,0.02)
ntheta=length(theta.vec)
einf.percent.vec=rep(NA,ntheta)
for(ii in c(1:ntheta))
{
 theta=theta.vec[ii]
 einf.percent.vec[ii]=100*gst1(r,na,iobs,zbdy,theta)[[4]]/ifix
}
plot(theta.vec,einf.percent.vec,ylim=c(0,100),pch=20,
     main="E(N) as percent of fixed sample size",ylab="E(N)",xlab="theta")

# Running the trial with observed group sizes

# Suppose we have
#
# Cumulative sample sizes on treatment and control

nt=c(20,40,60,82,95)
nc=c(20,40,60,82,95)

# and treatment effect estimates

theta.hat=c(0.10,0.06,0.21,0.31)

# Then observed information levels are

iobs=1/(sigma^2/nt+sigma^2/nc)
iobs

# and Z-values are

zobs=theta.hat*sqrt(iobs[1:4])
zobs

# Error is spent according to the observed information iobs

cume[2,1:na]=pmin(alpha*(iobs/imax)^rho,alpha)
cume[2,na]=alpha
cume[1,1:na]=pmin(beta*(iobs/imax)^rho,beta)

# Giving the boundary

out=ersp1b(r,na,iobs,delta,cume)
zbdy=out[[1]]
zbdy

plot(c(1:na),zbdy[1,],col="blue",ylim=c(-2,4),xlim=c(0.5,5.2),pch=20,
     main="Error spending boundary",ylab="Z",xlab="Analysis")
points(c(1:na),zbdy[2,],col="red",pch=20)

# Comment in ersp1b:
#
# With too little information: The lower and upper boundaries do not
# meet up. To use these in practice, take zbdy[na,2] as the final
# critical value, and ignore zbdy[1,na].

# Add the observed data to the plot

points(c(1:4),zobs,col="dark green",pch=8,cex=0.5)

# The efficacy boundary would have been crossed at analysis 4

# At that point the boundary would have been calculated as

rho=2
cume4=array(0,c(2,4))
cume4[2,1:4]=pmin(alpha*(iobs[1:4]/imax)^rho,alpha)
cume4[1,1:4]=pmin(beta*(iobs[1:4]/imax)^rho,beta)
cume4

iobs4=iobs[1:4]
out=ersp1b(r,4,iobs4,delta,cume4)
zbdy4=out[[1]]
zbdy4

plot(c(1:4),zbdy4[1,],col="blue",ylim=c(-2,4),xlim=c(0.5,5.2),pch=20,
     main="Error spending boundary",ylab="Z",xlab="Analysis")
points(c(1:4),zbdy4[2,],col="red",pch=20)
points(c(1:4),zobs[1:4],col="dark green",pch=8,cex=0.5)

# With a binding futility boundary
#
# Still with imax=74.39,alpha=0.025, beta=0.1, delta=0.4
#
# Cumulative error rates are as before:

rho=2
cume=array(0,c(2,na))
cume[2,1:na]=pmin(alpha*(iobs/imax)^rho,alpha)
cume[1,1:na]=pmin(beta*(iobs/imax)^rho,beta)
cume[2,na]=alpha

out=ersp1c(r,na,iobs,delta,cume)
zbdy_binding_futility=out[[1]]
zbdy_binding_futility

# Compare with the non-binding boundary

zbdy
