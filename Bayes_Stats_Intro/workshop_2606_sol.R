#### Example 1: Disease|Test

## test positive
prior = 0.01
likelihood = 0.99
marginal_likelihood = 0.99*prior + 0.05*(1-prior)

posterior = likelihood*prior/marginal_likelihood
show(posterior)

#### Example 2: Multiple updates
# What's the posterior if second test is negative, third is positive again?
# Use posterior of previous update as prior for next update!

## Second test negative
prior_2 = posterior
likelihood = 1-0.99
marginal_likelihood = (1-0.99)*prior_2 + (1-0.05)*(1-prior_2)

posterior_2 = likelihood*prior_2/marginal_likelihood
show(posterior_2)

## Third test positive
prior_3 = posterior_2
likelihood = 0.99
marginal_likelihood = 0.99*prior_3 + 0.05*(1-prior_3)

posterior_3 = likelihood*prior_3/marginal_likelihood
show(posterior_3)

rm(list=ls())

#### Example 3a: Order relevant?
# Reverse the order of test 2 and 3 - does that make a difference to the final posterior?
# (i.e. first two tests positive, then last negative)

## First test positive
prior = 0.01
likelihood = 0.99
marginal_likelihood = 0.99*prior + 0.05*(1-prior)

posterior = likelihood*prior/marginal_likelihood
show(posterior)

## Second test positive
prior_2 = posterior
likelihood = 0.99
marginal_likelihood = 0.99*prior_2 + 0.05*(1-prior_2)

posterior_2 = likelihood*prior_2/marginal_likelihood
show(posterior_2)

## Third test negative
prior_3 = posterior_2
likelihood = 1-0.99
marginal_likelihood = (1-0.99)*prior_3 + (1-0.05)*(1-prior_3)

posterior_3 = likelihood*prior_3/marginal_likelihood
show(posterior_3)

rm(list=ls())

#### Example 3b: Order relevant II?
# Reverse the order of test 1 and 2 - does that make a difference to the final posterior?
# (i.e. first test negative, then both positive)

## First test negative
prior = 0.01
likelihood = 1-0.99
marginal_likelihood = (1-0.99)*prior + (1-0.05)*(1-prior)

posterior = likelihood*prior/marginal_likelihood
show(posterior)

## Second test positive
prior_2 = posterior
likelihood = 0.99
marginal_likelihood = 0.99*prior_2 + 0.05*(1-prior_2)

posterior_2 = likelihood*prior_2/marginal_likelihood
show(posterior_2)

## Third test positive
prior_3 = posterior_2
likelihood = 0.99
marginal_likelihood = 0.99*prior_3 + 0.05*(1-prior_3)

posterior_3 = likelihood*prior_3/marginal_likelihood
show(posterior_3)

rm(list=ls())
#------------------------ Other examples DISCRETE --------------------------------
# adapted from Kruschke
# example: Is a coin fair?
# estimate: theta; if fair: theta = 0.5

graphics.off()
source("DBDA2E-utilities.R")
source("BernGrid.R")

# discrete example 1
# 1 flip, symmetrical prior
Theta = seq( 0 , 1 , length=5 )  # Sparse teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,0),rep(1,1))      # Single flip with 1 head

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=FALSE , showpD=FALSE )
# saveGraph(file="BernGridExample0",type="eps")
#------------------------------------------------------------------------------

# discrete example 2
# same as before, theta has more possible values
Theta = seq( 0 , 1 , length=11 )  # Sparse teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,0),rep(1,1))      # Single flip with 1 head

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=FALSE , showpD=FALSE )
saveGraph(file="BernGridExample1",type="eps")
#------------------------------------------------------------------------------

#------------------------ Other examples CONTINUOUS --------------------------------
# adapted from Kruschke

# example 3
# flat prior - has no effect!
Theta = seq( 0 , 1 , length=1001 ) # Fine teeth for Theta.
pTheta = rep(1,length(Theta))      # Uniform (horizontal) shape for pTheta.
pTheta = pTheta/sum(pTheta)        # Make pTheta sum to 1.0
Data = c(rep(0,0),rep(1,1))        # Single flip with 1 head

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=FALSE , showpD=FALSE )
saveGraph(file="BernGridExample2",type="eps")
#------------------------------------------------------------------------------

# example 4
# Prior allows for only two values!
Theta = seq( 0 , 1 , length=1001 ) # Fine teeth for Theta.
pTheta = rep(0,length(Theta))      # Only extremes are possible!
pTheta[2] = 1                      # Only extremes are possible!
pTheta[length(pTheta)-1] = 1       
pTheta = pTheta/sum(pTheta)        # Make pTheta sum to 1.0
Data = c(rep(0,0),rep(1,1))        # Single flip with 1 head

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=FALSE , showpD=FALSE )
saveGraph(file="BernGridExample3",type="eps")
#------------------------------------------------------------------------------

# example 5
# small data set - prior is very influential
Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,3),rep(1,1))      # 25% heads, N=4

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )
saveGraph(file="BernGridExample4",type="eps")
#------------------------------------------------------------------------------

# example 6
# small data set with highly precise prior - heavily influences posterior!
Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
pTheta = pTheta^10               # Sharpen pTheta !
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,3),rep(1,1))      # 25% heads, N=4

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )
saveGraph(file="BernGridExample5",type="eps")
#------------------------------------------------------------------------------

# example 7
# Completely uninformative prior - even in small data set posterior dominated by likelihood!
Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
pTheta = pTheta^0.1              # Flatten pTheta !
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,3),rep(1,1))      # 25% heads, N=4

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )
saveGraph(file="BernGridExample6",type="eps")
#------------------------------------------------------------------------------

# example 8
# larger sample - enough to overwhelm prior!
Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,30),rep(1,10))    # 25% heads, N=40

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )
saveGraph(file="BernGridExample7",type="eps")
#------------------------------------------------------------------------------

# example 9
# both precise prior and fairly large sample - posterior is a compromise
Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
pTheta = pTheta^10               # Sharpen pTheta !
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,30),rep(1,10))    # 25% heads, N=40

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )
saveGraph(file="BernGridExample8",type="eps")
#------------------------------------------------------------------------------

# example 10
# vague prior revisited.
Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
pTheta = pTheta^0.1              # Flatten pTheta !
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,30),rep(1,10))    # 25% heads, N=40

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )
saveGraph(file="BernGridExample9",type="eps")
#------------------------------------------------------------------------------

# example 11
# getting a bit silly here.
Theta = seq( 0 , 1 , length=1000 )  # Fine teeth for Theta.
# Two triangular peaks on a small non-zero floor:
pTheta = c( rep(1,200),seq(1,100,length=50),seq(100,1,length=50),rep(1,200) , 
            rep(1,200),seq(1,100,length=50),seq(100,1,length=50),rep(1,200) )
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,13),rep(1,14)) 

openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=FALSE , showpD=FALSE )
saveGraph(file="BernGridExample10",type="eps")
#------------------------------------------------------------------------------

graphics.off()
rm(list=ls())

#--------------------------- Example Bayes Factor --------------------------------------
# update.packages('Matrix')
# update.packages('coda')
# install.packages("BayesFactor")
library(BayesFactor)
data(sleep)

show(sleep)

summary(aov(extra ~ group + Error(ID/group), data = sleep))

anovaBF(extra ~ group + ID, data = sleep, whichRandom = "ID", progress=FALSE)

ttestBF(x = sleep$extra[sleep$group==1], y = sleep$extra[sleep$group==2], paired=TRUE)

#------------------------ Intention of Sampling --------------------------------
# Adapted from Kruschke

#### Exercise sampling intentions
# Assume we have a die and want to test whether probability that six-dotted face comes up is fair.
# Therefore, there are two possible outcomes: six-dots or not six-dots.
# If die is fair, then probability of six-dotted face is 1/6.

# a) Consider intention to stop after 45 throws and we throw 3 six-dots.
# Data:
N = 45 ; z = 3
theta = 1/6

# Consider the low tail because z/N = 3/45 is less than expected p=1/6:
lowTailZ = 0:z

# Cumulative low tail probability given by binomial distribution:
lowTailP = sum( choose(N,lowTailZ) * theta^lowTailZ * (1-theta)^(N-lowTailZ) )

# Two-tail probability:
TwoTailP = 2 * lowTailP
show(TwoTailP)

# b) Consider intention to stop after 3 six-dots, which takes 45 throws.
# Cumulative low tail probability is now given by negative binomial
complN = z:(N-1)

complP = sum( (z/complN) * choose(complN,z) * theta^z * (1-theta)^(complN-z) )
lowTailP = 1-complP

# Two-tail probability:
TwoTailP = 2 * lowTailP
show(TwoTailP)

#### Exercise 2 sampling intentions: extention to confidence intervals
# confidence interval as range of parameters that would not be rejected
# a) Consider intention to stop after 45 throws and we throw 3 six-dots.
N = 45 ; z = 3
theta = 1/6 
# z/N = 3/45 = 0.06666. 
# For candidate theta values from 0.170 to 0.190, which are greater than z/N observed, 
# compute the left-tail p value:
lowTailZ = 0:z
for ( theta in seq(0.170 , 0.190 , 0.001) ) {
  show( c(theta,2*sum( choose(N,lowTailZ) * theta^lowTailZ * (1-theta)^(N-lowTailZ) )))
}
# The columns are theta and p value, highest theta not rejected: 0.182

# For candidate theta values from 0.005 to 0.020, which are less than z/N observed, # compute the right-tail p value:
highTailZ = z:N
for ( theta in seq(0.005 , 0.020 , 0.001) ) {
  show( c(theta,2*sum( choose(N,highTailZ) * theta^highTailZ * (1-theta)^(N-highTailZ) )))
}
# The columns are theta and p value, lowest theta: 0.014

# Therefore CI = [0.014,0.182]

# b) Consider intention to stop after 3 six-dots, which takes 45 throws.
# For candidate theta values GREATER than z/N observed, compute the LEFT-tail p
# value:
complN = z:(N-1)
for ( theta in seq(0.150 , 0.160 , 0.001) ) {
  show( c(theta,2*(1-sum( (z/complN) * choose(complN,z) * theta^z * (1-theta)^(complN-z)))))
}
# highest theta not rejected: 0.154

# For candidate theta values LESS than z/N observed, compute the RIGHT-tail p
# value:
highTailN = z:N # Notice N not N-1
for ( theta in seq( 0.005 , 0.020 , 0.001) ) {
  show( c(theta,2*sum((z/highTailN) * choose(highTailN,z) * theta^z * (1-theta)^(highTailN-z))))
}
# lowest theta not rejected: 0.014

# Therefore CI = [0.014,0.154]

graphics.off()
rm(list=ls())

#---------------------- Difference Model comaprison/estimation ------------------------------
source("DBDA2E-utilities.R")
source("BernGrid.R")

#### Use previous example, compute Bayes Factor:
source("BernBeta.R")
z = 7 ; N = 24
theta = 0.5
pDgTheta = theta^z * (1-theta)^(N-z)
print( pDgTheta )
# Result is 5.96e-08

#### Use previous example, compute ROPE (=region of practical equivalence) and HDI (h)
# a = 1 ; b = 1
# a = 20 ; b = 5
# a = 5 ; b = 20
a = 2000 ; b = 2000
openGraph(width=5,height=7)
BernBeta( c(a,b) , c(rep(0,N-z),rep(1,z)) , ROPE=c(0.48,0.52) ,
          plotType="Bars" , showCentTend="Mode" , showHDI=TRUE , showpD=TRUE )
# Result is p(D)=6.02e-08, displayed in graph

# Parameter estimation seems more informative than model comparison because
# estimation gives explicit information about the magnitude of the parameter
# that describes the data, while the BF only indicates which prior distribution
# (model index) gets more credibility without saying anything about the
# parameter magnitude. 

graphics.off()
rm(list=ls())

#--------------------------- Proper Modelling MCMC --------------------------------------
### To run this program, install the Bayesian MCMC sampling program JAGS from http://mcmc-jags.sourceforge.net/
### Also, type install.packages("rjags")

# Get the functions loaded into R's working memory:
source("BEST.R")

# Specify data as vectors: 
y1 = c(101,100,102,104,102,97,105,105,98,101,100,123,105,103,100,95,102,106,
       109,102,82,102,100,102,102,101,102,102,103,103,97,97,103,101,97,104,
       96,103,124,101,101,100,101,101,104,100,101)
y2 = c(99,101,100,101,102,100,97,101,104,101,102,102,100,105,88,101,100,
       104,100,100,100,101,102,103,97,101,101,100,101,99,101,100,100,
       101,100,99,101,100,102,99,100,99)

# Run the Bayesian analysis:
mcmcChain = BESTmcmc( y1 , y2 ) 

postInfo = BESTplot( y1 , y2 , mcmcChain , pairsPlot=TRUE )
# Show detailed summary info on console:
show( postInfo ) 
# You can save the plot(s) using the pull-down menu in the R graphics window,
# or by using the following:
# saveGraph( file="BESTexample" , type="eps" )
# saveGraph( file="BESTexample" , type="jpeg" )

# Save the data and results for future use:
save( y1, y2, mcmcChain, postInfo, file="BESTexampleMCMC.Rdata" )
# To re-load the saved data and MCMC chain, type: load( "BESTexampleMCMC.Rdata" ) 

# compare to t-test:
t.test(y1,y2,paired=FALSE)

# compare to Bayes factor:
ttestBF(y1,y2,paired=FALSE)
#-------------------------------------------------------------------------------

#--------------------------- Even more Proper Modelling MCMC with pre-defined ROPEs --------------------------------------
#### Example 1 Alcohol preferences of sexually deprived flies by Shohat-Ophir et al. (2012)

# "One cohort, rejected-isolated, was subjected to courtship conditioning; they experiences 1h sessions of 
# sexual rejeciton by mated females, three times a day, for 4 days. [...] FLies in the mated-group cohort experienced
# 6h sessions of mating with multiple receptive virign females (ratio 1:5) for 4 days. Flies from each cohort were then tested 
# in a two-choice preference assay, in which they voluntarily choose to consume food with or without 15% ethanol supplementation."
# 
# a) Test for a difference in preference ratio between groups.

myDataFrame = read.csv( file="ShohatOphirKAMH2012dataReduced.csv" )
xName="Group"
yName="PreferenceIndex"
fileNameRoot="ShohatOphirKAMH2012data-PI-"
RopeMuDiff=c(-0.1,0.1) ; RopeSdDiff=c(-0.1,0.1) ; RopeEff=c(-0.1,0.1)

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom2grp-MrobustHet.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC(datFrm=myDataFrame , yName=yName , xName=xName ,
                   numSavedSteps=50000 , saveName=fileNameRoot)
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , RopeMuDiff=RopeMuDiff , 
                        RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xName=xName , 
          RopeMuDiff=RopeMuDiff , RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
          pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

# b) Test for a difference in absolute amount food consumed between groups.
myDataFrame = read.csv( file="ShohatOphirKAMH2012dataReduced.csv" )
xName="Group"
yName="GrandTotal"
fileNameRoot="ShohatOphirKAMH2012data-GT-"
RopeMuDiff=c(-0.1,0.1) ; RopeSdDiff=c(-0.1,0.1) ; RopeEff=c(-0.1,0.1)

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom2grp-MrobustHet.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=50000 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , RopeMuDiff=RopeMuDiff , 
                        RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC(mcmcCoda , datFrm=myDataFrame , yName=yName , xName=xName , 
         RopeMuDiff=RopeMuDiff , RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
         pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType)
#------------------------------------------------------------------------------- 

#--------------------------- Proper Modelling MCMC --------------------------------------
#### Example 2 Berger et al. (1988)

# a) Does restricting the diet change the lifespan and/or variance of lifespan of rats?
# Also check for outliers
myDataFrame = read.csv( file="RatLives.csv" )
xName="Group"
yName="DaysLive"
fileNameRoot = "RatLives-"
RopeMuDiff=c(-10,10) ; RopeSdDiff=c(-10,10) ; RopeEff=c(-0.1,0.1)

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom2grp-MrobustHet.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=50000 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , RopeMuDiff=RopeMuDiff , 
                        RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xName=xName , 
          RopeMuDiff=RopeMuDiff , RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
          pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

# b) Repeat analysis but account for scew - values are squared to get rid of skew to left
myDataFrame = read.csv( file="RatLives.csv" )
xName="Group"
myDataFrame = cbind( myDataFrame , DaysLiveSq = myDataFrame$DaysLive^2 )
yName="DaysLiveSq"
fileNameRoot = "RatLives-DaySq-"
RopeMuDiff=c(-100,100) ; RopeSdDiff=c(-100,100) ; RopeEff=c(-0.1,0.1)

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom2grp-MrobustHet.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=50000 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , RopeMuDiff=RopeMuDiff , 
                        RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xName=xName , 
          RopeMuDiff=RopeMuDiff , RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
          pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

#--------------------------- Bayesian modelling of belief updating --------------------------------------
# Setting:
# Cues are 85% valid
# There are visual and auditory cues
# given cues and outcome of a trial, subjects have to infer whether currently shape or tone matters

# Trial 1
# assume flat prior over vision/sound <=> p(vision,sound)=c(0.5,0.5)
# assume bad visual cue, good auditory cue
# assume outcome is good, thus likelihood p(outcome|vision) = 0.15, p(outcome|sound)=0.85

prior = c(0.5,0.5) # c(vision,sound)
likelihood = c(0.15,0.85)
evidence = sum(prior*likelihood)

posterior_1 = likelihood*prior/evidence
show(posterior_1)
# shift toward auditory

# Trial 2
# posterior of trial 1 is prior of trial 2
# assume bad visual and bad auditory cue and bad outcome
# thus, likelihood for both is 0.85

prior = posterior_1# c(vision,sound)
likelihood = c(0.85,0.85)
evidence = sum(prior*likelihood)

posterior_2 = likelihood*prior/evidence
show(posterior_2)
# nothing changes, because trial not informative & unsurprising

# Trial 3
# posterior of trial 2 is prior of trial 3
# assume good visual and good auditory cue and bad outcome
# thus, likelihood for both is 0.15

prior = posterior_2# c(vision,sound)
likelihood = c(0.15,0.15)
evidence = sum(prior*likelihood)

posterior_3 = likelihood*prior/evidence
show(posterior_3)
# nothing changes, because trial not informative

# Trial 4
# posterior of trial 3 is prior of trial 4
# assume good visual and bad auditory cue and bad outcome
# thus, likelihood p(outcome|vision) = 0.15, p(outcome|sound)=0.85 again

prior = posterior_3# c(vision,sound)
likelihood = c(0.15,0.85)
evidence = sum(prior*likelihood)

posterior_4 = likelihood*prior/evidence
show(posterior_4)
# credibility of auditory relevance increases

# Trial 5
# posterior of trial 4 is prior of trial 5
# assume good visual and bad auditory cue and good outcome
# thus, likelihood p(outcome|vision) = 0.85, p(outcome|sound)=0.15 (conflicting prior belief)

prior = posterior_4# c(vision,sound)
likelihood = c(0.85,0.15)
evidence = sum(prior*likelihood)

posterior_5 = likelihood*prior/evidence
show(posterior_5)
# credibility of auditory decreases - back to level before

# [...]
# imaging regressor: shift in beliefs

library(entropy)
-KL.plugin(posterior_1,posterior_2)
-KL.plugin(posterior_2,posterior_3)
KL.plugin(posterior_3,posterior_4)
KL.plugin(posterior_4,posterior_5)
