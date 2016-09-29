library(R2jags)
setwd("F:/research stuff/FS_PostDoc/outside_consult/JuliaBeetle_Chap4")  # Set workspace to location with 'Data_compiled' workspace.
#setwd("~/Documents/Rdata/JuliaBeetle_Chap4")
load("Data_compiled.RData")

# Extract and scale (z-score) covariate inputs #
rows <- c("QMD","EarlInf","MidInf","snag")
cols <- c("mean","sd")
zscore.factors <- matrix(NA,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols))
rm(rows,cols)
zscore.factors["QMD",] <- c(mean(Plot.data$QMD),sd(Plot.data$QMD))
zscore.factors["snag",] <- c(mean(Plot.data$snag),sd(Plot.data$snag))
zscore.factors["EarlInf",] <- c(mean(as.matrix(Plot.data[,c("EarlInf_2014","EarlInf_2014","EarlInf_2015","EarlInf_2016")])),
                                   sd(as.matrix(Plot.data[,c("EarlInf_2014","EarlInf_2014","EarlInf_2015","EarlInf_2016")])))
zscore.factors["MidInf",] <- c(mean(as.matrix(Plot.data[,c("MidInf_2014","MidInf_2014","MidInf_2015","MidInf_2016")])),
                            sd(as.matrix(Plot.data[,c("MidInf_2014","MidInf_2014","MidInf_2015","MidInf_2016")])))

snag <- (Plot.data$snag - zscore.factors["snag","mean"])/zscore.factors["snag","sd"]
QMD <- (Plot.data$QMD - zscore.factors["QMD","mean"])/zscore.factors["QMD","sd"]

EarlInf <- # Note this is a matrix because values change across years.
  as.matrix(cbind(Plot.data[,"EarlInf_2014"],Plot.data[,c("EarlInf_2014","EarlInf_2015","EarlInf_2016")]))
EarlInf <- (EarlInf - zscore.factors["EarlInf","mean"])/zscore.factors["EarlInf","sd"]

MidInf <- # Note this is a matrix because values change across years.
  as.matrix(cbind(Plot.data[,"MidInf_2014"],Plot.data[,c("MidInf_2014","MidInf_2015","MidInf_2016")]))
MidInf <- (MidInf - zscore.factors["MidInf","mean"])/zscore.factors["MidInf","sd"]

##---Analysis of occupancy (assumes closure between visits within a year)---##
Y <- apply(Y.arry,c(1,2,3),function(x) 1*(x>0)) # Collapse counts to detection-nondetection data
Y <- # Collapse further to represent the number of visits when species detected for each plot X year occassion.
  apply(Y,c(1,2),function(x) sum(x,na.rm=T))
Y[which(is.element(substr(plot,1,2),c("1a","2a","1b","2b"))),1] <- NA # Sites 1 and 2 not surveyed in 2013.

nplot <- length(plot)
nyear <- length(year)

Z.init <- apply(Y,c(1,2),function(x) sum(x>0)) # Initial values for occupancy state parameter.

# Assemble the data names list for JAGS #
data <- list("Y","nplot","nyear","n.visits","snag","QMD","EarlInf","MidInf")

# Assemble the initial values for JAGS.  Here z is the apparent occupancy state 

inits <- function(){list(Z=Z.init)}

# Assemble the parameters vector for JAGS (What we want to track).
parameters <- c("p","beta0","beta1","beta.EarlInf","beta.MidInf","beta.snag","beta.QMD","test")

sink("model.txt")
cat("

# The model.
model{

# Prior for detection probability
p ~ dunif(0,1) 
  
# Priors for occupancy model parameters.
beta0 ~ dnorm(0,0.01)T(-10,10) # mean probability of occupancy for sites previously unoccupied (= colonization or gamma)
beta1 ~ dnorm(0,0.01)T(-10,10) # Persistence offset (added to beta0 if occupied in previous year) 
# colonization probability (gamma in literature) = expit(beta0)
# persistence probability (phi in literature) = expit(beta0 + beta1)

beta.QMD ~ dnorm(0,0.01)T(-10,10)
beta.EarlInf ~ dnorm(0,0.01)T(-10,10)
beta.MidInf ~ dnorm(0,0.01)T(-10,10)
beta.snag ~ dnorm(0,0.01)T(-10,10)
  
#Model
psi0 ~ dunif(0,1) # Probability of occupancy in year prior to sampling   
for(j in 1:nplot){  # Loop for year 1 - uses z0 to calculate persistence offset 
  z0[j] ~ dbern(psi0) # Occupancy state in year prior to sampling
  ###Model for occupancy probability in year 1 ###
  logit(psi[j,1]) <- beta0 + beta1*z0[j] + beta.EarlInf*EarlInf[j,1] + beta.MidInf*MidInf[j,1] + 
    beta.snag*snag[j] + beta.QMD*QMD[j]
  ###Model for detection probability ###
  Z[j,1] ~ dbern(psi[j,1]) # Occupancy state
  pZ[j,1] <- p*Z[j,1]  # Unconditional probability of detection
  Y[j,1] ~ dbin(pZ[j,1],n.visits[1]) # Data model: Data ~ binomial(trials = number of visits)

  #____________Bayesian GOF_________________________________
  ynew[j,1] ~ dbin(pZ[j,1],n.visits[1])  #simulated new data y under model
    
  LLsim[j,1] <- (ynew[j,1]*log(p)+
    (n.visits[1]-ynew[j,1])*log(1-p))*Z[j,1]  #log-likelihood simulated data
  LLdata[j,1]<- (Y[j,1]*log(p)+
    (n.visits[1]-Y[j,1])*log(1-p))*Z[j,1]    #log-likelihood observed data
  #_________________________________________________________
  for(t in 2:nyear){  # Loop through remaining years 2 and up - uses Z[,,t-1] to calculate persistence offset
    ###Model for occupancy probability ###
    logit(psi[j,t]) <- beta0 + beta1*Z[j,(t-1)] + beta.EarlInf*EarlInf[j,t] + beta.MidInf*MidInf[j,t] + 
      beta.snag*snag[j] + beta.QMD*QMD[j]
    ###Model for detection probability ###
    Z[j,t] ~ dbern(psi[j,t]) # Occupancy state
    pZ[j,t] <- p*Z[j,t] # Unconditional probability of detection
    Y[j,t] ~ dbin(pZ[j,t],n.visits[t]) # Data model
    #____________Bayesian GOF_________________________________
    ynew[j,t] ~ dbin(pZ[j,t],n.visits[t])  #simulated new data y under model
    
    LLsim[j,t] <- (ynew[j,t]*log(p)+
      (n.visits[t]-ynew[j,t])*log(1-p))*Z[j,t]  #log-likelihood simulated data
    LLdata[j,t]<- (Y[j,t]*log(p)+
      (n.visits[t]-Y[j,t])*log(1-p))*Z[j,t]    #log-likelihood observed data
    #_________________________________________________________
    }
  }

#________Bayesian GOF___________
#deviance
dev_sim <- (-2)*sum(LLsim[,])
dev_data <- (-2)*sum(LLdata[,])

#test statistics should be ~0.5 if model fits
test<-step(dev_data-dev_sim)
#________________________________

}
",fill=TRUE)
sink()

# MCMC values.
nc <- 4
nb <- 5000
ni <- 25000
nt <- 1

# Send it all to JAGS and hope for the best!

# To help track time.
starttime <- Sys.time()

bugout <- jags(data, inits, parameters, model.file="model.txt", n.chains=nc, n.iter=ni,
n.burnin=nb, n.thin=nt)

endtime <- Sys.time()
runtime <- endtime - starttime
runtime  # Total time it took to run model.

#Check n.effectives and R.hats for parameters
length(which(bugout$BUGSoutput$summary[,"n.eff"]<100))/length(bugout$BUGSoutput$summary[,"n.eff"])
min(bugout$BUGSoutput$summary[,"n.eff"])
sort(bugout$BUGSoutput$summary[,"n.eff"])[1:50]
max(bugout$BUGSoutput$summary[,"Rhat"])
sort(bugout$BUGSoutput$summary[,"Rhat"],decreasing=T)[1:50]

#use bugout object to manipulate results in R environment

#save results
library(R.utils)
saveObject(bugout,"Mod_Pers_LogInfestSnagQMD")
