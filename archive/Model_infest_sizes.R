library(R2jags)
setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/JuliaH_dissert/Chap4")  # Set workspace to location with 'Data_compiled' workspace.
#setwd("~/Documents/Rdata/JuliaBeetle_Chap4")
load("Data_compiled.RData")

# Extract and scale (z-score) covariate inputs #
rows <- c("EarlInf_lt23","EarlInf_23to38","EarlInf_gt38")
cols <- c("mean","sd")
zscore.factors <- matrix(NA,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols))
rm(rows,cols)
zscore.factors["EarlInf_lt23",] <- c(mean(as.matrix(Plot.data[,c("EarlInf_lt23_2014","EarlInf_lt23_2014","EarlInf_lt23_2015","EarlInf_lt23_2016")])),
                                   sd(as.matrix(Plot.data[,c("EarlInf_lt23_2014","EarlInf_lt23_2014","EarlInf_lt23_2015","EarlInf_lt23_2016")])))
zscore.factors["EarlInf_23to38",] <- c(mean(as.matrix(Plot.data[,c("EarlInf_23to38_2014","EarlInf_23to38_2014","EarlInf_23to38_2015","EarlInf_23to38_2016")])),
                                     sd(as.matrix(Plot.data[,c("EarlInf_23to38_2014","EarlInf_23to38_2014","EarlInf_23to38_2015","EarlInf_23to38_2016")])))
zscore.factors["EarlInf_gt38",] <- c(mean(as.matrix(Plot.data[,c("EarlInf_gt38_2014","EarlInf_gt38_2014","EarlInf_gt38_2015","EarlInf_gt38_2016")])),
                                     sd(as.matrix(Plot.data[,c("EarlInf_gt38_2014","EarlInf_gt38_2014","EarlInf_gt38_2015","EarlInf_gt38_2016")])))

EarlInf_lt23.x <- as.matrix(cbind(Plot.data[,"EarlInf_lt23_2014"],Plot.data[,c("EarlInf_lt23_2014","EarlInf_lt23_2015",
                                                                               "EarlInf_lt23_2016")]))
EarlInf_lt23.z <- (EarlInf_lt23.x - zscore.factors["EarlInf_lt23","mean"])/zscore.factors["EarlInf_lt23","sd"]

EarlInf_23to38.x <- as.matrix(cbind(Plot.data[,"EarlInf_23to38_2014"],Plot.data[,c("EarlInf_23to38_2014","EarlInf_23to38_2015",
                                                                               "EarlInf_23to38_2016")]))
EarlInf_23to38.z <- (EarlInf_23to38.x - zscore.factors["EarlInf_23to38","mean"])/zscore.factors["EarlInf_23to38","sd"]

EarlInf_gt38.x <- as.matrix(cbind(Plot.data[,"EarlInf_gt38_2014"],Plot.data[,c("EarlInf_gt38_2014","EarlInf_gt38_2015",
                                                                               "EarlInf_gt38_2016")]))
EarlInf_gt38.z <- (EarlInf_gt38.x - zscore.factors["EarlInf_gt38","mean"])/zscore.factors["EarlInf_gt38","sd"]

dimnames(EarlInf_lt23.x)[[2]] <- dimnames(EarlInf_23to38.x)[[2]] <- dimnames(EarlInf_gt38.x)[[2]] <-
dimnames(EarlInf_lt23.z)[[2]] <- dimnames(EarlInf_23to38.z)[[2]] <- dimnames(EarlInf_gt38.z)[[2]] <-
  c("2013","2014","2015","2016")

# Correlation matrix #
cor(cbind(EarlInf_lt23=as.numeric(EarlInf_lt23.x),EarlInf_23to38=as.numeric(EarlInf_23to38.x),
          EarlInf_gt38=as.numeric(EarlInf_gt38.x)))

cor(cbind(EarlInf_lt23=as.numeric(EarlInf_lt23.x),EarlInf_gt23=as.numeric(EarlInf_23to38.x + EarlInf_gt38.x)))

###---Analysis of occupancy (assumes closure between visits within a year)---###
Y <- apply(Y.arry,c(1,2,3),function(x) 1*(x>0)) # Collapse counts to detection-nondetection data
Y <- # Collapse further to represent the number of visits when species detected for each plot X year occassion.
  apply(Y,c(1,2),function(x) sum(x,na.rm=T))
Y[which(is.element(substr(plot,1,2),c("1a","2a","1b","2b"))),1] <- NA # Sites 1 and 2 not surveyed in 2013.

# Remove observations for pointxyear occassions following clear-cut (optional) #
#Y[which(is.element(substr(plot,1,2),c("1a7","1a8","1a9"))),c(3,4)] <- NA
#Y[which(plot=="2a5"),4] <- NA

nplot <- length(plot)
nyear <- length(year)

Z.init <- apply(Y,c(1,2),function(x) sum(x>0)) # Initial values for occupancy state parameter.

# Assemble the data names list for JAGS #
data <- list("Y","nplot","nyear","n.visits","EarlInf_lt23.z","EarlInf_23to38.z","EarlInf_gt38.z")

# Assemble the initial values for JAGS.  Here z is the apparent occupancy state 

inits <- function(){list(Z=Z.init)}

# Assemble the parameters vector for JAGS (What we want to track).
parameters <- c("p","beta0","beta1","beta.lt23","beta.23to38","beta.gt38","test","psi.fs")

sink("model.txt")
cat("

# The model.
model{

# Prior for detection probability
p ~ dunif(0,1) 
  
# Priors for occupancy model parameters.
beta0 ~ dnorm(0,0.01)T(-30,30) # Mean logit colonization probability (previously unoccupied sites)
beta1 ~ dnorm(0,0.01)T(-30,30) # Mean offset for persistence probability (previously occupied sites)

beta.lt23 ~ dnorm(0,0.01)T(-30,30)
beta.23to38 ~ dnorm(0,0.01)T(-30,30)
beta.gt38 ~ dnorm(0,0.01)T(-30,30)

#Model
psi0 ~ dunif(0,1)
for(j in 1:nplot){  # Loop for year 1 - uses z0 to calculate persistence offset
  Z0[j] ~ dbern(psi0)
  ###Model for occupancy probability ###
  logit(psi[j,1]) <- beta0 + beta1*Z0[j] + beta.lt23*EarlInf_lt23.z[j,1] + beta.23to38*EarlInf_23to38.z[j,1] +
    beta.gt38*EarlInf_gt38.z[j,1]
  ###Model for detection probability ###
  Z[j,1] ~ dbern(psi[j,1]) # Occupancy state
  pZ[j,1] <- p*Z[j,1] # Unconditional probability of detection
  Y[j,1] ~ dbin(pZ[j,1],n.visits[1]) # Data model
  #____________Bayesian GOF_________________________________
  ynew[j,1] ~ dbin(pZ[j,1],n.visits[1])  #simulated new data y under model
    
  LLsim[j,1] <- (ynew[j,1]*log(p)+
    (n.visits[1]-ynew[j,1])*log(1-p))*Z[j,1]  #log-likelihood simulated data
  LLdata[j,1]<- (Y[j,1]*log(p)+
    (n.visits[1]-Y[j,1])*log(1-p))*Z[j,1]    #log-likelihood observed data
  #_________________________________________________________
  for(t in 2:nyear){  # Loop through remaining years 2 and up - uses Z[,,t-1] to calculate persistence offset
    ###Model for occupancy probability ###
    logit(psi[j,t]) <- beta0 + beta1*Z[j,(t-1)] + beta.lt23*EarlInf_lt23.z[j,t] + beta.23to38*EarlInf_23to38.z[j,t] +
      beta.gt38*EarlInf_gt38.z[j,t]
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

psi.fs <- sum(Z[,])/(nyear*nplot)
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
nb <- 10000
ni <- 1010000
nt <- 5

# Send it all to JAGS and hope for the best!

# To help track time.
starttime <- Sys.time()

bugout <- jags(data, inits, parameters, model.file="model.txt", n.chains=nc, n.iter=ni,
n.burnin=nb, n.thin=nt)
#ni <- 300000
#bugout <- update(bugout,n.iter = ni)

endtime <- Sys.time()
runtime <- endtime - starttime
runtime  # Total time it took to run model.

#use bugout object to manipulate results in R environment

#save results
library(R.utils)
saveObject(bugout,"Mod_Pers_EInf_X_size")
