library(R2jags)
setwd("F:/research stuff/FS_PostDoc/outside_consult/JuliaBeetle_Chap4")  # Set workspace to location with 'Data_compiled' workspace.
#setwd("~/Documents/Rdata/JuliaBeetle_Chap4")
load("Data_compiled.RData")

# Extract and scale (z-score) covariate inputs #
rows <- c("QMD","EarlInf","MidInf","snag")
cols <- c("mean","sd")
zscore.factors <- matrix(NA,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols))
rm(rows,cols)
zscore.factors["QMD",] <- c(mean(as.matrix(Plot.data[,c("QMD_2014","QMD_2014","QMD_2015","QMD_2016")])),
                                   sd(as.matrix(Plot.data[,c("QMD_2014","QMD_2014","QMD_2015","QMD_2016")])))
zscore.factors["snag",] <- c(mean(as.matrix(Plot.data[,c("snag_2014","snag_2014","snag_2015","snag_2016")])),
                                   sd(as.matrix(Plot.data[,c("snag_2014","snag_2014","snag_2015","snag_2016")])))
zscore.factors["EarlInf",] <- c(mean(as.matrix(Plot.data[,c("EarlInf_2014","EarlInf_2014","EarlInf_2015","EarlInf_2016")])),
                                   sd(as.matrix(Plot.data[,c("EarlInf_2014","EarlInf_2014","EarlInf_2015","EarlInf_2016")])))
zscore.factors["MidInf",] <- c(mean(as.matrix(Plot.data[,c("MidInf_2014","MidInf_2014","MidInf_2015","MidInf_2016")])),
                            sd(as.matrix(Plot.data[,c("MidInf_2014","MidInf_2014","MidInf_2015","MidInf_2016")])))

EarlInf.x <- as.matrix(cbind(Plot.data[,"EarlInf_2014"],Plot.data[,c("EarlInf_2014","EarlInf_2015","EarlInf_2016")]))
EarlInf.z <- (EarlInf.x - zscore.factors["EarlInf","mean"])/zscore.factors["EarlInf","sd"]

MidInf.x <- as.matrix(cbind(Plot.data[,"MidInf_2014"],Plot.data[,c("MidInf_2014","MidInf_2015","MidInf_2016")]))
MidInf.z <- (MidInf.x - zscore.factors["MidInf","mean"])/zscore.factors["MidInf","sd"]

snag.x <- as.matrix(cbind(Plot.data[,"snag_2014"],Plot.data[,c("snag_2014","snag_2015","snag_2016")]))
snag.z <- (snag.x - zscore.factors["snag","mean"])/zscore.factors["snag","sd"]

QMD.x <- as.matrix(cbind(Plot.data[,"QMD_2014"],Plot.data[,c("QMD_2014","QMD_2015","QMD_2016")]))
QMD.z <- (QMD.x - zscore.factors["QMD","mean"])/zscore.factors["QMD","sd"]

dimnames(EarlInf.x)[[2]] <- dimnames(MidInf.x)[[2]] <- dimnames(snag.x)[[2]] <- dimnames(QMD.x)[[2]] <-
dimnames(EarlInf.z)[[2]] <- dimnames(MidInf.z)[[2]] <- dimnames(snag.z)[[2]] <- dimnames(QMD.z)[[2]] <-
  c("2013","2014","2015","2016")

# Correlation matrix #
cor(cbind(EarlInf=as.numeric(EarlInf.x),MidInf=as.numeric(MidInf.x),snag=as.numeric(snag.x),
          QMD=as.numeric(QMD.x)))

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

#--- Occupancy model with persistence parameter (save finite-sample estimates for plotting)---#

# Bin assignments for finite-sample estimates #
  # Bins for early infestation (mean): 0 (0), 1-8 (3.08), 9-50 (21.08)
EarlInf0 <- (EarlInf.z==min(EarlInf.z))*1
EarlInf.low <- (EarlInf.z>min(EarlInf.z)&EarlInf.z<0.2)*1
EarlInf.high <- (EarlInf.z>0.2)*1

  # Bins for mid infestation (mean): 0-1 (0.29), 2-6 (3.13), 7-39 (14.19)
MidInf.low <- (MidInf.z<(-0.6))*1
MidInf.mod <- (MidInf.z>(-0.6)&MidInf.z<0.1)*1
MidInf.high <- (MidInf.z>0.1)*1

  # Bins for snag (mean): 0-1 (0.61), 2-4 (2.68), 5-32 (9.77)
snag.low <- (snag.z<(-0.65))*1
snag.mod <- (snag.z>(-0.65)&snag.z<0.1)
snag.high <- (snag.z>0.1)

  # Bins for QMD (mean): 32-520 (305), 520-810 (633), 810-2871 (1136)
QMD.low <- (QMD.z<(-0.39))
QMD.mod <- (QMD.z>(-0.39)&QMD.z<0.277)
QMD.high <- (QMD.z>0.277)

  # Exclude 2013 at transects 1a and 2a from finite-sample estimates
ind <- which(is.element(substr(plot,1,2),c("1a","2a")))
EarlInf0[ind,"2013"] <- EarlInf.low[ind,"2013"] <- EarlInf.high[ind,"2013"] <- MidInf.low[ind,"2013"] <-
  MidInf.mod[ind,"2013"] <- MidInf.high[ind,"2013"] <- snag.low[ind,"2013"] <- snag.mod[ind,"2013"] <-
  snag.high[ind,"2013"] <- QMD.low[ind,"2013"] <- QMD.mod[ind,"2013"] <- QMD.high[ind,"2013"] <- 0
rm(ind)

  # Remove values for clear-cut pointXyear occassions #
#ind <- which(is.element(substr(plot,1,2),c("1a7","1a8","1a9")))
#EarlInf0[ind,c("2015","2016")] <- EarlInf.low[ind,c("2015","2016")] <- EarlInf.high[ind,c("2015","2016")] <-
#  MidInf.low[ind,c("2015","2016")] <- MidInf.mod[ind,c("2015","2016")] <- MidInf.high[ind,c("2015","2016")] <-
#  snag.low[ind,c("2015","2016")] <- snag.mod[ind,c("2015","2016")] <- snag.high[ind,c("2015","2016")] <-
#  QMD.low[ind,c("2015","2016")] <- QMD.mod[ind,c("2015","2016")] <- QMD.high[ind,c("2015","2016")] <- 0

#ind <- which(plot=="2a5")
#EarlInf0[ind,"2016"] <- EarlInf.low[ind,"2016"] <- EarlInf.high[ind,"2016"] <-
#  MidInf.low[ind,"2016"] <- MidInf.mod[ind,"2016"] <- MidInf.high[ind,"2016"] <-
#  snag.low[ind,"2016"] <- snag.mod[ind,"2016"] <- snag.high[ind,"2016"] <-
#  QMD.low[ind,"2016"] <- QMD.mod[ind,"2016"] <- QMD.high[ind,"2016"] <- 0
#rm(ind)

# Assemble the data names list for JAGS #
data <- list("Y","nplot","nyear","n.visits","EarlInf.z","MidInf.z","snag.z","QMD.z","EarlInf0","EarlInf.low","EarlInf.high",
             "MidInf.low","MidInf.mod","MidInf.high","snag.low","snag.mod","snag.high","QMD.low","QMD.mod","QMD.high")

# Assemble the initial values for JAGS.  Here z is the apparent occupancy state 

inits <- function(){list(Z=Z.init)}

# Assemble the parameters vector for JAGS (What we want to track).
parameters <- c("p","beta0","beta1","beta.EarlInf","beta.MidInf","beta.snag","beta.QMD","test","psi.EI","psi.MI",
                "psi.snag","psi.QMD","psi.fs")

sink("model.txt")
cat("

# The model.
model{

# Prior for detection probability
p ~ dunif(0,1) 
  
# Priors for occupancy model parameters.
beta0 ~ dnorm(0,0.01)T(-30,30) # Mean logit colonization probability (previously unoccupied sites)
beta1 ~ dnorm(0,0.01)T(-30,30) # Mean offset for persistence probability (previously occupied sites)

beta.QMD ~ dnorm(0,0.01)T(-30,30)
beta.EarlInf ~ dnorm(0,0.01)T(-30,30)
beta.MidInf ~ dnorm(0,0.01)T(-30,30)
beta.snag ~ dnorm(0,0.01)T(-30,30)
  
#Model
psi0 ~ dunif(0,1)
for(j in 1:nplot){  # Loop for year 1 - uses z0 to calculate persistence offset
  Z0[j] ~ dbern(psi0)
  ###Model for occupancy probability ###
  logit(psi[j,1]) <- beta0 + beta1*Z0[j] + beta.EarlInf*EarlInf.z[j,1] + beta.MidInf*MidInf.z[j,1] + 
    beta.snag*snag.z[j,1] + beta.QMD*QMD.z[j,1]
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
    logit(psi[j,t]) <- beta0 + beta1*Z[j,(t-1)] + beta.EarlInf*EarlInf.z[j,t] + beta.MidInf*MidInf.z[j,t] + 
      beta.snag*snag.z[j,t] + beta.QMD*QMD.z[j,t]
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
    #___________Zs for finite-sample estimates________________
  for(t in 1:nyear){  # Loop through remaining years 2 and up - uses Z[,,t-1] to calculate persistence offset
    Z.EI0[j,t] <- Z[j,t]*EarlInf0[j,t]
    Z.EIlow[j,t] <- Z[j,t]*EarlInf.low[j,t]
    Z.EIhigh[j,t] <- Z[j,t]*EarlInf.high[j,t]
      
    Z.MIlow[j,t] <- Z[j,t]*MidInf.low[j,t]
    Z.MImod[j,t] <- Z[j,t]*MidInf.mod[j,t]
    Z.MIhigh[j,t] <- Z[j,t]*MidInf.high[j,t]

    Z.snaglo[j,t] <- Z[j,t]*snag.low[j,t]
    Z.snagmd[j,t] <- Z[j,t]*snag.mod[j,t]
    Z.snaghi[j,t] <- Z[j,t]*snag.high[j,t]
    
    Z.QMDlo[j,t] <- Z[j,t]*QMD.low[j,t]
    Z.QMDmd[j,t] <- Z[j,t]*QMD.mod[j,t]
    Z.QMDhi[j,t] <- Z[j,t]*QMD.high[j,t]
    #_________________________________________________________
    }
  }

# Derive finite-sample occupancy estimates #
psi.EI[1] <- sum(Z.EI0[,])/sum(EarlInf0[,])
psi.EI[2] <- sum(Z.EIlow[,])/sum(EarlInf.low[,])
psi.EI[3] <- sum(Z.EIhigh[,])/sum(EarlInf.high[,])

psi.MI[1] <- sum(Z.MIlow[,])/sum(MidInf.low[,])
psi.MI[2] <- sum(Z.MImod[,])/sum(MidInf.mod[,])
psi.MI[3] <- sum(Z.MIhigh[,])/sum(MidInf.high[,])

psi.snag[1] <- sum(Z.snaglo[,])/sum(snag.low[,])
psi.snag[2] <- sum(Z.snagmd[,])/sum(snag.mod[,])
psi.snag[3] <- sum(Z.snaghi[,])/sum(snag.high[,])

psi.QMD[1] <- sum(Z.QMDlo[,])/sum(QMD.low[,])
psi.QMD[2] <- sum(Z.QMDmd[,])/sum(QMD.mod[,])
psi.QMD[3] <- sum(Z.QMDhi[,])/sum(QMD.high[,])

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
saveObject(bugout,"Mod_Pers_EInfMInfSnagQMD")
#saveObject(bugout,"Mod_Pers_EInfMInfSnagQMD_noClrCut")
