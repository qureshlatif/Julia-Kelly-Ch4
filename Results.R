library(R2jags)
library(R.utils)
setwd("F:/research stuff/FS_PostDoc/outside_consult/JuliaBeetle_Chap4")  # Set workspace to location with 'Data_compiled' workspace.
#setwd("~/Documents/Rdata/JuliaBeetle_Chap4")
load("Data_compiled.RData")

##### Re-run covariate extraction, scaling, and compile detection data ######
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
out <- cor(cbind(EarlInf=as.numeric(EarlInf.x),MidInf=as.numeric(MidInf.x),snag=as.numeric(snag.x),
          QMD=as.numeric(QMD.x)))
write.csv(out,"Covariate_correlations.csv")
rm(out)

Y <- apply(Y.arry,c(1,2,3),function(x) 1*(x>0)) # Collapse counts to detection-nondetection data
Y <- # Collapse further to represent the number of visits when species detected for each plot X year occassion.
  apply(Y,c(1,2),function(x) sum(x,na.rm=T))
Y[which(is.element(substr(plot,1,2),c("1a","2a","1b","2b"))),1] <- NA # Sites 1 and 2 not surveyed in 2013.

nplot <- length(plot)
nyear <- length(year)

Z.init <- apply(Y,c(1,2),function(x) sum(x>0)) # Initial values for occupancy state parameter.

mod <- loadObject("Mod_Pers_EInfMInfSnagQMD")
write.csv(mod$BUGSoutput$summary,"Model_summary.csv") # Export summary of model parameters and output.

#mod <- loadObject("Mod_Pers_EInfMInfSnagQMD_noClrCut")
#write.csv(mod$BUGSoutput$summary,"Model_summary_noClrCut.csv") # Export summary of model parameters and output.

##### Plot modeled occupancy relationships #####
library(ggplot2)
library(cowplot)
expit <- function(x) exp(x)/(1+exp(x))

# Early infestation #
  # finite-sample estimates #
X <- c(0,3.08,21.08) # mean values for finite-sample bins
psi <- apply(mod$BUGSoutput$sims.list$psi.EI,2,median)
psi.lo <- apply(mod$BUGSoutput$sims.list$psi.EI,2,function(x) quantile(x,prob=0.025,type=8))
psi.hi <- apply(mod$BUGSoutput$sims.list$psi.EI,2,function(x) quantile(x,prob=0.975,type=8))
dat.est <- data.frame(cbind(X,psi,psi.lo,psi.hi))

X <- seq(min(EarlInf.x),max(EarlInf.x),length.out=20)
beta0 <- mod$BUGSoutput$sims.list$beta0
beta1 <- mod$BUGSoutput$sims.list$beta1
psifs <- mod$BUGSoutput$sims.list$psi.fs
betaX <- mod$BUGSoutput$sims.list$beta.EarlInf

Z <- (X - zscore.factors["EarlInf","mean"])/zscore.factors["EarlInf","sd"]
psi <- psi.lo <- psi.hi <- numeric(length=length(X))
for(i in 1:length(X)) {
  y <- expit(beta0 + beta1*psifs + betaX*Z[i])
  psi[i] <- median(y)
  psi.lo[i] <- quantile(y,prob=0.025,type=8)
  psi.hi[i] <- quantile(y,prob=0.975,type=8)
}
dat.prd <- data.frame(cbind(X,psi,psi.lo,psi.hi))

pEI <- ggplot(data = dat.prd,aes(x=X,y=psi)) +
  geom_line(size=1,linetype="solid") +
  geom_line(aes(y=psi.lo),size=1,linetype="dashed") +
  geom_line(aes(y=psi.hi),size=1,linetype="dashed") +
  geom_point(data=dat.est,aes(x=X,y=psi),size=5) + 
  geom_errorbar(data=dat.est,aes(x=X,ymin=psi.lo,ymax=psi.hi),size=1,width=1) +
  ylab(NULL) + xlab("Number Early Infested Trees") +
  scale_y_continuous(lim=c(0,1.05),breaks=c(0,0.25,0.5,0.75,1)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=25)) +
  guides(shape=FALSE,linetype=FALSE) +
  geom_text(aes(x=1,y=1.05),label = "A",size=10)

# Mid infestation #
  # finite-sample estimates #
X <- c(0.29,3.13,14.19) # mean values for finite-sample bins
psi <- apply(mod$BUGSoutput$sims.list$psi.MI,2,median)
psi.lo <- apply(mod$BUGSoutput$sims.list$psi.MI,2,function(x) quantile(x,prob=0.025,type=8))
psi.hi <- apply(mod$BUGSoutput$sims.list$psi.MI,2,function(x) quantile(x,prob=0.975,type=8))
dat.est <- data.frame(cbind(X,psi,psi.lo,psi.hi))

X <- seq(min(MidInf.x),max(MidInf.x),length.out=20)
beta0 <- mod$BUGSoutput$sims.list$beta0
beta1 <- mod$BUGSoutput$sims.list$beta1
psifs <- mod$BUGSoutput$sims.list$psi.fs
betaX <- mod$BUGSoutput$sims.list$beta.MidInf

Z <- (X - zscore.factors["MidInf","mean"])/zscore.factors["MidInf","sd"]
psi <- psi.lo <- psi.hi <- numeric(length=length(X))
for(i in 1:length(X)) {
  y <- expit(beta0 + beta1*psifs + betaX*Z[i])
  psi[i] <- median(y)
  psi.lo[i] <- quantile(y,prob=0.025,type=8)
  psi.hi[i] <- quantile(y,prob=0.975,type=8)
}
dat.prd <- data.frame(cbind(X,psi,psi.lo,psi.hi))

pMI <- ggplot(data = dat.prd,aes(x=X,y=psi)) +
  geom_line(size=1,linetype="solid") +
  geom_line(aes(y=psi.lo),size=1,linetype="dashed") +
  geom_line(aes(y=psi.hi),size=1,linetype="dashed") +
  geom_point(data=dat.est,aes(x=X,y=psi),size=5) + 
  geom_errorbar(data=dat.est,aes(x=X,ymin=psi.lo,ymax=psi.hi),size=1,width=1) +
  ylab(NULL) + xlab("Number Mid-infested Trees") +
  scale_y_continuous(lim=c(0,1.05),breaks=c(0,0.25,0.5,0.75,1)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=25)) +
  guides(shape=FALSE,linetype=FALSE) +
  geom_text(aes(x=1,y=1.05),label = "B",size=10)

# Snags #
  # finite-sample estimates #
X <- c(0.61,2.68,9.77) # mean values for finite-sample bins
psi <- apply(mod$BUGSoutput$sims.list$psi.snag,2,median)
psi.lo <- apply(mod$BUGSoutput$sims.list$psi.snag,2,function(x) quantile(x,prob=0.025,type=8))
psi.hi <- apply(mod$BUGSoutput$sims.list$psi.snag,2,function(x) quantile(x,prob=0.975,type=8))
dat.est <- data.frame(cbind(X,psi,psi.lo,psi.hi))

X <- seq(min(snag.x),max(snag.x),length.out=20)
beta0 <- mod$BUGSoutput$sims.list$beta0
beta1 <- mod$BUGSoutput$sims.list$beta1
psifs <- mod$BUGSoutput$sims.list$psi.fs
betaX <- mod$BUGSoutput$sims.list$beta.snag

Z <- (X - zscore.factors["snag","mean"])/zscore.factors["snag","sd"]
psi <- psi.lo <- psi.hi <- numeric(length=length(X))
for(i in 1:length(X)) {
  y <- expit(beta0 + beta1*psifs + betaX*Z[i])
  psi[i] <- median(y)
  psi.lo[i] <- quantile(y,prob=0.025,type=8)
  psi.hi[i] <- quantile(y,prob=0.975,type=8)
}
dat.prd <- data.frame(cbind(X,psi,psi.lo,psi.hi))

pSnag <- ggplot(data = dat.prd,aes(x=X,y=psi)) +
  geom_line(size=1,linetype="solid") +
  geom_line(aes(y=psi.lo),size=1,linetype="dashed") +
  geom_line(aes(y=psi.hi),size=1,linetype="dashed") +
  geom_point(data=dat.est,aes(x=X,y=psi),size=5) + 
  geom_errorbar(data=dat.est,aes(x=X,ymin=psi.lo,ymax=psi.hi),size=1,width=1) +
  ylab(NULL) + xlab("Number of Snags") +
  scale_y_continuous(lim=c(0,1.05),breaks=c(0,0.25,0.5,0.75,1)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=25)) +
  guides(shape=FALSE,linetype=FALSE) +
  geom_text(aes(x=1,y=1.05),label = "C",size=10)

# QMD #
  # finite-sample estimates #
X <- c(305,633,1136) # mean values for finite-sample bins
psi <- apply(mod$BUGSoutput$sims.list$psi.QMD,2,median)
psi.lo <- apply(mod$BUGSoutput$sims.list$psi.QMD,2,function(x) quantile(x,prob=0.025,type=8))
psi.hi <- apply(mod$BUGSoutput$sims.list$psi.QMD,2,function(x) quantile(x,prob=0.975,type=8))
dat.est <- data.frame(cbind(X,psi,psi.lo,psi.hi))

X <- seq(min(QMD.x),max(QMD.x),length.out=20)
beta0 <- mod$BUGSoutput$sims.list$beta0
beta1 <- mod$BUGSoutput$sims.list$beta1
psifs <- mod$BUGSoutput$sims.list$psi.fs
betaX <- mod$BUGSoutput$sims.list$beta.QMD

Z <- (X - zscore.factors["QMD","mean"])/zscore.factors["QMD","sd"]
psi <- psi.lo <- psi.hi <- numeric(length=length(X))
for(i in 1:length(X)) {
  y <- expit(beta0 + beta1*psifs + betaX*Z[i])
  psi[i] <- median(y)
  psi.lo[i] <- quantile(y,prob=0.025,type=8)
  psi.hi[i] <- quantile(y,prob=0.975,type=8)
}
dat.prd <- data.frame(cbind(X,psi,psi.lo,psi.hi))

pQMD <- ggplot(data = dat.prd,aes(x=X,y=psi)) +
  geom_line(size=1,linetype="solid") +
  geom_line(aes(y=psi.lo),size=1,linetype="dashed") +
  geom_line(aes(y=psi.hi),size=1,linetype="dashed") +
  geom_point(data=dat.est,aes(x=X,y=psi),size=5) + 
  geom_errorbar(data=dat.est,aes(x=X,ymin=psi.lo,ymax=psi.hi),size=1,width=50) +
  ylab(NULL) + xlab(expression("QMD ("~cm^2~")")) +
  scale_y_continuous(lim=c(0,1.05),breaks=c(0,0.25,0.5,0.75,1)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=25)) +
  guides(shape=FALSE,linetype=FALSE) +
  geom_text(aes(x=50,y=1.05),label = "D",size=10)

p <- ggdraw() +  
  draw_plot(pEI, x = 0.05, y = 0.5, width = .47, height = .5) + 
  draw_plot(pMI, x = 0.52, y = 0.5, width = .47, height = .5) + 
  draw_plot(pSnag, x = 0.05, y = 0, width = .47, height = .5) + 
  draw_plot(pQMD, x = 0.52, y = 0, width = .46, height = .5) + 
  draw_plot_label(label="Point occupancy",size=30,x=0,y=0.3,angle=90)

save_plot("Plots.jpeg", p, ncol = 3, nrow = 3, dpi=600) 


##### Model summaries and selection table (obsolete) #####

#models <- c("NoPers_EInfMInfSnagQMD","Pers_EInfMInfSnagQMD",
            "PersDyn_EInfMInfSnagQMD","IndYrs_EInfMInfSnagQMD")
#cols <- c("DIC","WAIC","pD")
#out <- data.frame(matrix(NA,nrow=length(models),ncol=length(cols),
                         dimnames=list(models,cols)),stringsAsFactors=F)

#for(m in 1:length(models)) {
#  md <- loadObject(paste("Mod_",models[m],sep=""))
#  sum <- md$BUGSoutput$summary
#  write.csv(sum[-which(substr(dimnames(sum)[[1]],1,3)=="psi"),],
#            paste("Summary_",models[m],".csv",sep=""))
#  out[m,"DIC"] <- round(md$BUGSoutput$DIC,digits=1)
#  psi <- md$BUGSoutput$sims.list$psi
#  p <- md$BUGSoutput$sims.list$p
#  ## Computed log pointwise predictive density (Gelman eq 3)
#  pr.y <- array(dim = dim(psi))
#  for(i in 1:dim(pr.y)[1]) {
#    psi.p <- psi[i,,]*p[i]
#    for(t in 1:nyear) pr.y[i,,t] <- dbinom(x=Y[,t],prob=psi.p[,t],size=n.visits[t])
#  }
#  pr.y.hat <- apply(pr.y,c(2,3),mean,na.rm=T)
#  clpd <- sum(log(pr.y.hat[is.finite(pr.y.hat)]))
#  
#  ## Estimated Effective Number of Parameters (Gelman eq 6)
#  sample.variance <- apply(log(pr.y),c(2,3),var,na.rm=T)
#  eff.par <- sum(sample.variance[is.finite(sample.variance)])
#  
#  out[m,"WAIC"] <- round(-2*(clpd - eff.par),digits=1)
#  out[m,"pD"] <- round(eff.par,digits=1)
#}
#out <- out[order(out[,"WAIC"]),]
#write.csv(out,"Model_selection.csv")

# Cleanup #
#rm(md,models,i,pr.y,pr.y.hat,clpd,t,eff.par,psi.p,p,psi,sum,m)

##### Plot relationships for selected model ######
mod <- loadObject("Mod_NoPers_EInfMInfSnagQMD")

