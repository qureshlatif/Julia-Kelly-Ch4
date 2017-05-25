library(R2jags)
library(R.utils)
setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/JuliaH_dissert/Chap4")  # Set workspace to location with 'Data_compiled' workspace.
#setwd("~/Documents/Rdata/JuliaBeetle_Chap4")
load("Data_compiled.RData")

##### Re-run covariate extraction, scaling, and compile detection data ######
rows <- c("QMD","EarlInf","MidInf","snag")
cols <- c("mean","sd")
zscore.factors <- matrix(NA,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols))
rm(rows,cols)
zscore.factors["QMD",] <- c(mean(as.matrix(Plot.data[,c("QMD_pien_2014","QMD_pien_2014","QMD_pien_2015","QMD_pien_2016")])),
                            sd(as.matrix(Plot.data[,c("QMD_pien_2014","QMD_pien_2014","QMD_pien_2015","QMD_pien_2016")])))
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

QMD.x <- as.matrix(cbind(Plot.data[,"QMD_pien_2014"],Plot.data[,c("QMD_pien_2014","QMD_pien_2015","QMD_pien_2016")]))
QMD.z <- (QMD.x - zscore.factors["QMD","mean"])/zscore.factors["QMD","sd"]

dimnames(EarlInf.x)[[2]] <- dimnames(MidInf.x)[[2]] <- dimnames(snag.x)[[2]] <- dimnames(QMD.x)[[2]] <-
  dimnames(EarlInf.z)[[2]] <- dimnames(MidInf.z)[[2]] <- dimnames(snag.z)[[2]] <- dimnames(QMD.z)[[2]] <-
  c("2013","2014","2015","2016")

# Correlation matrix #
#cor(cbind(EarlInf=as.numeric(EarlInf.x),MidInf=as.numeric(MidInf.x),snag=as.numeric(snag.x),
#          QMD=as.numeric(QMD.x)))

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

mod <- loadObject("Mod_Pers_EInfMInfSnagQMD")
#mod <- loadObject("Mod_Pers_EInfMInfSnagQMD_PriorTrunc")
write.csv(mod$BUGSoutput$summary,"Model_summary.csv") # Export summary of model parameters and output.

#mod <- loadObject("Mod_Pers_EInfMInfSnagQMD_noClrCut")
#write.csv(mod$BUGSoutput$summary,"Model_summary_noClrCut.csv") # Export summary of model parameters and output.

##### Plot modeled occupancy relationships #####
library(ggplot2)
library(cowplot)
expit <- function(x) exp(x)/(1+exp(x))

## Plot covariate parameter estimates ##

rows <- c("beta0","beta1","EarlInf","MidInf","Snag","QMD")
cols <- c("med","lo","hi")
dat <- data.frame(matrix(0,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols)))

dat$med <- apply(mod$BUGSoutput$sims.array[,,c("beta0","beta1","beta.EarlInf","beta.MidInf",
                                                   "beta.snag","beta.QMD")],3,median)
dat$lo <- apply(mod$BUGSoutput$sims.array[,,c("beta0","beta1","beta.EarlInf","beta.MidInf",
                                                   "beta.snag","beta.QMD")],3,function(x)
                                                     quantile(x,prob=0.025))
dat$hi <- apply(mod$BUGSoutput$sims.array[,,c("beta0","beta1","beta.EarlInf","beta.MidInf",
                                                   "beta.snag","beta.QMD")],3,function(x)
                                                     quantile(x,prob=0.975))
dat$X <- c(6.5,5.5,4:1)

p <- ggplot(data = dat,aes(x=X,y=med)) +
  geom_errorbar(aes(ymin=lo,ymax=hi),size=1,width=0) +
  geom_point(size=2.5) + 
  geom_hline(yintercept=0,linetype="dotted") +
  coord_flip() +
  scale_x_continuous(breaks=c(1:4,5.5,6.5),labels=c(expression(beta["QMD"]),
                                                     expression(beta["sng"]),
                                                     expression(beta["Minf"]),
                                                     expression(beta["Einf"]),
                                                     expression(beta[1]),
                                                     expression(beta[0]))) +
  ylab("Parameter estimate") + xlab("              Covariate effects                             Intercept terms") +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.text.y=element_text(size=25))

save_plot("Occupancy_parameters.jpeg", p, ncol = 3, nrow = 3, dpi=600) 

# Early infestation #
  # finite-sample estimates #
X <- c(0,3.08,21.08) # mean values for finite-sample bins
psi <- apply(mod$BUGSoutput$sims.list$psi.EI,2,median)
psi.lo <- apply(mod$BUGSoutput$sims.list$psi.EI,2,function(x) quantile(x,prob=0.025,type=8))
psi.hi <- apply(mod$BUGSoutput$sims.list$psi.EI,2,function(x) quantile(x,prob=0.975,type=8))
dat.est <- data.frame(cbind(X,psi,psi.lo,psi.hi))

  # Observed values (covariate value for point X year detections and non-detections) #
#det <- apply(Y,c(1,2),function(x) sum(x>0,na.rm=T)*1)
#det[is.na(Y)] <- NA
#detections <- data.frame(as.numeric(EarlInf.x[which(det==1)]))
#nondetects <- data.frame(as.numeric(EarlInf.x[which(det==0)]))
#names(detections) <- names(nondetects) <- "x"
#detections <- detections + runif(nrow(detections),-0.3,0.3)
#nondetects <- nondetects + runif(nrow(nondetects),-0.3,0.3)

  # predicted probabilities #
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

# With raw data instead of finite sample estimates #
#p <- ggplot(data = dat.prd,aes(x=X,y=psi)) +
#  geom_line(size=1,linetype="solid") +
#  geom_line(aes(y=psi.lo),size=1,linetype="dashed") +
#  geom_line(aes(y=psi.hi),size=1,linetype="dashed") +
#  geom_point(data=detections,aes(x=x,y=1.01),size=3,alpha=0.5) + 
#  geom_point(data=nondetects,aes(x=x,y=0),size=3,alpha=0.5) + 
#  ylab("Occupancy probability") + xlab("Number Early Infested Trees") +
#  scale_y_continuous(lim=c(0,1.05),breaks=c(0,0.25,0.5,0.75,1)) +
#  theme(axis.title.x=element_text(size=30)) +
#  theme(axis.title.y=element_text(size=30)) +
#  theme(axis.text.x=element_text(size=20)) +
#  theme(axis.text.y=element_text(size=25)) +
#  guides(shape=FALSE,linetype=FALSE)

# With finite-sample estimates #
p <- ggplot(data = dat.prd,aes(x=X,y=psi)) +
  geom_line(size=1,linetype="solid") +
  geom_line(aes(y=psi.lo),size=1,linetype="dashed") +
  geom_line(aes(y=psi.hi),size=1,linetype="dashed") +
  geom_point(data=dat.est,aes(x=X,y=psi),size=5) + 
  geom_errorbar(data=dat.est,aes(x=X,ymin=psi.lo,ymax=psi.hi),size=1,width=1) +
  ylab("Point occupancy") + xlab("Number Early Infested Trees") +
  scale_y_continuous(lim=c(0,1.05),breaks=c(0,0.25,0.5,0.75,1)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=25)) +
  guides(shape=FALSE,linetype=FALSE)

save_plot("EarlInf_relation.jpeg", p, ncol = 3, nrow = 3, dpi=600) 
