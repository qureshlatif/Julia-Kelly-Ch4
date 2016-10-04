#setwd("C:/Users/qlatif/Desktop/JuliaBeetle_Chap4") # Change to environment
setwd("F:/research stuff/FS_PostDoc/outside_consult/JuliaBeetle_Chap4") # Change to environment

## Load and clean data ##
Detection.data <- read.csv("BirdsAll72Plots.csv",header=T,stringsAsFactors=F)
Detection.data <- Detection.data[-which(Detection.data$point==88),]

  # Cleanup plot IDs 
Detection.data <- Detection.data[-which(Detection.data$plot=="11a"),]
Detection.data <- Detection.data[-which(Detection.data$plot=="11b"),]
Detection.data <- Detection.data[-which(Detection.data$plot=="4a"),]
Detection.data <- Detection.data[-which(Detection.data$plot=="4b"),]
Detection.data <- Detection.data[-which(Detection.data$plot=="5a"),]
Detection.data <- Detection.data[-which(Detection.data$plot=="5b"),]
Detection.data <- Detection.data[-which(Detection.data$plot==""),]
Detection.data <- Detection.data[-which(Detection.data$plot=="2a"),]
Detection.data$plot[which(Detection.data$plot=="11a 1")] <- "11a1"
Detection.data$plot[which(Detection.data$plot=="11a 7")] <- "11a7"
Detection.data$plot[which(Detection.data$plot=="11a 8")] <- "11a8"
Detection.data$plot[which(Detection.data$plot=="11a 9")] <- "11a9"
plot <- sort(unique(Detection.data$plot))

  # Keep only three-toed detections
Detection.data <- Detection.data[which(Detection.data$species=="attw"),]

  # Check distribution of distances #
sort(unique(Detection.data$distance),na.last=T)
#hist(Detection.data$distance,breaks=seq(0,200,by=25))
Detection.data <- Detection.data[which(Detection.data$distance<=75),] # Distance truncation

  # Check minute values
#sort(unique(Detection.data$minute),na.last=T)
#tapply(Detection.data$minute,Detection.data$minute,length)
#Detection.data <- # One detection with missing minute. Removed for now. Need to add back in if we just analyze raw counts.
#  Detection.data[-is.na(Detection.data$minute),] 

  # Check years
sort(unique(Detection.data$year))

  # Check visit numbers
sort(unique(Detection.data$survey_round))

  # Check level of clustering #
sort(unique(Detection.data$cluster_code),na.last=T) # None
sort(unique(Detection.data$cluster_size),na.last=T)
tapply(Detection.data$cluster_size,Detection.data$cluster_size,length)
for(i in which(Detection.data$cluster_size=="2")) # Replicates record for clusters of 2 individuals.
  Detection.data <- rbind(Detection.data,Detection.data[i,])

  # Find maximum number of individuals #
#occassionID <- paste(Detection.data$plot,"_",Detection.data$year,"_",Detection.data$survey_round,"_",
#                     Detection.data$minute,sep="")
#max.ind <- max(tapply(occassionID,occassionID,length)) # Maximum number of individuals per minute
occassionID <- paste(Detection.data$plot,"_",Detection.data$year,"_",Detection.data$survey_round,sep="")
max(tapply(occassionID,occassionID,length)) # Maximum number of individuals per visit
hist(tapply(occassionID,occassionID,length))

##---Compile detection data array (N-mixture model - 2 levels)---##

year <- sort(unique(Detection.data$year))
n.visits <- c(1,3,3,3) # 1 visit in 2013, then 3 visits in each other year
Y.arry <- Dist.arry <- array(0,c(length(plot),length(year),max(n.visits)))
for(j in 1:dim(Y.arry)[1]) for(t in 1:dim(Y.arry)[2]) for(v in 1:dim(Y.arry)[3]) {
  obs <- Detection.data[which(Detection.data$plot==plot[j]&Detection.data$year==year[t]&
                                Detection.data$survey_round==v),]
  if(nrow(obs)>0) Y.arry[j,t,v] <- nrow(obs)
  }
Y.arry[which(is.element(substr(plot,1,2),c("1a","2a","1b","2b"))),1,] <- NA # Sites 1 and 2 not surveyed in 2013.
Y.arry[,1,2:3] <- NA # Remove repeat visits from 2013

##---Compile plot-level covariate matrix---##
Plot.data <- read.csv("PlotLevel_Data.csv",header=T,stringsAsFactors=F)
sort(unique(Plot.data$plot),na.last=T)
Plot.data <- Plot.data[order(Plot.data$plot),] # Order plots to match detection data.

Plot.data$treat_cont[which(Plot.data$treat_cont=="t ")] <- "t"
unique(Plot.data$treat_cont)

Plot.data <- Plot.data[,c("plot","treat_cont")]

# Calculate plot-level values for tree data #
Tree.data <- read.csv("TreesAll72Plots.csv",header=T,stringsAsFactors=F)
Tree.data <- Tree.data[-which(Tree.data$site==""),]
sum(!is.element(Tree.data$plot,Plot.data$plot))
sum(!is.element(Plot.data$plot,Tree.data$plot))

Tree.data$core_outsideplot[which(Tree.data$core_outsideplot=="yes")] <- "y"
sort(unique(Tree.data$core_outsideplot))

  # Check tree species 
Tree.data$species[which(Tree.data$species=="aba")] <- "abla"
Tree.data$species[which(Tree.data$species=="abla ")] <- "abla"
Tree.data$species[which(Tree.data$species=="albla")] <- "abla"
Tree.data$species[which(Tree.data$species=="albla")] <- "abla"
Tree.data$species[which(Tree.data$species=="piea")] <- "pien"
Tree.data$species[which(Tree.data$species=="pien ")] <- "pien"
Tree.data$species[which(Tree.data$species=="ppien")] <- "pien"
Tree.data <- Tree.data[-which(Tree.data$species==""),] # Check with Julia
Tree.data <- Tree.data[-which(Tree.data$species=="na"),] # Check with Julia
sort(unique(Tree.data$species))

  # Check dbh columns
sort(unique(Tree.data$snag_dbh),na.last=T)
sort(unique(Tree.data$stump_dbh),na.last=T)
sort(unique(Tree.data$tree_dbh),na.last=T)

  # Check that there is only one dbh value per tree and look for missing DBH
sum(!is.na(Tree.data$snag_dbh))+sum(!is.na(Tree.data$stump_dbh))+sum(!is.na(Tree.data$tree_dbh)) # 9 trees missing DBH
sum((!is.na(Tree.data$snag_dbh))+(!is.na(Tree.data$stump_dbh))+(!is.na(Tree.data$tree_dbh))) # Should be same as prev row
  # Only one value per tree found. 9 missing values found.

  # Combine dbh columns and create indicator of tree type
Tree.data$type <- ""
Tree.data$type[which(!is.na(Tree.data$snag_dbh))] <- "snag"
Tree.data$type[which(!is.na(Tree.data$tree_dbh))] <- "tree"
Tree.data$type[which(!is.na(Tree.data$stump_dbh))] <- "stump"
Tree.data$type[which(Tree.data$type=="")] <- "tree"
Tree.data$dbh <- 0
for(i in 1:nrow(Tree.data)) Tree.data$dbh[i] <-
  sum(Tree.data$snag_dbh[i],Tree.data$tree_dbh[i],Tree.data$stump_dbh[i],na.rm=T)

  # Fix additional dbh errors
Tree.data$dbh[which(Tree.data$plot=="5b2"&Tree.data$dbh==3)] <- 13
Tree.data$dbh[which(Tree.data$plot=="5b4"&Tree.data$dbh==3)] <- 34
Tree.data$dbh[which(Tree.data$plot=="5b5"&Tree.data$dbh==1)] <- 11
Tree.data$dbh[which(Tree.data$plot=="5b9"&Tree.data$dbh==3)] <- 31
Tree.data$dbh[which(Tree.data$dbh==0)] <- NA
sort(unique(Tree.data$dbh),na.last=T)

  # Check status
Tree.data$status[which(Tree.data$status=="Br")] <- "br"
Tree.data$status[which(Tree.data$status=="g")] <- "gr"
Tree.data$status[which(Tree.data$status=="grr")] <- "gr"
Tree.data$status[which(Tree.data$status=="i")] <- "li"
Tree.data$status[which(Tree.data$status=="l")] <- "li"
Tree.data$status[which(Tree.data$status=="lo")] <- "li"
Tree.data$status[which(Tree.data$status=="te")] <- "tw"
Tree.data$status[which(Tree.data$status=="r")] <- "br"
Tree.data$status[which(Tree.data$status=="w")] <- "tw"
Tree.data <- # Remove this one tree because not found in raw data
  Tree.data[-which(Tree.data$plot=="1a5"&Tree.data$species=="abla"&Tree.data$tree_dbh==13),]
Tree.data$status[which(is.element(Tree.data$status,c("br","sn","de")))] <- "snag"
Tree.data$status[which(is.element(Tree.data$status,c("nd","tw")))] <- "mid Inf"
Tree.data$status[which(is.element(Tree.data$status,c("gr","ye")))] <- "early Inf"
Tree.data$status[which(Tree.data$status=="li")] <- "healthy"
sort(unique(Tree.data$status),na.last=T)
#u = stump

  # Fill in 3 missing dbh values
mn <- mean(Tree.data$dbh[which(Tree.data$plot=="5b5"&Tree.data$species=="pien")],na.rm=T)
Tree.data$dbh[which(Tree.data$plot=="5b5"&is.na(Tree.data$dbh))] <- mn

mn <- mean(Tree.data$dbh[which(Tree.data$plot=="11a2"&Tree.data$species=="pien")],na.rm=T)
Tree.data$dbh[which(Tree.data$plot=="11a2"&is.na(Tree.data$dbh))] <- mn

mn <- mean(Tree.data$dbh[which(Tree.data$plot=="4a6"&Tree.data$species=="abla")],na.rm=T)
Tree.data$dbh[which(Tree.data$plot=="4a6"&is.na(Tree.data$dbh))] <- mn

rm(mn)

  # Compile tree summary values into Plot-level table #
Plot.data$QMD_2014 <- Plot.data$QMD_2015 <- Plot.data$QMD_2016 <-
  Plot.data$EarlInf_2016 <- Plot.data$EarlInf_2015 <- Plot.data$EarlInf_2014 <- 
  Plot.data$MidInf_2016 <- Plot.data$MidInf_2015 <- Plot.data$MidInf_2014 <- 0
  Plot.data$snag_2016 <- Plot.data$snag_2015 <- Plot.data$snag_2014 <- 0
for(i in 1:nrow(Plot.data)) {
  obs2014 <- Tree.data[which(Tree.data$plot==Plot.data$plot[i]&Tree.data$year==2014&Tree.data$core_outsideplot!="y"),]
  obs2015 <- Tree.data[which(Tree.data$plot==Plot.data$plot[i]&Tree.data$year==2015&Tree.data$core_outsideplot!="y"),]
  obs2016 <- Tree.data[which(Tree.data$plot==Plot.data$plot[i]&Tree.data$year==2016&Tree.data$core_outsideplot!="y"),]
  Plot.data$EarlInf_2014[i] <- sum(obs2014$status=="early Inf")
  Plot.data$EarlInf_2015[i] <- sum(obs2015$status=="early Inf")
  Plot.data$EarlInf_2016[i] <- sum(obs2016$status=="early Inf")
  Plot.data$MidInf_2014[i] <- sum(obs2014$status=="mid Inf")
  Plot.data$MidInf_2015[i] <- sum(obs2015$status=="mid Inf")
  Plot.data$MidInf_2016[i] <- sum(obs2016$status=="mid Inf")
  Plot.data$snag_2014[i] <- sum(obs2014$status=="snag")
  Plot.data$snag_2015[i] <- sum(obs2015$status=="snag")
  Plot.data$snag_2016[i] <- sum(obs2016$status=="snag")
  Plot.data$QMD_2014[i] <- mean(obs2014$dbh[which(obs2014$type!="stump"&obs2014$status!="mid Inf"&
                                                    obs2014$status!="snag")]^2,na.rm=T)
  Plot.data$QMD_2015[i] <- mean(obs2015$dbh[which(obs2015$type!="stump"&obs2015$status!="mid Inf"&
                                                    obs2015$status!="snag")]^2,na.rm=T)
  Plot.data$QMD_2016[i] <- mean(obs2016$dbh[which(obs2016$type!="stump"&obs2016$status!="mid Inf"&
                                                    obs2016$status!="snag")]^2,na.rm=T)
}

#cor(rbind(as.matrix(Plot.data[,c("QMD","snag","EarlInf_2014","MidInf_2014")]), # Inspect correlation plot if desired.
#          as.matrix(Plot.data[,c("QMD","snag","EarlInf_2015","MidInf_2015")]),make.row.names = F))

# Variance inflation factors #
vars <- c("QMD_2015","EarlInf_2015","MidInf_2015","snag_2015")
VIFs <- matrix(NA,nrow=length(vars),ncol=1,dimnames=list(vars,"VIF"))
for(i in 1:length(vars)) {
  m <- lm(as.formula(paste(vars[i],"~",paste(vars[-i],collapse="+",sep=""),sep="")),data=Plot.data)
  VIFs[i] <- 1/(1 - summary(m)$r.squared)
}
VIFs

# Cleanup
rm(i,j,m,t,v,vars,obs2014,obs2015,obs2016,obs,occassionID)

save.image("Data_compiled.RData")
