# all necker results & analyses figures

# means for response latency
# for additional graphs, etc. use "necker_latency1.R"

# x <- df1
# x <- x[which(x$istest==0),]
# x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
# x$time2 <- round(x$ctime/0.05)*0.05
# unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
# unik <- unlist(unik, recursive = T, use.names=F)
# tmp  <- seq_along(x$time2)[unik]
# x <- x[tmp,]
# rownames(x) <- 1:nrow(x)
# trains <- subset(x, istest==0)
# trains$mismatch <- !trains$mismatch
# 
# # cut off first five seconds
# trains<-trains[-which(trains$time2<5),]
# 
# # now generate the timestamp vector for each trial
# setwd("/Users/trbrooks/documents/past semesters/spring 2016/necker/data directory/necker_round1/data")
# dir()->file.names
# grep("csv", file.names)->temp
# file.names[temp]->pydata1
# 
# # ** ALL filenames in the pydata vector should start with a subject number!!**
# pydata1
# 
# # remove un-numbered cases
# pydata1<-pydata1[grep("^[0-9].*", pydata1)]
# 
# # remove deleted cases
# exclude<-c("001_","002_","004_", "005", "008_", "010_")
# exclude<-grep(paste(exclude,collapse="|"),pydata1,value=TRUE)
# pydata1<-pydata1[! pydata1 %in% exclude]
# 
# # get only the training trials
# exclude<-c("test")
# exclude<-grep(paste(exclude,collapse="|"),pydata1,value=TRUE)
# pydata<-pydata1[! pydata1 %in% exclude]
# 
# doubler <- function(x) {
#   strsplit(as.character((ts1[,x])[1]),",")[[1]]->temp
#   as.numeric(gsub("\\[|\\]", temp, replacement=""))->out
#   options(scipen=999)
#   return (out)
# }
# for (i in 1:length(pydata)){
#   options(stringsAsFactors=FALSE)
#   
#   read.csv(pydata[i], header=TRUE)->ts1
#   
#   charact.Er(11)[1]->participant
#   length(doubler(3)[which(doubler(3)<150 & doubler(3)>=5)])->stamps
#   charact.Er(5)[1]->iscong
#   
#   temp<-c(participant,iscong,stamps)
#   
#   ifelse(i==1, num_swi<-temp, rbind(num_swi,temp)->num_swi)
#   print(i)
# }
# 
# setEPS()
# postscript("barplots.eps")
# plot.new()
# 
# num_swi<-as.data.frame(apply(num_swi, 2, as.numeric))
# colnames(num_swi) <- c("participant","iscong","number")
# 
# 
# temp  <- aggregate(trains$mismatch, by=list(trains$participant, trains$iscong), FUN=sum, simplify=T)
# colnames(temp) <- c("participant", "iscong", "sum")
# switches_num2 <- as.data.frame(merge(temp, num_swi))
# 
# ##
# 
# library(plotrix)
# switches_num2$number<-as.integer(switches_num2$number)
# switches_num2$mean_ms <- (switches_num2$sum*20)/switches_num2$number
# means <- aggregate(switches_num2$mean_ms, by=list(switches_num2$iscong), FUN=mean, simplify=T)
# stds<-aggregate(switches_num2$mean_ms, by=list(switches_num2$iscong), FUN=std.error)
# 
# #aggregate(switches_num2$number, list(switches_num2$iscong), mean)
# 
# t.test(switches_num2$mean_ms[switches_num2$iscong==1],
#        switches_num2$mean_ms[switches_num2$iscong==0], paired=T)
# 
# par(mfrow=c(1,2))
# plot2a <-barplot(means[,2], xpd = FALSE, col=c("white", "gray"), ylab="milliseconds",
#                 ylim=c(0,max(means[,2])+max(stds[,2])))
# segments(plot2a, means[,2] - stds[,2], plot2a, means[,2] + stds[,2], lwd=2)
# abline(h=0,lty=1, lwd=2)
# par(cex.main=1)
# #title("b.")
# par(cex.main=.75)
# title("   Compatible          Incompatible", line=-25)
# 
# 
# # Mean number of switches and SD for test phase.
# 
# x <- df1
# x <- x[which(x$istest==1),]
# x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
# x$time2 <- round(x$ctime/0.05)*0.05
# unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
# unik <- unlist(unik, recursive = T, use.names=F)
# tmp  <- seq_along(x$time2)[unik]
# x <- x[tmp,]
# rownames(x) <- 1:nrow(x)
# 
# 
# switchr <- function(inn){
#   lag<-c(inn[-1],inn[length(inn)])
#   return(as.numeric(inn!=lag))
# }
# 
# library(reshape)
# x$swi<-unlist(by(x$cstate, x$trial, FUN=switchr))
# num.swi<-aggregate(x$swi, by=list(x$istest, x$iscong, x$participant), FUN=sum)
# num.swi<-subset(num.swi, Group.1==1)
# num.swi$Group.1 <- NULL
# names(num.swi) <- c("iscong","participant","number") 
# 
# mean(num.swi$number)
# sd(num.swi$number)
# 
# # Number of switches. Test Phase Congruent condition descriptives
# mean(num.swi[which(num.swi$iscong==1),]$number)
# sd(num.swi[which(num.swi$iscong==1),]$number)
# 
# mean(num.swi[which(num.swi$iscong==0),]$number)
# sd(num.swi[which(num.swi$iscong==0),]$number)
# 
# temp  <- num.swi
# means <- aggregate(temp[,3], by=list(temp[,1]), FUN=mean, simplify=T)
# means <- means[order(c(1,0)),]
# stds  <- aggregate(temp[,3], by=list(temp[,1]), FUN=std.error, simplify=T)
# stds <- stds[order(c(1,0)),]
# 
# means
# stds
# 
# t.test(num.swi[which(num.swi$iscong==1),]$number,
#        num.swi[which(num.swi$iscong==0),]$number,
#        paired = T)
# 
# 
# plot2b <-barplot(means[,2], xpd = FALSE, col=c("white", "gray"), ylab="switches",
#                 ylim=c(0,max(means[,2])+max(stds[,2])))
# segments(plot2b, means[,2] - stds[,2], plot2b, means[,2] + stds[,2], lwd=2)
# abline(h=0,lty=1, lwd=2)
# par(cex.main=1)
# #title("c.")
# par(cex.main=.75)
# title(" Compatible          Incompatible", line=-25)
# dev.off()
# 

## PREFERENCE (total time in each state)
x <- df1
x <- x[which(x$istest==1),]
x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
x$time2 <- round(x$ctime/0.05)*0.05
unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
unik <- unlist(unik, recursive = T, use.names=F)
tmp  <- seq_along(x$time2)[unik]
x <- x[tmp,]
rownames(x) <- 1:nrow(x)
x <- subset(x, istest==1)
x <- x[which(x$time2>=5),]
head(x)

pref<-function(inn){
  temp<-sum(inn)-(length(inn)-sum(inn))
  return(temp/20)
}

prefs<-aggregate(x$cstate, by=list(x$participant, x$iscong), FUN=pref, simplify=T)
aggregate(prefs$x, by=list(prefs$Group.2), FUN=mean)
t.test(prefs[1:13,3],prefs[14:26,3],paired=T)
t.test(prefs$x)

# Descriptives for total orientation preference, by condition
# 1 = congruent, 0 = incongruent
mean(prefs$x)
sd(prefs$x)
aggregate(prefs$x, by=list(prefs$Group.2), FUN=mean)
aggregate(prefs$x, by=list(prefs$Group.2), FUN=sd)

# Random walk plot for switches (Figure 4a)
# setEPS()
# postscript("switchwalk_plot.eps")
# plot.new()
# x <- df1
# x <- x[which(x$istest==1),]
# x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
# x$time2 <- round(x$ctime/0.05)*0.05
# unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
# unik <- unlist(unik, recursive = T, use.names=F)
# tmp  <- seq_along(x$time2)[unik]
# x <- x[tmp,]
# rownames(x) <- 1:nrow(x)
# #x <- subset(x, time2>=5)
# x <- subset(x, istest==1)
# 
# # create a "random walk" version of the switching vector
# switchr <- function(inn){
#   lag<-c(inn[-1],inn[length(inn)])
#   return(as.numeric(inn!=lag))
# }
# 
# library(reshape)
# x$swi<-unlist(by(x$cstate, x$trial, FUN=switchr))
# x$rand <- x$swi
# temp <- by(x$rand, x$trial, FUN = cumsum, simplify = T)
# x$walk <- unlist(temp, recursive = T, use.names=F)
# 
# x$mod_res <- residuals(gcmt2)
# # aggregate(x$walk, by=list(x$iscong), FUN = sum, simplify = T)
# # x$walk<-unlist(by(x$walk, x$trial, detrend))
# # 4A: Random walk plot for preference
# 
# library(plotrix)
# 
# temp <- subset(x, iscong=="0")
# #temp$walk <- unlist(by(temp$walk, temp$trial, detrend))
# temp1 <- aggregate(temp$walk, list(time = temp$time2), mean)
# # temp2 <- aggregate(temp$walk, list(time = temp$time2), std.error)
# temp2 <- aggregate(temp$mod_res, list(time = temp$time2), mean)
# itraj <- temp1[-length(temp1[,1]),]
# istde <- temp2[-length(temp2[,1]),]
# 
# temp <- subset(x, iscong=="1")
# #temp$walk <- unlist(by(temp$walk, temp$trial, detrend))
# temp1 <- aggregate(temp$walk, list(time = temp$time2), mean)
# #temp2 <- aggregate(temp$walk, list(time = temp$time2), std.error)
# temp2 <- aggregate(temp$mod_res, list(time = temp$time2), mean)
# #temp2$x <- temp2$x^2
# ctraj <- temp1[-length(temp1[,1]),]
# cstde <- temp2[-length(temp2[,1]),]
# 
# plot.new()
# par(mfrow=c(1,1))
# 
# lims <- c(0,30)
# plot(itraj$time,rep(0,length(itraj$time)),ylim=lims,type='l',
#      ylab="Detrended number of perceptual switches", xlab="time (seconds)")
# title("a.")
# polygon(c(itraj$time,rev(itraj$time)),c(itraj$x+istde$x,
#                                         rev(itraj$x-istde$x)),col=gray(0.8),border=F)
# polygon(c(itraj$time,rev(itraj$time)),c(ctraj$x+cstde$x,
#                                         rev(ctraj$x-cstde$x)),col=gray(0.8),border=F)
# 
# lines(itraj$time, y=itraj$x, lty=2)
# lines(ctraj$time, y=ctraj$x, cex = .8, lty=1)
# 
# legend(x=100,y=10,legend=c("Incongruent   ", "Congruent", "Model SS Residuals"),
#        lwd=c(1,1,10),col=c("black","black","gray"),lty=c(2,1,1),cex=0.7)
# dev.off()






# Figure 4B: random walk of preference

setEPS()
postscript("prefwalk_plot.eps")
plot.new()
library(plotrix)
#par(mfrow=c(1,2))
x <- df3
x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
x$time2 <- round(x$ctime/0.05)*0.05
unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
unik <- unlist(unik, recursive = T, use.names=F)
tmp  <- seq_along(x$time2)[unik]
x <- x[tmp,]
rownames(x) <- 1:nrow(x)
x <- subset(x, istest==1)

# create a "random walk" version of the state vector
# and rescale to milliseconds
x$rand <- (x$cstate-0.5)/10
temp <- by(x$rand, x$trial, FUN = cumsum, simplify = T)
x$walk <- unlist(temp, recursive = T, use.names=F)

temp <- subset(x, iscong==0)
temp1 <- aggregate(temp$walk, list(time = temp$time2), mean)
temp2 <- aggregate(temp$walk, list(time = temp$time2), std.error)
itraj <- temp1[-length(temp1[,1]),]
istde <- temp2[-length(temp2[,1]),]

temp <- subset(x, iscong==1)
temp1 <- aggregate(temp$walk, list(time = temp$time2), mean)
temp2 <-aggregate(temp$walk, list(time = temp$time2), std.error)
ctraj <- temp1[-length(temp1[,1]),]
cstde <- temp2[-length(temp2[,1]),]


#istde<-aggregate(x$res, list(time = x$time2), mean)[-1,]
#istde<-x$res[which(x$iscong=="0")]
#cstde<-x$res[which(x$iscong=="1")]

par(mfrow=c(1,1))

#lims <- c(-max(c(ctraj$x,itraj$x)), max(c(ctraj$x,itraj$x)))
lims <- c(-20,30)
plot(itraj$time,rep(0,length(itraj$time)),type='l',ylab="Mean Cumulative Percept Score", xlab="time (seconds)", ylim=lims)
polygon(c(itraj$time,rev(itraj$time)),c(itraj$x+istde$x,
                                        rev(itraj$x-istde$x)),col=gray(0.8),border=F)
polygon(c(itraj$time,rev(itraj$time)),c(ctraj$x+cstde$x,
                                        rev(ctraj$x-cstde$x)),col=gray(0.8),border=F)
lines(itraj, lty=2)
lines(ctraj, cex = .8, lty=1)

legend(x=90,y=-10,legend=c("Incompatible", "Compatible", "Standard Error"),
       lwd=c(1,1,10),col=c("black","black","gray"),lty=c(2,1,1),cex=0.7)
# downcube
cube <- c(5,300,200,16,350,0,250,11)
cube <- c(4,21,17,16,22,0,18,11)
#cube[c(2,3,5,7)] <-cube[c(2,3,5,7)]-10
cube[c(1,4,6,8)] <-cube[c(1,4,6,8)]+40
segments(cube[1],cube[2],x1=cube[1],y1=cube[3])
segments(cube[4],cube[2],x1=cube[4],y1=cube[3])
segments(cube[6],cube[5],x1=cube[1],y1=cube[2])
segments(cube[1],cube[3],x1=cube[6],y1=cube[7])
segments(cube[1],cube[3],x1=cube[4],y1=cube[3])
segments(cube[1],cube[2],x1=cube[4],y1=cube[2])
segments(cube[4],cube[2],x1=cube[8],y1=cube[5])
segments(cube[6],cube[5],x1=cube[8],y1=cube[5])
segments(cube[6],cube[7],x1=cube[6],y1=cube[5])
#upcube
cube[c(2,3,5,7)] <-cube[c(2,3,5,7)]-8
#cube[c(1,4,6,8)] <-cube[c(1,4,6,8)]+40
segments(cube[6],-cube[2],x1=cube[6],y1=-cube[3])
segments(cube[8],-cube[2],x1=cube[8],y1=-cube[3])
segments(cube[1],-cube[5],x1=cube[6],y1=-cube[2])
segments(cube[8],-cube[2],x1=cube[4],y1=-cube[5])
segments(cube[6],-cube[3],x1=cube[8],y1=-cube[3])
segments(cube[6],-cube[2],x1=cube[8],y1=-cube[2])
segments(cube[1],-cube[5],x1=cube[4],y1=-cube[5])
segments(cube[4],-cube[7],x1=cube[4],y1=-cube[5])
segments(cube[4],-cube[7],x1=cube[8],y1=-cube[3])
dev.off()


# Figure 5A

x <- subset(df1, istest==1)
x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
x$time2 <- round(x$ctime/0.05)*0.05
unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
unik <- unlist(unik, recursive = T, use.names=F)
tmp  <- seq_along(x$time2)[unik]
x <- x[tmp,]
rownames(x) <- 1:nrow(x)

mark.alp <- function(x){
  temp  <- x+1
  lag   <- c(1,temp)-1.5
  unik  <- (temp*head(lag,-1))[-1]
  alpha <- length(which(unik==-1))/(length(which(unik==1))+length(which(unik==0.5)))
  #beta  <- length(which(unik==0.5))/length(which(unik==-0.5))
  return(alpha)
}
mark.bet <- function(x){
  temp  <- x+1
  lag   <- c(1,temp)-1.5
  unik  <- (temp*head(lag,-1))[-1]
  #alpha <- length(which(unik==-1))/length(which(unik==1))
  beta  <- length(which(unik==0.5))/(length(which(unik==-0.5))+length(which(unik==-1)))
  return(beta)
}

temp  <- x
marko <- aggregate(temp$cstate, by=list(temp$participant, temp$iscong), FUN=mark.alp)
marko$beta <- aggregate(temp$cstate, by=list(temp$participant, temp$iscong), FUN=mark.bet)$x
names(marko) <- c("part","iscong","alpha","beta")
marko_all <- marko

mean(marko_all$alpha)
mean(marko_all$beta)
sd(marko_all$alpha)
sd(marko_all$beta)

t.test(marko_all$alpha[0:13],marko_all$beta[0:13],paired=T)
t.test(marko_all$alpha[14:26],marko_all$beta[14:26],paired=T)

marko_all$bpref<-marko_all$beta-marko_all$alpha
mean(marko_all$bpref[1:13])
mean(marko_all$bpref[14:26])
sd(marko_all$bpref[1:13])
sd(marko_all$bpref[14:26])

temp <- marko_all[,1:3]
temp <- cbind(temp, rep(1, nrow(temp)))
names(temp) <- c("part", "cong", "value", "is.alpha")
temp <- rbind(temp, setNames(cbind(marko_all[,1:2], marko_all[,4],rep(0, nrow(temp))),names(temp)))
means <- as.data.frame(aggregate(temp$value, by=list(temp$cong,temp$is.alpha), FUN=mean, simplify=T))
names(means) <- c("cong", "is.alpha", "value")
means[,2] <- as.character(means[,2])
means <- means[order(-rank(means$cong),-rank(means$is.alpha)),]
library(plotrix)
std <- aggregate(temp$value, by=list(temp$cong,temp$is.alpha), FUN=std.error, simplify=T)
library(lattice)

segments(0.2,1,0.2,0)

segments(c(1,2,3,4), means[,3]-std[,3],
         c(1,2,3,4),means[,3] + std[,3], lwd=2)
segments(0.2,means[2,3]-std[2,3],0.2,means[2,3]+std[2,3])

# means<- means[rev(rownames(means)),]
std  <- std[rev(rownames(std)),]
std  <- std$x
std  <- std[c(1,3,2,4)]

counts <- matrix(means$value,ncol=2,byrow=T)
counts <- t(counts)

# Markov plot

setEPS()
postscript("markov_plot.eps")
plot.new()
par(mfrow=c(1,2))
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
#       widths=c(3,1), heights=c(1,2))
par(cex.main=1)
plot5A <- barplot(counts, xlab="", ylab="Probability",
                 col=c("white", "gray93", "gray77", "gray66"), beside=T,
                 ylim=c(0,max(counts)+max(std)))
segments(plot5A, counts[] - std, plot5A, counts[] + std, lwd=2)
segments(0.8,0,6.2,0, lwd=1.5)
par(cex.main=1)
title(expression(alpha), line = -17, adj=0.1)
title(expression(beta), line = -17, adj=0.3)
title(expression(alpha), line = -17, adj=0.7)
title(expression(beta), line = -17, adj=0.9)
par(cex.main=.75)
title("Compatible", line=-27, adj=0.1)
title("Incompatible", line=-27, adj=0.9)

# Plot 5B
mean(marko_all$alpha)
std.error(marko_all$alpha)
mean(marko_all$beta)
std.error(marko_all$beta)


marko_all$diff <- marko$alpha - marko_all$beta
t.test(marko_all$diff)

# show means/sds:
# alpha-congruent | beta-congruent | alpha-incongruent | beta-incongruent
# mean(marko_all[14:26,3])
# std.error(marko_all[14:26,3])
# mean(marko_all[14:26,4])
# std.error(marko_all[14:26,4])
# mean(marko_all[1:13,3])
# std.error(marko_all[1:13,3])
# mean(marko_all[1:13,4])
# std.error(marko_all[1:13,4])

t.test(marko_all[1:13,3],marko_all[1:13,4], paired = T)

# means to plot
ma <- c(mean(marko_all$alpha[14:26]),mean(marko_all$alpha[0:13]))
mb <- c(mean(marko_all$beta[14:26]),mean(marko_all$beta[0:13]))

dif<- c(mean((marko_all$beta-marko_all$alpha)[14:26]),
        mean((marko_all$beta-marko_all$alpha)[0:13]))

# t.test((marko_all$beta-marko_all$alpha)[14:26],(marko_all$beta-marko_all$alpha)[0:13])

library(plotrix)

sa <- c(std.error(marko_all$alpha[14:26]),std.error(marko_all$alpha[0:13]))

sb <- c(std.error(marko_all$beta[14:26]),
        std.error(marko_all$beta[0:13]))

sdif<- c(std.error((marko_all$alpha-marko_all$beta)[14:26]),
         std.error((marko_all$alpha-marko_all$beta)[0:13]))

plot5B <-barplot(dif, xpd = FALSE, col=c("white", "gray66"), 
                 ylab=expression(paste(beta, " - ",alpha)),
                 ylim=c(0,max(c(dif+sdif,dif-sdif))))
segments(plot5B, dif - sdif, plot5B, dif + sdif, lwd=2)
segments(0,0,2.5,0, lwd=1.5)
par(cex.main=1)
#title("b.")
par(cex.main=.75)
title("Compatible", line=-27, adj=0.1)
title("Incompatible", line=-27, adj=0.85)
dev.off()


# ANOVA and MANOVA
{
marko <- aggregate(x$cstate, by=list(x$participant, x$iscong), FUN=mark.alp)
names(marko) <- c("part","iscong","value")
marko$isalp <- rep("1",nrow(marko))

marko2 <- aggregate(x$cstate, by=list(x$participant, x$iscong), FUN=mark.bet)
names(marko2) <- c("part","iscong","value")
marko2$isalp <- rep("0",nrow(marko2))
markova <- rbind(marko,marko2)

anova1 <- aov(value~iscong*isalp, data=markova)
summary(anova1)

anova2 <- aov(value~(iscong*isalp)+Error(part), data=markova)
summary(anova2)


library(lme4)
mod1 <- manova(cbind(alpha,beta) ~ iscong, data=marko_all)
summary(mod1, test="Wilks")
summary(mod1)
}

           
           
# GCM for Necker

setwd("/Users/trbrooks/Documents/Past Semesters/Spring 2016/necker/datadirectory/necker_round1")


x <- subset(df1, istest==1)
x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
x$time2 <- round(x$ctime/0.05)*0.05
unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
unik <- unlist(unik, recursive = T, use.names=F)
tmp  <- seq_along(x$time2)[unik]
x <- x[tmp,]
rownames(x) <- 1:nrow(x)

# create walk variable
temp <- by(x$change, x$trial, FUN = cumsum, simplify = F)
temp <- unlist(temp, recursive = T, use.names=F)
x$cumu <- temp
x$rand <- (x$cstate-0.5)*2
temp <- by(x$rand, x$trial, FUN = cumsum, simplify = T)
x$walk <- unlist(temp, recursive = T, use.names=F)


# create the "skeleton" of a person-level format file.

# int <- by(x, list(x$participant, x$iscong), 
#    function(data) coefficients(lm(walk ~ time2, data = data))[[1]])
# as.numeric((names(int)))->id
# int <- unlist(int)
# cbind(id,int)->tmp.int
# colnames(tmp.int)<-c("int0","int1")
# summary(tmp.int)
# 

#x2<-x[order(x$participant),]
# print(x2)
# x2[1,]
# summary(x2)
# rate <- by(x2, list(x2$participant, x2$iscong), 
#           function(data) coefficients(lm(x$walk ~ x$time2, data = data))[[2]])
# rate <- unlist(rate)
# names(rate) <- NULL
# cbind(id,rate)->tmp.rate
# colnames(tmp.rate)<-c("rate0","rate1")
# summary(tmp.rate) 
# 
# rsq <- by(x2, list(x2$participant, x2$iscong), 
#           function(data) summary(lm(x$walk ~ x$time2, data = data))$r.squared)
# rate <- unlist(rsq)
# names(rsq) <- NULL
# cbind(id, rsq)->tmp.rsq
# colnames(tmp.rsq)<-c("rsq0","rsq1")
# summary(tmp.rsq)
# 
# 
# cbind(tmp.int,tmp.rate,tmp.rsq)->tmp
# 
# library(ggplot2)
# summary(tmp)
# histogram(tmp[,1])
# tmp

x$timeq <- x$time2^2
x$time_cu <- x$time2^3

# Build the gcm models
library(lme4)
gcm1<-lmer(walk ~ time2 * iscong + (1 + time2|participant), data=x)
summary(gcm1)
gcm2<-lmer(walk ~ time2 * iscong + timeq + (1 + time2|participant), data=x)
summary(gcm2)
gcm3<-lmer(walk ~ iscong:time2 + time2 + timeq + timeq:iscong + (1 + time2|participant), data=x)
summary(gcm3)
gcm4 <- lmer(walk ~ -1 + iscong:time2 + time2 + timeq + timeq:iscong + (1 + time2|participant), data=x)
summary(gcm4)


gcm5<-lmer(walk ~ time2 + timeq + iscong:time2 + (1 + time2|participant), data=x)
summary(gcm5)
gcm6<-lmer(walk ~ -1 + time2 + timeq + iscong:time2 + (1 + time2|participant), data=x)
summary(gcm6)


anova(gcm1, gcm2)
anova(gcm2, gcm3)
anova(gcm4, gcm5)
anova(gcm3, gcm5)
anova(gcm3,gcm4)
anova(gcm4,gcm6)

