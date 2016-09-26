# all necker results & analyses figures

# means for response latency
# for additional graphs, etc. use "necker_latency1.R"

x <- df1
x <- x[which(x$istest==0),]
x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
x$time2 <- round(x$ctime/0.05)*0.05
unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
unik <- unlist(unik, recursive = T, use.names=F)
tmp  <- seq_along(x$time2)[unik]
x <- x[tmp,]
rownames(x) <- 1:nrow(x)

tests <- subset(x, istest==0)

x$mismatch <- !x$mismatch

temp  <- aggregate(x$mismatch, by=list(x$participant, x$iscong), FUN=sum, simplify=T)
means <- aggregate(temp[,3], by=list(temp[,2]), FUN=mean, simplify=T)
stds  <- aggregate(temp[,3], by=list(temp[,2]), FUN=std.error, simplify=T)

t.test(temp$x[1:13],temp$x[14:26], paired=T)

par(mfrow=c(1,2))
plot2a <-barplot(means[,2]*20, xpd = FALSE, col=c("lightsteelblue2", "lightpink1"), ylab="milliseconds",
                ylim=c(0,max(means[,2]*20)+max(stds[,2]*20)))
segments(plot2a, means[,2]*20 - stds[,2]*20, plot2a, means[,2]*20 + stds[,2]*20, lwd=2)
abline(h=0,lty=1, lwd=2)
par(cex.main=1)
title("a.")
par(cex.main=.75)
title("     Congruent          Incongruent", line=-17)


# Mean number of switches and SE for test phase.

x <- df1
x <- x[which(x$istest==1),]
x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
x$time2 <- round(x$ctime/0.05)*0.05
unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
unik <- unlist(unik, recursive = T, use.names=F)
tmp  <- seq_along(x$time2)[unik]
x <- x[tmp,]
rownames(x) <- 1:nrow(x)


switchr <- function(inn){
  lag<-c(inn[-1],inn[length(inn)])
  return(as.numeric(inn!=lag))
}

library(reshape)
x$swi<-unlist(by(x$cstate, x$trial, FUN=switchr))
num.swi<-aggregate(x$swi, by=list(x$istest, x$iscong, x$participant), FUN=sum)
num.swi<-subset(num.swi, Group.1==1)
num.swi$Group.1 <- NULL
names(num.swi) <- c("iscong","participant","number") 

mean(num.swi$number)
std.error(num.swi$number)

mean(num.swi[which(num.swi$iscong==1),]$number)
std.error(num.swi[which(num.swi$iscong==1),]$number)

mean(num.swi[which(num.swi$iscong==0),]$number)
std.error(num.swi[which(num.swi$iscong==0),]$number)

temp  <- num.swi
means <- aggregate(temp[,3], by=list(temp[,1]), FUN=mean, simplify=T)
means <- means[order(c(1,0)),]
stds  <- aggregate(temp[,3], by=list(temp[,1]), FUN=std.error, simplify=T)
stds <- stds[order(c(1,0)),]

means
stds

#t.test(num.swi[which(num.swi$iscong==1),]$number,
#       num.swi[which(num.swi$iscong==0),]$number,
#       paired = T)


plot2b <-barplot(means[,2], xpd = FALSE, col=c("lightsteelblue2", "lightpink1"), ylab="switches",
                ylim=c(0,max(means[,2])+max(stds[,2])))
segments(plot2b, means[,2] - stds[,2], plot2b, means[,2] + stds[,2], lwd=2)
abline(h=0,lty=1, lwd=2)
par(cex.main=1)
title("b.")
par(cex.main=.75)
title(" Congruent           Incongruent", line=-17)


## PREFERENCE (total time in each state)

# Random walk plot for switches (Figure 4a)

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

# create a "random walk" version of the switching vector
switchr <- function(inn){
  lag<-c(inn[-1],inn[length(inn)])
  return(as.numeric(inn!=lag))
}

library(reshape)
x$swi<-unlist(by(x$cstate, x$trial, FUN=switchr))
x$rand <- x$swi
temp <- by(x$rand, x$trial, FUN = cumsum, simplify = T)
x$walk <- unlist(temp, recursive = T, use.names=F)
# aggregate(x$walk, by=list(x$iscong), FUN = sum, simplify = T)


# Random walk plot for preference

library(plotrix)

temp <- subset(x, iscong=="0")
temp1 <- aggregate(temp$walk, list(time = temp$time2), mean)
temp2 <- aggregate(temp$walk, list(time = temp$time2), std.error)
itraj <- temp1[-length(temp1[,1]),]
istde <- temp2[-length(temp2[,1]),]

temp <- subset(x, iscong=="1")
temp1 <- aggregate(temp$walk, list(time = temp$time2), mean)
temp2 <- aggregate(temp$walk, list(time = temp$time2), std.error)
ctraj <- temp1[-length(temp1[,1]),]
cstde <- temp2[-length(temp2[,1]),]

plot.new()
par(mfrow=c(1,1))

lims <- c(0,30)
plot(itraj$time,rep(0,length(itraj$time)),ylim=lims,type='l',
     ylab="Mean Running Total of Perceptual Switches", xlab="time (seconds)")
title("a.")
polygon(c(itraj$time,rev(itraj$time)),c(itraj$x+istde$x,
                                        rev(itraj$x-istde$x)),col=gray(0.8),border=F)
polygon(c(itraj$time,rev(itraj$time)),c(ctraj$x+cstde$x,
                                        rev(ctraj$x-cstde$x)),col=gray(0.8),border=F)
lines(itraj$time, y=itraj$x, col="dark red")
lines(ctraj$time, y=ctraj$x, cex = .8, col = "dark blue")

legend(x=100,y=6,legend=c("Incongruent   ", "Congruent", "Standard Error"),
       lwd=c(1,1,10),col=c("dark red","blue","gray"),cex=0.7)

# Figure 4B: random walk of preference

x <- df1
x$cstate[which(x$iscong==0)] <- abs(x$cstate[which(x$iscong==0)]-1)
x$time2 <- round(x$ctime/0.05)*0.05
unik <- by(x$time2, x$trial, FUN = function(x) !duplicated(x))
unik <- unlist(unik, recursive = T, use.names=F)
tmp  <- seq_along(x$time2)[unik]
x <- x[tmp,]
rownames(x) <- 1:nrow(x)
x <- subset(x, istest==1)

# create a "random walk" version of the state vector
x$rand <- (x$cstate-0.5)*2
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

plot.new()
par(mfrow=c(1,1))

lims <- c(-max(c(ctraj$x,itraj$x)), max(c(ctraj$x,itraj$x)))
plot(itraj$time,rep(0,length(itraj$time)),ylim=lims,type='l',
     ylab="Preference for Down Position", xlab="time (seconds)")
polygon(c(itraj$time,rev(itraj$time)),c(itraj$x+istde$x,
                                        rev(itraj$x-istde$x)),col=gray(0.8),border=F)
polygon(c(itraj$time,rev(itraj$time)),c(ctraj$x+cstde$x,
                                        rev(ctraj$x-cstde$x)),col=gray(0.8),border=F)
lines(itraj, col="dark red")
lines(ctraj, cex = .8, col = "dark blue")

# lo1 <- loess(itraj$time2~itraj$x)
# lo2 <-

legend(x=100,y=-170,legend=c("Incongruent", "Congruent", "Standard Error"),
       lwd=c(1,1,10),col=c("dark red","blue","gray"),cex=0.7)
# downcube
#cube <- c(0,5,11,16,3000,5000,4000,6000)
segments(5,300,x1=5,y1=200)
segments(16,300,x1=16,y1=200)
segments(0,350,x1=5,y1=300)
segments(5,200,x1=0,y1=250)
segments(5,200,x1=16,y1=200)
segments(5,300,x1=16,y1=300)
segments(16,300,x1=11,y1=350)
segments(0,350,x1=11,y1=350)
segments(0,250,x1=0,y1=350)
#upcube
segments(0,-300,x1=0,y1=-200)
segments(11,-300,x1=11,y1=-200)
segments(5,-350,x1=0,y1=-300)
segments(11,-300,x1=16,y1=-350)
segments(0,-200,x1=11,y1=-200)
segments(0,-300,x1=11,y1=-300)
segments(5,-350,x1=16,y1=-350)
segments(16,-250,x1=16,y1=-350)
segments(16,-250,x1=11,y1=-200)
title("b.")

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

# t.test(marko_all$alpha[0:13],marko_all$beta[0:13],paired=T)
# t.test(marko_all$alpha[14:26],marko_all$beta[14:26],paired=T)

temp <- marko_all[,1:3]
temp <- cbind(temp, rep(1, nrow(temp)))
names(temp) <- c("part", "iscong", "value", "is.alpha")
temp <- rbind(temp, setNames(cbind(marko_all[,1:2], marko_all[,4],rep(0, nrow(temp))),names(temp)))
means <- aggregate(temp$value, by=list(temp$iscong,temp$is.alpha), FUN=mean, simplify=T)
names(means) <- c("iscong", "is.alpha", "value")
means[,2] <- as.character(means[,2])

std <- aggregate(temp$value, by=list(temp$iscong,temp$is.alpha), FUN=std.error, simplify=T)
library(lattice)
segments(0.2,1,0.2,0)

segments(c(1,2,3,4), means[,3]-std[,3],
         c(1,2,3,4),means[,3] + std[,3], lwd=2)
segments(0.2,means[2,3]-std[2,3],0.2,means[2,3]+std[2,3])

means<- means[rev(rownames(means)),]
std  <- std[rev(rownames(std)),]
std  <- std$x
std  <- std[c(1,3,2,4)]

plot.new()
par(mfrow=c(1,2))

counts <- matrix(means$value,ncol=2,byrow=T)
counts <- t(counts)
par(mfrow=c(1,2))
par(cex.main=1)

plot5A <- barplot(counts, main="a.", xlab="", 
                 col=c("lightsteelblue2", "lightsteelblue3", "lightpink2", "lightpink3"), beside=T,
                 ylim=c(0,max(counts)+max(std)))
segments(plot1, counts[] - std, plot1, counts[] + std, lwd=2)
segments(0.8,0,6.2,0, lwd=1.5)
par(cex.main=1)
title("α     β           α     β", line = -14, adj=0.4)
par(cex.main=.75)
title("Congruent  Incongruent", line=-17)

# Plot 5B

marko_all$diff <- marko$alpha - marko_all$beta
# t.test(marko_all$diff)

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

plot5B <-barplot(dif, xpd = FALSE, col=c("lightsteelblue2", "lightpink1"), ylab="α - β",
                 ylim=c(0,max(c(dif+sdif,dif-sdif))))
segments(plot5B, dif - sdif, plot3, dif + sdif, lwd=2)
segments(0,0,2.5,0, lwd=1.5)
par(cex.main=1)
title("b.")
par(cex.main=.75)
title("Congruent  Incongruent", line=-17)



