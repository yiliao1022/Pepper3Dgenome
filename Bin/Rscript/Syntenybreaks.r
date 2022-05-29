

jfunc <- function(x) { s <- sd(x); m <- mean(x); return(c(m+c(-1, 0, 1)*s)) }

#################For tomato scale #############################################
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/Synteny/SL4")

dat1 <- matrix(as.numeric(scan(file = 'HQ.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs1 = rowSums(dat1)
tmp1 <- apply(X=dat1, MARGIN = 2, function(x) x/rs1) 
ans1 <- apply(X = tmp1, MARGIN = 2, FUN = jfunc)

dat2 <- matrix(as.numeric(scan(file = 'RH89A.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs2 = rowSums(dat2)
tmp2 <- apply(X=dat2, MARGIN = 2, function(x) x/rs2) 
ans2 <- apply(X = tmp2, MARGIN = 2, FUN = jfunc)

dat3 <- matrix(as.numeric(scan(file = 'RH89B.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs3 = rowSums(dat3)
tmp3 <- apply(X=dat3, MARGIN = 2, function(x) x/rs3) 
ans3 <- apply(X= tmp3, MARGIN = 2, FUN = jfunc)

dat4 <- matrix(as.numeric(scan(file = 'CA59.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs4 = rowSums(dat4)
tmp4 <- apply(X=dat4, MARGIN = 2, function(x) x/rs4) 
ans4 <- apply(X = tmp4, MARGIN = 2, FUN = jfunc)


HQ <- as.numeric(scan(file = 'SL4.HQ.scale.bed'))
HQp <- HQ/sum(HQ)
RH89A <- as.numeric(scan(file = 'SL4.RH89A.scale.bed'))
RH89Ap <- RH89A/sum(RH89A)
RH89B <- as.numeric(scan(file = 'SL4.RH89B.scale.bed'))
RH89Bp <- RH89B/sum(RH89B)
SL4<- as.numeric(scan(file = 'SL4.CA59.scale.bed'))
SL4p <- SL4/sum(SL4)

plot(xlab="Solanum lycopersicum (Tomato) TADs",ylab="Synteny breakpoints [%]",-1,-1,ylim = range(HQp), xlim = c(1,20)) ##
polygon(y = t(ans1[c(1,3),]), x = c(1:20, 20:1), col = 'snow1', border = NA)
polygon(y = t(ans2[c(1,3),]), x = c(1:20, 20:1), col = 'snow2', border = NA)
polygon(y = t(ans3[c(1,3),]), x = c(1:20, 20:1), col = 'snow3', border = NA)
polygon(y = t(ans4[c(1,3),]), x = c(1:20, 20:1), col = 'snow4', border = NA)

lines(ans1[2,],col='chocolate',lty="dashed")
lines(ans2[2,],col='steelblue',lty="dashed")
lines(ans3[2,],col='seagreen4',lty="dashed")
lines(ans4[2,],col='maroon',lty="dashed")

lines(HQp, col = 'chocolate', lwd=2)
lines(RH89Ap, col = 'steelblue', lwd=2)
lines(RH89Bp, col = 'seagreen4', lwd=2)
lines(SL4p, col = 'maroon', lwd=2)

abline(v=5.5,lty="dashed")
abline(v=15.5,lty="dashed")


legend(4.5,0.06, legend=c("Eggplant (HQ)", "Potato (RH89A)", "Potato (RH89B)", "Pepper (CA59)"),
       col=c("chocolate", "steelblue","seagreen4","maroon"), lty=1:2, cex=0.8)


###########################For tomato real ##############################
#################For tomato scale #############################################
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/Synteny/SL4")

dat1 <- matrix(as.numeric(scan(file = 'HQ.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs1 = rowSums(dat1)
tmp1 <- apply(X=dat1, MARGIN = 2, function(x) x/rs1) 
ans1 <- apply(X = tmp1, MARGIN = 2, FUN = jfunc)

dat2 <- matrix(as.numeric(scan(file = 'RH89A.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs2 = rowSums(dat2)
tmp2 <- apply(X=dat2, MARGIN = 2, function(x) x/rs2) 
ans2 <- apply(X = tmp2, MARGIN = 2, FUN = jfunc)

dat3 <- matrix(as.numeric(scan(file = 'RH89B.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs3 = rowSums(dat3)
tmp3 <- apply(X=dat3, MARGIN = 2, function(x) x/rs3) 
ans3 <- apply(X= tmp3, MARGIN = 2, FUN = jfunc)

dat4 <- matrix(as.numeric(scan(file = 'CA59.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs4 = rowSums(dat4)
tmp4 <- apply(X=dat4, MARGIN = 2, function(x) x/rs4) 
ans4 <- apply(X = tmp4, MARGIN = 2, FUN = jfunc)


HQ <- as.numeric(scan(file = 'SL4.HQ.real.bed'))
HQp <- HQ/sum(HQ)
RH89A <- as.numeric(scan(file = 'SL4.RH89A.real.bed'))
RH89Ap <- RH89A/sum(RH89A)
RH89B <- as.numeric(scan(file = 'SL4.RH89B.real.bed'))
RH89Bp <- RH89B/sum(RH89B)
SL4<- as.numeric(scan(file = 'SL4.CA59.real.bed'))
SL4p <- SL4/sum(SL4)

plot(xlab="Solanum lycopersicum (Tomato) TADs",ylab="Synteny breakpoints [%]",-1,-1,ylim = range(HQp), xlim = c(1,20)) ##
polygon(y = t(ans1[c(1,3),]), x = c(1:20, 20:1), col = 'snow1', border = NA)
polygon(y = t(ans2[c(1,3),]), x = c(1:20, 20:1), col = 'snow2', border = NA)
polygon(y = t(ans3[c(1,3),]), x = c(1:20, 20:1), col = 'snow3', border = NA)
polygon(y = t(ans4[c(1,3),]), x = c(1:20, 20:1), col = 'snow4', border = NA)

lines(ans1[2,],col='chocolate',lty="dashed")
lines(ans2[2,],col='steelblue',lty="dashed")
lines(ans3[2,],col='seagreen4',lty="dashed")
lines(ans4[2,],col='maroon',lty="dashed")

lines(HQp, col = 'chocolate', lwd=2)
lines(RH89Ap, col = 'steelblue', lwd=2)
lines(RH89Bp, col = 'seagreen4', lwd=2)
lines(SL4p, col = 'maroon', lwd=2)

abline(v=5.5,lty="dashed")
abline(v=15.5,lty="dashed")


legend(4.5,0.06, legend=c("Eggplant (HQ)", "Potato (RH89A)", "Potato (RH89B)", "Pepper (CA59)"),
       col=c("chocolate", "steelblue","seagreen4","maroon"), lty=1:2, cex=0.8)



#################For pepper scale #############################################

setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/Synteny/CA59")

dat1 <- matrix(as.numeric(scan(file = 'HQ.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs1 = rowSums(dat1)
tmp1 <- apply(X=dat1, MARGIN = 2, function(x) x/rs1) 
ans1 <- apply(X = tmp1, MARGIN = 2, FUN = jfunc)

dat2 <- matrix(as.numeric(scan(file = 'RH89A.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs2 = rowSums(dat2)
tmp2 <- apply(X=dat2, MARGIN = 2, function(x) x/rs2) 
ans2 <- apply(X = tmp2, MARGIN = 2, FUN = jfunc)

dat3 <- matrix(as.numeric(scan(file = 'RH89B.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs3 = rowSums(dat3)
tmp3 <- apply(X=dat3, MARGIN = 2, function(x) x/rs3) 
ans3 <- apply(X= tmp3, MARGIN = 2, FUN = jfunc)

dat4 <- matrix(as.numeric(scan(file = 'SL4.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs4 = rowSums(dat4)
tmp4 <- apply(X=dat4, MARGIN = 2, function(x) x/rs4) 
ans4 <- apply(X = tmp4, MARGIN = 2, FUN = jfunc)


HQ <- as.numeric(scan(file = 'CA59.HQ.scale.bed'))
HQp <- HQ/sum(HQ)
RH89A <- as.numeric(scan(file = 'CA59.RH89A.scale.bed'))
RH89Ap <- RH89A/sum(RH89A)
RH89B <- as.numeric(scan(file = 'CA59.RH89B.scale.bed'))
RH89Bp <- RH89B/sum(RH89B)
CA59<- as.numeric(scan(file = 'CA59.SL4.scale.bed'))
CA59p <- CA59/sum(CA59)

plot(xlab="Capsicum annuum (Pepper) TADs",ylab="Synteny breakpoints [%]",-1,-1,ylim = range(HQp), xlim = c(1,20)) ##
polygon(y = t(ans1[c(1,3),]), x = c(1:20, 20:1), col = 'snow1', border = NA)
polygon(y = t(ans2[c(1,3),]), x = c(1:20, 20:1), col = 'snow2', border = NA)
polygon(y = t(ans3[c(1,3),]), x = c(1:20, 20:1), col = 'snow3', border = NA)
polygon(y = t(ans4[c(1,3),]), x = c(1:20, 20:1), col = 'snow4', border = NA)

lines(ans1[2,],col='chocolate',lty="dashed")
lines(ans2[2,],col='steelblue',lty="dashed")
lines(ans3[2,],col='seagreen4',lty="dashed")
lines(ans4[2,],col='maroon',lty="dashed")

lines(HQp, col = 'chocolate', lwd=2)
lines(RH89Ap, col = 'steelblue', lwd=2)
lines(RH89Bp, col = 'seagreen4', lwd=2)
lines(CA59p, col = 'maroon', lwd=2)

abline(v=5.5,lty="dashed")
abline(v=16,lty="dashed")

legend(6.5,0.08, legend=c("Eggplant (HQ)", "Potato (RH89A)", "Potato (RH89B)", "tomato (SL4)"),
       col=c("chocolate", "steelblue","seagreen4","maroon"), lty=1:2, cex=0.8)



#################For pepper real #############################################

setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/Synteny/CA59")

dat1 <- matrix(as.numeric(scan(file = 'HQ.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs1 = rowSums(dat1)
tmp1 <- apply(X=dat1, MARGIN = 2, function(x) x/rs1) 
ans1 <- apply(X = tmp1, MARGIN = 2, FUN = jfunc)

dat2 <- matrix(as.numeric(scan(file = 'RH89A.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs2 = rowSums(dat2)
tmp2 <- apply(X=dat2, MARGIN = 2, function(x) x/rs2) 
ans2 <- apply(X = tmp2, MARGIN = 2, FUN = jfunc)

dat3 <- matrix(as.numeric(scan(file = 'RH89B.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs3 = rowSums(dat3)
tmp3 <- apply(X=dat3, MARGIN = 2, function(x) x/rs3) 
ans3 <- apply(X= tmp3, MARGIN = 2, FUN = jfunc)

dat4 <- matrix(as.numeric(scan(file = 'SL4.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs4 = rowSums(dat4)
tmp4 <- apply(X=dat4, MARGIN = 2, function(x) x/rs4) 
ans4 <- apply(X = tmp4, MARGIN = 2, FUN = jfunc)


HQ <- as.numeric(scan(file = 'CA59.HQ.real.bed'))
HQp <- HQ/sum(HQ)
RH89A <- as.numeric(scan(file = 'CA59.RH89A.real.bed'))
RH89Ap <- RH89A/sum(RH89A)
RH89B <- as.numeric(scan(file = 'CA59.RH89B.real.bed'))
RH89Bp <- RH89B/sum(RH89B)
CA59<- as.numeric(scan(file = 'CA59.SL4.real.bed'))
CA59p <- CA59/sum(CA59)

plot(xlab="Capsicum annuum (Pepper) TADs",ylab="Synteny breakpoints [%]",-1,-1,ylim = range(HQp), xlim = c(1,20)) ##
polygon(y = t(ans1[c(1,3),]), x = c(1:20, 20:1), col = 'snow1', border = NA)
polygon(y = t(ans2[c(1,3),]), x = c(1:20, 20:1), col = 'snow2', border = NA)
polygon(y = t(ans3[c(1,3),]), x = c(1:20, 20:1), col = 'snow3', border = NA)
polygon(y = t(ans4[c(1,3),]), x = c(1:20, 20:1), col = 'snow4', border = NA)

lines(ans1[2,],col='chocolate',lty="dashed")
lines(ans2[2,],col='steelblue',lty="dashed")
lines(ans3[2,],col='seagreen4',lty="dashed")
lines(ans4[2,],col='maroon',lty="dashed")

lines(HQp, col = 'chocolate', lwd=2)
lines(RH89Ap, col = 'steelblue', lwd=2)
lines(RH89Bp, col = 'seagreen4', lwd=2)
lines(CA59p, col = 'maroon', lwd=2)

abline(v=5.5,lty="dashed")
abline(v=16,lty="dashed")

legend(6.5,0.08, legend=c("Eggplant (HQ)", "Potato (RH89A)", "Potato (RH89B)", "tomato (SL4)"),
       col=c("chocolate", "steelblue","seagreen4","maroon"), lty=1:2, cex=0.8)


#################For Potato scale #############################################
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/Synteny/RH89")

dat1 <- matrix(as.numeric(scan(file = 'HQ.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs1 = rowSums(dat1)
tmp1 <- apply(X=dat1, MARGIN = 2, function(x) x/rs1) 
ans1 <- apply(X = tmp1, MARGIN = 2, FUN = jfunc)

dat2 <- matrix(as.numeric(scan(file = 'SL4.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs2 = rowSums(dat2)
tmp2 <- apply(X=dat2, MARGIN = 2, function(x) x/rs2) 
ans2 <- apply(X = tmp2, MARGIN = 2, FUN = jfunc)

dat3 <- matrix(as.numeric(scan(file = 'CA59.chain.filter.tnet.synnet.bed.bed.background.scale.bed')), nc = 20)
rs3 = rowSums(dat3)
tmp3 <- apply(X=dat3, MARGIN = 2, function(x) x/rs3) 
ans3 <- apply(X = tmp3, MARGIN = 2, FUN = jfunc)

HQ <- as.numeric(scan(file = 'RH89A.HQ.scale.bed'))
HQp <- HQ/sum(HQ)
SL4 <- as.numeric(scan(file = 'RH89A.SL4.scale.bed'))
SL4p <- SL4/sum(SL4)
CA59<- as.numeric(scan(file = 'RH89A.CA59.scale.bed'))
CA59p <- CA59/sum(CA59)

plot(xlab="Solanum tuberosum (Potato) TADs",ylab="Synteny breakpoints [%]",-1,-1,ylim = c(0.03,0.08), xlim = c(1,20)) ##
polygon(y = t(ans1[c(1,3),]), x = c(1:20, 20:1), col = 'snow1', border = NA)
polygon(y = t(ans2[c(1,3),]), x = c(1:20, 20:1), col = 'snow2', border = NA)
polygon(y = t(ans3[c(1,3),]), x = c(1:20, 20:1), col = 'snow3', border = NA)

lines(ans1[2,],col='maroon',lty="dashed")
lines(ans2[2,],col='steelblue',lty="dashed")
lines(ans3[2,],col='darkorchid4',lty="dashed")

lines(HQp, col = 'maroon', lwd=2)
lines(SL4p, col = 'steelblue', lwd=2)
lines(CA59p, col = 'darkorchid4', lwd=2)

abline(v=5.5,lty="dashed")
abline(v=15.5,lty="dashed")

legend(6.5,0.08, legend=c("Eggplant (HQ)", "Tomato (SL4)", "Pepper (CA59)"),
       col=c("maroon", "steelblue","darkorchid4"), lty=1:2, cex=0.8)


#################  For potato real
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/Synteny/RH89")

dat1 <- matrix(as.numeric(scan(file = 'HQ.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs1 = rowSums(dat1)
tmp1 <- apply(X=dat1, MARGIN = 2, function(x) x/rs1) 
ans1 <- apply(X = tmp1, MARGIN = 2, FUN = jfunc)

dat2 <- matrix(as.numeric(scan(file = 'SL4.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs2 = rowSums(dat2)
tmp2 <- apply(X=dat2, MARGIN = 2, function(x) x/rs2) 
ans2 <- apply(X = tmp2, MARGIN = 2, FUN = jfunc)

dat3 <- matrix(as.numeric(scan(file = 'CA59.chain.filter.tnet.synnet.bed.bed.background.bed')), nc = 20)
rs3 = rowSums(dat3)
tmp3 <- apply(X=dat3, MARGIN = 2, function(x) x/rs3) 
ans3 <- apply(X = tmp3, MARGIN = 2, FUN = jfunc)

HQ <- as.numeric(scan(file = 'RH89A.HQ.real.bed'))
HQp <- HQ/sum(HQ)
SL4 <- as.numeric(scan(file = 'RH89A.SL4.real.bed'))
SL4p <- SL4/sum(SL4)
CA59<- as.numeric(scan(file = 'RH89A.CA59.real.bed'))
CA59p <- CA59/sum(CA59)

plot(xlab="Solanum tuberosum (Potato) TADs",ylab="Synteny breakpoints [%]",-1,-1,ylim = c(0.03,0.08), xlim = c(1,20)) ##
polygon(y = t(ans1[c(1,3),]), x = c(1:20, 20:1), col = 'snow1', border = NA)
polygon(y = t(ans2[c(1,3),]), x = c(1:20, 20:1), col = 'snow2', border = NA)
polygon(y = t(ans3[c(1,3),]), x = c(1:20, 20:1), col = 'snow3', border = NA)

lines(ans1[2,],col='maroon',lty="dashed")
lines(ans2[2,],col='steelblue',lty="dashed")
lines(ans3[2,],col='darkorchid4',lty="dashed")

lines(HQp, col = 'maroon', lwd=2)
lines(SL4p, col = 'steelblue', lwd=2)
lines(CA59p, col = 'darkorchid4', lwd=2)

abline(v=5.5,lty="dashed")
abline(v=15.5,lty="dashed")

legend(6.5,0.08, legend=c("Eggplant (HQ)", "Tomato (SL4)", "Pepper (CA59)"),
       col=c("maroon", "steelblue","darkorchid4"), lty=1:2, cex=0.8)


