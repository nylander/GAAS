#!/usr/bin/Rscript

args <- commandArgs(TRUE)
cat("Il y a",length(args),"arguments:",paste(args,collapse=" - "),"\n");

if (length(args) != 3) {
cat("Il faut 3 arguments: Tab1, tab2, resuname","\n");
quit
}

SUP2a=read.table(args[1], sep=" ")
SSO2a=read.table(args[2], sep=" ")
objp<-sapply(1:nrow(SUP2a),function(x) SUP2a[x,2]>0)
obju<-sapply(1:nrow(SSO2a),function(x) SSO2a[x,2]>0)
SUP2=SUP2a[objp,]
SSO2=SSO2a[obju,]

tiff(filename=args[3], width = 800, height = 800, units = "px", pointsize = 26, compression="lzw")
par(lwd=1) 
plot(SUP2[,2]~SUP2[,1], xlim=c(0,1),ylim=c(0,1), col=2, pch = 0, xlab="synonymous mutation rate", ylab="ω")
regi<-lm(SUP2[,2]~SUP2[,1])
abline(regi, col=2)
par(new=TRUE, lwd=1, pch = 2)
plot(SSO2[,2]~SSO2[,1], xlim=c(0,1),ylim=c(0,1), col=3, pch = 2, xlab="", ylab="")
rego<-lm(SSO2[,2]~SSO2[,1])
abline(rego, col=3)
abline(v = 170)


#title("mean W value according to \nthe total number of mutation")
legend("topright", col=(2:3), pch=(0:1), c("UPint", "SOall"))
dev.off()

