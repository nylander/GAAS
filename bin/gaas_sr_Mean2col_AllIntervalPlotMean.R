#!/anaconda3/envs/agat/bin/Rscript

args <- commandArgs(TRUE)
cat("Il y a",length(args),"arguments:",paste(args,collapse=" - "),"\n");

if (length(args) != 3) {
cat("Il faut 3 arguments: Tab1, tab2, resuname","\n");
quit
}

SUP2a=as.matrix(read.table(args[1], sep="\t", he=T))
SSO2a=as.matrix(read.table(args[2], sep="\t", he=T))
objp<-sapply(1:nrow(SUP2a),function(x) SUP2a[x,2]>0)
obju<-sapply(1:nrow(SSO2a),function(x) SSO2a[x,2]>0)
SUP2=SUP2a[objp,]
SSO2=SSO2a[obju,]

tiff(filename=args[3], width = 800, height = 800, units = "px", pointsize = 26, compression="lzw")
#pdf(args[3])
plot(SUP2, xlim=c(0,1),ylim=c(0,1.2), col=2, xlab="Number of synonymous mutation", ylab="W")
#regi<-lm(SUP2[,2]~SUP2[,1])
#abline(regi, col=2)
par(new=TRUE)
plot(SSO2, xlim=c(0,1),ylim=c(0,1.2), col=3, xlab="", ylab="")
#rego<-lm(SSO2[,2]~SSO2[,1])
#abline(rego, col=3)
title("mean W value according to \nthe total number of mutation")
legend("topright", col=(2:3), lty=1, c(args[1], args[2]))
dev.off()

