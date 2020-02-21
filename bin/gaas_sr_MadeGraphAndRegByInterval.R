#!/anaconda3/envs/agat/bin/Rscript

args <- commandArgs(TRUE)
cat("Il y a",length(args),"arguments:",paste(args,collapse=" - "),"\n");

if (length(args) != 5) {
cat("Il faut 5 arguments: Tab1, tab2, interval, fin, resuname","\n");
quit
}

count=as.numeric(args[3])
maxi=as.numeric(args[4])
interval=count
pdf(args[5])
while ( count < maxi) {
cat("while")
if (count == interval) {
SUP2a=as.matrix(read.table(args[1]))
SSO2a=as.matrix(read.table(args[2])) 
SUP2=SUP2a[SUP2a[,2]>0 & SUP2a[,1]>=0 & SUP2a[,1]<=count , ]
SSO2=SSO2a[SSO2a[,2]>0 & SSO2a[,1]>=0 & SSO2a[,1]<=count , ]
print(SUP2)
if(nrow(SUP2) != 0){
	cat("if1\n")
	plot(SUP2, xlim=c(0,500),ylim=c(0,1.2), col=2)
	cat("if1\n")
	regi<-lm(SUP2[,2]~SUP2[,1])
	abline(regi, col=2)
	if(nrow(SSO2) != 0){
		cat("if2\n")
		par(new=TRUE)
		plot(SSO2, xlim=c(0,500),ylim=c(0,1.2), col=3)
		rego<-lm(SSO2[,2]~SSO2[,1])
		abline(rego, col=3)
		}
	}
else { 	
	cat("else1\n")
	if(nrow(SSO2) != 0){
		cat("if3\n")
		plot(SSO2, xlim=c(0,500),ylim=c(0,1.2), col=3)
		rego<-lm(SSO2[,2]~SSO2[,1])
		abline(rego, col=3)
		}
}
count<-count+interval
cat("EndFirst1\n")
}
else {
cat("Else")
counttmp<-count+interval
SUP2new=SUP2a[SUP2a[,2]>0 & SUP2a[,1]>=count & SUP2a[,1]<=counttmp , ]
SSO2new=SSO2a[SSO2a[,2]>0 & SSO2a[,1]>=count & SSO2a[,1]<=counttmp , ]
if(nrow(SUP2new) != 0){
#	par(new=TRUE)
	plot(SUP2new, xlim=c(0,500),ylim=c(0,1.2), col=2)
	regi<-lm(SUP2new[,2]~SUP2new[,1])
	abline(regi, col=2)
	
	if(nrow(SSO2new) != 0){
		par(new=TRUE)
		plot(SSO2new, xlim=c(0,500),ylim=c(0,1.2), col=3)
		rego<-lm(SSO2new[,2]~SSO2new[,1])
		abline(rego, col=3)
		}
	}
	
else{
	if(nrow(SSO2new) != 0){
		plot(SSO2new, xlim=c(0,500),ylim=c(0,1.2), col=3)
		regu<-lm(SSO2new[,2]~SSO2new[,1])
		abline(regu, col=3)
		}
	}
	
count<-count+interval
}
}
title("mean W value according to the total number of mutation")
legend("topright", col=(2:3), lty=1, c(args[1], args[2]))
dev.off()

