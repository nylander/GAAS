#!/anaconda3/envs/agat/bin/Rscript

# To prepare data do:
#awk '{print $18}' maker_a1_p0_c0.gtf | sed s/[\",\;]//g >> AED_a1p0c0.csv

#Use this programm
# Just type the file containing row of aed value as argument

args <- commandArgs(TRUE)
cat("You gave",length(args),"arguments:",paste(args,collapse=" - "),"\n");

if (length(args) == 0) {
    cat("You have to give a least one file","\n");
    quit();
}
# create an empty vector
listFileName <- c();

cat("My file=",length(args),"\n");
pdf("outputPlot.pdf")
for (i in 1:length(args)) {
	

	if(i != 1){
		par(new=TRUE);
	}

	myfileName <- "";
	myfileName <- args[i] ;
	#take only file name
	myfileName <- basename(myfileName);
	#remove extensiom
	myfileName <- sub("^([^.]*).*", "\\1", myfileName);
#	cat("My file=",myfileName,"\n");

    listFileName<-c(listFileName,myfileName)
	
	myfileData=as.matrix(read.table(args[i]));
	plot(density(myfileData),xlim=c(0,1),ylim=c(0,5), col=i, xlab="", ylab="", main="")
}

legend("topright", col=(1:length(args)), lty=1, c(listFileName))
dev.off()




