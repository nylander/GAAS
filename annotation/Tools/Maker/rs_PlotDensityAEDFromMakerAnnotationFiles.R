#!/sw/R/3.0.2/bin/Rscript

# usage: script.R gtf file1.gtf
# You can use one or many gff or gtf files (should be well formatted)
# /!\ Do not forget to load R

# Jacques Dainat 04/2014 - bils.se

# get arguments
args <- commandArgs(TRUE)
cat("You gave",length(args),"arguments:",paste(args,collapse=" - "),"\n");

if (length(args) < 1) {
    cat("You have to give a least one file. Do nont forget to load R module.","\n");
    quit();
}
# create an empty vector
listlegend <- c();

nbFile=length(args);
cat("Number file(s) =",nbFile,"\n");
pdf("outputPlot.pdf")

for (i in 1:nbFile) {
    
    cat("file studied =",args[i],"\n");
    
    #If more than one file to plot
    if(i != 1){
        par(new=TRUE);
    }
    myfileName <- "";
    myfileName <- args[i] ;
    #take only file name
    myfileName <- basename(myfileName);
    #get extension
    myExt <- strsplit(myfileName, "\\.")[[1]];
    myExt <- myExt[2];
    
    #remove extensiom
    myfileName <- sub("^([^.]*).*", "\\1", myfileName);
    
    # call bash command
    if(myExt == "gtf"){
        command=paste('awk  \'!($12 in arr){cpt++;arr[$12]++;arrScore[cpt]=$18}END{for (x in arrScore) if(arrScore[x] != "") print arrScore[x]}\' ',args[i],' | sed s/[\\\",\\;]//g')
    }else if (myExt == "gff"){
        command=paste('awk \'{if($3=="mRNA") {split($9,a,";"); if(!(a[3] in arr)){ cpt++;arr[a[3]]++;arrScore[cpt]=a[4]}}} END {for (x in arrScore) if(arrScore[x] != "") print arrScore[x]}\' ',args[i],' | sed s/[\\\",\\;AED\\_\\=]//g')
    }
    else {
        cat("The extension ",myExt,"is not recognized by the program. Use only gtf or/and gff extension !\n");quit();
    }
    myData <- system(command, intern=TRUE)
    myGoodData=as.matrix(as.numeric(myData))
    
    legendInfo=paste(myfileName,"(",length(myGoodData),"mRNAs )")
    listlegend<-c(listlegend,legendInfo)
    
    #make plot
    if (nbFile == 1){
        plot(density(myGoodData),xlim=c(0,1), col=i, xlab="AED score", main="")
    }
    else{plot(density(myGoodData),xlim=c(0,1),ylim=c(0,18), col=i, xlab="AED score", main="")} # You can modify the value 18 according to the value in your graph
}
# Add Title
title(main="AED distribution")

#Add Legend
legend("topright", col=(1:nbFile), lty=1, c(listlegend))

#END
dev.off()
