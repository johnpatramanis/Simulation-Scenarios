#plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))
#class(density)



mydirectory=getwd()
setwd(mydirectory)

#Requires that you run the python file 
ComusMinMax <- read.table("COMUS_MIN_MAX",sep='\t')


############################################################
############################################################  
#COMUStats input

ComusFile <- "COMUSTATS_OUT"
con  <- file(ComusFile, open = "r")
ComusList <- list()
j=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  
  ComusList[j] <- strsplit(line,'\t')
  j <- j + 1
}

close(con)  

comusdataframe=data.frame(matrix(as.numeric(unlist(ComusList[-1])),nrow=length(ComusList)-1,byrow=T))

j=1
newcomusdistributions=list()
for (k in 1:ncol(comusdataframe)){
  
mycomusdataframe=na.omit(comusdataframe[,k])
par_min=ComusMinMax[k,1]
par_max=ComusMinMax[k,2]
  
newcomusdistributions[j]=list(density(mycomusdataframe,n=11,from=par_min,to=par_max)$y)
j=j+1
  
  
}



for (c in 1:length(newcomusdistributions)){
  
  for (w in 1:length(newcomusdistributions[[c]])){
    
  cat(paste(newcomusdistributions[[c]][w],'\t'),file='SAMPLE_FOR_ABC',append=TRUE,sep='\t')  
  }
}

cat('',file='SAMPLE_FOR_ABC',append=TRUE,sep='\n')
