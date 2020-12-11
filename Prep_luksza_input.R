args = commandArgs(trailingOnly=TRUE)
if (length(args)<2){
    stop("Two arguments must be supplied, directory and sample")
} 
#### Make Luksza Input Files ####
#path = args[2]
setwd(args[1])
i = args[2]
# Read in data
valid_data <- as.matrix(read.table(paste("Validated_",i,".txt", sep=""), header=TRUE))
data_MT <- as.matrix(read.table(paste(i, "_netCTL.out", sep=""), header=TRUE))
data_WT <- as.matrix(read.table(paste(i, ".WTmatch_netMHC.out", sep=""), header=TRUE))
data_MT_MHC <- as.matrix(read.table(paste(i, "_netMHC.out", sep=""), header=TRUE))
# Match data
data_MT <- data_MT[!duplicated(data_MT[,1]),]
total_data <- matrix(0, nrow=length(valid_data[,1]),ncol=7)

for (i in 1:length(valid_data[,1])){
  total_data[i,1] = valid_data[i,5]
  total_data[i,2] = valid_data[i,2]
  total_data[i,5] = "MUT"
  for (y in 1:length(data_WT[,1])){
    if (data_WT[y,2] == total_data[i,2] && data_WT[y,3] == total_data[i,1]){
      total_data[i,3] = data_WT[y,1]
      total_data[i,6] = data_WT[y,6]
    }
  }
  for (y in 1:length(data_MT_MHC[,1])){
    if (data_MT_MHC[y,2] == total_data[i,2] && data_MT_MHC[y,1] == total_data[i,1]){
      total_data[i,4] = data_MT_MHC[y,1]
      total_data[i,7] = data_MT_MHC[y,5]
    }
  }
}

colnames(total_data) <- cbind("genename", "id", "WT_peptide", "MT_peptide", "mut", "WT_MHC", "MT_MHC")
total_data <- as.data.frame(total_data)
total_data <- total_data[which(total_data$WT_peptide != 0),]


# Print data to file
options(max.print=1000000)
sink(file="All_lukszaprogram.in")
print(total_data)
sink(NULL)
