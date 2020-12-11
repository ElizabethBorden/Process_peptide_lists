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
data_MT_stab <- as.matrix(read.table(paste(i, "_netMHCstab.out", sep=""), header=TRUE))
data_WT <- as.matrix(read.table(paste(i, ".WTmatch_netMHC.out", sep=""), header=TRUE))
data_MT_MHC <- as.matrix(read.table(paste(i, "_netMHC.out", sep=""), header=TRUE))
sequence_alignment <- as.matrix(read.table(paste(i, ".epitopes.annotated.tsv", sep=""), 
                                               header=TRUE, row.names = NULL))
TCR_binding <- as.matrix(read.table("neoantigen_fitness_All.txt", header=TRUE))
# Match data
data_MT <- data_MT[!duplicated(data_MT[,1]),]
total_data <- matrix(0, nrow=length(valid_data[,1]),ncol=16)

for (i in 1:length(valid_data[,1])){
  total_data[i,1] = valid_data[i,5]
  total_data[i,16] = valid_data[i,9]
  total_data[i,2] = valid_data[i,2]
  for (y in 1:length(data_WT[,1])){
    if (data_WT[y,2] == total_data[i,2] && data_WT[y,3] == total_data[i,1]){
      total_data[i,3] = data_WT[y,1]
      total_data[i,5] = data_WT[y,6]
    }
  }
  for (y in 1:length(data_MT_MHC[,1])){
    if (data_MT_MHC[y,2] == total_data[i,2] && data_MT_MHC[y,1] == total_data[i,1]){
      total_data[i,4] = data_MT_MHC[y,1]
      total_data[i,6] = data_MT_MHC[y,5]
    }
  }
  for (y in 1:length(data_MT[,1])){
    if (data_MT[y,2] == total_data[i,2] && data_MT[y,1] == total_data[i,1]){
      total_data[i,7] = data_MT[y,4]
      total_data[i,8] = data_MT[y,5]
    }
  }
  for (y in 1:length(sequence_alignment[,1])){
    if (sequence_alignment[y,1] == total_data[i,1]){
      total_data[i,9] = sequence_alignment[y,4]
    }
  }
  for (y in 1:length(TCR_binding[,1])){
    if (TCR_binding[y,3] == total_data[i,2] && TCR_binding[y,2] == total_data[i,1]){
      total_data[i,10] = TCR_binding[y,7]
      total_data[i,11] = TCR_binding[y,8]
      total_data[i,12] = TCR_binding[y,9]
      total_data[i,13] = TCR_binding[y,10]
      total_data[i,14] = TCR_binding[y,6]
    }
  }
  for (y in 1:length(data_MT_stab[,1])){
    if (data_MT_stab[y,2] == total_data[i,2] && data_MT_stab[y,1] == total_data[i,1]){
      total_data[i,15] = data_MT_stab[y,3]
    }
  }
}

total_data <- as.data.frame(total_data)
total_data <- total_data[which(total_data$WT_peptide != 0),]


# Print data to file
options(max.print=1000000)
sink(file="Carreno_All_final_terms.in")
print(total_data)
sink(NULL)
