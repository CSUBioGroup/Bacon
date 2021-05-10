#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Estimate density
#===============
library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1)    #not to overload your computer
registerDoParallel(cl)

args = commandArgs(trailingOnly=TRUE)
loop_file<-read.table(args[1],sep="\t")
valid_pet<-read.table(args[2],sep="\t")
outdir<-paste(args[3],"loop_ES_rate",sep="/")

print(loop_file)

loop_file[,"sd"]<-(loop_file$V3-loop_file$V2)+(loop_file$V6-loop_file$V5)
loop_file[,"s1"]<-loop_file$V2-loop_file$sd
loop_file[,"e1"]<-loop_file$V3+loop_file$sd
loop_file[,"s2"]<-loop_file$V5-loop_file$sd
loop_file[,"e2"]<-loop_file$V6+loop_file$sd
loop_file[,"ES"]<-0

loop_final<-loop_file[1,]
loop_final<-loop_file[-1,]
chr_list<-unique(loop_file$V1)

print(length(chr_list))

loop_final <- foreach(j=1:length(chr_list), .combine=rbind) %dopar% {
  print("entering loop")
  loop_tmp<-loop_file[loop_file$V1==as.character(chr_list[j]),]
  valid_tmp<-valid_pet[valid_pet$V1==as.character(chr_list[j]),]
  
  for(i in 1:nrow(loop_tmp)){
    tmp<-valid_tmp[valid_tmp$V1==as.character(loop_tmp[i,1]) & valid_tmp$V4==as.character(loop_tmp[i,4]) &  valid_tmp$V2>=as.character(loop_tmp[i,"s1"]) & valid_tmp$V3<=as.character(loop_tmp[i,"e1"]) & valid_tmp$V5>=as.character(loop_tmp[i,"s2"]) & valid_tmp$V6<=as.character(loop_tmp[i,"e2"]),]    
    loop_tmp[i,"ES"]<-nrow(tmp)
    print(paste("ES:",nrow(tmp)))
  }
  loop_tmp[,"ES_rate"]<-(loop_tmp$V7/loop_tmp$ES)*5
  write.table(loop_tmp,paste(outdir,as.character(chr_list[j]),sep="_"),sep="\t",quote = F,row.names = F,col.names = F)
  loop_tmp
}

stopCluster(cl)

write.table(loop_final,outdir,sep="\t",quote = F,row.names = F,col.names = F)
print("Running successfully!")

