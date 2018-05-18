
#=============== Set the stage ========
get_conc_grad <- function(c)
{
  step = c/7
  v=c(8:1)
  v=v*step
  v = c-v
  v = round(c(v[2:8],c),5)
  return(v)
}

#=============== Benchmarking CBs ========

file_id = paste0(Load_dir,"Benchmarking_fitness_updated.txt")
Benchmarking_fit = read.table(file_id,sep="\t",header = F,row.names = 1)
Benchmarking_fit = as.data.frame(t(Benchmarking_fit))

file_id = paste0(Load_dir,"Benchmarking_database_updated.txt")
Benchmarking_database = read.table(file_id,sep="\t",header = T)

horiz_drugs = as.character(get_feature(names(Benchmarking_fit),feature = "Horizontal",Benchmarking_database))
vert_drugs = as.character(get_feature(names(Benchmarking_fit),feature = "Vertical",Benchmarking_database))
horiz_startCs = as.numeric(get_feature(names(Benchmarking_fit),feature = "Starting_Horizontal",Benchmarking_database))
vert_startCs = as.numeric(get_feature(names(Benchmarking_fit),feature = "Starting_Vertical",Benchmarking_database))
colors = as.character(get_feature(names(Benchmarking_fit),feature = "Color",Benchmarking_database))
colors = paste0("#",colors)
bugs = as.character(get_feature(names(Benchmarking_fit),feature = "Bug",Benchmarking_database))
replicates = as.character(get_feature(names(Benchmarking_fit),feature = "Replicate",Benchmarking_database))

pdf(paste0(Out_dir,"Benchmarking_CBs.pdf"))
for(c in 1:length(Benchmarking_fit))
{
  data = Benchmarking_fit[[c]]
  cb_id = names(Benchmarking_fit)[c]
  data_mat = as.data.frame(matrix(data,ncol=8))
  row.names(data_mat) = rev(as.character(get_conc_grad(vert_startCs[c])))
  names(data_mat) = as.character(get_conc_grad(horiz_startCs[c]))
  
  grad_color = colorRampPalette(c("white",colors[c]))
  
  heatmap.2(as.matrix(data_mat),trace="n",Rowv=F,Colv=F,dendrogram="none",col=grad_color,breaks=seq(0,1,by=0.05),
            mar = c(6,6),xlab = horiz_drugs[c],ylab = vert_drugs[c],na.color = "grey",
            main=names(Benchmarking_fit)[c],key=F)
}

dev.off()

#=============== Clinical isolates CBs ========

file_id = paste0(Load_dir,"ClinicIsolates_fitness.txt")
ClinicIsolates_fit = read.table(file_id,sep="\t",header = F,row.names = 1)
ClinicIsolates_fit = as.data.frame(t(ClinicIsolates_fit))

file_id = paste0(Load_dir,"ClinicIsolates_database.txt")
ClinicIsolates_database = read.table(file_id,sep="\t",header = T)

horiz_drugs = as.character(get_feature(names(ClinicIsolates_fit),feature = "Horiz_drug",ClinicIsolates_database))
vert_drugs = as.character(get_feature(names(ClinicIsolates_fit),feature = "Vert_drug",ClinicIsolates_database))
horiz_startCs = as.numeric(get_feature(names(ClinicIsolates_fit),feature = "Horiz_startC",ClinicIsolates_database))
vert_startCs = as.numeric(get_feature(names(ClinicIsolates_fit),feature = "Vert_startC",ClinicIsolates_database))
colors = as.character(get_feature(names(ClinicIsolates_fit),feature = "Color",ClinicIsolates_database))
colors = paste0("#",colors)
bugs = as.character(get_feature(names(ClinicIsolates_fit),feature = "bug",ClinicIsolates_database))
replicates = as.character(get_feature(names(ClinicIsolates_fit),feature = "Replicate",ClinicIsolates_database))

pdf(paste0(Out_dir,"ClinicIsolates_CBs.pdf"))
for(c in 1:length(ClinicIsolates_fit))
{
  data = ClinicIsolates_fit[[c]]
  cb_id = names(ClinicIsolates_fit)[c]
  data_mat = as.data.frame(matrix(data,ncol=8))
  row.names(data_mat) = rev(as.character(get_conc_grad(vert_startCs[c])))
  names(data_mat) = as.character(get_conc_grad(horiz_startCs[c]))
  
  grad_color = colorRampPalette(c("white",colors[c]))
  
  heatmap.2(as.matrix(data_mat),trace="n",Rowv=F,Colv=F,dendrogram="none",col=grad_color,breaks=seq(0,1,by=0.05),
            mar = c(6,6),xlab = horiz_drugs[c],ylab = vert_drugs[c],na.color = "grey",
            main=names(ClinicIsolates_fit)[c],key=F)
}

#********* plot replicate correlation *********

colors = as.character(get_feature(names(ClinicIsolates_fit),feature = "Color",ClinicIsolates_database))
colors = paste0("#",colors)
colors = unique(colors)

rep1 = ClinicIsolates_fit[seq(1,length(ClinicIsolates_fit),by=2)]
rep2 = ClinicIsolates_fit[seq(2,length(ClinicIsolates_fit),by=2)]

rep1_EC = rep1[grep("EC",names(rep1))]
rep2_EC = rep2[grep("EC",names(rep2))]

rep1_KP = rep1[grep("KP",names(rep1))]
rep2_KP = rep2[grep("KP",names(rep2))]

plot(as.vector(as.matrix(rep1_EC)),as.vector(as.matrix(rep2_EC)),ylab="replicate 2",xlab="replicate 1",
     cex=0.5,pch=19,col=colors[1],ylim=c(0,1.2),xlim=c(0,1.2))
n = length(as.vector(as.matrix(rep1_EC)))
r = round(cor(as.vector(as.matrix(rep1_EC)),as.vector(as.matrix(rep2_EC)),use="complete.obs"),3)
text(x=0.2,y=1,paste0("n=",n),adj=0)
text(x=0.2,y=0.9,paste0("R=",r),adj=0)

plot(as.vector(as.matrix(rep1_KP)),as.vector(as.matrix(rep2_KP)),ylab="replicate 2",xlab="replicate 1",
     cex=0.5,pch=19,col=colors[2],ylim=c(0,1.2),xlim=c(0,1.2))
n = length(as.vector(as.matrix(rep1_KP)))
r = round(cor(as.vector(as.matrix(rep1_KP)),as.vector(as.matrix(rep2_KP)),use="complete.obs"),3)
text(x=0.2,y=1,paste0("n=",n),adj=0)
text(x=0.2,y=0.9,paste0("R=",r),adj=0)

dev.off()

#Clean up
rm(horiz_drugs,vert_drugs,horiz_startCs,vert_startCs,colors,bugs,replicates,
   file_id,data,cb_id,data_mat,grad_color,n,r,rep1,rep2,rep1_EC,rep1_KP,rep2_EC,rep2_KP)



