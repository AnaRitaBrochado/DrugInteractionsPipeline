
load(paste0(Load_dir,"Brochado2018.RData"))

P_vals_dir = paste0(here_path,"Sensitivity Analysis/P-vals/")
Out_dir = paste0(here_path,"Sensitivity Analysis/OutData/")

#=============== Set parameters ===========

exp_fit_threshold = 0.2
range = "both" #both, complete or reduced
ROC_name = "ROC_p_val.txt"

multpTest_cutoff = 0.05
weak_interac_cutoff = 1
interac_cutoff = 0
opposite_control = F

#multpTest_cutoffs = seq(1e-04,1,by=0.01) #High performance computing recommended
multpTest_cutoffs = c(0.05,0.06,1)
for(i in 1:length(multpTest_cutoffs))
{
  roc = Sensitivity_aux(multpTest_cutoff,interac_cutoff,weak_interac_cutoff,Benchmarking_file = paste0(Load_dir,"BenchmarkingResults.txt"),
                        opposite_control=F, #both, complete or reduced
                        range = "both",exp_fit_threshold = 0.2)
  if(i==1)
  {ROC = roc} else
  {ROC = as.data.frame(rbind(ROC,roc))}
  
} # multpTest_cutoffs
row.names(ROC) = multpTest_cutoffs

#file_id = paste0(Out_dir,ROC_name)
#write.table(ROC_all,file_id,sep="\t",row.names=T,quote=FALSE,col.names = T)
#file_id = paste0(Out_dir,"ROC_p_val.txt")
#ROC = read.table(file_id,header=T, sep="\t", na.strings="NA", dec=".", strip.white=T,row.names=1)
ROC = as.data.frame(rbind(rep(0,length(ROC)),ROC))

# ========== Fix p-val, add 1sided interac ===========
multpTest_cutoff = 0.05
roc_1side = Sensitivity_aux(multpTest_cutoff,interac_cutoff,weak_interac_cutoff,Benchmarking_file = paste0(Load_dir,"BenchmarkingResults.txt"),
                      opposite_control=T, #both, complete or reduced
                      range = "both",exp_fit_threshold = 0.2)

# ========== Add interac_cutoff = 0.1 ===========
interac_cutoff = 0.1
roc_interac01 = Sensitivity_aux(multpTest_cutoff,interac_cutoff,weak_interac_cutoff,Benchmarking_file = paste0(Load_dir,"BenchmarkingResults.txt"),
                      opposite_control=T, #both, complete or reduced
                      range = "both",exp_fit_threshold = 0.2)
# ========== Add interac_cutoff = 0.2 ===========
interac_cutoff = 0.2
roc_interac02 = Sensitivity_aux(multpTest_cutoff,interac_cutoff,weak_interac_cutoff,Benchmarking_file = paste0(Load_dir,"BenchmarkingResults.txt"),
                                opposite_control=T, #both, complete or reduced
                                range = "both",exp_fit_threshold = 0.2)
# ========== Add weak_interac_cutoff = 0 ===========
interac_cutoff = 0.1
weak_interac_cutoff = 0
roc_weak_interac0 = Sensitivity_aux(multpTest_cutoff,interac_cutoff,weak_interac_cutoff,Benchmarking_file = paste0(Load_dir,"BenchmarkingResults.txt"),
                                opposite_control=T, #both, complete or reduced
                                range = "both",exp_fit_threshold = 0.2)
# ========== Add weak_interac_cutoff = 0.06 ===========
weak_interac_cutoff = 0.06
roc_weak_interac006 = Sensitivity_aux(multpTest_cutoff,interac_cutoff,weak_interac_cutoff,Benchmarking_file = paste0(Load_dir,"BenchmarkingResults.txt"),
                                    opposite_control=T, #both, complete or reduced
                                    range = "both",exp_fit_threshold = 0.2)


ROC_all = as.data.frame(rbind(ROC[grep(0.05,row.names(ROC)),],roc_1side,roc_interac01,roc_interac02,roc_weak_interac0,roc_weak_interac006))
row.names(ROC_all) = c("p-val=0.05","+1-sided interac","+interac_cutoff=0.1","or +interac_cutoff=0.2","+weakInterac=0.06","+weakInterac=0")
# ========== Plot ===========
pdf(paste0(Out_dir,"ROC_final.pdf"),useDingbats=FALSE)

colors = c("gray80","#2E97D0","#2E97D0","#2E97D0","#FFA52D","#FF9200")
shapes = c(4,2,17,6,19,1)
plot(ROC[[2]],ROC[[1]],ylim=c(0,1),xlim=c(0,1),type="l",lwd=2,col=colors[1],ylab ="TPR",xlab="FPR",main="ROC curve",cex=1.5)
points(c(0,1),c(0,1),type="l",col="grey",lty=2)
for(i in 1:length(ROC_all[[1]]))
{points(ROC_all[i,2],ROC_all[i,1],pch=shapes[i],col=colors[i],cex=1.5)}
legend("bottomright",legend = c("p-val=0.05","+1-sided interac","+interac_cutoff=0.1","or +interac_cutoff=0.2","+weakInterac=0.06","+weakInterac=0"),
       col=colors,pch=shapes,lwd=2,bty="n")

dev.off()

