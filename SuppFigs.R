
#==== ED Fig 2 & 3 ========
#========== Plot ditsributions of Q1 and Q3 for all bugs =========
par(mfrow=c(2,2),mar=c(3,3,3,3))
for(sp in 1:length(Multi_Species))
{
  species = Multi_Species[sp]
  bugs = Species_strains[[match(species,names(Species_strains))]]
  
  bug1 = bugs[1]
  bug2 = bugs[2]
  
  bug1_q = Interaction_Qs[[match(bug1,names(Interaction_Qs))]][c(11,1,2)]
  bug2_q = Interaction_Qs[[match(bug2,names(Interaction_Qs))]][c(11,1,2)]
  
  plot(density(as.numeric(bug1_q[[1]])),main=paste0(species," median"),xlim=c(-0.8,0.8),ylim=c(0,20))
  points(density(as.numeric(bug2_q[[1]])),type="l",lty=2)
  points(density(as.numeric(bug1_q[[2]])),type="l",col=col_synergy,lwd=2)
  points(density(as.numeric(bug2_q[[2]])),type="l",col=col_synergy,lwd=2,lty=2)
  points(density(as.numeric(bug1_q[[3]])),type="l",col=col_antagonism,lwd=2)
  points(density(as.numeric(bug2_q[[3]])),type="l",col=col_antagonism,lwd=2,lty=2)
  legend(x=0.1,y=20,legend = bugs,lty=c(1,2),bty="n",cex=0.9,x.intersp=0.5,y.intersp=0.8)
  
  text(x=-0.8,y=18,adj=0,cex=0.9,paste0(bug1," medians"))
  text(x=-0.8,y=16.5,adj=0,cex=0.8,round(median(as.numeric(bug1_q[[2]])),3),col=col_synergy)
  text(x=-0.8,y=15,adj=0,cex=0.8,round(median(as.numeric(bug1_q[[1]])),4))
  text(x=-0.8,y=13.5,adj=0,cex=0.8,round(median(as.numeric(bug1_q[[3]])),3),col=col_antagonism)
  
  text(x=-0.8,y=11.5,adj=0,cex=0.9,paste0(bug2," medians"))
  text(x=-0.8,y=10,adj=0,cex=0.8,round(median(as.numeric(bug2_q[[2]])),3),col=col_synergy)
  text(x=-0.8,y=8.5,adj=0,cex=0.8,round(median(as.numeric(bug2_q[[1]])),4))
  text(x=-0.8,y=7,adj=0,cex=0.8,round(median(as.numeric(bug2_q[[3]])),3),col=col_antagonism)
}

#==================== Replicate correlation ====================

for(b in 1:length(Fitness))
{
  bug = names(Fitness[b])
  fitness = Fitness[[b]]
  
  control_data = fitness[grep("Control",names(fitness))]
  batch_info = All_bugs_batch_info[grep(bug,All_bugs_batch_info[[1]]),]
  batch_info = batch_info[match(names(control_data),as.character(batch_info[[2]])),]
  batch_info = as.character(batch_info[[3]])
  
  unique_batch = unique(batch_info)
  Rs = c()
  for(bt in 1:length(unique_batch))
  {
    batch = unique_batch[bt]
    batch_data = control_data[grep_exact(batch,batch_info)]
    
    cor_matrix = cor(batch_data,use="complete.obs")
    cor_matrix = cor_matrix[lower.tri(cor_matrix)]
    Rs = c(Rs,cor_matrix)
  }
  
  if(b==1)
  {Control_corrs = list()}
  Control_corrs[[length(Control_corrs)+1]] = Rs
  
}
names(Control_corrs) = names(Fitness)

#Correlate wells with plates across plate
for(b in 1:length(Fitness))
{
  bug = names(Fitness[b])
  fitness = Fitness[[b]]
  
  fitness_single = fitness[grep("SingleDrug",row.names(fitness)),]
  fitness = fitness[-grep("LB",row.names(fitness)),]
  fitness = fitness[-grep("SingleDrug",row.names(fitness)),]
  
  rep2_data = fitness[grep(" rep 2",row.names(fitness)),]
  rep1_data = fitness[-grep(" rep 2",row.names(fitness)),]
  
  row.names(rep2_data) = substr(row.names(rep2_data),1,nchar(row.names(rep2_data))-6)
  rep2_data = rep2_data[match(row.names(rep1_data),row.names(rep2_data)),]
  
  rep1_data = as.data.frame(rbind(rep1_data,fitness_single[1:3,]))
  rep2_data = as.data.frame(rbind(rep2_data,fitness_single[4:6,]))
  
  rep1_data = as.data.frame(t(rep1_data))
  rep2_data = as.data.frame(t(rep2_data))
  
  Rs = c()
  names_Rs = c()
  medians = c()
  for(w in 1:length(rep1_data))
  {
    flag1 = 0
    flag2 = 0
    if(length(rep1_data[[w]]) > length(grep(T,is.na(rep1_data[[w]]))))
    {flag1=1}
    if(length(rep2_data[[w]]) > length(grep(T,is.na(rep2_data[[w]]))))
    {flag2=1}
    
    if(flag1==1 && flag2==1)
    {
      r = try(cor(rep1_data[[w]],rep2_data[[w]],use="complete.obs"))
      Rs = c(Rs,r)
      names_Rs = c(names_Rs,names(rep1_data)[w])
      
      rep1_med = median(rep1_data[[w]],na.rm=T)
      rep2_med = median(rep2_data[[w]],na.rm=T)
      med = min(rep1_med,rep2_med)
      medians = c(medians,med)
    }
  }
  
  Rs = as.data.frame(cbind(Rs,medians))
  row.names(Rs) = names_Rs
  rm(names_Rs)
  
  if(b==1)
  {Wells_corrs = list()}
  Wells_corrs[[length(Wells_corrs)+1]] = Rs
  
}
names(Wells_corrs) = names(Fitness)
for(b in 1:length(Wells_corrs))
{
  Rs = Wells_corrs[[b]]
  
  if(length(grep("Error",Rs[[1]]))>0)
  {
    names_Rs = row.names(Rs)
    names_Rs = names_Rs[-grep("Error",Rs[[1]])]
    medians = as.character(Rs[[2]])
    medians = medians[-grep("Error",Rs[[1]])]
    medians = as.numeric(medians)
    rs = as.character(Rs[[1]])
    rs = rs[-grep("Error",Rs[[1]])]
    rs = as.numeric(as.matrix(rs))
    Rs=as.data.frame(cbind(rs,medians))
    row.names(Rs) = names_Rs
    rm(names_Rs)
  }
  
  Rs = Rs[-grep(T,Rs[[2]]<0.1),] 
  Wells_corrs[[b]]=Rs
}

#Plot R versus median growth
par(mfrow=c(3,2),mar = c(4,4,4,4))
for(b in 1:length(Wells_corrs))
{
  bug = names(Wells_corrs)[b]
  bug_color = unlist(Bug_colors)[b]
  Rs = Wells_corrs[[b]]
  plot(Rs[[2]],Rs[[1]],xlab="Median growth well",ylab="Pearson correlation",ylim=c(0,1),xlim=c(0,1),main=bug,
       cex=0.7,pch=19,col=bug_color)
  text(x=0.0,y=0.2,paste0("n=",length(Rs[[2]])," wells"),adj=0)
  points(c(-0.5,1.5),c(median(Rs[[1]]),median(Rs[[1]])),type="l",lty=2)
}

#Plotting
par(mfrow=c(1,1))
plot_data=list()
for(b in 1:length(Control_corrs))
{
  plot_data[[length(plot_data)+1]] = Control_corrs[[b]]
  plot_data[[length(plot_data)+1]] = Wells_corrs[[b]][[1]]
}
here_border = rep(unlist(Bug_colors),each=2)
here_colors = here_border
here_colors[seq(1,length(here_colors),by=2)]="white"
boxplot(plot_data,ylim=c(0,1),las=2,border=here_border,pch=19,cex=0.5,
        col = here_colors,main="Replicate correlation",names=rep(names(Control_corrs),each=2),
        ylab = "Pearson correlation")

ns_open = unlist(lapply(Control_corrs,FUN=length))
text(x=rep(5,6),y=rev(seq(0.1,0.6,by=0.1)),paste0("nr correlations open boxes =",ns_open),col=here_colors[c(2,4,6,8,10,12)])

#to add to Results_list
mean_corr_all_replicates = mean(unlist(plot_data))
fraction_noisy_wells = na_well_count/well_count

Results_list[[length(Results_list)+1]] = mean_corr_all_replicates
names(Results_list)[length(Results_list)] = "mean_corr_all_replicates"
Results_list[[length(Results_list)+1]] = fraction_noisy_wells
names(Results_list)[length(Results_list)] = "fraction_noisy_wells"

#=========== Plot robust fitness againts 12 single drug wells ===========
for(b in 1:length(Fitness))
{
  bug = names(Fitness[b])
  fitness = Fitness[[b]]
  bug_color = unlist(Bug_colors)[b]
  
  single_drug_wells = fitness[grep("Single",row.names(fitness)),]
  median_single_wells = sapply(single_drug_wells,FUN=median,na.rm=T)
  
  robust_fitness = All_bugs_batch_info[grep(bug,All_bugs_batch_info[[1]]),]
  robust_fitness = robust_fitness[match(names(single_drug_wells),as.character(robust_fitness[[2]])),]
  
  if(b==1)
  {
    plot(median_single_wells,robust_fitness[[4]],ylab = "Robust fitness estimator",xlab = "median 12 single drug wells",
         main="Comparing query drugs finess estimators",
         cex=0.5,pch=19,col=bug_color)
  } else
  {points(median_single_wells,robust_fitness[[4]],cex=0.5,pch=19,col=bug_color)}
  
  query = as.data.frame(cbind(median_single_wells,robust_fitness[[4]]))
  row.names(query) = c()
  
  if(b==1)
  {All_query = query} else
  {All_query = as.data.frame(rbind(All_query,query))}
  
}

#to add to Results_list
robust_fingle_fitness_correlation = cor(All_query[[1]],All_query[[2]])
Results_list[[length(Results_list)+1]] = robust_fingle_fitness_correlation
names(Results_list)[length(Results_list)] = "robust_fingle_fitness_correlation"

#================= Plot all Bliss scores distributions per bug =============
par(mfrow=c(3,3))
for(b in 1:length(Interac_distrib))
{
  bug = names(Interac_distrib)[b]
  color = unlist(Bug_colors)[b]
  plot(density(unlist(Interac_distrib[[b]]),na.rm=T),main="Interactions score",col=color,lwd=2)
  text(x=-0.9,y=6,paste0("n=",length(na.omit(unlist(Interac_distrib[[b]])))),adj=0)
}




#==================== Expected fitness vs Bliss scores ========
par(mfrow=c(1,1),mar=c(6,6,6,6))
b=1
bug = names(Interac_distrib)[b]
bug_interac_distrib = Interac_distrib[[b]]
bug_expect_fitness = Exp_fit_dist[[b]]

plot(unlist(bug_expect_fitness),unlist(bug_interac_distrib),cex=0.2,,pch=19,col="grey",
     xlab="Expected fitness",ylab="Bliss scores",main="Blind Spots")
points(unlist(bug_expect_fitness)[unlist(bug_interac_distrib)>0.2],unlist(bug_interac_distrib)[unlist(bug_interac_distrib)>0.2],
       col=col_antagonism,cex=0.2,pch=19)
points(unlist(bug_expect_fitness)[unlist(bug_interac_distrib)<(-0.2)],unlist(bug_interac_distrib)[unlist(bug_interac_distrib)<(-0.2)],
       col=col_synergy,cex=0.2,pch=19)
points(c(-2,2),c(0,0),type="l")
points(c(0.2,0.2),c(-1.5,0),type="l")
points(c(0.8,0.8),c(0,1.5),type="l")
text(x=0.1,y=-1,"Synergy\nblind spot",srt=90,adj=0)
text(x=0.8,y=0.2,"Antagonism\nblind spot",srt=90,adj=0)

#==== Supplementary Fig 3 ========
#======= plot nr interactions per drugs per bug =========

#nr interactions per drug. Weak & strong interactions taken into account

for(sp in 1:length(Multi_Species))
{
  species = Multi_Species[sp]
  sp_hits = Spieces_conservation[[sp]]
  bugs = Species_strains[[match(species,names(Species_strains))]]
  
  bug1_hits = sp_hits[1]
  bug1_hits = sp_hits[grep(F,is.na(bug1_hits[[1]])),]
  weak = sp_hits[grep_exact("weak_conserved",sp_hits[[3]]),]
  bug1_hits = as.data.frame(rbind(bug1_hits,weak[grep(T,is.na(weak[[1]])),]))
  
  bug2_hits = sp_hits[2]
  bug2_hits = sp_hits[grep(F,is.na(bug2_hits[[1]])),]
  weak = sp_hits[grep_exact("weak_conserved",sp_hits[[3]]),]
  bug2_hits = as.data.frame(rbind(bug2_hits,weak[grep(T,is.na(weak[[2]])),]))
  
  if(sp==1)
  {ns = c()}
  ns = c(ns,length(bug1_hits[[1]]),length(bug2_hits[[1]]))
  
  #Bug1
  #sp_hits = sp_hits[c(grep_exact("conserved",sp_hits[[3]]),grep_exact("weak_conserved",sp_hits[[3]])),]
  sp_hits = bug1_hits
  interaction = sp_hits[[1]]
  interaction[is.na(interaction)] = sp_hits[is.na(interaction),2]
  interaction = as.character(as.vector(interaction))
  drug1 = unlist(strsplit(row.names(sp_hits),split="_"))[seq(1,2*length(row.names(sp_hits)),by=2)]
  drug2 = unlist(strsplit(row.names(sp_hits),split="_"))[seq(2,2*length(row.names(sp_hits)),by=2)]
  Hits = as.data.frame(cbind(drug1,interaction,drug2,row.names(sp_hits)))
  names(Hits) = c("DrugA","Interaction","DrugB","DrugA_B")
  
  #Count interaction per drug
  drugs = unique(c(as.character(Hits[[1]]),as.character(Hits[[3]])))
  drug_counts = sort(unlist(lapply(sapply(drugs,FUN=grep,x = as.character(Hits[[4]])),FUN=length)))
  
  all_drugs = as.character(Attr_table[[1]])
  bug_counts = rep(0,length(all_drugs))
  bug_counts[match(names(drug_counts),all_drugs)]=drug_counts
  
  bug_counts = as.data.frame(bug_counts)
  row.names(bug_counts) = all_drugs
  names(bug_counts) = bugs[1]
  
  if(sp==1)
  {Drug_counts = bug_counts} else
  {Drug_counts = as.data.frame(cbind(Drug_counts,bug_counts))}
  
  #Bug2
  #sp_hits = sp_hits[c(grep_exact("conserved",sp_hits[[3]]),grep_exact("weak_conserved",sp_hits[[3]])),]
  sp_hits = bug2_hits
  interaction = sp_hits[[1]]
  interaction[is.na(interaction)] = sp_hits[is.na(interaction),2]
  interaction = as.character(as.vector(interaction))
  drug1 = unlist(strsplit(row.names(sp_hits),split="_"))[seq(1,2*length(row.names(sp_hits)),by=2)]
  drug2 = unlist(strsplit(row.names(sp_hits),split="_"))[seq(2,2*length(row.names(sp_hits)),by=2)]
  Hits = as.data.frame(cbind(drug1,interaction,drug2,row.names(sp_hits)))
  names(Hits) = c("DrugA","Interaction","DrugB","DrugA_B")
  
  #Count interaction per drug
  drugs = unique(c(as.character(Hits[[1]]),as.character(Hits[[3]])))
  drug_counts = sort(unlist(lapply(sapply(drugs,FUN=grep,x = as.character(Hits[[4]])),FUN=length)))
  
  all_drugs = as.character(Attr_table[[1]])
  bug_counts = rep(0,length(all_drugs))
  bug_counts[match(names(drug_counts),all_drugs)]=drug_counts
  
  bug_counts = as.data.frame(bug_counts)
  row.names(bug_counts) = all_drugs
  names(bug_counts) = bugs[2]
  
  Drug_counts = as.data.frame(cbind(Drug_counts,bug_counts))
} 

#Clear from drugs not queried for specific species

Drug_counts[match("Flucytosine",row.names(Drug_counts)),1:4]=NA
Drug_counts[match(c("PMS","Pseudomonic acid","Pyocyanin"),row.names(Drug_counts)),5:6] = NA

par(mfrow=c(1,1),mar=c(6,6,6,6))
for(b in 1:length(Drug_counts))
{
  if(b==1)
  {
    plot(density(Drug_counts[[b]],na.rm=T),col=unlist(Bug_colors)[b],ylim=c(0,0.12),xlim=c(0,40),lwd=2,
         main = "Interactions per drug, strong & weak")
  } else
  {points(density(Drug_counts[[b]],na.rm=T),col=unlist(Bug_colors)[b],type="l",lwd=2)}
  #hist(Drug_counts[[b]],breaks = 10)
}
legend("topright",legend=ns,lwd=2,col=unlist(Bug_colors))
med_interactions = sapply(Drug_counts,FUN=median,na.rm=T)
text(x=rep(15,6),y=rev(seq(0.06,0.11,by=0.01)),paste0("median int=",med_interactions),col=unlist(Bug_colors))

Drug_counts = as.data.frame(t(Drug_counts))
here_col= colorRampPalette(c("black","#558BC8"))
heatmap.2(as.matrix(Drug_counts),trace="none",col=here_col,key = T,Rowv = T,Colv=F,dendrogram="row",
          main="Strong & weak interactions",na.color = "grey")
rm(here_col)

#============= Test for higher nr synergies/drug in membrane drugs =========
Absolute_interactions_syn = Absolute_interactions[grep("Syn",Absolute_interactions$direction),]
combs = as.character(Absolute_interactions_syn[[1]])
drug1 = unlist(strsplit(combs,split="_"))[seq(1,2*length(combs),by=2)]
drug2 = unlist(strsplit(combs,split="_"))[seq(2,2*length(combs),by=2)]
syn_drug = table(c(drug1,drug2))

memb_drugs = as.character(Attr_table[[1]])[grep("membrane",as.character(Attr_table$Master_Class))]
memb_syn_drug = syn_drug[match(memb_drugs,names(syn_drug))]
syn_drug = syn_drug[-match(names(memb_syn_drug),names(syn_drug))]

memb_col = as.character(Attr_table[match("membrane",as.character(Attr_table$Master_Class)),4])
plot(density(syn_drug),ylab="nr synergies per drug",lwd=2, main="Membrane participate in more synergies than the rest",xlim=c(0,120))
points(density(memb_syn_drug),type="l",lwd=2,col=memb_col)

test = wilcox.test(syn_drug,memb_syn_drug)
text(x=60,y=0.015,paste0("Wilcoxon p-val=",test$p.value),adj=0)
legend("bottomright",legend=c("no membrane drugs","membrane drugs"),lwd=2,col=c("black",memb_col))

rm(Absolute_interactions_syn,combs,drug1,drug12,drug2,syn_drug,master_c_drugs,memb_col)
#======= Plot nr interactions versus lowest fitness in screen ===========
file_id = paste0(Load_dir,"Lowest_single_fitness.txt")
Bugs_single_min = read.table(file_id,header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
names(Bugs_single_min) = substr(names(Bugs_single_min),1,(nchar(names(Bugs_single_min))-4))

Total_ints_counts = c()
for(b in 1:length(Bugs_single_min))
{
  bug = names(Bugs_single_min)[b]
  fitness = Bugs_single_min[[b]]
  bug_color = unlist(Bug_colors)[b]
  
  bug_ints = as.vector(Absolute_interactions[grep(bug,as.vector(Absolute_interactions$bug)),1])
  ints_counts = c()
  for(d in 1:length(Bugs_single_min[[1]]))
  {
    drug = row.names(Bugs_single_min)[d]
    count = length(grep(drug,bug_ints))
    
    if(length(count)==0)
    {count=NA}
    
    ints_counts = c(ints_counts,count) 
  }
  
  if(b==1)
  {plot(fitness,ints_counts,pch=19,col=bug_color,ylab="Nr interactions",xlab="Lowest single fitness in screen",ylim=c(0,40))} else
  {points(fitness,ints_counts,pch=19,col=bug_color)}
  
  Total_ints_counts = c(Total_ints_counts,ints_counts)
}
r = cor(as.vector(as.matrix(Bugs_single_min)),Total_ints_counts,use="complete.obs")
a = cor.test(as.vector(as.matrix(Bugs_single_min)),Total_ints_counts)
p_val = a$p.value
text(x=0.8,y=35,paste0("R=",round(r,2),adj=0))
text(x=0.8,y=30,paste0("p-val=",p_val,adj=0))


#==== Supplementary Fig 4 ========
#================ Benchmarking plots ========
par(mfrow=c(1,1),mar=c(6,6,6,6))

file_id = paste0(Load_dir,"BenchmarkingResults.txt")
Benchmarking = read.table(file_id,sep="\t",header = T)

Screen = rep(NA,length(Benchmarking[[1]]))
Conservation = rep(NA,length(Benchmarking[[1]]))
for(c in 1:length(Benchmarking[[1]]))
{
  comb = as.character(Benchmarking$Combination[c])
  bug = as.character(Benchmarking$Strain[c])
  rev = as.character(Benchmarking$Reverse[c])
  
  sp_hits = Absolute_interactions[grep(bug,as.character(Absolute_interactions$bug)),]
  
  if(length(grep(comb,sp_hits$combination))>0)
  {Screen[c] = as.character(sp_hits[grep(comb,sp_hits$combination),2])}
  
  if(length(grep(rev,sp_hits$combination))>0)
  {Screen[c] = as.character(sp_hits[grep(rev,sp_hits$combination),2])}
  
  if(length(grep(comb,sp_hits$combination))>0)
  {Conservation[c] = as.character(sp_hits[grep(comb,sp_hits$combination),3])}
  
  if(length(grep(rev,sp_hits$combination))>0)
  {Conservation[c] = as.character(sp_hits[grep(rev,sp_hits$combination),3])}
}

Screen[grep("Antag",Screen)] = "A"
Screen[grep("Syn",Screen)] = "S"
Screen[grep(T,is.na(Screen))] = "N"

Benchmarking = as.data.frame(cbind(Benchmarking,Screen,Conservation))
Benchmarking = Benchmarking[c(1,2,3,6,4,5)]

Hit_classif = rep(NA,length(Benchmarking[[1]]))
dummy = paste0(as.character(Benchmarking$Benchmarking),as.character(Benchmarking$Screen))
Hit_classif[grep("NN",dummy)]="TN"
Hit_classif[grep("AA",dummy)]="TP"
Hit_classif[grep("SS",dummy)]="TP"
Hit_classif[grep("NS",dummy)]="FP"
Hit_classif[grep("NA",dummy)]="FP"
Hit_classif[grep("SN",dummy)]="FN"
Hit_classif[grep("AN",dummy)]="FN"

Benchmarking = as.data.frame(cbind(Benchmarking,Hit_classif))
Counts = table(Hit_classif)
Precision = round(Counts[4]/(Counts[2]+Counts[4]),2)
Recall = round(Counts[4]/(Counts[4]+Counts[1]),2)

Screen_total_S = as.numeric(Abs_counts_summary$Synergy)[7]
Screen_total_A = as.numeric(Abs_counts_summary$Antagonism)[7]
Screen_total_N = as.numeric(Abs_counts_summary$Total_probed_ns)[7] - as.numeric(Abs_counts_summary$Synergy)[7] - as.numeric(Abs_counts_summary$Antagonism)[7]

pie(c(Screen_total_A,Screen_total_N,Screen_total_S),col =c(col_antagonism,"white",col_synergy),labels = c("A","N","S"))
pie(table(Benchmarking$Benchmarking),col =c(col_antagonism,"white",col_synergy))

pie(table(Hit_classif))
text(x=0.6,y=1,paste0("recall=",Recall),adj=0)
text(x=0.6,y=-0.8,paste0("precision=",Precision),adj=0)
text(x=0.6,y=-1,paste0("n=",length(Hit_classif)),adj=0)

here_colors = unlist(Bug_colors[c(1,3,2)])
barplot(table(as.character(Benchmarking$Strain)),col=here_colors,ylim=c(0,60),ylab="nr benchmarked combinations",las=2)

benchmar_conservation = as.data.frame(table(as.character(Benchmarking$Conservation)))
row.names(benchmar_conservation) = as.character(benchmar_conservation[[1]])
benchmar_conservation = benchmar_conservation[-c(1)]
benchmar_conservation[3,1] = benchmar_conservation[3,1]+benchmar_conservation[4,1]
benchmar_conservation = as.data.frame(t(benchmar_conservation))
benchmar_conservation = benchmar_conservation[-c(4)]
benchmar_conservation = as.data.frame(cbind(benchmar_conservation,c(length(Benchmarking[[1]])-sum(benchmar_conservation[1,]))))
names(benchmar_conservation)[5] = "Neutral"
benchmar_conservation = benchmar_conservation[c(5,1,4,3,2)]

barplot(t(benchmar_conservation),beside = F,col=c("white","dodgerblue4","lightskyblue","gray50","gray80"),xlim = c(0,250), xlab="# interactions",
        horiz=T,main="Classification of benchmarking interactions\naccording to conservation")

benchmar_conservation_weakCons = table(Benchmarking[grep("weak_conserved",Benchmarking$Conservation),grep("Hit_classif",names(Benchmarking))])

barplot((benchmar_conservation_weakCons),beside = F,col=c("white","dodgerblue4","lightskyblue","gray50","gray80"),xlim = c(0,250), xlab="# interactions",
        horiz=T,main="TP vs FN within weak conserved interactions")
text(x=100,y=4, benchmar_conservation_weakCons[[4]])
text(x=100,y=2, benchmar_conservation_weakCons[[2]])
