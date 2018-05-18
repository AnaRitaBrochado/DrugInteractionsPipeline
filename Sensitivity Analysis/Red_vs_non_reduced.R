
load(paste0(Load_dir,"Brochado2018.RData"))

P_vals_dir = paste0(here_path,"Sensitivity Analysis/P-vals/")
Out_dir = paste0(here_path,"Sensitivity Analysis/OutData/")
#====== Iterate over exp_fit_threshold =====

multpTest_cutoff = 0.05
interac_cutoff = 0.1
weak_interac_cutoff = 0.06
range = "both" #"complete or reduced
ROC_name = "ROC_ExpFitSens.txt"

exp_fit_thresholds = c(0,0.1,0.2,0.3,0.4,0.5)
for(i in 1:length(exp_fit_thresholds))
{
  exp_fit_threshold=exp_fit_thresholds[i] 
  
  #=============== Get strong synergies & antagonisms =========missing plots check line 945-87===========
  exp_fit_threshold_syn = exp_fit_threshold
  exp_fit_threshold_antag = 1- exp_fit_threshold
  
  #Obtain Q1, Q3 & p-vals
  for(b in 1:length(Interac_distrib))
  {
    bug = names(Interac_distrib)[b]
    bug_Interac_distrib = Interac_distrib[[b]]
    #bug_pvals = P_values[[b]]
    #names(bug_pvals) = c("CompleteConcRange","AntagConcRange","SynConcRange","AntagPot","SynPot")
    
    file_id= paste0(P_vals_dir,"P-vals_",bug,"_",exp_fit_threshold,".txt")
    bug_pvals = read.table(file_id,header=T, sep="\t", na.strings="NA", dec=".", strip.white=T)
    bug_pvals = bug_pvals[c(1,3,5,2,4)]
    names(bug_pvals) = c("CompleteConcRange","AntagConcRange","SynConcRange","AntagPot","SynPot")
    
    combs = names(bug_Interac_distrib)
    bug_pvals = bug_pvals[match(combs,row.names(bug_pvals)),]
    
    bug_Q1 = as.vector(unlist(lapply(bug_Interac_distrib,FUN=quantile,probs=c(0.25))))
    bug_Q3 = as.vector(unlist(lapply(bug_Interac_distrib,FUN=quantile,probs=c(0.75))))
    bug_int_score = bug_Q3
    bug_int_score[grep(T,abs(bug_Q1) > abs(bug_Q3))] = bug_Q1[grep(T,abs(bug_Q1) > abs(bug_Q3))]
    bug_median = as.vector(unlist(lapply(bug_Interac_distrib,FUN=quantile,probs=c(0.5))))
    
    #Reduce to relevant conc
    bug_exp_fit = Exp_fit_dist[[b]]
    
    bug_Interac_distrib_Syn = list() #Only for drug concenentrations relevant for Synergy
    bug_Interac_distrib_Antag = list() #Only for drug concenentrations relevant for Antagonism
    for(c in 1:length(bug_Interac_distrib))
    {
      comb_exp_fit = bug_exp_fit[[c]]
      comb_exp_fit[comb_exp_fit>1] = 1
      comb_Interac_distrib = bug_Interac_distrib[[c]]
      
      comb_Interac_distrib_syn = comb_Interac_distrib[comb_exp_fit>exp_fit_threshold_syn]
      comb_Interac_distrib_antag = comb_Interac_distrib[comb_exp_fit<exp_fit_threshold_antag]
      
      bug_Interac_distrib_Syn[[length(bug_Interac_distrib_Syn)+1]] = comb_Interac_distrib_syn
      bug_Interac_distrib_Antag[[length(bug_Interac_distrib_Antag)+1]] =  comb_Interac_distrib_antag
    }
    names(bug_Interac_distrib_Syn) = combs
    names(bug_Interac_distrib_Antag) = combs
    bug_Q1_reduced = as.vector(unlist(lapply(bug_Interac_distrib_Syn,FUN=quantile,probs=c(0.25),na.rm=T)))
    bug_Q3_reduced = as.vector(unlist(lapply(bug_Interac_distrib_Antag,FUN=quantile,probs=c(0.75),na.rm=T)))
    
    bug_interactions = as.data.frame(cbind(bug_Q1,bug_Q3,bug_Q1_reduced,bug_Q3_reduced,bug_pvals,bug_int_score,bug_median))
    row.names(bug_interactions) = combs
    
    if(b==1)
    {Interaction_Qs = list()}
    Interaction_Qs[[length(Interaction_Qs)+1]] = bug_interactions
    
  }
  names(Interaction_Qs) = names(Interac_distrib)
  #clean up
  rm(bug_Interac_distrib,bug_exp_fit,bug_Q1,bug_Q3,bug_pvals,bug_Q1_reduced,bug_Q3_reduced,c,b)
  rm(comb_exp_fit,comb_Interac_distrib,comb_Interac_distrib_syn,comb_Interac_distrib_antag,bug_Interac_distrib_Syn,
     bug_Interac_distrib_Antag,bug_int_score)
  
  #Correct for multiple testing: fdr & classify interactions
  for(b in 1:length(Interaction_Qs))
  {
    bug = names(Interaction_Qs)[b]
    bug_interactions = Interaction_Qs[[b]]
    
    #Correct for multiple testing
    bug_interactions[5:7] = sapply(bug_interactions[5:7],FUN=p.adjust,method = "BH")
    
    #Assign interactions
    bug_interactions_sign = bug_interactions[bug_interactions$CompleteConcRange<=multpTest_cutoff,]
    #bug_interactions_sign = bug_interactions[bug_interactions$AntagConcRange<=multpTest_cutoff,]
    
    bug_interactions_sign = as.data.frame(rbind(bug_interactions_sign,
                                                bug_interactions[bug_interactions$SynConcRange<=multpTest_cutoff,]))
    bug_interactions_sign = as.data.frame(rbind(bug_interactions_sign,
                                                bug_interactions[bug_interactions$AntagConcRange<=multpTest_cutoff,]))
    
    bug_interactions_sign = bug_interactions_sign[grep(F,duplicated(bug_interactions_sign)),]
    bug_interactions_sign = bug_interactions_sign[!is.na(bug_interactions_sign[[1]]),]
    
    #synergy
    syn = bug_interactions_sign[bug_interactions_sign$bug_Q1 <= (-interac_cutoff),]
    syn = syn[syn$CompleteConcRange <= multpTest_cutoff,]
    syn_red =  bug_interactions_sign[bug_interactions_sign$bug_Q1_reduced <= (-interac_cutoff),]
    syn_red = syn_red[syn_red$SynConcRange <= multpTest_cutoff,]
    
    complt_specific = grep(F,row.names(syn) %in% row.names(syn_red))
    hit_classif = rep("both",length(syn[[1]]))
    hit_classif[complt_specific] = "complete_specific"
    syn = as.data.frame(cbind(syn,hit_classif))
    
    reduce_specific = grep(F,row.names(syn_red) %in% row.names(syn))
    hit_classif = rep("both",length(syn_red[[1]]))
    hit_classif[reduce_specific] = "reduced_specific"
    syn_red = as.data.frame(cbind(syn_red,hit_classif))
    
    syn = as.data.frame(rbind(syn,syn_red))
    syn = syn[grep(F,duplicated(syn)),]
    syn = syn[!is.na(syn[[1]]),]
    syn = as.data.frame(cbind(syn,rep("Synergy",length(syn[[1]]))))
    names(syn)[length(syn)] = c("direction")
    rm(syn_red)
    
    #antagonism
    antag = bug_interactions_sign[bug_interactions_sign$bug_Q3 >= (interac_cutoff),]
    antag = antag[antag$CompleteConcRange <= multpTest_cutoff,]
    antag_red = bug_interactions_sign[bug_interactions_sign$bug_Q3_reduced >= (interac_cutoff),]
    antag_red = antag_red[antag_red$AntagConcRange <= multpTest_cutoff,]
    
    complt_specific = grep(F,row.names(antag) %in% row.names(antag_red))
    hit_classif = rep("both",length(antag[[1]]))
    hit_classif[complt_specific] = "complete_specific"
    antag = as.data.frame(cbind(antag,hit_classif))
    
    reduce_specific = grep(F,row.names(antag_red) %in% row.names(antag))
    hit_classif = rep("both",length(antag_red[[1]]))
    hit_classif[reduce_specific] = "reduced_specific"
    antag_red = as.data.frame(cbind(antag_red,hit_classif))
    
    antag = as.data.frame(rbind(antag,antag_red))
    antag = antag[grep(F,duplicated(antag)),]
    antag = antag[!is.na(antag[[1]]),]
    antag = as.data.frame(cbind(antag,rep("Antagonism",length(antag[[1]]))))
    names(antag)[length(antag)] = c("direction")
    rm(antag_red)
    
    #Check for combinations with both antag & syn and remove from both
    opposite_int = row.names(antag)[grep(T,row.names(antag) %in% row.names(syn))]
    if(length(opposite_int)>0)
    {
      antag = antag[-match(opposite_int,row.names(antag)),]
      syn = syn[-match(opposite_int,row.names(syn)),]
    }
    
    #Create a single table with interactions per bug
    bug_hits = as.data.frame(rbind(syn,antag))
    if(b==1)
    {Bugs_hits = list()}
    Bugs_hits[[length(Bugs_hits)+1]] = bug_hits
    
  }
  names(Bugs_hits) = names(Interaction_Qs)
  
  #clean up
  rm(bug_hits,bug_interactions,bug_interactions_sign,b,syn,antag,opposite_int)
  
  #============== Restrict to reduced/complete range ===========
  
  for(b in 1:length(Bugs_hits))
  {
    bug_hits = Bugs_hits[[b]]
    
    hit_classif = as.character(bug_hits$hit_classif)
    range_classification_counts = as.data.frame(t(rep(0,3)))
    names(range_classification_counts) = c("both","complete_specific","reduced_specific")
    
    range_classification_counts[match(names(table(hit_classif)),names(range_classification_counts))] = table(hit_classif)
    
    if(range == "complete" && length(grep("reduced_specific",hit_classif))>0)
    {bug_hits = bug_hits[-grep("reduced_specific",hit_classif),]}
    if(range == "reduced" && length(grep("complete_specific",hit_classif))>0)
    {bug_hits = bug_hits[-grep("complete_specific",hit_classif),]}
    
    Bugs_hits[[b]] = bug_hits
    
    if(b==1)
    {Range_classification_counts = list()}
    Range_classification_counts[[length(Range_classification_counts)+1]] = range_classification_counts
    
  }
  names(Range_classification_counts) = names(Bugs_hits)
  #=============== Count detectable synergy & antagonism ====================
  
  for(b in 1:length(P_values))
  {
    bug = names(P_values)[b]
    bug_pvals = P_values[[b]]
    
    min_wells = min(bug_pvals[grep(T,bug_pvals$Antag_p_val <= multpTest_cutoff),4])
    possible_antagonsim = length(bug_pvals[(bug_pvals$Pos_pot >= min_wells),1])
    min_wells = min(bug_pvals[grep(T,bug_pvals$Syn_p_val <= multpTest_cutoff),5])
    possible_synergy = length(bug_pvals[(bug_pvals$Neg_pot >= min_wells),1])
    
    detectable_syn_antag = as.data.frame(c(possible_synergy,possible_antagonsim))
    row.names(detectable_syn_antag) = c("detectableSyn","detectableAntag")
    names(detectable_syn_antag) = bug
    
    if(b==1)
    {DetectSynAntag = detectable_syn_antag} else
    {DetectSynAntag = as.data.frame(cbind(DetectSynAntag,detectable_syn_antag))}
  }
  DetectSynAntag = as.data.frame(t(DetectSynAntag))
  #clean up
  rm(b,bug_pvals,min_wells,possible_antagonsim,possible_synergy,detectable_syn_antag)
  
  #=============== Conservation within species ================
  for(sp in 1:length(Multi_Species))
  {
    species = Multi_Species[sp]
    bugs = Species_strains[[match(species,names(Species_strains))]]
    bug1 = bugs[1]
    bug2 = bugs[2] 
    
    bug1_interac = Interaction_Qs[[match(bug1,names(Interaction_Qs))]]
    bug2_interac = Interaction_Qs[[match(bug2,names(Interaction_Qs))]]
    
    bug1_expFit=Exp_fit_dist[[match(bug1,names(Exp_fit_dist))]]
    bug2_expFit=Exp_fit_dist[[match(bug2,names(Exp_fit_dist))]]
    
    bug1_hits = Bugs_hits[[grep(bug1,names(Bugs_hits))]]
    bug2_hits = Bugs_hits[[grep(bug2,names(Bugs_hits))]]
    all_hits = unique(c(row.names(bug1_hits),row.names(bug2_hits)))
    
    # bug1_hits_restruct = bug1_hits[c(grep("Syn",bug1_hits$direction),grep("Antag",bug1_hits$direction)),]
    # interac_score = c(bug1_hits[grep("Syn",bug1_hits$direction),match("bug_Q1",names(bug1_hits))],
    #                   bug1_hits[grep("Antag",bug1_hits$direction),match("bug_Q3",names(bug1_hits))])
    # bug1_hits_restruct = as.data.frame(cbind(bug1_hits_restruct,interac_score))  
    # bug1_hits_restruct = bug1_hits_restruct[match(c("direction","interac_score"),names(bug1_hits_restruct))]
    # 
    # bug2_hits_restruct = bug2_hits[c(grep("Syn",bug2_hits$direction),grep("Antag",bug2_hits$direction)),]
    # interac_score = c(bug2_hits[grep("Syn",bug2_hits$direction),match("bug_Q1",names(bug2_hits))],
    #                   bug2_hits[grep("Antag",bug2_hits$direction),match("bug_Q3",names(bug2_hits))])
    # bug2_hits_restruct = as.data.frame(cbind(bug2_hits_restruct,interac_score))  
    # bug2_hits_restruct = bug2_hits_restruct[match(c("direction","interac_score"),names(bug2_hits_restruct))]
    
    #new - took reduced Q1 and Q3 - only relevant expfit wells
    bug1_hits_restruct = bug1_hits[c(grep("Syn",bug1_hits$direction),grep("Antag",bug1_hits$direction)),]
    interac_score = c(bug1_hits[grep("Syn",bug1_hits$direction),match("bug_Q1_reduced",names(bug1_hits))],
                      bug1_hits[grep("Antag",bug1_hits$direction),match("bug_Q3_reduced",names(bug1_hits))])
    bug1_hits_restruct = as.data.frame(cbind(bug1_hits_restruct,interac_score))  
    bug1_hits_restruct = bug1_hits_restruct[match(c("direction","interac_score"),names(bug1_hits_restruct))]
    
    bug2_hits_restruct = bug2_hits[c(grep("Syn",bug2_hits$direction),grep("Antag",bug2_hits$direction)),]
    interac_score = c(bug2_hits[grep("Syn",bug2_hits$direction),match("bug_Q1_reduced",names(bug2_hits))],
                      bug2_hits[grep("Antag",bug2_hits$direction),match("bug_Q3_reduced",names(bug2_hits))])
    bug2_hits_restruct = as.data.frame(cbind(bug2_hits_restruct,interac_score))  
    bug2_hits_restruct = bug2_hits_restruct[match(c("direction","interac_score"),names(bug2_hits_restruct))]
    
    dir_bug1 = rep(NA,length(all_hits))
    dir_bug2 = rep(NA,length(all_hits))
    int_score_bug1 = rep(NA,length(all_hits))
    int_score_bug2 = rep(NA,length(all_hits))
    
    dir_bug1[match(row.names(bug1_hits_restruct),all_hits)] = as.character(bug1_hits_restruct[[1]])
    dir_bug2[match(row.names(bug2_hits_restruct),all_hits)] = as.character(bug2_hits_restruct[[1]])
    int_score_bug1[match(row.names(bug1_hits_restruct),all_hits)] = bug1_hits_restruct[[2]]
    int_score_bug2[match(row.names(bug2_hits_restruct),all_hits)] = bug2_hits_restruct[[2]]
    
    species_conserved = as.data.frame(cbind(dir_bug1,dir_bug2,int_score_bug1,int_score_bug2))
    row.names(species_conserved) = all_hits
    
    #Assess conservation
    conserved = rep("not_conserved",length(all_hits)) 
    conserved[grep(T,dir_bug1 == dir_bug2)] = "conserved"
    conserved[grep(paste0("Antagonism","Synergy"),paste0(dir_bug1,dir_bug2))] = "specific"
    conserved[grep(paste0("Synergy","Antagonism"),paste0(dir_bug1,dir_bug2))] = "specific"
    
    #Find integraction scores for weak interactions
    #bug1
    strong_bug1 = row.names(species_conserved)[grep(T,is.na(species_conserved$dir_bug2))]
    strong_bug1_bug2 = bug2_interac[match(strong_bug1,row.names(bug2_interac)),1:2]
    strong_bug1_interac_bug2 = strong_bug1_bug2[[2]]
    strong_bug1_interac_bug2[grep(T,abs(strong_bug1_bug2[[1]]) > abs(strong_bug1_bug2[[2]]))] = strong_bug1_bug2[grep(T,abs(strong_bug1_bug2[[1]]) > abs(strong_bug1_bug2[[2]])),1]
    
    int_score_bug2 = as.numeric(as.vector(species_conserved[[4]]))
    int_score_bug2[match(strong_bug1,row.names(species_conserved))] = strong_bug1_interac_bug2
    species_conserved[4] = int_score_bug2
    rm(strong_bug1,strong_bug1_bug2,strong_bug1_interac_bug2)
    
    #bug2
    strong_bug2 = row.names(species_conserved)[grep(T,is.na(species_conserved$dir_bug1))]
    strong_bug2_bug1 = bug1_interac[match(strong_bug2,row.names(bug1_interac)),1:2]
    strong_bug2_interac_bug1 = strong_bug2_bug1[[2]]
    strong_bug2_interac_bug1[grep(T,abs(strong_bug2_bug1[[1]]) > abs(strong_bug2_bug1[[2]]))] = strong_bug2_bug1[grep(T,abs(strong_bug2_bug1[[1]]) > abs(strong_bug2_bug1[[2]])),1]
    
    int_score_bug1 = as.numeric(as.vector(species_conserved[[3]]))
    int_score_bug1[match(strong_bug2,row.names(species_conserved))] = strong_bug2_interac_bug1
    species_conserved[3] = int_score_bug1
    rm(strong_bug2,strong_bug2_bug1,strong_bug2_interac_bug1)
    
    #Synergy
    strong_bug1_syn = species_conserved[grep("Syn",species_conserved$dir_bug1),]
    weak_conserved_syn_b2 = strong_bug1_syn[grep(T,strong_bug1_syn$int_score_bug2 <= (-weak_interac_cutoff)),]
    weak_conserved_syn_b2 = weak_conserved_syn_b2[grep(T,is.na(weak_conserved_syn_b2$dir_bug2)),]
    strong_bug2_syn = species_conserved[grep("Syn",species_conserved$dir_bug2),]
    weak_conserved_syn_b1 = strong_bug2_syn[grep(T,strong_bug2_syn$int_score_bug1 <= (-weak_interac_cutoff)),]
    weak_conserved_syn_b1 = weak_conserved_syn_b1[grep(T,is.na(weak_conserved_syn_b1$dir_bug1)),]
    
    weak_conserved_syn = c(row.names(weak_conserved_syn_b1),row.names(weak_conserved_syn_b2))
    conserved[match(weak_conserved_syn,all_hits)] = "weak_conserved"
    rm(strong_bug1_syn,weak_conserved_syn_b2,strong_bug2_syn,weak_conserved_syn_b1,weak_conserved_syn)
    
    #Antagonism
    strong_bug1_antag = species_conserved[grep("Antag",species_conserved$dir_bug1),]
    weak_conserved_antag_b2 = strong_bug1_antag[grep(T,strong_bug1_antag$int_score_bug2 >= (weak_interac_cutoff)),]
    weak_conserved_antag_b2 = weak_conserved_antag_b2[grep(T,is.na(weak_conserved_antag_b2$dir_bug2)),]
    strong_bug2_antag = species_conserved[grep("Antag",species_conserved$dir_bug2),]
    weak_conserved_antag_b1 = strong_bug2_antag[grep(T,strong_bug2_antag$int_score_bug1 >= (weak_interac_cutoff)),]
    weak_conserved_antag_b1 = weak_conserved_antag_b1[grep(T,is.na(weak_conserved_antag_b1$dir_bug1)),]
    
    weak_conserved_antag = c(row.names(weak_conserved_antag_b1),row.names(weak_conserved_antag_b2))
    conserved[match(weak_conserved_antag,all_hits)] = "weak_conserved"
    rm(strong_bug1_antag,weak_conserved_antag_b2,strong_bug2_antag,weak_conserved_antag_b1,weak_conserved_antag)
    
    #Find non-comparable interactions
    not_conserved_combs = all_hits[grep("not_conserved",conserved)]
    sp_p_vals = c()
    for(c in 1:length(bug1_expFit))
    {
      comb = names(bug1_expFit)[c]
      comb_expfit_bug1 = bug1_expFit[[c]]
      comb_expfit_bug2 = bug2_expFit[[match(comb,names(bug2_expFit))]]
      
      WT = wilcox.test(comb_expfit_bug1,comb_expfit_bug2)
      p_val = WT$p.value
      sp_p_vals = c(sp_p_vals,p_val)
      
    }
    sp_p_vals_BH = p.adjust(sp_p_vals,"BH")
    
    not_conserved_pvals = as.data.frame(sp_p_vals_BH[match(not_conserved_combs,names(bug1_expFit))])
    row.names(not_conserved_pvals) = not_conserved_combs
    names(not_conserved_pvals) = "Hits_pvals_BH"
    
    not_comparable = row.names(not_conserved_pvals)[grep(T,not_conserved_pvals[[1]] <= 0.05)]
    conserved[match(not_comparable,all_hits)] = "not_comparable"
    
    #Assign specific interactions
    conserved[grep(paste0("Antagonism","Synergy"),paste0(dir_bug1,dir_bug2))] = "specific"
    conserved[grep(paste0("Synergy","Antagonism"),paste0(dir_bug1,dir_bug2))] = "specific"
    
    species_conserved = as.data.frame(cbind(species_conserved,conserved))
    species_conserved = species_conserved[c(1,2,5,3,4)]
    
    if(sp==1)
    {Spieces_conservation = list()}
    Spieces_conservation[[length(Spieces_conservation)+1]] = species_conserved
    
    file_id = paste0(Out_dir,species,"_speciesInterac.txt")
    write.table(species_conserved,file_id,sep="\t",row.names=T,quote=FALSE,col.names = T)
  }
  names(Spieces_conservation) = Multi_Species
  
  #Count conservation and plot
  for(sp in 1:length(Multi_Species))
  {
    species = Multi_Species[sp]
    bugs = Species_strains[[match(species,names(Species_strains))]]
    
    sp_conser_table = Spieces_conservation[[sp]]
    
    sp_counts = c(length(grep_exact("conserved",sp_conser_table[[3]])),
                  length(grep_exact("weak_conserved",sp_conser_table[[3]])),
                  length(grep_exact("specific",sp_conser_table[[3]])),
                  length(grep_exact("not_conserved",sp_conser_table[[3]])),
                  length(grep_exact("not_comparable",sp_conser_table[[3]]))
    )
    
    sp_counts = as.data.frame(t(sp_counts))
    names(sp_counts) = c("Conserved","Weak conserved","Specific","Not conserved","Not comparable")
    row.names(sp_counts) = species
    
    if(sp==1)
    {Strain_comp_counts = sp_counts} else
    {Strain_comp_counts = as.data.frame(rbind(Strain_comp_counts,sp_counts))}
    
    #plot
    int_score_bug1 = sp_conser_table$int_score_bug1
    int_score_bug2 = sp_conser_table$int_score_bug2
    conserved = sp_conser_table$conserved
    
    # plot(int_score_bug1,int_score_bug2,cex=0.7,pch=19,col="grey",
    #      main=paste0("Conserved interactions ",species),ylim=c(-1,1),xlim=c(-1,1),ylab=bugs[[2]],xlab=bugs[[1]])
    # points(c(-10,10),c(0,0),type="l")
    # points(c(0,0),c(-10,10),type="l")
    # points(c(-10,10),c(-10,10),type="l",lty=2)
    # points(int_score_bug1[grep_exact("weak_conserved",conserved)],int_score_bug2[grep_exact("weak_conserved",conserved)],cex=0.7,pch=19,col="lightskyblue")
    # points(int_score_bug1[grep_exact("conserved",conserved)],int_score_bug2[grep_exact("conserved",conserved)],
    #        cex=0.7,pch=19,col="dodgerblue4")
    # points(int_score_bug1[grep_exact("specific",conserved)],int_score_bug2[grep_exact("specific",conserved)],cex=0.7,pch=19,col="red3")
    # r = round(cor(unlist(int_score_bug1),unlist(int_score_bug2)),2)
    # text(x=-1,y=0.9,adj=0,paste0("R=",r))
    # text(x=-1,y=0.8,adj=0,paste0("n=",length(conserved)))
    # 
    # legend("bottomright",legend=c("Conserved","Weak conserved","Strain specific","Not conserved"),pch=19,col=c("dodgerblue4","lightskyblue","red3","grey"),bty="n")
    
  }
  
  # barplot(sapply(as.data.frame(t(Strain_comp_counts)),FUN=sum)
  #         ,ylim=c(0,700),main="Strain conservation",col="grey",
  #         lwd=1:2, angle=45, density=seq(10,10,5))
  # barplot(t(Strain_comp_counts[1:4]),col=c("dodgerblue4","lightskyblue","red3","grey82"),add=T)
  # legend("topright",legend=c("Conserved","Weak conserved","Strain specific","Not conserved","Not comparable"),
  #        fill=c("dodgerblue4","lightskyblue","red3","grey","grey"),bty="n",cex=0.8)
  
  #======== Absolute counts for all bugs ===========
  flag=0
  for(sp in 1:length(Multi_Species))
  {
    species = Multi_Species[sp]
    sp_hits = Spieces_conservation[[sp]]
    
    keys = c(paste0("Antagonism","_","conserved"),paste0("Synergy","_","conserved"),
             paste0("Antagonism","_","weak_conserved"),paste0("Synergy","_","weak_conserved"),
             paste0("Antagonism","_","not_conserved"),paste0("Synergy","_","not_conserved"),
             paste0("Antagonism","_","specific"),paste0("Synergy","_","specific"),
             paste0("Antagonism","_","not_comparable"),paste0("Synergy","_","not_comparable")
    )
    
    #for b=1
    bug = unlist(Species_strains[[sp]])[1]
    temp = paste0(sp_hits[[1]],"_",sp_hits[[3]])
    
    absolute_counts = c()
    for(k in 1:length(keys))
    {
      key = keys[k]
      absolute_counts = c(absolute_counts,length(grep_exact(key,temp)))
    }
    
    absolute_counts = as.data.frame((absolute_counts))
    names(absolute_counts)=bug
    row.names(absolute_counts) = keys
    
    if(sp == 1 && flag==0)
    {Abs_counts = absolute_counts} else
    {Abs_counts = as.data.frame(cbind(Abs_counts,absolute_counts))}
    flag=1
    
    #for b=2
    bug = unlist(Species_strains[[sp]])[2]
    temp = paste0(sp_hits[[2]],"_",sp_hits[[3]])
    
    absolute_counts = c()
    for(k in 1:length(keys))
    {
      key = keys[k]
      absolute_counts = c(absolute_counts,length(grep_exact(key,temp)))
    }
    
    absolute_counts = as.data.frame((absolute_counts))
    names(absolute_counts)=bug
    row.names(absolute_counts) = keys
    
    Abs_counts = as.data.frame(cbind(Abs_counts,absolute_counts))
    
  }
  rm(flag,keys,bug,temp,sp,b,absolute_counts,species,sp_hits)
  
  Abs_counts[grep("Antagonism_weak_conserved",row.names(Abs_counts)),1:2] = sum(Abs_counts[grep("Antagonism_weak_conserved",row.names(Abs_counts)),1:2])
  Abs_counts[grep("Antagonism_weak_conserved",row.names(Abs_counts)),3:4] = sum(Abs_counts[grep("Antagonism_weak_conserved",row.names(Abs_counts)),3:4])
  Abs_counts[grep("Antagonism_weak_conserved",row.names(Abs_counts)),5:6] = sum(Abs_counts[grep("Antagonism_weak_conserved",row.names(Abs_counts)),5:6])
  Abs_counts[grep("Synergy_weak_conserved",row.names(Abs_counts)),1:2] = sum(Abs_counts[grep("Synergy_weak_conserved",row.names(Abs_counts)),1:2])
  Abs_counts[grep("Synergy_weak_conserved",row.names(Abs_counts)),3:4] = sum(Abs_counts[grep("Synergy_weak_conserved",row.names(Abs_counts)),3:4])
  Abs_counts[grep("Synergy_weak_conserved",row.names(Abs_counts)),5:6] = sum(Abs_counts[grep("Synergy_weak_conserved",row.names(Abs_counts)),5:6])
  
  #*********************Count by interactions sign (Anatagonism & Synergy)
  
  Total_probed_ns = unlist(lapply(Interac_distrib,FUN=length))
  hit_rates = sapply(Abs_counts,FUN=sum)/Total_probed_ns
  
  hit_rates_species = as.data.frame(c(sum(sapply(Abs_counts,FUN=sum)[1:2])/sum(Total_probed_ns[1:2]), #EC
                                      sum(sapply(Abs_counts,FUN=sum)[3:4])/sum(Total_probed_ns[3:4]), #ST
                                      sum(sapply(Abs_counts,FUN=sum)[5:6])/sum(Total_probed_ns[5:6]))) #PA
  row.names(hit_rates_species) = Multi_Species
  names(hit_rates_species) = "hit_rates_species"
  
  Total_syn = sum(as.matrix(as.vector(Abs_counts[grep("Syn",row.names(Abs_counts)),])))
  Total_antag = sum(as.matrix(as.vector(Abs_counts[grep("Antag",row.names(Abs_counts)),])))
  Total_antag_int = Total_syn+Total_antag
  
  #plot Total Syn and Antag over possible
  normalized_Syn_fraction = round(Total_syn/sum(DetectSynAntag[[1]]),3)
  normalized_Antag_fraction = round(Total_antag/sum(DetectSynAntag[[2]]),3)
  # barplot(c(normalized_Antag_fraction,normalized_Syn_fraction),
  #         beside=T,col=c(col_antagonism,col_synergy),ylim=c(0,0.16),las=2,
  #         main="Prevalent antagonism: Weak & Strong")
  # text(x=1.1,y=0.13,paste0("syn frac=",normalized_Syn_fraction),adj=0,col=col_synergy)
  # text(x=1.1,y=0.145,paste0("antag frac=",normalized_Antag_fraction),adj=0,col=col_antagonism)
  
  #Output tables
  Abs_counts_summary = as.data.frame(cbind(Total_probed_ns,sapply(Abs_counts,FUN=sum),hit_rates,
                                           sapply(Abs_counts[grep("Syn",row.names(Abs_counts)),],FUN=sum),
                                           sapply(Abs_counts[grep("Antag",row.names(Abs_counts)),],FUN=sum)))
  names(Abs_counts_summary)[c(2,4,5)] = c("Strong&weak_interactions","Synergy","Antagonism")
  Abs_counts_summary = as.data.frame(rbind(Abs_counts_summary,sapply(Abs_counts_summary,FUN=sum)))
  row.names(Abs_counts_summary)[7] = "Totals"
  Abs_counts_summary[7,3] = Abs_counts_summary[7,2]/Abs_counts_summary[7,1]
  
  DetectSynAntag = as.data.frame(cbind(DetectSynAntag,DetectSynAntag/Abs_counts_summary[1:6,1]))
  names(DetectSynAntag)[c(3,4)] = paste0("fraction_",names(DetectSynAntag)[c(3,4)])
  DetectSynAntag = as.data.frame(rbind(DetectSynAntag,sapply(DetectSynAntag,FUN=sum)))
  row.names(DetectSynAntag)[7] = "Totals"
  DetectSynAntag[7,3] = DetectSynAntag[7,1]/Abs_counts_summary[7,1]
  DetectSynAntag[7,4] = DetectSynAntag[7,2]/Abs_counts_summary[7,1]
  
  #***********************************Count by conservation
  Not_conserved = Abs_counts[grep("not_conserved",row.names(Abs_counts)),]
  Not_conserved = as.data.frame(rbind(Not_conserved,Abs_counts[grep("_specific",row.names(Abs_counts)),]))
  Abs_counts_blind = Abs_counts[c(1,3,5,7,9),]+Abs_counts[c(2,4,6,8,10),]
  Abs_counts_blind = as.data.frame(rbind(Abs_counts_blind,
                                         sapply(Abs_counts_blind,FUN=sum)))
  row.names(Abs_counts_blind) = c("conserved","weak_conserved","not_conserved","specific","not_comparable","total")
  
  #nr of strong interactions for at least 1 strain. Fig 3a&b, ED2
  ns_strong1strain = c(length(Spieces_conservation[[1]][,1]),length(Spieces_conservation[[1]][,1]),
                       length(Spieces_conservation[[2]][,1]),length(Spieces_conservation[[2]][,1]),
                       length(Spieces_conservation[[3]][,1]),length(Spieces_conservation[[3]][,1]))
  
  #Fractions of withihn species conserved/non-conserved interactions
  fraction_all_conserved = (Abs_counts_blind[1,]+Abs_counts_blind[2,])/ns_strong1strain
  fraction_all_conserved_corrected = (Abs_counts_blind[1,]+Abs_counts_blind[2,])/(ns_strong1strain-Abs_counts_blind[5,])
  fraction_all_non_conserved = (Abs_counts_blind[3,]+Abs_counts_blind[4,]+Abs_counts_blind[5,])/ns_strong1strain
  fraction_all_non_conserved_corrected = (Abs_counts_blind[3,]+Abs_counts_blind[4,])/(ns_strong1strain-Abs_counts_blind[5,])
  Conservation_fractions_strain = as.data.frame(rbind(fraction_all_conserved,fraction_all_conserved_corrected,
                                                      fraction_all_non_conserved,fraction_all_non_conserved_corrected,
                                                      ns_strong1strain))
  row.names(Conservation_fractions_strain) = c("Conserved/ns_strong1strain","Conserved/ns_strong1strain corrected",
                                               "Non-onserved/ns_strong1strain","Non-onserved/ns_strong1strain corrected",
                                               "ns_strong1strain")
  
  #Get all interactions
  flag=0
  for(sp in 1:length(Multi_Species))
  {
    species = Multi_Species[sp]
    sp_hits = Spieces_conservation[[sp]]
    
    #for b=1
    bug = unlist(Species_strains[[sp]])[1]
    bug1_int = sp_hits[grep(F,is.na(sp_hits[[1]])),]
    bug1_weak = sp_hits[grep_exact("weak_conserved",sp_hits[[3]]),]
    bug1_weak = bug1_weak[grep(T,is.na(bug1_weak[[1]])),]
    bug1_weak[1] = bug1_weak[2]
    bug1_int = as.data.frame(cbind(row.names(bug1_int),bug1_int))
    bug1_weak = as.data.frame(cbind(row.names(bug1_weak),bug1_weak))
    row.names(bug1_int) = c()
    row.names(bug1_weak) = c()
    names(bug1_int)[1] = "combination"
    names(bug1_weak)[1] = "combination"
    bug1_int = as.data.frame(rbind(bug1_int,bug1_weak))
    bug1_int = as.data.frame(cbind(bug1_int[c(1,2,4)],rep(bug,length(bug1_int[[1]]))))
    names(bug1_int)[c(2,4)] = c("direction","bug")
    
    if(flag==0)
    {
      Absolute_interactions = bug1_int
      flag=1
    } else
    {Absolute_interactions = as.data.frame(rbind(Absolute_interactions,bug1_int))}
    
    #for b=2
    bug = unlist(Species_strains[[sp]])[2]
    bug1_int = sp_hits[grep(F,is.na(sp_hits[[2]])),]
    bug1_weak = sp_hits[grep_exact("weak_conserved",sp_hits[[3]]),]
    bug1_weak = bug1_weak[grep(T,is.na(bug1_weak[[2]])),]
    bug1_weak[2] = bug1_weak[1]
    bug1_int = as.data.frame(cbind(row.names(bug1_int),bug1_int))
    bug1_weak = as.data.frame(cbind(row.names(bug1_weak),bug1_weak))
    row.names(bug1_int) = c()
    row.names(bug1_weak) = c()
    names(bug1_int)[1] = "combination"
    names(bug1_weak)[1] = "combination"
    bug1_int = as.data.frame(rbind(bug1_int,bug1_weak))
    bug1_int = as.data.frame(cbind(bug1_int[c(1,3,4)],rep(bug,length(bug1_int[[1]]))))
    names(bug1_int)[c(2,4)] = c("direction","bug")
    
    rm(bug1_weak)
    Absolute_interactions = as.data.frame(rbind(Absolute_interactions,bug1_int))
  }
  rm(flag)
  
  #*************add MAster class info
  drug_ints = as.vector(Absolute_interactions[[1]])
  
  drug1 = unlist(strsplit(drug_ints,split="_"))[seq(1,2*length(drug_ints),by=2)]
  drug2 = unlist(strsplit(drug_ints,split="_"))[seq(2,2*length(drug_ints),by=2)]
  
  class1 = get_feature(drug1,feature = "Master_Class",Attr_table)
  class2 = get_feature(drug2,feature = "Master_Class",Attr_table)
  
  ints_classes = paste0(class1,"_",class2)
  Absolute_interactions = as.data.frame(cbind(Absolute_interactions,ints_classes))
  rm(ints_classes,drug_ints,drug1,drug2,class1,class2)
  
  #*************add class info
  drug_ints = as.vector(Absolute_interactions[[1]])
  
  drug1 = unlist(strsplit(drug_ints,split="_"))[seq(1,2*length(drug_ints),by=2)]
  drug2 = unlist(strsplit(drug_ints,split="_"))[seq(2,2*length(drug_ints),by=2)]
  
  class1 = get_feature(drug1,feature = "Class",Attr_table)
  class2 = get_feature(drug2,feature = "Class",Attr_table)
  
  ints_classes = paste0(class1,"_",class2)
  Absolute_interactions = as.data.frame(cbind(Absolute_interactions,ints_classes))
  rm(ints_classes,drug_ints,drug1,drug2,class1,class2)
  
  names(Absolute_interactions) = c("combination","direction","conserved","bug","Master Class","Class")
  
  #========== Create Results_list =========
  
  nr_drugs_1int_per_bug = c()
  for(b in 1:length(unlist(Species_strains)))
  {
    bug = unlist(Species_strains)[b]
    bug_ints = Absolute_interactions[grep(bug,Absolute_interactions$bug),]
    
    if(length(row.names(bug_ints))>0)
    {
      combs = as.character(bug_ints[[1]])
      drug1 = unlist(strsplit(combs,split="_"))[seq(1,2*length(combs),by=2)]
      drug2 = unlist(strsplit(combs,split="_"))[seq(2,2*length(combs),by=2)]
      drugs12 = c(drug1,drug2)
    } else
    {drugs12 = c()}
    
    nr_drugs_1int_per_bug = c(nr_drugs_1int_per_bug,length(table(drugs12)))
  }
  
  
  Results_list = list(Abs_counts, #Total, Syn and Antag counts
                      Abs_counts_summary,
                      hit_rates_species,
                      DetectSynAntag, #Detectable Syn & Antagonism
                      Abs_counts_blind, #Conservation_within species
                      Conservation_fractions_strain, #Fractions of withihn species conserved/non-conserved interactions
                      nr_drugs_1int_per_bug,
                      Range_classification_counts,
                      Absolute_interactions #Interactions table
  )
  
  names(Results_list) = c("Abs_counts","Abs_counts_summary","hit_rates_species",
                          "Detectable_Syn_&_Antagonism",
                          "Conservation_within_species",
                          "Conservation_fractions_strain",
                          "nr_drugs_1int_per_bug",
                          "Range_classification_counts",
                          "Interactions_table")
  
  #============ Benchmark =========
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
  
  if(length(grep("AS",dummy))>0)
  {Hit_classif[grep("AS",dummy)]="FP"}
  if(length(grep("SA",dummy))>0)
  {Hit_classif[grep("SA",dummy)]="FP"}
  
  Benchmarking = as.data.frame(cbind(Benchmarking,Hit_classif))
  counts = table(Hit_classif)
  
  counts = as.data.frame(t(as.data.frame(t(counts)[1,])))
  Counts = as.data.frame(matrix(0,ncol = 4,nrow = 1))
  row.names(Counts)=NULL
  names(Counts) = c("TP","TN","FP","FN")
  
  Counts[match(names(counts),names(Counts))]=counts
  
  #Precision = round(Counts[4]/(Counts[2]+Counts[4]),2)
  #Recall = round(Counts[4]/(Counts[4]+Counts[1]),2)
  #FPR = round(Counts[2]/(Counts[2]+Counts[3]),2)
  #TPR = Recall 
  
  Precision = round(Counts$TP/(Counts$TP+Counts$FP),2)
  Recall = round(Counts$TP/(Counts$TP+Counts$FN),2)
  FPR = round(Counts$FP/(Counts$FP+Counts$TN),2)
  TPR = Recall
  
  nr_interac = Results_list$Abs_counts_summary
  nr_interac = nr_interac$`Strong&weak_interactions`
  nr_interac = nr_interac[length(nr_interac)]
  
  roc = as.data.frame(t(c(TPR,FPR,Precision,nr_interac)))
  names(roc) = c("TPR","FPR","Precision","total_ints")
  row.names(roc) = interac_cutoff
  
  if(i==1)
  {ROC = roc} else
  {ROC = as.data.frame(rbind(ROC,roc))}
  
  Range_classification_counts = as.data.frame(Range_classification_counts)
  Range_classification_counts = as.data.frame(cbind(t(Range_classification_counts[grep("both",names(Range_classification_counts))]),
                                                    t(Range_classification_counts[grep("complete_specific",names(Range_classification_counts))]),
                                                    t(Range_classification_counts[grep("reduced_specific",names(Range_classification_counts))])))
  row.names(Range_classification_counts) = names(Bugs_hits)
  names(Range_classification_counts) = c("Common","Compelte ExpFit range","Reduced ExpFit range")
  
  if(i==1)
  {Range_counts = list()}
  Range_counts[[length(Range_counts)+1]] = Range_classification_counts
  
} # exp_fit_thresholds = c(0,0.2,0.3,0.4,0.5)

names(Range_counts) = exp_fit_thresholds
row.names(ROC) = exp_fit_thresholds

#====== Plot =====

for(i in 1:length(exp_fit_thresholds))
{
  Range_classification_counts = Range_counts[[i]]
  cutoff_counts = as.data.frame(sapply(Range_classification_counts,FUN=sum))
  names(cutoff_counts) = exp_fit_thresholds[i]
  
  if(i==1)
  {Range_counts_plot = cutoff_counts} else
  {Range_counts_plot = as.data.frame(cbind(Range_counts_plot,cutoff_counts))}
}

pdf(paste0(Out_dir,"ExpFit_win_analysis.pdf"))
barplot(ROC$total_ints,ylim=c(0,3000),names.arg = exp_fit_thresholds,col = "white",
        ylab="#interactions",xlab="ExpectedFitness cutoff",
        main="Expected fitness window analysis")
barplot(as.matrix(Range_counts_plot),ylim=c(0,3000),add=T,axes = F, axisnames = F,legend.text = row.names(Range_counts_plot))
dev.off()

file_id = paste0(Out_dir,ROC_name)
write.table(ROC,file_id,sep="\t",row.names=T,quote=FALSE,col.names = T)





