
#=========== Set stage =======

#Need to run "Brochado2017.R" until line 807 before running this script

EC_hits = Spieces_conservation[[1]]
ST_hits = Spieces_conservation[[2]]
PA_hits = Spieces_conservation[[3]]

#Get commonly probed combinations
common_count = rep(0,length(Unique_interacs))
for(b in 1:length(Interac_distrib))
{
  bug_interacs = names(Interac_distrib[[b]])
  common_count[match(bug_interacs,Unique_interacs)] = common_count[match(bug_interacs,Unique_interacs)]+1
}
common_labels = rep(NA,length(common_count))
common_labels[grep(6,common_count)] = "ECSTPA"
common_labels[grep(4,common_count)] = "ECST"
common_labels[grep(2,common_count)] = "PA"

Probed_combs = as.data.frame(common_labels)
row.names(Probed_combs) = Unique_interacs

# ******************** Remove strain specific hits
EC_hits = EC_hits[-grep("specific",EC_hits[[3]]),]
ST_hits = ST_hits[-grep("specific",ST_hits[[3]]),]
PA_hits = PA_hits[-grep("specific",PA_hits[[3]]),]

all_hits = unique(c(row.names(EC_hits),row.names(ST_hits),row.names(PA_hits)))
#sanity check for rev
drug1=unlist(strsplit(all_hits,split="_"))[seq(1,length(unlist(strsplit(all_hits,split="_"))),by=2)]
drug2=unlist(strsplit(all_hits,split="_"))[seq(2,length(unlist(strsplit(all_hits,split="_"))),by=2)]
if(length(all_hits[all_hits %in% paste0(drug2,"_",drug1)]))
{
  for(h in 1:length(all_hits[all_hits %in% paste0(drug2,"_",drug1)]))
  {
    hit = all_hits[all_hits %in% paste0(drug2,"_",drug1)][h]
    rev_hit = paste0(strsplit(hit,split="_")[[1]][2],"_",strsplit(hit,split="_")[[1]][1])
    
    if(length(grep(rev_hit,row.names(EC_hits)))>0)
    {row.names(EC_hits)[grep(rev_hit,row.names(EC_hits))] = hit}
    
    if(length(grep(rev_hit,row.names(ST_hits)))>0)
    {row.names(ST_hits)[grep(rev_hit,row.names(ST_hits))] = hit}
    
    if(length(grep(rev_hit,row.names(PA_hits)))>0)
    {row.names(PA_hits)[grep(rev_hit,row.names(PA_hits))] = hit}
    
  }
}
all_hits = unique(c(row.names(EC_hits),row.names(ST_hits),row.names(PA_hits)))

all_hits = all_hits[all_hits %in% row.names(Probed_combs)[grep("ECSTPA",Probed_combs[[1]])]]

flag = c()
for(h in 1:length(all_hits))
{
  hit = all_hits[h]
  
  hit_EC = EC_hits[grep(hit,row.names(EC_hits)),]
  hit_ST = ST_hits[grep(hit,row.names(ST_hits)),]
  hit_PA = PA_hits[grep(hit,row.names(PA_hits)),]
  
  flag1 = "EC"
  flag2 = "ST"
  flag3 = "PA"
    
  #if(is.null(hit_EC[1,1][[1]]))
  if(length(row.names(hit_EC))==0)
  {flag1 = ""}
  
  #if(is.null(hit_ST[1,1][[1]]))
  if(length(row.names(hit_ST))==0)
  {flag2 = ""}
  
  #if(is.null(hit_PA[1,1][[1]]))
  if(length(row.names(hit_PA))==0)
  {flag3 = ""}
  
  flag = c(flag,paste0(flag1,flag2,flag3))  
}

EC_hits[grep(T,is.na(unlist(EC_hits[[1]]))),1] = unlist(EC_hits[[2]])[grep(T,is.na(unlist(EC_hits[[1]])))]
ST_hits[grep(T,is.na(unlist(ST_hits[[1]]))),1] = unlist(ST_hits[[2]])[grep(T,is.na(unlist(ST_hits[[1]])))]
PA_hits[grep(T,is.na(unlist(PA_hits[[1]]))),1] = unlist(PA_hits[[2]])[grep(T,is.na(unlist(PA_hits[[1]])))]

#=========== Count overlaps and conflicts =======

conflicts = c()
direction = c()
for(h in 1:length(all_hits))
{
  hit = all_hits[h]
  flag1 = flag[h]  
  flag_conf = "no_conflict"
  dir = "non_determined"
  
  if(flag1=="ECSTPA")
  {
    hit_EC = EC_hits[grep(hit,row.names(EC_hits)),1][[1]]
    hit_ST = ST_hits[grep(hit,row.names(ST_hits)),1][[1]]
    hit_PA = PA_hits[grep(hit,row.names(PA_hits)),1][[1]]
    dir = as.character(hit_PA)
    
    if(hit_EC!=hit_ST || hit_EC!=hit_PA || hit_PA!=hit_ST)
    {flag_conf = "conflict"}
  }
  
  if(flag1=="ECST")
  {
    hit_EC = EC_hits[grep(hit,row.names(EC_hits)),1][[1]]
    hit_ST = ST_hits[grep(hit,row.names(ST_hits)),1][[1]]
    dir = as.character(hit_EC)
    
    if(hit_EC!=hit_ST)
    {flag_conf = "conflict"}
  }
  
  if(flag1=="ECPA")
  {
    hit_EC = EC_hits[grep(hit,row.names(EC_hits)),1][[1]]
    hit_PA = PA_hits[grep(hit,row.names(PA_hits)),1][[1]]
    dir = as.character(hit_EC)
    
    if(hit_EC!=hit_PA)
    {flag_conf = "conflict"}
  }
  
  if(flag1=="STPA")
  {
    hit_PA = PA_hits[grep(hit,row.names(PA_hits)),1][[1]]
    hit_ST = ST_hits[grep(hit,row.names(ST_hits)),1][[1]]
    dir = as.character(hit_ST)
    
    if(hit_PA!=hit_ST)
    {flag_conf = "conflict"}
  }
  
  conflicts = c(conflicts,flag_conf) 
  direction = c(direction,dir)
}

Overlaps = as.data.frame(cbind(all_hits,flag,direction,conflicts))
Conflicts = Overlaps[grep_exact("conflict",conflicts),]
Conflicts = Conflicts[-c(3)]
Overlaps = Overlaps[-grep_exact("conflict",conflicts),]

#Add all interactions to Overlaps
for(f in 1:length(unique(as.character(Overlaps[[2]]))))
{
  flag1 = unique(as.character(Overlaps[[2]]))[f]
  
  hits = as.character(Overlaps[grep_exact(flag1,Overlaps[[2]]),1])
  if(flag1=="EC")
  {    
    dirs = as.character(unlist(EC_hits[match(hits,row.names(EC_hits)),1]))
    Overlaps[match(hits,Overlaps[[1]]),3] = dirs
  }
  
  if(flag1=="ST")
  {    
    dirs = as.character(unlist(ST_hits[match(hits,row.names(ST_hits)),1]))
    Overlaps[match(hits,Overlaps[[1]]),3] = dirs
  }
  
  if(flag1=="PA")
  {    
    dirs = as.character(unlist(PA_hits[match(hits,row.names(PA_hits)),1]))
    Overlaps[match(hits,Overlaps[[1]]),3] = dirs
  }
    
}

#Add all interactions to Conflicts
Conflicts = as.data.frame(cbind(Conflicts,rep("EC_hit",length(Conflicts[[1]])),rep("ST_hit",length(Conflicts[[1]])),rep("PA_hit",length(Conflicts[[1]]))))
names(Conflicts)[c(4,5,6)] = c("EC_hit","ST_hit","PA_hit")

EC_conf = EC_hits[match(as.character(Conflicts[[1]]),row.names(EC_hits)),1]
EC_conf[grep(T,lapply(EC_conf,FUN=is.null))] = NA
Conflicts[4] = unlist(EC_conf)

ST_conf = ST_hits[match(as.character(Conflicts[[1]]),row.names(ST_hits)),1]
ST_conf[grep(T,lapply(ST_conf,FUN=is.null))] = NA
Conflicts[5] = unlist(ST_conf)

PA_conf = PA_hits[match(as.character(Conflicts[[1]]),row.names(PA_hits)),1]
PA_conf[grep(T,lapply(PA_conf,FUN=is.null))] = NA
Conflicts[6] = unlist(PA_conf)

# ==================== Get the classes & write files ===========
Overlaps = as.data.frame(cbind(Overlaps,get_class_comb(as.character(Overlaps[[1]]),Attr_table)))
names(Overlaps)[length(Overlaps)] = c("class-class")

Conflicts = as.data.frame(cbind(Conflicts,get_class_comb(as.character(Conflicts[[1]]),Attr_table)))
names(Conflicts)[length(Conflicts)] = c("class-class")

file_id = paste0(Out_dir,"Species Conflicts.txt")
write.table(Conflicts,file_id,sep="\t",row.names=T,quote=FALSE,col.names = T)
file_id = paste0(Out_dir,"Species Overlaps.txt")
write.table(Overlaps,file_id,sep="\t",row.names=T,quote=FALSE,col.names = T)

#=========== Plot Venn diagram =======

require(VennDiagram)


EC_ST_PA = length(grep("ECSTPA",Overlaps[[2]]))
EC_ST = length(grep("ECST",Overlaps[[2]])) #+ EC_ST_PA
ST_PA = length(grep("STPA",Overlaps[[2]])) #+ EC_ST_PA
EC_PA = length(grep("ECPA",Overlaps[[2]])) + EC_ST_PA
EC = length(grep("EC",Overlaps[[2]])) #- EC_ST_PA #+ EC_ST + EC_PA
ST = length(grep("ST",Overlaps[[2]])) #- EC_ST_PA #+ EC_ST + ST_PA
PA = length(grep("PA",Overlaps[[2]])) #- EC_ST_PA #+ EC_PA + ST_PA

#Set dummy par
plot(c(0,0),c(0,0),col="white",frame=F,axes=F,ylab="",xlab="")

draw.triple.venn(EC, ST, PA, EC_ST, ST_PA, EC_PA, EC_ST_PA,
                 category = c("E. coli", "Salmonella", "Pseudomonas"),
                 fill = c("orange","red3","darkgreen"), lwd=c(2.5,2.5,2.5),
                 col=c(c("orange","red3","darkgreen")),alpha=c(0.5,0.5,0.5),ext.percent = rep(0.05,3),margin=c(0.2),
                 cat.dist = rep(0.07,3),euler.d = T,scaled=T)
#clean up
rm(all_hits,conflicts,direction,flag,flag_conf,flag2,flag1,flag3,common_labels,common_count,
   EC_conf,ST_conf,PA_conf,EC_hits,ST_hits,PA_hits,EC,ST,PA,ST_PA,EC_ST,EC_PA,EC_ST_PA,h)

#========= Create file for conserved species network Cytoscape ======

Conserved2sp_interac = Overlaps[grep(T,nchar(as.character(Overlaps$flag))>2),]
Conserved2sp_interac = Conserved2sp_interac[c(1:3)]

new_flag = rep("",length(Conserved2sp_interac[[1]]))
new_flag[grep(T,nchar(as.character(Conserved2sp_interac$flag))==6)]="fully"
new_flag[grep(T,nchar(as.character(Conserved2sp_interac$flag))==4)]="partially"

combs = as.character(Conserved2sp_interac[[1]])
drug1 = unlist(strsplit(combs,split="_"))[seq(1,2*length(combs),by=2)]
drug2 = unlist(strsplit(combs,split="_"))[seq(2,2*length(combs),by=2)]
drug1 = get_feature(drug1,feature = "Abbreviation",Attr_table)
drug2 = get_feature(drug2,feature = "Abbreviation",Attr_table)

Conserved2sp_interac = as.data.frame(cbind(drug1,as.character(Conserved2sp_interac$direction),drug2,new_flag))
names(Conserved2sp_interac)[c(2,4)] = c("interaction","conservation")

file_id = paste0(Out_dir,"SpeciesNetwork_Cytoscape.txt")
write.table(Conserved2sp_interac,file_id,sep="\t",row.names=F,quote=FALSE,col.names = T)

#clean up
rm(combs,drug1,drug2,new_flag)


