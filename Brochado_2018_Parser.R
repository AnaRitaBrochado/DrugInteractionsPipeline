
my_list_na_rm <- function(x)
{
  for (i in 1:length(x)) 
  {attr(x[[i]], "na.action") <- NULL}
  
  return(x)
}

#Assign file names
file_id_Bliss = "All_strains_Bliss_scores.txt"
file_id_Pvals = "All_strains_P-vals.txt"
file_id_ExpFit = "All_strains_ExpectedFitness.txt"
file_id_Fit = "All_strains_doubleFitness.txt"
file_id_ArrayFit = "ArraySingleFitness.txt"
file_id_QueryFit = "QuerySingleFitness.txt"

#Rename files
if(!file.exists(paste0(Load_dir,file_id_Bliss)))
{file.rename(paste0(Load_dir,"Supplementary File 1.txt"), paste0(Load_dir,file_id_Bliss))}
if(!file.exists(paste0(Load_dir,file_id_Pvals)))
{file.rename(paste0(Load_dir,"Supplementary File 1.txt"), paste0(Load_dir,file_id_Pvals))}
if(!file.exists(paste0(Load_dir,file_id_ExpFit)))
{file.rename(paste0(Load_dir,"Supplementary File 1.txt"), paste0(Load_dir,file_id_ExpFit))}
if(!file.exists(paste0(Load_dir,file_id_Fit)))
{file.rename(paste0(Load_dir,"Supplementary File 1.txt"), paste0(Load_dir,file_id_Fit))}
if(!file.exists(paste0(Load_dir,file_id_ArrayFit)))
{file.rename(paste0(Load_dir,"Supplementary File 1.txt"), paste0(Load_dir,file_id_ArrayFit))}
if(!file.exists(paste0(Load_dir,file_id_QueryFit)))
{file.rename(paste0(Load_dir,"Supplementary File 1.txt"), paste0(Load_dir,file_id_QueryFit))}
   
#read Bliss scores
Bliss_scores = read.table(paste0(Load_dir,file_id_Bliss),header=T, sep="\t", na.strings="NA", dec=".", strip.white=T)
Interac_distrib = list()
for(b in 1:length(All_Bugs))
{
  bug = All_Bugs[b]
  bug_data = Bliss_scores[grep(bug,Bliss_scores$bug),]
  combs = as.character(bug_data$combination)
  bug_data = as.data.frame(t(bug_data[-c(1,2)]))
  names(bug_data) = combs
  
  bug_data = as.list(bug_data)
  bug_data = lapply(bug_data,FUN=na.omit)
  bug_data = my_list_na_rm(bug_data)
  
  Interac_distrib[[length(Interac_distrib)+1]] = bug_data
}
names(Interac_distrib) = All_Bugs  

#read ExpectedFitness
ExpFit_table = read.table(paste0(Load_dir,file_id_ExpFit),header=T, sep="\t", na.strings="NA", dec=".", strip.white=T)
Exp_fit_dist = list()
for(b in 1:length(All_Bugs))
{
  bug = All_Bugs[b]
  bug_data = ExpFit_table[grep(bug,ExpFit_table$bug),]
  combs = as.character(bug_data$combination)
  bug_data = as.data.frame(t(bug_data[-c(1,2)]))
  names(bug_data) = combs
  
  bug_data = as.list(bug_data)
  bug_data = lapply(bug_data,FUN=na.omit)
  bug_data = my_list_na_rm(bug_data)
  
  Exp_fit_dist[[length(Exp_fit_dist)+1]] = bug_data
}
names(Exp_fit_dist) = All_Bugs

#read P-values
P_val_table = read.table(paste0(Load_dir,file_id_Pvals),header=T, sep="\t", na.strings="NA", dec=".", strip.white=T)
P_values = list()
for(b in 1:length(All_Bugs))
{
  bug = All_Bugs[b]
  bug_data = P_val_table[grep(bug,P_val_table$bug),]
  combs = as.character(bug_data$combination)
  bug_data = as.data.frame(bug_data[-c(1,2)])
  row.names(bug_data) = combs
  
  P_values[[length(P_values)+1]] = bug_data
}
names(P_values) = All_Bugs  

#read Fitness
Fit_table = read.table(paste0(Load_dir,file_id_Fit),header=T, sep="\t", na.strings="NA", dec=".", strip.white=T)
Fitness = list()
well_count = 0
na_well_count = 0
for(b in 1:length(All_Bugs))
{
  bug = All_Bugs[b]
  bug_data = Fit_table[grep(bug,Fit_table$Bug),]
  array_drugs = as.character(as.matrix(bug_data[1,]))
  query_drugs = as.character(bug_data[[2]])
  names(bug_data) = array_drugs
  row.names(bug_data) = query_drugs
  bug_data = bug_data[-c(1,2)]
  bug_data = bug_data[-1,]
  
  fitness = as.data.frame(t(bug_data))
  
  fitness_2 = as.data.frame(matrix(as.numeric(as.matrix(fitness)),ncol=length(fitness)))
  names(fitness_2) = names(fitness)
  row.names(fitness_2) = row.names(fitness)
  fitness = fitness_2
  
  well_count = well_count + length(as.vector(as.matrix(fitness)))
  na_well_count = na_well_count + length(grep(T,is.na(fitness)))
  
  Fitness[[length(Fitness)+1]] = fitness
  
}
names(Fitness) = All_Bugs

#read Query fitness
All_bugs_batch_info = read.table(paste0(Load_dir,file_id_QueryFit),sep="\t",header = T)

#read Array fitness
All_bugs_arrayFit = read.table(paste0(Load_dir,file_id_ArrayFit),sep="\t",header = T)

rm(Bliss_scores,ExpFit_table,P_val_table,Fit_table,combs,my_list_na_rm,bug_data,fitness,fitness_2)