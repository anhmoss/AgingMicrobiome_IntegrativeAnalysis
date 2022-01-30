
#Step 1: reads in counts table, labeled by cohort name
  #input: one file input (filepaths and cohort names as columns)
  #output: list of OTU counts table for all 11 datasets

raw_q2021_filespaths = read.table(myFilePath, header = TRUE, sep = "\t")
raw_q2021_files = list()

for(i in 1:nrow(raw_q2021_filespaths)) {
  raw_q2021_files[[i]] = read.table(raw_q2021_filespaths[i,1], header = TRUE, sep = "\t")
}

names(raw_q2021_files) = raw_q2021_filespaths[,2]


#Step 2: lognormalize counts for each dataset, perform stats test (taxa ~ age), permanova test, alpha diversity analysis at the genus level
  #input: counts table for each dataset (ie each element from the list)
  #output: stats calculations (non-parametric and parametric), permanova results per cohort, alpha diversity results and plots (richness and shannon diversity)
  
stats_fdrs_list = list()
lm_results_list = list()
coh_l = vector()


for(i in 1:11){
  raw_file = raw_q2021_files[[i]]
  taxa_table = raw_file[,grep("g__",colnames(raw_file))]
  taxa_table = taxa_table[rowSums(taxa_table)!=0,]
  age_meta = raw_file$Age[rowSums(taxa_table)!=0]
  
  lognorm_file = lognorm_function(taxa_table)
  names_short=sapply(strsplit(colnames(lognorm_file),"g__"),"[[",2)
  colnames(lognorm_file)=names_short
  lognorm_file=lognorm_file[,colSums(lognorm_file)!=0]

  coh_l[i]=raw_file[1,1]
  
  stat_output = stat_correlation_function(counts=lognorm_file, age = age_meta, corrMethod = "kendall")
  stats_fdrs_list[[i]] = stat_output
  
  lm_output = stat_simpleLM_function(counts=lognorm_file, age = age_meta)
  lm_results_list[[i]] = lm_output

}

## Step 2: Check for taxa that are common across all datasets, then check for the significant (FDR 5%) among the common taxa

list1 = rownames(stats_fdrs_list[[1]])
list2 = rownames(stats_fdrs_list[[2]])
list3 = rownames(stats_fdrs_list[[3]])
list4 = rownames(stats_fdrs_list[[4]])
list5 = rownames(stats_fdrs_list[[5]])
list6 = rownames(stats_fdrs_list[[6]])
list7 = rownames(stats_fdrs_list[[7]])
list8 = rownames(stats_fdrs_list[[8]])
list9 = rownames(stats_fdrs_list[[9]])
list10 = rownames(stats_fdrs_list[[10]])
list11 = rownames(stats_fdrs_list[[11]])

all_list = list(list1, list2, list3, list4, list5, list6,
                list7, list8, list9, list10, list11)

reduced_list = Reduce(intersect, all_list) #86 total
reduced_list_narrow = reduced_list[-(grep("UCG|uncultured", reduced_list))] #70 total, exclues genera that are uncultured or unclassified

#pulls stats results for only taxa that are common in all 11 cohorts (reduce_list_narrow)
commonGenus_narrow_list = list()

for(i in 1:11){
  commonGenus_narrow_list[[i]] = stats_fdrs_list[[i]][rownames(stats_fdrs_list[[i]]) %in% reduced_list_narrow,]
}


#among the common taxa, pull the ones that have fdr adjusted pvals  < 0.05 for each cohort
commonAndSignificant_list = list()
names(commonAndSignificant_list) = coh_l

for(i in 1:11){
  commonAndSignificant_list[[i]] = commonGenus_narrow_list[[i]][as.numeric(commonGenus_narrow_list[[i]][,3]) < 0.05,,drop=F]

}

#Bifidobacterium: morgan, gloor, goodrich, agp
#baxter cohort: had 3 that were significant at 5% fdr threhold....
cs_taxa_morgan = rownames(commonAndSignificant_list$Morgan)
cs_taxa_baxter = rownames(commonAndSignificant_list$Baxter)
cs_taxa_goodrich = rownames(commonAndSignificant_list$Goodrich)
cs_taxa_gloor = rownames(commonAndSignificant_list$Gloor)
cs_taxa_agp = rownames(commonAndSignificant_list$AGP)
cs_taxa_nogbcn0 = rownames(commonAndSignificant_list$Nog_BCN0)

cs_taxa_allList = list(cs_taxa_morgan, cs_taxa_baxter, cs_taxa_goodrich, cs_taxa_gloor,
                       cs_taxa_agp, cs_taxa_nogbcn0)
names(cs_taxa_allList) = c("Morgan", "Baxter", "Goodrich", "Gloor", "AGP", "NogBCN0")

#checking if taxa in baxter match w any other ones
cs_taxa_baxter %in% cs_taxa_goodrich # 0 matches
cs_taxa_baxter %in% cs_taxa_gloor # all 3 are in gloor
cs_taxa_baxter %in% cs_taxa_agp #all 3 are in agp
cs_taxa_baxter %in% cs_taxa_nogbcn0 # 0 matches

# " " goodrich
sum(cs_taxa_goodrich %in% cs_taxa_gloor) #6 matches
cs_taxa_goodrich[cs_taxa_goodrich %in% cs_taxa_gloor] #: Parabacteroides, Haemophilus, Haemophilus
cs_taxa_goodrich[cs_taxa_goodrich %in% cs_taxa_agp]  #3 matches 
sum(cs_tax__goodrich %in% cs_taxa_nogbcn0) #1 match : Parabacteroides

#gloor
sum(cs_taxa_gloor %in% cs_taxa_agp) #17 matches
cs_taxa_gloor[cs_taxa_gloor %in% cs_taxa_agp] %in% cs_taxa_nogbcn0]
sum(cs_taxa_gloor %in% cs_taxa_nogbcn0) #9 matches
cs_taxa_gloor[cs_taxa_gloor %in% cs_taxa_agp]
#agp
sum(cs_taxa_agp %in% cs_taxa_nogbcn0) # 3 matches: "Bacteroides" "Howardella"  "Akkermansia"
#cs_taxa_agp[cs_taxa_agp %in% cs_taxa_nogbcn0]


#check
"Bifidobacterium" %in% cs_taxa_morgan #true
"Bifidobacterium" %in% cs_taxa_baxter
"Bifidobacterium" %in% cs_taxa_goodrich #true
"Bifidobacterium" %in% cs_taxa_gloor #true
"Bifidobacterium" %in% cs_taxa_agp #true
"Bifidobacterium" %in% cs_taxa_nogbcn0




