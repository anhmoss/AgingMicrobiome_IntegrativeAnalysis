

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
reduced_list_narrow = reduced_list[-(grep("UCG|uncultured", reduced_list))] #70 total


testmerge_lists = list()

commonGenus_narrow_list = list()

for(i in 1:11){
  commonGenus_narrow_list[[i]] = stats_fdrs_list[[i]][rownames(stats_fdrs_list[[i]]) %in% reduced_list_narrow,]
}

commonGenus_narrow_list[[1]][(commonGenus_narrow_list[[1]][,3] < 0.05),,drop=F]

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


countsMatched = NULL
cs_taxa_matched_lists = list()

for(i in 1:5){
  for(j in (i+1):6){
  countsMatched[i] = sum(cs_taxa_allList[[i]] %in% cs_taxa_allList[[j]])
  cs_taxa_matched_lists[[j]] =cs_taxa_allList[[i]][cs_taxa_allList[[i]] %in% cs_taxa_allList[[j]]]
  }
}

for(i in 1:5){
  for(k in (i+1):6){
  # for(j in 1:length(cs_taxa_matched_lists[[i]])){
    if(cs_taxa_allList[[i]] %in% cs_taxa_allList[[j]]){
      
    }
  }
}

length(cs_taxa_allList[[5]])
