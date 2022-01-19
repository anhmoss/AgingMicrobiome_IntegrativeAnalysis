##parsing sex metadata for all cohorts except AGP

#check metadata gloor

gloor_a = read.table("/Users/anhil/Desktop/gloor_A.txt", header = TRUE, sep = "\t")
gloor_b = read.table("/Users/anhil/Desktop/gloor_B.txt", header = TRUE, sep = "\t")

check_gloor_ab = merge(gloor_a, gloor_b, by="SampleName")

dim(check_gloor_ab)

all.equal(sort(check_gloor_ab$Run), sort(rawCounts_gloor$SampleID))
all.equal(sort(check_gloor_ab$Age), sort(rawCounts_gloor$Age))

write.table(check_gloor_ab, "gloor_metadata_ageSex.txt", quote = F, row.names = F)

## morgan
morgan_A = read.table("/Users/anhil/Desktop/morgan_A.txt",header = TRUE, sep = "\t")
morgan_B = read.table("/Users/anhil/Desktop/morgan_B.txt",header = TRUE, sep = "\t")

check_morgan_ab = merge(morgan_A, morgan_B, by.x = "sample", by.y="ID")
write.table(check_morgan_ab, "/Users/anhil/Desktop/cohorts_metadata/morgan_ageSex.txt",quote = F, row.names = F)

all.equal(sort(check_morgan_ab$Run), sort(rawCounts_morgan$SampleID))

## zeller germany and france
zeller_all = read.table("/Users/anhil/Desktop/cohorts_metadata/zeller_A.txt", header = TRUE, sep = "\t")
zellerfrance_A = read.table("/Users/anhil/Desktop/cohorts_metadata/zeller_france_16s_meta.tsv", header = TRUE, sep = "\t")
zellergermany_A = read.table("/Users/anhil/Desktop/cohorts_metadata/zeller_germany_16s_metadata.tsv", header = TRUE, sep = "\t")
zeller_zg_B = read.table("/Users/anhil/Desktop/cohorts_metadata/zeller_B.txt", header = TRUE, sep = "\t")

check_zf = merge(zellerfrance_A, zeller_all, by.x="runID", by.y="Run")
check_zg = merge(zellergermany_A, zeller_all, by.x="runID", by.y="Run") ## come back to this...
#maybe pull id's from our current metadata, and merge it with zeller_all...
check_zf$AGE=NULL #removing duplcate column
check_zf$geographic_location_.countryand.orsea = NULL
names(check_zf) = c("Run", "Age", "Country", "Diagnosis", "sex")
write.table(check_zf, "/Users/anhil/Desktop/cohorts_metadata/zellerfrance_ageSex.txt", quote = F, row.names = F)

## for zellergermany metadata
myzg_sampleId = rawCounts_zellergermany$SampleID
my_zg_check = zeller_zg_B[zeller_zg_B$Run %in% myzg_sampleId,]
all.equal(sort(my_zg_check$Run), sort(rawCounts_zellergermany$SampleID))
names(my_zg_check) = c("Run", "Age", "AssayType", "sampleLocation", "sex", "country", "status", "sampleType")

write.table(my_zg_check, "/Users/anhil/Desktop/cohorts_metadata/zellergermany_ageSex.txt", quote = F, row.names = F)


## goodrich
goodrich_A = read.table("/Users/anhil/Desktop/cohorts_metadata/goodrich1_SraRunTable(3).txt",header = TRUE, sep = ",")
goodrich_B = read.table("/Users/anhil/Desktop/cohorts_metadata/goodrich2_SraRunTable(3).txt",header = TRUE, sep = ",")

goodrich_ab = rbind(goodrich_A, goodrich_B)

goodrich_ab_df = data.frame(goodrich_ab$Run, goodrich_ab$individualid, goodrich_ab$AGE, goodrich_ab$sex)
colnames(goodrich_ab)

length(unique(goodrich_ab_df$goodrich_ab.individualid)) #951 unqiue ids vs 1046
length(duplicated(goodrich_ab_df$goodrich_ab.individualid))

myGoodrich_sampleID = rawCounts_goodrich$SampleID
myGoodrich_df = goodrich_ab_df[goodrich_ab_df$goodrich_ab.Run %in% myGoodrich_sampleID,]

#check
all.equal(sort(rawCounts_goodrich$SampleID), sort(myGoodrich_df$goodrich_ab.Run))
write.table(myGoodrich_df,"/Users/anhil/Desktop/cohorts_metadata/goodrich_ageSex.txt" , quote = F, row.names = F)

# agp
agp_fullMetadata = read.table("/Users/anhil/Desktop/cohorts_metadata/agp_SraRunTable.txt", header = TRUE, sep = "\t")
grep("male",agp_fullMetadata)


#baxter
baxter_fullSRATable = read.table("/Users/anhil/Desktop/cohorts_metadata/baxter_SraRunTable.txt", header = TRUE, sep = "\t")
myBaxter_sampleIDs = rawCounts_baxter$SampleID
myBaxter_check_df = baxter_fullSRATable[baxter_fullSRATable$Run %in% myBaxter_sampleIDs,]

all.equal(myBaxter_check_df$Run, myBaxter_sampleIDs)
myBaxter_smaller_df = data.frame(myBaxter_check_df$Run, myBaxter_check_df$AGE,
                                 myBaxter_check_df$gender, myBaxter_check_df$Diagnosis)
colnames(myBaxter_smaller_df) = c("Run", "Age", "Sex", "Diagnosis")
write.table(myBaxter_smaller_df, "/Users/anhil/Desktop/cohorts_metadata/baxter_ageSex.txt", quote = F, row.names = F)
