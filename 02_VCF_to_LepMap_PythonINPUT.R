#Rishi De-Kayne 2017/18 for De-Kayne & Feulner 2018
#vcf to lepmap script - rishi - starts by converting vcf to 'startlines' and 'genotypes' files
#these then need to be inputted to python script to carry out conversion to lepmap format

# FIRST- produce the sex_info_file which is produced manually the format is as follows:
# column1: family number (can be all 1 if you do not have separate families)
# column2: individual name 
# column3: ID of father (can be blank for both parents if you only have F1s)
# column4: ID of mother (can be blank for both parents if you only have F1s)
# column5: sex 0 = unknown, 1 = male and 2 = female
# column6: phenotype, can be 0

#then vcf filtering for one loci/rad locus as in joinmap prep script
require(stringr)
library(plyr)
library(dplyr)
#load in vcf file NOTE THIS FILE IS NAMED 'FileS1' ON THE JOURNAL WEBSITE 
vcf<-read.table("ConcatParentsQ30Thresh20HaploCalled_Total.g5mac3.biallelic.recode.vcf")

####### we don't want any snps from rad loci with multiple (>3 hits)
# make df a with counts occurances of each consensus number (i.e. leaves uniques - 13699)
a <- count(vcf, vars = vcf$V1)
# leave entries where freq is <3 (leaves 12394)
b <- subset(a, (a$n <3))
#now filter vcf based on b (leaves 15175)
filteredVCF <- subset(vcf, vcf$V1 %in% b$vars)
####### of those loci left we only want one per rad read i.e. those with 2 snps per locus need to be filtered (leaves 12394)
uniquevcf <- filteredVCF %>% distinct(filteredVCF$V1, .keep_all = TRUE)
newvcf <- uniquevcf[,1:167]

startlines <- newvcf[,1:9]
marker <- startlines[,1:1]
#create the top line of final file which will become a marker name list
marker <- as.data.frame(marker)
marker <- t(marker)
#and the next lot of columns containing genotypes
genotype <- newvcf[,10:length(newvcf)]
#extract the pattern to remove the many digits following x/y:0 ***
newgenotype <- genotype
for (i in 1:length(genotype)){
  newgenotype[,i] <- as.vector(str_extract(genotype[,i],pattern = "./."))
}
#transform the matrix to fit the lepmap format with individuals as rows and markers as columns
trans_newgenotype <- as.data.frame(t(newgenotype), na.rm = TRUE)

#now to bind the marker matrix with the genotype info
total <- rbind(marker, trans_newgenotype)
total2 <- as.data.frame(t(total))
x1 <- subset(total2, !(total2$V166 == "0/0" & total2$V167 == "1/1"))
x2 <- subset(x1, !(x1$V166 == "0/0" & x1$V167 == "0/0"))
x3 <- subset(x2, !(x2$V166 == "1/1" & x2$V167 == "0/0"))
x4 <- subset(x3, !(x3$V166 == "1/1" & x3$V167 == "1/1"))
x5 <- subset(x4, !(x4$V166 == "./."))
x6 <- subset(x5, !(x5$V167 == "./."))

x6$missing <- rowSums(x6 == "./.")
#missing data allowed: 20% i.e. can be max numb individuals * 0.2 ./. before we reject
rejectloc <- 156*0.2
x6filt <- subset(x6, x6$missing < rejectloc)
x6_2 <- x6filt[,1:((ncol(x6filt))-1)]
total3 <- as.data.frame(t(x6_2))


#lepmap_vcf <- merge(startlines, newgenotype)
#save each of the elements as a csv for pandas
write.csv(startlines, file = "startlines_158TEST.csv")
write.csv(total3, file = "One_perLocus_trans_newgenotype_TEST158_NO_HOMOZYGOUS.csv")
