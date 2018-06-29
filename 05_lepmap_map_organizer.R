#Rishi De-Kayne 2018 for De-Kayne & Feuler 2018
#this script allows the analysis of linkage map data resulting from file 04_LepMapCommands.txt
#input files are: 
#     1. the map data here called LOD16_SA_noho_data.csv
#     2. the input file One_perLocus_trans_newgenotype_TEST158_NO_HOMOZYGOUS.csv
#     3. a csv. of the sam file from mapping our parental denovo assembly to the Atlantic Salmon genome: Stampy.07.csv stored on dryad at ***
#     4. the info file with lengths of salmon chromosomes for plotting SalmonChromosomes.csv
#     5. the info file with salmon chromsome arm information inc. lenghts and LORe/AORe assignment: SalmonChromosomeArms.csv
#########################################
#########################################
#     Species average maps
#######################################
#######################################
#first made .txt a .csv in calc to separate cols properly
SAmap <- read.csv(file = "LOD16_SA_noho_data.csv", header = TRUE)
SAmap$LG <- gsub(".SA", "SA", SAmap$LG)
SAmap$LG <- gsub(".txt", "", SAmap$LG)
SAmap$LG <- as.character(SAmap$LG)
for (i in 1:length(SAmap$LG)){
  if ((nchar(SAmap$LG)[i]) < 4){
    SAmap$LG[i] <- gsub("SA", "SA0", SAmap$LG[i])
  }
}
SAmap$LG <- as.factor(SAmap$LG)
#can order to check marker numbers from lepmap do not exceed number of markers in starting file
orderedSAmap<- SAmap[order(SAmap$id),]
#load in the loci/genotype file used before for python script
loci <- read.csv("One_perLocus_trans_newgenotype_TEST158_NO_HOMOZYGOUS.csv")
#transpose so can match rows in id columns
locit <- as.data.frame(t(loci))
#remove marker row
locit <- locit[2:nrow(locit),]
#locit$V1

#do loop to extract the correct id for each given lepmap id
for (i in 1:length(SAmap$id)){
  SAmap$realID[i] <- as.character(locit$V1[as.numeric(SAmap$id[i])]) 
}

#now to get map length
lengths <- tapply(SAmap$cpos, SAmap$LG, max)
lengths
sum(lengths)
t <- table(SAmap$LG)

#########################################
#########################################
#     Species specific maps
#######################################
#######################################
#first made .txt a .csv in calc to separate cols properly
SSmap <- read.csv(file = "LOD16_SS_noho_data.csv", header = TRUE)
SSmap$LG <- as.character(SSmap$LG)
SSmap$LG <- gsub(".SS", "SS", SSmap$LG)
SSmap$LG <- gsub(".txt", "", SSmap$LG)
for (i in 1:length(SAmap$LG)){
  if ((nchar(SSmap$LG)[i]) < 4){
    SSmap$LG[i] <- gsub("SS", "SS0", SSmap$LG[i])
  }
}
SSmap$LG <- as.factor(SSmap$LG)

#can order to check marker numbers from lepmap do not exceed number of markers in starting file
orderedSSmap<- SSmap[order(SSmap$id),]
#load in the loci/genotype file used before for python script
#loci <- read.csv("/home/rishi/Dropbox/Rishi_Linux/lepmap3/PrepInput/One_perLocus_trans_newgenotype_TEST158.csv")
#transpose so can match rows in id columns
locit <- as.data.frame(t(loci))
#remove marker row
locit <- locit[2:nrow(locit),]
#locit$V1

#do loop to extract the correct id for each given lepmap id
for (i in 1:length(SSmap$id)){
  SSmap$realID[i] <- as.character(locit$V1[as.numeric(SSmap$id[i])]) 
}

#now to get map length
Male_lengths <- tapply(SSmap$mpos, SSmap$LG, max)
Male_lengths
sum(Male_lengths)

Fem_lengths <- tapply(SSmap$fpos, SSmap$LG, max)
Fem_lengths
sum(Fem_lengths)


#now want a way of taking the unique positions for SA, female map and male map - remove redundancy
Avemap <- as.data.frame(SAmap$realID)
Avemap$pos <- SAmap$cpos
Avemap$LG <- SAmap$LG
str(Avemap)

AveFilt <- as.data.frame(vector())
for (i in 1:length(levels(Avemap$LG))){
  a <- subset(Avemap, Avemap$LG == as.character((levels(Avemap$LG))[i]))
  b <- as.data.frame(a[!duplicated(a$pos), ])
  AveFilt <- rbind(AveFilt, b)
}

Malemap <- as.data.frame(SSmap$realID)
Malemap$pos <- SSmap$mpos
Malemap$LG <- SSmap$LG
str(Malemap)
MaleFilt <- as.data.frame(vector())
for (i in 1:length(levels(Malemap$LG))){
  a <- subset(Malemap, Malemap$LG == as.character((levels(Malemap$LG))[i]))
  b <- as.data.frame(a[!duplicated(a$pos), ])
  MaleFilt <- rbind(MaleFilt, b)
}


Femmap <- as.data.frame(SSmap$realID)
Femmap$pos <- SSmap$fpos
Femmap$LG <- SSmap$LG
str(Femmap)
FemFilt <- as.data.frame(vector())
for (i in 1:length(levels(Femmap$LG))){
  a <- subset(Femmap, Femmap$LG == as.character((levels(Femmap$LG))[i]))
  b <- as.data.frame(a[!duplicated(a$pos), ])
  FemFilt <- rbind(FemFilt, b)
}

#see if hetero loci present in M/F maps
heteros01 <- subset(locit, locit$V158 == "0/1")
heteros0101 <- subset(heteros01, heteros01$V159 == "0/1")
heterolistSS <- SSmap$realID [SSmap$realID %in% heteros0101$V1]
heterolistSA <- SAmap$realID [SAmap$realID %in% heteros0101$V1]
sames <- SAmap$realID[SAmap$realID %in% SSmap$realID]

write.csv(SAmap, file = "SAmap_forSYNTENY_IN_LOG16_renamed_noho.csv")


Lepmapsummary <- as.data.frame(levels(SAmap$LG))
Lepmapsummary$markernumber <- t
Lepmapsummary$lglengths <- lengths
Lepmapsummary$markersPerCm <- (Lepmapsummary$lglengths/Lepmapsummary$markernumber)
Lepmapsummary$femalelength <- Fem_lengths
Lepmapsummary$fem_markersPerCm <- (Lepmapsummary$femalelength/Lepmapsummary$markernumber)
Lepmapsummary$malelength <- Male_lengths
Lepmapsummary$male_markersPerCm <- (Lepmapsummary$malelength/Lepmapsummary$markernumber)
colnames(Lepmapsummary)[1] <- "LG"
write.csv(Lepmapsummary, file = "MapSummariesforPaper_rename_noho.csv", row.names = FALSE)


#now we want plots for each LG showing both sa and ss markers and colour by which markers are male informative, female informative, heterozygous
#get vectors like heteros0101 but for male/female
#FEMALE
mum01all <- subset(locit, locit$V158 == "0/1")
feminformativeloci <- subset(mum01all, mum01all$V159 == "0/0" | mum01all$V159 == "1/1")
#to make this a vector of those present in the SSmap
feminformativelociSS <- SSmap$realID [SSmap$realID %in% feminformativeloci$V1]
#MALE
dad01all <- subset(locit, locit$V159 == "0/1")
maleinformativeloci <- subset(dad01all, dad01all$V158 == "0/0" | dad01all$V158 == "1/1")
#to make this a vector of those present in the SSmap
maleinformativelociSS <- SSmap$realID [SSmap$realID %in% maleinformativeloci$V1]

#now want to join SS and SA maps so we have all positions
SAhalf <- SAmap
colnames(SAhalf)[3] <- "mpos"
colnames(SAhalf)[4] <- "fpos"
SShalf <- SSmap
SShalf$X <- SShalf$realID
colnames(SShalf)[5] <- "realID"
SShalf <- SShalf[,1:5]

#now add colour column, SA all in forest green, SS heterozygotes in goldenrod, female informative in deeppink2, male informative in darkturquoise
SAhalf$col <- "forestgreen"
for (i in 1:nrow(SShalf)){
  tempy <- subset(SShalf, SShalf$realID == SShalf$realID[i])
  SShalf$col[i] <- as.character("black")
  if (length(samey <- tempy$realID [tempy$realID %in% heteros0101$V1]) > 0){
    SShalf$col[i] <- as.character("goldenrod")}
  if (length(samey <- tempy$realID [tempy$realID %in% feminformativeloci$V1]) > 0){
    SShalf$col[i] <- as.character("deeppink2")}
  if (length(samey <- tempy$realID [tempy$realID %in% maleinformativeloci$V1]) > 0){
    SShalf$col[i] <- as.character("darkturquoise")}
}
FullMarkerset <- rbind(SAhalf, SShalf)
FullMarkerset$LG <- gsub("SA", "SS", FullMarkerset$LG)
FullMarkerset$LG <- as.factor(FullMarkerset$LG)
#nrow(subset(SSmap, SSmap$col == "goldenrod"))
black <- subset(FullMarkerset, FullMarkerset$col == "black")
#tiff("trial.pic",height=30,width=30,units="in",res=300,compression="lzw")
#par(mfrow = c(8, 5))

for (i in 1:40){
  grouping <- subset(FullMarkerset, FullMarkerset$LG == as.character(levels(FullMarkerset$LG)[i]))
  name <- as.character(levels(FullMarkerset$LG)[i])
  filename <- paste(name, ".tiff", sep = "")
  tiff(filename,height=6,width=6,units="in",res=300,compression="lzw")
  plot(grouping$mpos, grouping$fpos, main = name, xlab = "male position (cM)", ylab = "female position (cM)", col = grouping$col, pch = 16, cex = 1)
  legend("topleft", c("species averaged", "F info", "M info", "heterozygous"), fill = c("forestgreen", "deeppink2", "darkturquoise", "goldenrod"), cex = 0.5)
  dev.off()
}

#dev.off()
samfile <- read.csv("Stampy.07.csv")
samfile <- samfile[- grep("@", samfile$X.HD),]
test <- samfile[1:10,]
SyntMap <- SAmap
#first make sure we reduce the set to mapped markers:
MappedLociList <- as.data.frame(samfile$X.HD [samfile$X.HD %in% SyntMap$realID ])
MappedLociSam <- samfile[which( samfile$X.HD %in% MappedLociList$`samfile$X.HD[samfile$X.HD %in% SyntMap$realID]`),]
#now test one row
reduced <- subset(samfile, samfile$X.HD == as.character(SyntMap$realID[1]))

#now for the whole dataset, goes line by line matching the ID with that in the sam file and then outputting salm chrom, pos and qual of mapping
for (i in 1:length(SyntMap$realID)){
  reduced <- subset(samfile, samfile$X.HD == as.character(SyntMap$realID[i]))
  SyntMap$SalmChrom[i] <- as.character(reduced$SO.unsorted)
  SyntMap$SalmPos[i] <- as.character(reduced$X)
  SyntMap$qual[i] <- as.character(reduced$X.1)
}
#now get rid of non ssa lines
SyntMap <- SyntMap[grep("NC", SyntMap$SalmChrom),]
library(stringr)
#leave just SSa part of chromosome and make into Z
SalmChromsdf <- as.data.frame(str_split_fixed(SyntMap$SalmChrom, "_", 3))
SyntMap$SalmChrom <- SalmChromsdf$V3
SyntMap$SalmChrom <- gsub("ssa", "Z", SyntMap$SalmChrom)

#convert bp to -cm for plotting comparison to work:
#sum(lengths) gives length of SA map and since we think Salmon genome is 2966890203 bp then
# map length/bp positions give a cm correction factor for the salmon loci so they are similarly spaced
corFactor <- (sum(lengths)/2966890203)
SyntMap$cmpos <- ((corFactor*as.integer(SyntMap$SalmPos))-(2*(corFactor*as.integer(SyntMap$SalmPos))))

#add yval for plotting later
SyntMap$yval <- 0.5

#subset by quality
qual30 <- subset(SyntMap, as.integer(SyntMap$qual) > 30)

#now to find equivalent whitefish LG for each salmon choromosome  and put in new table consensus df
qual30$SalmChrom <- as.factor(qual30$SalmChrom)
Salmconsensusdf <- as.data.frame(levels(qual30$SalmChrom))
for (i in 1:length(levels(qual30$SalmChrom))){
  subby <- subset(qual30, qual30$SalmChrom == (levels(qual30$SalmChrom)[i]))
  t <- table(subby$LG)
  Salmconsensusdf$equivalent[i] <- names(t)[which.max(t)]
}
#now to find equivalent Salm chromosome for each whiteifish LG  and put in new table consensus df
qual30$LG <- as.factor(qual30$LG)
WFconsensusdf <- as.data.frame(levels(qual30$LG))
for (i in 1:length(levels(qual30$LG))){
  subby <- subset(qual30, qual30$LG == (levels(qual30$LG)[i]))
  t <- table(subby$SalmChrom)
  WFconsensusdf$equivalent[i] <- names(t)[which.max(t)]
}

#now want to make new column linking left and right columns consensus dfs 
WFconsensusdf$`levels(qual30$LG)` <- as.factor(WFconsensusdf$`levels(qual30$LG)`)
for (i in 1:length(qual30$LG)){
  temp <- subset(WFconsensusdf, as.character(WFconsensusdf$`levels(qual30$LG)`) == qual30$LG[i])
  qual30$equivalentSalm[i] <- temp$equivalent
} 

#now plot this info
dev.off()

for (i in 1:40){
  W <- subset(qual30, qual30$LG == (levels(qual30$LG)[i]))
  plot(W$SalmChrom, main = as.character(levels(qual30$LG)[i]))
}
for (i in 1:29){
  W <- subset(qual30, qual30$SalmChrom == (levels(qual30$SalmChrom)[i]))
  plot(W$LG, main = as.character(levels(qual30$SalmChrom)[i]))
}

#####################################################
#####################################################
#recombination study:
#first filter qual30 file so that only loci where the salmon match is the majority are retained:
recomb_df <- subset(qual30, qual30$SalmChrom == qual30$equivalentSalm)
#get positive salm pos again:
recomb_df$pos_salmpos <- ((recomb_df$cmpos)-(2*recomb_df$cmpos))
#then for each locus do wf_cm/salm_cm
recomb_df$recfrac <- recomb_df$cpos/recomb_df$pos_salmpos

dev.off()
for (i in 1:40){
  V <- subset(recomb_df, recomb_df$LG == (levels(recomb_df$LG)[i]))
  title <- paste("WhitefishLG-", as.character((levels(recomb_df$LG)[i])), "  SalmonChrom-", as.character(V$equivalentSalm), "  #markers", as.character(nrow(V)), sep = "")
  plot(V$pos_salmpos, V$recfrac, type = "p", col = "black", pch = 16, cex = 1.5, main = title[1])
}

for (i in 1:40){
  V <- subset(recomb_df, recomb_df$LG == (levels(recomb_df$LG)[i]))
  title <- paste("WhitefishLG-", as.character((levels(recomb_df$LG)[i])), "  SalmonChrom-", as.character(V$equivalentSalm), "  #markers", as.character(nrow(V)), sep = "")
  plot(V$cpos, V$SalmPos, type = "p", col = "black", pch = 16, cex = 1.5, main = title[1])
}

Lepmapsummary$equivalent <- WFconsensusdf$equivalent
#now sort LepmapSummary by equivalent
ordLepmap <- Lepmapsummary[order(Lepmapsummary$equivalent),]
ordLepmap$W <- "W"
#asign new order for synteny mapping
ordLepmap$rank <- c(1:40)
#make this column with a W for circlize
ordLepmap$SyntLG <- paste(ordLepmap$W, ordLepmap$rank, sep = "")
ordLepmap <- subset( ordLepmap, select = -c(W, rank) )
for (row in 1:nrow(ordLepmap)){
  x <- nchar(ordLepmap$SyntLG)[row]
  if (x == "2"){
    ordLepmap$SyntLG[row] <- gsub("W", "W0", ordLepmap$SyntLG[row])} 
}
#we know that we want to swap W01 and W03 to make the synt plot nicer, same with W19 and W20
for (row in 1:nrow(ordLepmap)){
  x <- ordLepmap$SyntLG[row]
  if (x == "W01"){
    ordLepmap$SyntLG[row] <- gsub("W01", "W02", ordLepmap$SyntLG[row])} 
  if (x == "W03"){
    ordLepmap$SyntLG[row] <- gsub("W03", "W01", ordLepmap$SyntLG[row])}
  if (x == "W14"){
    ordLepmap$SyntLG[row] <- gsub("W14", "W15", ordLepmap$SyntLG[row])}
  if (x == "W15"){
    ordLepmap$SyntLG[row] <- gsub("W15", "W14", ordLepmap$SyntLG[row])}
  if (x == "W02"){
    ordLepmap$SyntLG[row] <- gsub("W02", "W03", ordLepmap$SyntLG[row])}
  if (x == "W19"){
    ordLepmap$SyntLG[row] <- gsub("W19", "W20", ordLepmap$SyntLG[row])}
  if (x == "W20"){
    ordLepmap$SyntLG[row] <- gsub("W20", "W19", ordLepmap$SyntLG[row])}
  if (x == "W30"){
    ordLepmap$SyntLG[row] <- gsub("W30", "W31", ordLepmap$SyntLG[row])}
  if (x == "W31"){
    ordLepmap$SyntLG[row] <- gsub("W31", "W30", ordLepmap$SyntLG[row])}
}

#now reorder them
Lepmapsummary <- ordLepmap[order(ordLepmap$LG),]
#and ammend the summary file
write.csv(Lepmapsummary, file = "MapSummariesforPaper_rename_noho_Swap.csv", row.names = FALSE)

#now for synt mapping we want our qual30 dataset to have the equivalent "W" LG naming and the same with SyntMap
for (i in 1:nrow(qual30)){
  equi <- subset(Lepmapsummary, Lepmapsummary$LG == as.character(qual30$LG)[i])
  qual30$SyntLG[i] <- equi$SyntLG
}

for (i in 1:nrow(SyntMap)){
  equi <- subset(Lepmapsummary, Lepmapsummary$LG == as.character(SyntMap$LG)[i])
  SyntMap$SyntLG[i] <- equi$SyntLG
}


###############################################################################################
###############################################################################################
#           Synteny mapping 
###############################################################################################
###############################################################################################
#for synt mapping want specific format:
#need to make dat of same format here: columsn are: locus. pos(in cm), numb, group W01-W40..., and y value = 0.5
#this now fits the 'dat' format in the previous synt mapping script
SyntDat_Map_syntstructure <-(as.data.frame(SyntMap$realID))
SyntDat_Map_syntstructure$locus <- SyntDat_Map_syntstructure$`SyntMap$realID`
SyntDat_Map_syntstructure$pos <- SyntMap$cpos
SyntDat_Map_syntstructure$numb <- as.character(1:length(SyntDat_Map_syntstructure$pos))
SyntDat_Map_syntstructure$group <- SyntMap$SyntLG
SyntDat_Map_syntstructure$yval <- 0.5
SyntDat_Map_syntstructure <- SyntDat_Map_syntstructure[, 2:length(SyntDat_Map_syntstructure)]
#save as 'dat' to test previous script
dat <- SyntDat_Map_syntstructure
#get min and max for each syntLG
dat$group <- as.factor(dat$group)
#get mins
dats_min <- as.data.frame(levels(dat$group))
dats_min$group <- dats_min$`levels(dat$group)`
dats_min$pos <- 0
dats_min <- dats_min[,2:3]
#and maxs
dats_max <- as.data.frame(levels(dat$group))
dats_max$group <- dats_max$`levels(dat$group)`
for (i in 1:length(levels(dat$group))){
  maxs <- subset(SyntDat_Map_syntstructure, SyntDat_Map_syntstructure$group == as.character(levels(dat$group)[i]))
  dats_max$pos[i] <- max(maxs[,2])
}
dats_max <- dats_max[,2:3]
#now merge these dats_min and dats_max
dats <- rbind(dats_min, dats_max)
dats$numb <- c(1:80)
dats$locus <- as.factor(c(1:80))
dats$yval <- 0.5

#actual salmon data for chromosomes - almon chromosomes are the wrong way round
#salmon chromosome data taken from: https://www.ncbi.nlm.nih.gov/nuccore/NC_027307.1
# bp is then converted using a bp-cm conversion based on the map length (in excel)
salm <- read.csv("SalmonChromosomes.csv", sep = ",")
salm
salm$pos <- salm$length_BP*corFactor
salm <- salm[,4:7]
salm$yval <- 0.5
salm2 <- salm
salm2$pos <- (salm2$pos - (2*salm2$pos))
salm2$group <- gsub("SS", "", salm2$group)
salm2$salm <- "S"
salm2$NewGroup <- paste(salm2$salm, salm2$group, sep = "")
salm2$group <- salm2$NewGroup
for (row in 1:nrow(salm2)){
  x <- nchar(salm2$group)[row]
  if (x == "2"){
    salm2$group[row] <- gsub("S", "Z0", salm2$group[row])} else {
      salm2$group[row] <- gsub("S", "Z", salm2$group[row])
    }
}


salm_total <- salm2

#this leaves us with our dat data frame
# here dat$groupchar gives the LG (W01-W42)and dat$pos gives x location and yval gives 0.5
# and in salmon we now have salm_total$group with the Z01-Z29 - REMEMBER these are inveted for plotting!!
# we also have salm_total$pos with the starts and ends of chromosomes
dat <- dats
#and merge with my whitefish dat dataframe after cuttign both to the same size
dat <- dat[,1:5]
salm_total <- salm_total[,1:5]
#make sure dfs match up with type
str(dat)
str(salm_total)
#convert dat $ numb to int from char
dat$numb <- as.integer(dat$numb)
complete <- rbind(dat, salm_total)
complete$group <- as.factor(complete$group)
complete <- complete[order(complete$group),]


chromarmtestreal <- read.csv("SalmonChromosomeArms.csv", header = T)
chromarmtestreal$StartPos <- ((as.integer(chromarmtestreal$StartBP) - (2*as.integer(chromarmtestreal$StartBP)))*corFactor)
chromarmtestreal$EndPos <- ((as.integer(chromarmtestreal$EndBP) - (2*as.integer(chromarmtestreal$EndBP)))*corFactor)


########################################################################################################################
########################################################################################################################
#                           and plot everything - preliminary plot
########################################################################################################################
########################################################################################################################
#now try circos plot, track outlines are given using 'complete' - dots can be added from this - then want to add links
library(circlize)
circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = complete$group, x = complete$pos)
circos.track(factors = complete$group, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
#circos.trackPoints(complete$group, complete$pos, complete$yval, col = "black", pch = 16, cex = 0.5)
#circos.trackPoints(qual30$SyntLG, qual30$cpos, qual30$yval, col = "black", pch = 16, cex = 0.5)
#circos.trackPoints(salm2$group, salm2$pos, salm2$yval, col = "black", pch = 16, cex = 0.5)



vects <- c("red", "orange", "green", "blue", "purple")

vec <- c(add_transparency("red", transparency = 0.8), 
         add_transparency("orange", transparency = 0.8), 
         add_transparency("green", transparency = 0.8), 
         add_transparency("blue", transparency = 0.8), 
         add_transparency("purple", transparency = 0.8))
colvec <- rep(vec, 8)
colvects <- rep(vects, 8)

#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(qual30))){
  a <- qual30$SyntLG[row]
  b <- qual30$cpos[row]
  e <- qual30$SalmChrom[row]
  d <- qual30$cmpos[row]
  fact <- as.factor(qual30$SyntLG)[row]
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.5, col = colvec[fact])
  circos.lines(c(as.integer(b)), 1, as.character(a), col = colvects[fact], type = 'h')
  circos.lines(c(as.integer(d)), 1, as.character(e), col = colvects[fact], type = 'h')
}
#circos.trackPoints(qual30$SalmChrom, qual30$cmpos, qual30$yval, col = colvec[as.factor(qual30$SyntLG)], pch = 16, cex = 0.5)
#circos.trackPoints(qual30$SyntLG, qual30$cpos, qual30$yval, col = colvec[as.factor(qual30$SyntLG)], pch = 16, cex = 0.5)

#now alter salm_total to change naming of salmon chromosomes and then do the same for qual30
complete <- rbind(dat, salm_total)
complete$group <- as.factor(complete$group)
complete <- complete[order(complete$group),]

#first change numbering
for (i in 1:length(salm_total$locus)){
  salm_total$NewSalmnumb[i] <- (gsub("Z", "", salm_total$group[i]))
  salm_total$NewSalmnumb[i] <- as.integer(salm_total$NewSalmnumb)[i]
  salm_total$NewSalmnumb[i] <- as.character(30-as.integer(salm_total$NewSalmnumb)[i])
  salm_total$Z[i] <- "Z"
  salm_total$NewSalm[i] <- paste(salm_total$Z[i], salm_total$NewSalmnumb[i], sep = "")
}

#add the 0s
for (row in 1:nrow(salm_total)){
  x <- nchar(salm_total$NewSalm)[row]
  if (x == "2"){
    salm_total$NewSalm[row] <- gsub("Z", "Z0", salm_total$NewSalm[row])}
}


subs <- subset(salm_total, salm_total$group == qual30$SalmChrom[100])

for (i in 1:length(qual30$realID)){
  subs <- subset(salm_total, salm_total$group == qual30$SalmChrom[i])
  qual30$NewSalm[i] <- subs$NewSalm[1]
}

newsalm <- salm_total
newsalm$group <- newsalm$NewSalm
newsalm <- newsalm[,1:5]
newsalm <- newsalm[order(newsalm$group),]
complete <- rbind(dat, newsalm)
complete$group <- as.factor(complete$group)
complete <- complete[order(complete$group),]

########################################################################################################################
########################################################################################################################
#now try circos plot, track outlines are given using 'complete' - dots can be added from this - then want to add links
library(circlize)
circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = complete$group, x = complete$pos)
circos.track(factors = complete$group, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })

vects <- c("red", "orange", "green", "blue", "purple")

vec <- c(add_transparency("red", transparency = 0.8), 
         add_transparency("orange", transparency = 0.8), 
         add_transparency("green", transparency = 0.8), 
         add_transparency("blue", transparency = 0.8), 
         add_transparency("purple", transparency = 0.8))
colvec <- rep(vec, 8)
colvects <- rep(vects, 8)

#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(qual30))){
  a <- qual30$SyntLG[row]
  b <- qual30$cpos[row]
  e <- qual30$NewSalm[row]
  d <- qual30$cmpos[row]
  fact <- as.factor(qual30$SyntLG)[row]
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = colvec[fact])
  circos.lines(c(as.integer(b)), 1, as.character(a), col = colvects[fact], type = 'h')
  circos.lines(c(as.integer(d)), 1, as.character(e), col = colvects[fact], type = 'h')
}


########################################################################################################################
########################################################################################################################
#                           and plot everything with a gap
########################################################################################################################
########################################################################################################################
library(circlize)
circos.clear()
#tiff("synteny_figure_reorderedWFLGs.tiff",height=8,width=8,units="in",res=300,compression="lzw")
# Customize plot
par(fig=c(0,1,0,1),mar=c(1,1,1,1))

#make a vector of gaps
gapsnorm <- c(1)
gapswide <- c(1)
wfgaps <- rep(gapsnorm, 39)
salmgaps <- rep(gapswide, 28)
startgap <- c(9)
endgap <- c(9)
gapsall <- c(wfgaps, startgap, salmgaps, endgap)

Wnamevector <- c("W01", "W02", "W03", "W04", "W05", "W06", "W07", "W08", "W09", "W10", 
                 "W11", "W12", "W13", "W14", "W15", "W16", "W17", "W18", "W19", "W20", 
                 "W21", "W22", "W23", "W24", "W25", "W26", "W27", "W28", "W29", "W30", 
                 "W31", "W32", "W33", "W34", "W35", "W36", "W37", "W38", "W39", "W40")
Snamevector <- c("S29", "S28", "S27", "S26", "S25", "S24", "S23", "S22", "S21", "S20", 
                 "S19", "S18", "S17", "S16", "S15", "S14", "S13", "S12", "S11", "S10", 
                 "S09", "S08", "S07", "S06", "S05", "S04", "S03", "S02", "S01")

#create vector of segment names
namevector <- c(Wnamevector, Snamevector)
#then set circos plot parameters
circos.par("track.height" = 0.08 , start.degree=95, cell.padding = c(0.005, 0, 0.005, 0), gap.degree = gapsall)
circos.initialize(factors = complete$group, x = complete$pos)

circos.track(factors = complete$group, ylim = c(0,0.1), track.height = uh(1, "mm"),
             panel.fun = function(x, y) {
             }, bg.border = NA)
circos.track(factors = complete$group, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 2.2, 
                           namevector[CELL_META$sector.numeric.index], cex = 0.5)
               circos.axis(labels.cex = 0.5)
             })


# for (i in 1:69){
#   namey <- as.character(levels(complete$group)[i])
#   circos.update(sector.index = as.character(namey), track.index = 1, 
#                 bg.col = "white", bg.border = "white")
# }

for (i in 1:(nrow(chromarmtestreal))){
  ab <- as.integer(chromarmtestreal$StartPos[i])
  ba <- as.integer(chromarmtestreal$EndPos[i])
  circos.rect(ab, 0, ba, 0.1, as.character(chromarmtestreal$SALMName[i]), track.index = 1, col = as.character(chromarmtestreal$Colo[i]))
}
for (i in 1:(nrow(chromarmtestreal))){
  ab <- as.integer(chromarmtestreal$StartPos[i])
  ba <- as.integer(chromarmtestreal$EndPos[i])
  circos.text(((ab+ba)/2), 0.2, as.character(chromarmtestreal$SalmArm)[i], as.character(chromarmtestreal$SALMName[i]), track.index = 1, cex = 0.5)
}


vects <- c("red", "orange", "green", "blue", "purple")

vec <- c(add_transparency("red", transparency = 0.8), 
         add_transparency("orange", transparency = 0.8), 
         add_transparency("green", transparency = 0.8), 
         add_transparency("blue", transparency = 0.8), 
         add_transparency("purple", transparency = 0.8))
colvec <- rep(vec, 8)
colvects <- rep(vects, 8)

#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(qual30))){
  a <- qual30$SyntLG[row]
  b <- qual30$cpos[row]
  e <- qual30$NewSalm[row]
  d <- qual30$cmpos[row]
  fact <- as.factor(qual30$SyntLG)[row]
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = colvec[fact])
  circos.lines(c(as.integer(b)), 1, as.character(a), col = colvects[fact], type = 'h')
  circos.lines(c(as.integer(d)), 1, as.character(e), col = colvects[fact], type = 'h')
}

#dev.off()


#circos.trackPoints(complete$group, complete$pos, complete$yval, col = "black", pch = 16, cex = 0.5)
circos.trackPoints(dat$group, dat$pos, dat$yval, col = "black", pch = 20, cex = 0.5)
# circos.trackPoints(newsalmtotal$group, newsalmtotal$pos, newsalmtotal$yval, col = "black", pch = 16, cex = 0.5)  


#loop through linkage data frame and extract info to make links
for (row in 1:nrow(link_total30)){
  a <- link_total30$LGchar[row]
  b <- link_total30$pos[row]
  e <- link_total30$salmChrom[row]
  d <- link_total30$cmpos[row]
  fact <- as.factor(link_total30$LGchar)[row]
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "darkgrey")
  circos.trackPoints(as.character(e), as.integer(d), 0.5, col = "black", pch = 20, cex = 0.5)
  
}

for (row in 1:nrow(link_total30)){
  if (link_total30$LGchar[row] == "W13"){
    if(link_total30$salmChrom[row] == "Z13"){
      a <- link_total30$LGchar[row]
      b <- link_total30$pos[row]
      e <- link_total30$salmChrom[row]
      d <- link_total30$cmpos[row]
      fact <- as.factor(link_total30$LGchar)[row]
      circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d),(as.integer(d) + (100*0.00000175))), h=0.9, col = "red")
      circos.trackPoints(as.character(e), as.integer(d), 0.5, col = "red", pch = 20, cex = 0.5)
      circos.trackPoints(as.character(a), as.integer(b), 0.5, col = "red", pch = 20, cex = 0.5)
    }
  }
}

for (row in 1:nrow(link_total30)){
  if (link_total30$LGchar[row] == "W03"){
    a <- link_total30$LGchar[row]
    b <- link_total30$pos[row]
    e <- link_total30$salmChrom[row]
    d <- link_total30$cmpos[row]
    fact <- as.factor(link_total30$LGchar)[row]
    circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d),(as.integer(d) + (100*0.00000175))), h=0.9, col = "green")
    circos.trackPoints(as.character(e), as.integer(d), 0.5, col = "green", pch = 20, cex = 0.5)
    circos.trackPoints(as.character(a), as.integer(b), 0.5, col = "green", pch = 20, cex = 0.5)
  }
}

for (row in 1:nrow(link_total30)){
  if (link_total30$LGchar[row] == "W19"){
    a <- link_total30$LGchar[row]
    b <- link_total30$pos[row]
    e <- link_total30$salmChrom[row]
    d <- link_total30$cmpos[row]
    fact <- as.factor(link_total30$LGchar)[row]
    circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d),(as.integer(d) + (100*0.00000175))), h=0.9, col = "purple")
    circos.trackPoints(as.character(e), as.integer(d), 0.5, col = "purple", pch = 20, cex = 0.5)
    circos.trackPoints(as.character(a), as.integer(b), 0.5, col = "purple", pch = 20, cex = 0.5)
  }
}

for (row in 1:nrow(link_total30)){
  if (link_total30$LGchar[row] == "W09"){
    a <- link_total30$LGchar[row]
    b <- link_total30$pos[row]
    e <- link_total30$salmChrom[row]
    d <- link_total30$cmpos[row]
    fact <- as.factor(link_total30$LGchar)[row]
    circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d),(as.integer(d) + (100*0.00000175))), h=0.9, col = "blue")
    circos.trackPoints(as.character(e), as.integer(d), 0.5, col = "blue", pch = 20, cex = 0.5)      
    circos.trackPoints(as.character(a), as.integer(b), 0.5, col = "blue", pch = 20, cex = 0.5)
  }
}
#circos.trackPoints(link_total30$salmChrom, link_total30$cmpos, link_total30$yval, col = "black", pch = 16, cex = 0.5)
#circos.trackPoints(link_total30$LGchar, link_total30$pos, link_total30$yval, col = "black", pch = 16, cex = 0.5)
#dev.off()  

################################################################################################################################################################
################################################################################################################################################################
#
#                                               FIGURE 2 FOR PUBLICATION:
################################################################################################################################################################
################################################################################################################################################################
#no axies

library(circlize)
circos.clear()
#name of file/parameters
tiff("Figure2.tiff",height=8,width=8,units="in",res=300,compression="lzw")
# Customize plot parameters
par(fig=c(0,1,0,1),mar=c(1,1,1,1))
#create gaps vector
gapsnorm <- c(1)
gapswide <- c(1)
wfgaps <- rep(gapsnorm, 39)
salmgaps <- rep(gapswide, 28)
startgap <- c(9)
endgap <- c(9)
gapsall <- c(wfgaps, startgap, salmgaps, endgap)

#create names vector
Wnamevector <- c("W01", "W02", "W03", "W04", "W05", "W06", "W07", "W08", "W09", "W10", 
                 "W11", "W12", "W13", "W14", "W15", "W16", "W17", "W18", "W19", "W20", 
                 "W21", "W22", "W23", "W24", "W25", "W26", "W27", "W28", "W29", "W30", 
                 "W31", "W32", "W33", "W34", "W35", "W36", "W37", "W38", "W39", "W40")
Snamevector <- c("Ssa29", "Ssa28", "Ssa27", "Ssa26", "Ssa25", "Ssa24", "Ssa23", "Ssa22", "Ssa21", "Ssa20", 
                 "Ssa19", "Ssa18", "Ssa17", "Ssa16", "Ssa15", "Ssa14", "Ssa13", "Ssa12", "Ssa11", "Ssa10", 
                 "Ssa09", "Ssa08", "Ssa07", "Ssa06", "Ssa05", "Ssa04", "Ssa03", "Ssa02", "Ssa01")

namevector <- c(Wnamevector, Snamevector)

#then set circos plot parameters including the rotation, gaps between segments and gaps between tracks
circos.par("track.height" = 0.08 , start.degree=95, cell.padding = c(0.005, 0, 0.005, 0), gap.degree = gapsall, track.margin = c(0.0045, 0.0045))
#initialize
circos.initialize(factors = complete$group, x = complete$pos)

#initialize first track - chromosomes/arms
circos.track(factors = complete$group, ylim = c(0,0.1), track.height = uh(1.5, "mm"),
             panel.fun = function(x, y) {
             }, bg.border = NA)
#initialize second track colour bands
circos.track(factors = complete$group, ylim = c(0,0.1), track.height = uh(1, "mm"),
             panel.fun = function(x, y) {
             })

#recolour all bands with custom colours
circos.update(sector.index = "W01", track.index = 2, 
              bg.col = "firebrick", bg.border = "firebrick")
circos.update(sector.index = "W02", track.index = 2, 
              bg.col = "firebrick", bg.border = "firebrick")
circos.update(sector.index = "W03", track.index = 2, 
              bg.col = "firebrick", bg.border = "firebrick")
circos.update(sector.index = "Z29", track.index = 2, 
              bg.col = "firebrick", bg.border = "firebrick")
circos.update(sector.index = "W04", track.index = 2, 
              bg.col = "deeppink", bg.border = "deeppink")
circos.update(sector.index = "W05", track.index = 2, 
              bg.col = "deeppink", bg.border = "deeppink")
circos.update(sector.index = "Z28", track.index = 2, 
              bg.col = "black", bg.border = "black")
circos.update(sector.index = "Z27", track.index = 2, 
              bg.col = "deeppink", bg.border = "deeppink")
circos.update(sector.index = "W06", track.index = 2, 
              bg.col = "purple", bg.border = "purple")
circos.update(sector.index = "W07", track.index = 2, 
              bg.col = "purple", bg.border = "purple")
circos.update(sector.index = "Z26", track.index = 2, 
              bg.col = "purple", bg.border = "purple")
circos.update(sector.index = "W08", track.index = 2, 
              bg.col = "darkblue", bg.border = "darkblue")
circos.update(sector.index = "Z25", track.index = 2, 
              bg.col = "darkblue", bg.border = "darkblue")
circos.update(sector.index = "W09", track.index = 2, 
              bg.col = "royalblue", bg.border = "royalblue")
circos.update(sector.index = "Z24", track.index = 2, 
              bg.col = "royalblue", bg.border = "royalblue")
circos.update(sector.index = "W10", track.index = 2, 
              bg.col = "darkturquoise", bg.border = "darkturquoise")
circos.update(sector.index = "Z23", track.index = 2, 
              bg.col = "darkturquoise", bg.border = "darkturquoise")
circos.update(sector.index = "Z22", track.index = 2, 
              bg.col = "white", bg.border = "white")
circos.update(sector.index = "W11", track.index = 2, 
              bg.col = "green", bg.border = "green")
circos.update(sector.index = "W12", track.index = 2, 
              bg.col = "green", bg.border = "green")
circos.update(sector.index = "W13", track.index = 2, 
              bg.col = "green", bg.border = "green")
circos.update(sector.index = "Z21", track.index = 2, 
              bg.col = "green", bg.border = "green")
circos.update(sector.index = "W14", track.index = 2, 
              bg.col = "forestgreen", bg.border = "forestgreen")
circos.update(sector.index = "W15", track.index = 2, 
              bg.col = "forestgreen", bg.border = "forestgreen")
circos.update(sector.index = "Z20", track.index = 2, 
              bg.col = "forestgreen", bg.border = "forestgreen")
circos.update(sector.index = "W16", track.index = 2, 
              bg.col = "goldenrod", bg.border = "goldenrod")
circos.update(sector.index = "W17", track.index = 2, 
              bg.col = "goldenrod", bg.border = "goldenrod")
circos.update(sector.index = "Z19", track.index = 2, 
              bg.col = "goldenrod", bg.border = "goldenrod")
circos.update(sector.index = "W18", track.index = 2, 
              bg.col = "brown", bg.border = "brown")
circos.update(sector.index = "Z18", track.index = 2, 
              bg.col = "brown", bg.border = "brown")
circos.update(sector.index = "W19", track.index = 2, 
              bg.col = "darkorange3", bg.border = "darkorange3")
circos.update(sector.index = "W20", track.index = 2, 
              bg.col = "darkorange3", bg.border = "darkorange3")
circos.update(sector.index = "Z17", track.index = 2, 
              bg.col = "darkorange3", bg.border = "darkorange3")
circos.update(sector.index = "W21", track.index = 2, 
              bg.col = "orangered", bg.border = "orangered")
circos.update(sector.index = "Z16", track.index = 2, 
              bg.col = "orangered", bg.border = "orangered")
circos.update(sector.index = "W22", track.index = 2, 
              bg.col = "tan1", bg.border = "tan1")
circos.update(sector.index = "W23", track.index = 2, 
              bg.col = "tan1", bg.border = "tan1")
circos.update(sector.index = "Z15", track.index = 2, 
              bg.col = "tan1", bg.border = "tan1")
circos.update(sector.index = "W24", track.index = 2, 
              bg.col = "gold", bg.border = "gold")
circos.update(sector.index = "W25", track.index = 2, 
              bg.col = "gold", bg.border = "gold")
circos.update(sector.index = "Z14", track.index = 2, 
              bg.col = "gold", bg.border = "gold")
circos.update(sector.index = "W26", track.index = 2, 
              bg.col = "yellow", bg.border = "yellow")
circos.update(sector.index = "Z13", track.index = 2, 
              bg.col = "yellow", bg.border = "yellow")
circos.update(sector.index = "W27", track.index = 2, 
              bg.col = "yellowgreen", bg.border = "yellowgreen")
circos.update(sector.index = "W28", track.index = 2, 
              bg.col = "yellowgreen", bg.border = "yellowgreen")
circos.update(sector.index = "Z12", track.index = 2, 
              bg.col = "yellowgreen", bg.border = "yellowgreen")
circos.update(sector.index = "W29", track.index = 2, 
              bg.col = "peachpuff", bg.border = "peachpuff")
circos.update(sector.index = "Z11", track.index = 2, 
              bg.col = "peachpuff", bg.border = "peachpuff")
circos.update(sector.index = "W30", track.index = 2, 
              bg.col = "pink", bg.border = "pink")
circos.update(sector.index = "W31", track.index = 2, 
              bg.col = "pink", bg.border = "pink")
circos.update(sector.index = "Z10", track.index = 2, 
              bg.col = "pink", bg.border = "pink")
circos.update(sector.index = "W32", track.index = 2, 
              bg.col = "palevioletred4", bg.border = "palevioletred4")
circos.update(sector.index = "Z09", track.index = 2, 
              bg.col = "palevioletred4", bg.border = "palevioletred4")
circos.update(sector.index = "W33", track.index = 2, 
              bg.col = "thistle1", bg.border = "thistle1")
circos.update(sector.index = "Z08", track.index = 2, 
              bg.col = "thistle1", bg.border = "thistle1")
circos.update(sector.index = "W34", track.index = 2, 
              bg.col = "violet", bg.border = "violet")
circos.update(sector.index = "Z07", track.index = 2, 
              bg.col = "violet", bg.border = "violet")
circos.update(sector.index = "W35", track.index = 2, 
              bg.col = "mediumorchid", bg.border = "mediumorchid")
circos.update(sector.index = "Z06", track.index = 2, 
              bg.col = "mediumorchid", bg.border = "mediumorchid")
circos.update(sector.index = "W36", track.index = 2, 
              bg.col = "midnightblue", bg.border = "midnightblue")
circos.update(sector.index = "Z05", track.index = 2, 
              bg.col = "midnightblue", bg.border = "midnightblue")
circos.update(sector.index = "W37", track.index = 2, 
              bg.col = "lightslateblue", bg.border = "lightslateblue")
circos.update(sector.index = "Z03", track.index = 2, 
              bg.col = "lightslateblue", bg.border = "lightslateblue")
circos.update(sector.index = "Z04", track.index = 2, 
              bg.col = "black", bg.border = "black")
circos.update(sector.index = "W38", track.index = 2, 
              bg.col = "lightsteelblue", bg.border = "lightsteelblue")
circos.update(sector.index = "W39", track.index = 2, 
              bg.col = "lightsteelblue", bg.border = "lightsteelblue")
circos.update(sector.index = "Z02", track.index = 2, 
              bg.col = "lightsteelblue", bg.border = "lightsteelblue")
circos.update(sector.index = "W40", track.index = 2, 
              bg.col = "seashell3", bg.border = "seashell3")
circos.update(sector.index = "Z01", track.index = 2, 
              bg.col = "seashell3", bg.border = "seashell3")
# circos.track(factors = complete$group, ylim = c(0,1),
#              panel.fun = function(x, y) {
#                circos.text(CELL_META$xcenter, 1.4 + uy(5, "mm"), 
#                            namevector[CELL_META$sector.numeric.index], cex = 0.35)
#              })
#create a letter to search for
Wlet <- "W"

#initialize final track with the track labels based on the presence of W in the name of the segment i.e whitefish labels in closer, Salmon out further
circos.track(factors = complete$group, ylim = c(0,1),
             panel.fun = function(x, y) {
               if(grepl(Wlet, namevector[CELL_META$sector.numeric.index])){
                 circos.text(CELL_META$xcenter, 0.9 + uy(5, "mm"), 
                             namevector[CELL_META$sector.numeric.index], cex = 0.37)
               } else {
                 circos.text(CELL_META$xcenter, 1.45 + uy(5, "mm"), 
                             namevector[CELL_META$sector.numeric.index], cex = 0.37)
               }
             })


#add chromosome arm data
for (i in 1:(nrow(chromarmtestreal))){
  ab <- as.integer(chromarmtestreal$StartPos[i])
  ba <- as.integer(chromarmtestreal$EndPos[i])
  circos.rect(ab, 0, ba, 0.1, as.character(chromarmtestreal$SALMName[i]), track.index = 1, col = as.character(chromarmtestreal$Colo[i]))
}
#add chromosome arm labels
for (i in 1:(nrow(chromarmtestreal))){
  ab <- as.integer(chromarmtestreal$StartPos[i])
  ba <- as.integer(chromarmtestreal$EndPos[i])
  circos.text(((ab+ba)/2), 0.3, as.character(chromarmtestreal$Letter)[i], as.character(chromarmtestreal$SALMName[i]), track.index = 1, cex = 0.3)
}

#make links for LORe regions (make sure you run this after the end of the script)
#include to produce publication figure
# Lores <- subset(chromarmtestreal, chromarmtestreal$ORE == "L")
# from1 <- subset(Lores, Lores$SalmArm == "2p")
# to1 <- subset(Lores, Lores$SalmArm == "5q")
# circos.link(from1$SALMName, c(from1$StartPos, from1$EndPos), to1$SALMName, c(to1$StartPos, to1$EndPos), h=1, col = add_transparency("white", transparency = 0.75), border = "black")
# from1 <- subset(Lores, Lores$SalmArm == "2q")
# to1 <- subset(Lores, Lores$SalmArm == "12qa")
# circos.link(from1$SALMName, c(from1$StartPos, from1$EndPos), to1$SALMName, c(to1$StartPos, to1$EndPos), h=1, col = add_transparency("white", transparency = 0.75), border = "black")
# from1 <- subset(Lores, Lores$SalmArm == "3q")
# to1 <- subset(Lores, Lores$SalmArm == "6p")
# circos.link(from1$SALMName, c(from1$StartPos, from1$EndPos), to1$SALMName, c(to1$StartPos, to1$EndPos), h=1, col = add_transparency("white", transparency = 0.75), border = "black")
# from1 <- subset(Lores, Lores$SalmArm == "4p")
# to1 <- subset(Lores, Lores$SalmArm == "8q")
# circos.link(from1$SALMName, c(from1$StartPos, from1$EndPos), to1$SALMName, c(to1$StartPos, to1$EndPos), h=1, col = add_transparency("white", transparency = 0.75), border = "black")
# from1 <- subset(Lores, Lores$SalmArm == "7q")
# to1 <- subset(Lores, Lores$SalmArm == "17qb")
# circos.link(from1$SALMName, c(from1$StartPos, from1$EndPos), to1$SALMName, c(to1$StartPos, to1$EndPos), h=1, col = add_transparency("white", transparency = 0.75), border = "black")
# from1 <- subset(Lores, Lores$SalmArm == "11qa")
# to1 <- subset(Lores, Lores$SalmArm == "26")
# circos.link(from1$SALMName, c(from1$StartPos, from1$EndPos), to1$SALMName, c(to1$StartPos, to1$EndPos), h=1, col = add_transparency("white", transparency = 0.75), border = "black")
# from1 <- subset(Lores, Lores$SalmArm == "16qb")
# to1 <- subset(Lores, Lores$SalmArm == "17qa")
# circos.link(from1$SALMName, c(from1$StartPos, from1$EndPos), to1$SALMName, c(to1$StartPos, to1$EndPos), h=1, col = add_transparency("white", transparency = 0.75), border = "black")

#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(qual30))){
  a <- qual30$SyntLG[row]
  b <- qual30$cpos[row]
  e <- qual30$NewSalm[row]
  d <- qual30$cmpos[row]
  fact <- as.factor(qual30$SyntLG)[row]
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("darkgrey", transparency = 0.35))
  circos.lines(c(as.integer(b)), 1, as.character(a), col = "darkgrey", type = 'h')
  circos.lines(c(as.integer(d)), 1, as.character(e), col = "darkgrey", type = 'h')
}
# 
# #do links for W1,2,3
# for (row in 1:(nrow(qual30))){
#   a <- qual30$SyntLG[row]
#   b <- qual30$cpos[row]
#   e <- qual30$NewSalm[row]
#   d <- qual30$cmpos[row]
#   fact <- as.factor(qual30$SyntLG)[row]
#   if (as.character(a) == "W01"){
#     if (as.character(e) == "Z29"){
#       circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("red", transparency = 0.65))
#     } else {
#       circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("blue", transparency = 0.65))
#     }
#   } else if (as.character(a) == "W02"){
#     if (as.character(e) == "Z29"){
#       circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("red", transparency = 0.65))
#     } else {
#       circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("blue", transparency = 0.65))
#     }
#   } else if (as.character(a) == "W03"){
#     if (as.character(e) == "Z29"){
#       circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("red", transparency = 0.65))
#     } else {
#       circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("blue", transparency = 0.65))
#     }
#   }
# }

# #add links from qual 30 dataset
# for (row in 1:(nrow(qual30))){
#   a <- qual30$SyntLG[row]
#   b <- qual30$cpos[row]
#   e <- qual30$NewSalm[row]
#   d <- qual30$cmpos[row]
#   samey1 <- as.character(qual30$SalmChrom)[row]
#   samey2 <- as.character(qual30$equivalentSalm)[row]
#   fact <- as.factor(qual30$SyntLG)[row]
#   if (samey1 == samey2){
#     circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("darkgrey", transparency = 0.65)) 
#   } else {
#     circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("black", transparency = 0.65))
#   }
# }



#now add the lines within track three to identify same-same links or same-different
for (row in 1:(nrow(qual30))){
  a <- qual30$SyntLG[row]
  b <- qual30$cpos[row]
  e <- qual30$NewSalm[row]
  d <- qual30$cmpos[row]
  samey1 <- as.character(qual30$SalmChrom)[row]
  samey2 <- as.character(qual30$equivalentSalm)[row]
  fact <- as.factor(qual30$SyntLG)[row]
  if (samey1 == samey2){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "darkgrey", lwd = 0.5, type = 'h') 
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "darkgrey", lwd = 0.5, type = 'h')
    
  } else {
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "black", lwd = 0.5, type = 'h') 
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "black", lwd = 0.5, type = 'h')
  }
}

dev.off()

# ################################################################################
# ################################################################################
# #   RECOLOUR - this was a new plot with different colours
# ################################################################################
# ################################################################################
# library(circlize)
# circos.clear()
# tiff("synteny_figure_reorderedWFLGs_greyRecol2.tiff",height=8,width=8,units="in",res=300,compression="lzw")
# # Customize plot
# par(fig=c(0,1,0,1),mar=c(1,1,1,1))
# gapsnorm <- c(1)
# gapswide <- c(1)
# wfgaps <- rep(gapsnorm, 39)
# salmgaps <- rep(gapswide, 28)
# startgap <- c(9)
# endgap <- c(9)
# gapsall <- c(wfgaps, startgap, salmgaps, endgap)
# 
# Wnamevector <- c("W01", "W02", "W03", "W04", "W05", "W06", "W07", "W08", "W09", "W10", 
#                  "W11", "W12", "W13", "W14", "W15", "W16", "W17", "W18", "W19", "W20", 
#                  "W21", "W22", "W23", "W24", "W25", "W26", "W27", "W28", "W29", "W30", 
#                  "W31", "W32", "W33", "W34", "W35", "W36", "W37", "W38", "W39", "W40")
# Snamevector <- c("S29", "S28", "S27", "S26", "S25", "S24", "S23", "S22", "S21", "S20", 
#                  "S19", "S18", "S17", "S16", "S15", "S14", "S13", "S12", "S11", "S10", 
#                  "S09", "S08", "S07", "S06", "S05", "S04", "S03", "S02", "S01")
# 
# namevector <- c(Wnamevector, Snamevector)
# #then set circos plot parameters
# circos.par("track.height" = 0.08 , start.degree=95, cell.padding = c(0.005, 0, 0.005, 0), gap.degree = gapsall)
# circos.initialize(factors = complete$group, x = complete$pos)
# circos.track(factors = complete$group, ylim = c(0,0.1), track.height = uh(1, "mm"),
#              panel.fun = function(x, y) {
#              })
# circos.update(sector.index = "W01", track.index = 1, 
#               bg.col = "orangered", bg.border = "orangered")
# circos.update(sector.index = "W02", track.index = 1, 
#               bg.col = "orangered", bg.border = "orangered")
# circos.update(sector.index = "W03", track.index = 1, 
#               bg.col = "orangered", bg.border = "orangered")
# circos.update(sector.index = "Z29", track.index = 1, 
#               bg.col = "orangered", bg.border = "orangered")		  
# circos.update(sector.index = "W04", track.index = 1, 
#               bg.col = "darkorange3", bg.border = "darkorange3")
# circos.update(sector.index = "W05", track.index = 1, 
#               bg.col = "darkorange3", bg.border = "darkorange3")
# circos.update(sector.index = "Z27", track.index = 1, 
#               bg.col = "darkorange3", bg.border = "darkorange3")
# circos.update(sector.index = "Z28", track.index = 1, 
#               bg.col = "black", bg.border = "black")		  
# circos.update(sector.index = "W06", track.index = 1, 
#               bg.col = "tan1", bg.border = "tan1")
# circos.update(sector.index = "W07", track.index = 1, 
#               bg.col = "tan1", bg.border = "tan1")
# circos.update(sector.index = "Z26", track.index = 1, 
#               bg.col = "tan1", bg.border = "tan1")		  
# circos.update(sector.index = "W08", track.index = 1, 
#               bg.col = "gold", bg.border = "gold")
# circos.update(sector.index = "Z25", track.index = 1, 
#               bg.col = "gold", bg.border = "gold")		  
# circos.update(sector.index = "W09", track.index = 1, 
#               bg.col = "yellow", bg.border = "yellow")
# circos.update(sector.index = "Z24", track.index = 1, 
#               bg.col = "yellow", bg.border = "yellow")
# circos.update(sector.index = "W10", track.index = 1, 
#               bg.col = "yellowgreen", bg.border = "yellowgreen")
# circos.update(sector.index = "Z23", track.index = 1, 
#               bg.col = "yellowgreen", bg.border = "yellowgreen")
# circos.update(sector.index = "Z22", track.index = 1, 
#               bg.col = "white", bg.border = "white")
# circos.update(sector.index = "W11", track.index = 1, 
#               bg.col = "forestgreen", bg.border = "forestgreen")
# circos.update(sector.index = "W12", track.index = 1, 
#               bg.col = "forestgreen", bg.border = "forestgreen")
# circos.update(sector.index = "W13", track.index = 1, 
#               bg.col = "forestgreen", bg.border = "forestgreen")
# circos.update(sector.index = "Z21", track.index = 1, 
#               bg.col = "forestgreen", bg.border = "forestgreen")
# circos.update(sector.index = "W14", track.index = 1, 
#               bg.col = "darkturquoise", bg.border = "darkturquoise")
# circos.update(sector.index = "W15", track.index = 1, 
#               bg.col = "darkturquoise", bg.border = "darkturquoise")
# circos.update(sector.index = "Z20", track.index = 1, 
#               bg.col = "darkturquoise", bg.border = "darkturquoise")
# circos.update(sector.index = "W16", track.index = 1, 
#               bg.col = "royalblue", bg.border = "royalblue")
# circos.update(sector.index = "W17", track.index = 1, 
#               bg.col = "royalblue", bg.border = "royalblue")
# circos.update(sector.index = "Z19", track.index = 1, 
#               bg.col = "royalblue", bg.border = "royalblue")
# circos.update(sector.index = "W18", track.index = 1, 
#               bg.col = "lightsteelblue", bg.border = "lightsteelblue")
# circos.update(sector.index = "Z18", track.index = 1, 
#               bg.col = "lightsteelblue", bg.border = "lightsteelblue")
# circos.update(sector.index = "W19", track.index = 1, 
#               bg.col = "pink", bg.border = "pink")
# circos.update(sector.index = "W20", track.index = 1, 
#               bg.col = "pink", bg.border = "pink")
# circos.update(sector.index = "Z17", track.index = 1, 
#               bg.col = "pink", bg.border = "pink")
# circos.update(sector.index = "W21", track.index = 1, 
#               bg.col = "deeppink", bg.border = "deeppink")
# circos.update(sector.index = "Z16", track.index = 1, 
#               bg.col = "deeppink", bg.border = "deeppink")
# circos.update(sector.index = "W22", track.index = 1, 
#               bg.col = "brown", bg.border = "brown")
# circos.update(sector.index = "W23", track.index = 1, 
#               bg.col = "brown", bg.border = "brown")
# circos.update(sector.index = "Z15", track.index = 1, 
#               bg.col = "brown", bg.border = "brown")
# circos.update(sector.index = "W24", track.index = 1, 
#               bg.col = "orangered", bg.border = "orangered")
# circos.update(sector.index = "W25", track.index = 1, 
#               bg.col = "orangered", bg.border = "orangered")
# circos.update(sector.index = "Z14", track.index = 1, 
#               bg.col = "orangered", bg.border = "orangered")
# circos.update(sector.index = "W26", track.index = 1, 
#               bg.col = "darkorange3", bg.border = "darkorange3")
# circos.update(sector.index = "Z13", track.index = 1, 
#               bg.col = "darkorange3", bg.border = "darkorange3")
# circos.update(sector.index = "W27", track.index = 1, 
#               bg.col = "tan1", bg.border = "tan1")
# circos.update(sector.index = "W28", track.index = 1, 
#               bg.col = "tan1", bg.border = "tan1")
# circos.update(sector.index = "Z12", track.index = 1, 
#               bg.col = "tan1", bg.border = "tan1")
# circos.update(sector.index = "W29", track.index = 1, 
#               bg.col = "gold", bg.border = "gold")
# circos.update(sector.index = "Z11", track.index = 1, 
#               bg.col = "gold", bg.border = "gold")
# circos.update(sector.index = "W30", track.index = 1, 
#               bg.col = "yellow", bg.border = "yellow")
# circos.update(sector.index = "W31", track.index = 1, 
#               bg.col = "yellow", bg.border = "yellow")
# circos.update(sector.index = "Z10", track.index = 1, 
#               bg.col = "yellow", bg.border = "yellow")
# circos.update(sector.index = "W32", track.index = 1, 
#               bg.col = "yellowgreen", bg.border = "yellowgreen")
# circos.update(sector.index = "Z09", track.index = 1, 
#               bg.col = "yellowgreen", bg.border = "yellowgreen")
# circos.update(sector.index = "W33", track.index = 1, 
#               bg.col = "forestgreen", bg.border = "forestgreen")
# circos.update(sector.index = "Z08", track.index = 1, 
#               bg.col = "forestgreen", bg.border = "forestgreen")
# circos.update(sector.index = "W34", track.index = 1, 
#               bg.col = "darkturquoise", bg.border = "darkturquoise")
# circos.update(sector.index = "Z07", track.index = 1, 
#               bg.col = "darkturquoise", bg.border = "darkturquoise")
# circos.update(sector.index = "W35", track.index = 1, 
#               bg.col = "royalblue", bg.border = "royalblue")
# circos.update(sector.index = "Z06", track.index = 1, 
#               bg.col = "royalblue", bg.border = "royalblue")
# circos.update(sector.index = "W36", track.index = 1, 
#               bg.col = "lightsteelblue", bg.border = "lightsteelblue")
# circos.update(sector.index = "Z05", track.index = 1, 
#               bg.col = "lightsteelblue", bg.border = "lightsteelblue")
# circos.update(sector.index = "W37", track.index = 1, 
#               bg.col = "pink", bg.border = "pink")
# circos.update(sector.index = "Z03", track.index = 1, 
#               bg.col = "pink", bg.border = "pink")
# circos.update(sector.index = "Z04", track.index = 1, 
#               bg.col = "black", bg.border = "black")
# circos.update(sector.index = "W38", track.index = 1, 
#               bg.col = "deeppink", bg.border = "deeppink")
# circos.update(sector.index = "W39", track.index = 1, 
#               bg.col = "deeppink", bg.border = "deeppink")
# circos.update(sector.index = "Z02", track.index = 1, 
#               bg.col = "deeppink", bg.border = "deeppink")
# circos.update(sector.index = "W40", track.index = 1, 
#               bg.col = "brown", bg.border = "brown")
# circos.update(sector.index = "Z01", track.index = 1, 
#               bg.col = "brown", bg.border = "brown")
# 
# circos.track(factors = complete$group, ylim = c(0,1),
#              panel.fun = function(x, y) {
#                circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
#                            namevector[CELL_META$sector.numeric.index], cex = 0.5)
#              })
# 
# #loop through linkage data frame and extract info to make links
# for (row in 1:(nrow(qual30))){
#   a <- qual30$SyntLG[row]
#   b <- qual30$cpos[row]
#   e <- qual30$NewSalm[row]
#   d <- qual30$cmpos[row]
#   fact <- as.factor(qual30$SyntLG)[row]
#   circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d) + (100*corFactor))), h=0.9, col = add_transparency("darkgrey", transparency = 0.65))
#   circos.lines(c(as.integer(b)), 1, as.character(a), col = "darkgrey", type = 'h')
#   circos.lines(c(as.integer(d)), 1, as.character(e), col = "darkgrey", type = 'h')
# }
# 
# dev.off()
# 
# W02groupss01 <- subset(qual30, qual30$SyntLG == "W02" & qual30$SalmChrom == "Z01")
# W02groupss19 <- subset(qual30, qual30$SyntLG == "W02" & qual30$SalmChrom == "Z19")

#analysis of distribution of links AORe/LORe for revised manuscript
#the  essential df is qual30 
for (i in 1:length(qual30$qual)){
  correSalm <- as.character(qual30$NewSalm)[i]
  origlim <- as.integer(qual30$SalmPos)[i]
  inset <- subset(chromarmtestreal, as.character(chromarmtestreal$SALMName) == correSalm)
  if(origlim < as.integer(inset$EndBP)[1]){
    qual30$ARM[i] <- as.character(inset$SalmArm)[1]
  } else if(origlim < as.integer(inset$EndBP)[2]){
    qual30$ARM[i] <- as.character(inset$SalmArm)[2]
  } else if(origlim < as.integer(inset$EndBP)[3]){
    qual30$ARM[i] <- as.character(inset$SalmArm)[3]
  }
}

for (i in 1:length(qual30$SalmChrom)){
  if (as.character(qual30$SalmChrom)[i] == as.character(qual30$equivalentSalm)[i]){
    qual30$SameDiff[i] <- "grey"
  } else {
    qual30$SameDiff[i] <- "black"
  }
}

library(plyr)
TT <- count(qual30, qual30$ARM)
for (i in 1:length(TT$`qual30$ARM`)){
  arm <- as.character(TT$`qual30$ARM`)[i]
  armrow <- subset(chromarmtestreal, chromarmtestreal$SalmArm == arm)
  TT$armlength[i] <- armrow$LengthBP
  TT$ORE[i] <- as.character(armrow$ORE)
}

# add 8q data:
colnames(TT)
name8 <- "8q"
n8 <- 0
armlength8df <- subset(chromarmtestreal, chromarmtestreal$SalmArm == "8q")
armlength8 <- armlength8df$LengthBP
ORE8 <- "L"
nrow(TT)
row8q <- c(name8, n8, armlength8, ORE8)
TT <- rbind(TT, row8q)
#expected number of qual 30 loci per arm?
TotalSyntLoci <- nrow(qual30)
#expected number per bp?
LociPerBP <- TotalSyntLoci/(sum(as.numeric(TT$armlength)))

TT$armlength <- as.integer(TT$armlength)
TT$n <- as.integer(TT$n)

TT$expectedLoci <- LociPerBP*TT$armlength
TT$RATIO <- TT$expectedLoci/TT$n
dev.off()
plot(as.factor(TT$`qual30$ARM`), TT$RATIO)


#now order by ratio
TT$`qual30$ARM` <- as.factor(TT$`qual30$ARM`)
TT$`qual30$ARM` <- reorder(TT$`qual30$ARM`, TT$RATIO)
plot(TT$`qual30$ARM`, TT$RATIO, pch = 16, cex.axis = 0.5, type = "p", col = "red")

#now to get grey/black ratio for each chromosome arm:
for (i in 1:length(TT$armlength)){
  m <- as.character(TT$`qual30$ARM`)[i]
  TT$black[i] <- nrow(subset(qual30, as.character(qual30$ARM) == m & as.character(qual30$SameDiff) == "black"))
  TT$grey[i] <- nrow(subset(qual30, as.character(qual30$ARM) == m & as.character(qual30$SameDiff) == "grey"))
}

TT$blackgreyratio <- TT$black/TT$grey
greyblack <- as.data.frame(TT[,6:7])
greyblack$arm <- TT$`qual30$ARM`
greyblack$RATIO <- TT$RATIO


TT$ORE <- as.factor(TT$ORE)
TT$ORE <- as.character(TT$ORE)

#data is in the wrong format for a boxplot so:
plotDF_A <- subset(TT, TT$ORE == "A")
plotDF_P <- subset(TT, TT$ORE == "P")
plotDF_L <- subset(TT, TT$ORE == "L")

A_RATIO <- plotDF_A$RATIO
P_RATIO <- plotDF_P$RATIO
L_RATIO <- plotDF_L$RATIO

#ttest showing difference between AORe and LORe regions
wilcox_Res <- wilcox.test(A_RATIO, L_RATIO)
wilcox_Res

numbA <- nrow(plotDF_A)
numP <- nrow(plotDF_P)
numL <- nrow(plotDF_L)

#figure for paper showing the number of 
boxplot(A_RATIO, P_RATIO, L_RATIO, ylab = "Expected/Observed Number Of Mappings", xlab = "Chromosome Arm Ohnologue Resolution Model", names = c((paste("AORe", " N =", as.character(numbA))), (paste("Partial LORe", " N =", as.character(numP))), (paste("LORe", " N =", as.character(numL)))))

########################################################################################################################################
########################################################################################################################################
#
#                         FIGURE 3 In publication
########################################################################################################################################
########################################################################################################################################
boxplot(A_RATIO, L_RATIO, ylab = "Expected/Observed Number Of Mappings", xlab = "Chromosome Arm Ohnologue Resolution Model", names = c((paste("AORe", " N =", as.character(numbA))), (paste("LORe", " N =", as.character(numL)))))
abline(h = 1, col = "grey", lty = 2)
points(1.45, 10, pch = "***", cex = 1.5)
points(1.5, 10, pch = "***", cex = 1.5)
points(1.55, 10, pch = "***", cex = 1.5)

for (i in 1:length(P_RATIO)){
  points(1.5, P_RATIO[i])
}

# hist(as.integer(samfile$X.1))
# hist(as.integer(SyntMap$qual))
# subsetz28<- subset(SyntMap, SyntMap$SalmChrom == "Z28")
# subsetz29<- subset(SyntMap, SyntMap$SalmChrom == "Z29")
# hist(as.integer(subsetz28$qual))
# 
# hist(as.integer(subsetz29$qual))

write.csv(TT, "LOReAOReDF.csv")

# #checking some chromosomea arms
# ssa07q <- subset(qual30, qual30$ARM == "7q")
# ssa03q <- subset(qual30, qual30$ARM == "3q")
# 
# ssa01p <- subset(qual30, qual30$ARM == "1p")
# ssa01qa <- subset(qual30, qual30$ARM == "1qa")
# ssa01qb <- subset(qual30, qual30$ARM == "1qb")
# ssa19qa <- subset(qual30, qual30$ARM == "19qa")
# ssa19qb <- subset(qual30, qual30$ARM == "19qb")
# 
# W07 <- subset(qual30, qual30$SyntLG == "W07")
# W17 <- subset(qual30, qual30$SyntLG == "W17")
# W38 <- subset(qual30, qual30$SyntLG == "W38")
# W39 <- subset(qual30, qual30$SyntLG == "W39")
# ssa14p <- subset(qual30, qual30$ARM == "14qb")
