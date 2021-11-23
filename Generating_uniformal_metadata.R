setwd("~/Desktop/Projects_2018/Baby_pilot/01.Data_analysis/03.SCRIPTS/Clean_code/01.Metadata_sorting/")

##### LIBRARIES #####
library(stringr)
library(dplyr)
##### FUNCTIONS #####

##### INPUT DATA #####

metadata_VLP2 <- read.table("../../Analysis_for_paper/INPUT/summary_samples_metadata_TP_data.txt", sep='\t', header=T, row.names = 1)
metadata_VLP2$B_or_M <- if_else(metadata_VLP2$STATUS=="Baby", "B", "M")
metadata_VLP2$merging_id <- paste(metadata_VLP2$B_or_M, metadata_VLP2$LLNEXT_ID, metadata_VLP2$TIMEPOINT, sep='_')
metadata_VLP2$LABEL_ID <- row.names(metadata_VLP2)
trishla_metadata1 <- read.table("../../Analysis_for_paper/INPUT/phenotypes_complete_pilot_feb_2021.txt", sep='\t', header=T, row.names = 1)
trishla_metadata1$Family_ID <- str_pad(trishla_metadata1$Pair_ID, 2, pad = "0")
trishla_metadata1$Family_ID <- str_pad(trishla_metadata1$Family_ID, 3, pad = "F")
trishla_metadata1$pseudo_id <- str_pad(trishla_metadata1$Pair_ID, 2, pad = "0")
trishla_metadata1$B_or_M <- if_else(trishla_metadata1$Type=="Baby", "B", "M")
trishla_metadata1$merging_id <- paste(trishla_metadata1$B_or_M, trishla_metadata1$PSEUDOID_number, trishla_metadata1$Timepoint, sep = '_')

trishla_metadata1$Timepoint_continuous <- NA
trishla_metadata1[!is.na(trishla_metadata1$Timepoint) & trishla_metadata1$Timepoint=="P3",]$Timepoint_continuous <- 0
trishla_metadata1[!is.na(trishla_metadata1$Timepoint) & trishla_metadata1$Timepoint=="P7",]$Timepoint_continuous <- 4
trishla_metadata1[!is.na(trishla_metadata1$Timepoint) & trishla_metadata1$Timepoint=="B",]$Timepoint_continuous <- 6
trishla_metadata1[!is.na(trishla_metadata1$Timepoint) & trishla_metadata1$Timepoint=="M1",]$Timepoint_continuous <- 7
trishla_metadata1[!is.na(trishla_metadata1$Timepoint) & trishla_metadata1$Timepoint=="M2",]$Timepoint_continuous <- 8
trishla_metadata1[!is.na(trishla_metadata1$Timepoint) & trishla_metadata1$Timepoint=="M3",]$Timepoint_continuous <- 9

tmp_metadata_merged <- merge(trishla_metadata1[,c(1,3,4,11:13,18:25, 39,42,43)], metadata_VLP2, by="merging_id", all.y = T)
all.equal(tmp_metadata_merged$Timepoint, tmp_metadata_merged$TIMEPOINT)
tmp_metadata_merged$TIMEPOINT <- NULL
all.equal(tmp_metadata_merged$Type, tmp_metadata_merged$STATUS)
tmp_metadata_merged$STATUS <- NULL
all.equal(tmp_metadata_merged$PSEUDOID_number, tmp_metadata_merged$LLNEXT_ID)
tmp_metadata_merged$LLNEXT_ID <- NULL
tmp_metadata_merged$KINSHIP <- NULL
tmp_metadata_merged[tmp_metadata_merged$PSEUDOID_number=="44421250",40] <- "B2"
tmp_metadata_merged[tmp_metadata_merged$PSEUDOID_number=="44424334",40] <- "B2"
tmp_metadata_merged[tmp_metadata_merged$B_or_M=="B",40] <- "B1"
tmp_metadata_merged$new_sample_id <- paste(tmp_metadata_merged$B_or_M, tmp_metadata_merged$Family_ID, tmp_metadata_merged$Timepoint, sep = '_')
tmp_metadata_merged$new_individual_id <- paste(tmp_metadata_merged$B_or_M, tmp_metadata_merged$Family_ID, sep = '_')

virome_metadata <- tmp_metadata_merged
row.names(virome_metadata) <- virome_metadata$LABEL_ID


microbiome_metadata <- subset(trishla_metadata1, trishla_metadata1$Excluded_anlysis=="No")

a <- as.data.frame(matrix(NA, nrow = length(unique(microbiome_metadata$Family_ID)), ncol=2))
a$V1 <- unique(microbiome_metadata$Family_ID)
for (i in 1:30) {
  a[i,2] <- length(unique(microbiome_metadata[microbiome_metadata$Family_ID==a[i,1],]$PSEUDOID_number))
}

microbiome_metadata[microbiome_metadata$PSEUDOID_number=="44421250",41] <- "B2"
microbiome_metadata[microbiome_metadata$PSEUDOID_number=="44424334",41] <- "B2"
microbiome_metadata[microbiome_metadata$B_or_M=="B",41] <- "B1"
microbiome_metadata$new_sample_id <- paste(microbiome_metadata$B_or_M, microbiome_metadata$Family_ID, microbiome_metadata$Timepoint, sep = '_')
microbiome_metadata$new_individual_id <- paste(microbiome_metadata$B_or_M, microbiome_metadata$Family_ID, sep = '_')

b <- as.data.frame(matrix(NA, nrow = length(unique(microbiome_metadata[microbiome_metadata$Type=="Baby",]$PSEUDOID_number)), ncol=2))
b$V1 <- unique(microbiome_metadata[microbiome_metadata$Type=="Baby",]$PSEUDOID_number)
for (i in 1:32) {
  b[i,2] <- length(unique(microbiome_metadata[microbiome_metadata$PSEUDOID_number==b[i,1],]$Timepoint))
}
sum(b$V2>2) #26

sum((microbiome_metadata$Type=="Mother")==T) #93
sum((microbiome_metadata$Type=="Baby")==T) #90

sum((virome_metadata$Type=="Mother" & virome_metadata$Timepoint=="B")==T)
sum((virome_metadata$Type=="Baby" & virome_metadata$Timepoint=="M3")==T) 
for (i in 1:30) {
  a[i,3] <- length(unique(virome_metadata[virome_metadata$Family_ID==a[i,1],]$PSEUDOID_number))
}
sum(a$V3>1)
for (i in 1:32) {
  b[i,3] <- length(unique(virome_metadata[virome_metadata$PSEUDOID_number==b[i,1],]$Timepoint))
}
sum(b$V3>1)

colnames(a) <- c("Family_ID", "N_family_members_microbiome", "N_family_members_virome")
colnames(b) <- c("Baby_ID", "N_available_tp_microbiome", "N_available_tp_virome")


write.table(virome_metadata, "../INPUT/Virome_metadata_phenos.txt", sep='\t', quote=F)
write.table(microbiome_metadata, "../INPUT/Microbiome_metadata_phenos.txt", sep='\t', quote=F)
