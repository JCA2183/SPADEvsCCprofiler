SPADE_R1_10min <- as.data.table(read.csv('SPADE_R1_10min_protein_features.csv'))
SPADE_R2_10min <- as.data.table(read.csv('SPADE_proteinfeatures_R2_10min.csv'))
SPADE_R3_10min <- as.data.table(read.csv('SPADE_proteinfeatures_R3_10min.csv'))
SPADE_R1_30min <- as.data.table(read.csv('SPADE_proteinfeatures_R1_30.csv'))
SPADE_R2_30min <- as.data.table(read.csv('SPADE_proteinfeatures_R2_30min.csv'))
SPADE_R3_30min <- as.data.table(read.csv('SPADE_proteinfeatures_R3_30min.csv'))

CCprofiler_R1_10min <- as.data.table(read.csv('proteinfeaturesR1_10min.csv'))
CCprofiler_R2_10min <- as.data.table(read.csv('proteinfeaturesR2_10min.csv'))
CCprofiler_R3_10min <- as.data.table(read.csv('proteinfeaturesR3_10min.csv'))
CCprofiler_R1_30min <- as.data.table(read.csv('proteinfeaturesR1_30min.csv'))
CCprofiler_R2_30min <- as.data.table(read.csv('proteinfeaturesR2_30min.csv'))
CCprofiler_R3_30min <- as.data.table(read.csv('proteinfeaturesR3_30min.csv'))

CCprofiler_files <- ls(pattern="CCprofiler_R\\d_\\w+")

filtered_CCprofiler <- list()
patron= "DECOY_\\w+"
for (i in CCprofiler_files){
  var <- get("i")
  filtered_CCprofiler[[i]] <- var[!grepl(patron, var$protein_id),]
}



