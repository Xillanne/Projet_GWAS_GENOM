library(dplyr)
#### Code pour ajouter au fichier .sample les colonnes phénotypes et covariate avant SNPTEST####
## Répértoire où se trouve les fichiers
setwd("home/Alix/Etudes/M2/GENOM/Projet")
## Nom de l'output
output<-'chr5.sample'

# Lis les trois fichiers
sample <- read.table('phased_chr5.sample',header=TRUE)
phenotype <- read.table('phenotype_plink.txt',header=TRUE)
covar <- read.table('covar_plink.txt',header=TRUE)
# Homogénise le nom des colonnes
colnames(covar)[colnames(covar) == "gid"] <- "ID_1"
colnames(covar)[colnames(covar) == "gid.1"] <- "ID_2"
colnames(phenotype)[colnames(phenotype) == "ID"] <- "ID_1"
colnames(phenotype)[colnames(phenotype) == "FID"] <- "ID_2"

# Merge les trois Df, mets NA si individus par présents dans les autres
merged_df <- sample %>%
  merge(covar, by = c("ID_1", "ID_2","sex"), all = TRUE) %>%
    merge(phenotype, by = c("ID_1", "ID_2"), all = TRUE)
    
# Garde uniquement les individus qui étaient dans sample
merged_final <- sample %>%
  merge(merged_df)

# Nomme les colonne selon la norme de Impute2
merged_final[1,]['age']='C'
merged_final[1,]['BMI']='C'
merged_final[1,]['ApoA1']='P'
merged_final[1,]['HDL_CHOL']='P'

# Ecris les résultats
write.table(merged_final,output,quote=FALSE,row.names = FALSE)
