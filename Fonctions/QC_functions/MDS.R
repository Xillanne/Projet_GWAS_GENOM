### SCRIPT faisant partie des QC fonction ###
# Trace le plotMDS

#setwd("~/Etudes/M2/GENOM/Projet")
# Charger les résultats MDS générés par PLINK
mds_results <- read.table("mds_results.mds", header = TRUE)

mean1 <- mean(mds_results$C1)
std1 <- sd(mds_results$C1)

mean2 <- mean(mds_results$C2)
std2 <- sd(mds_results$C2)

# Filtrer les points en fonction de la moyenne et de l'écart type
to_keep_3 <- mds_results[abs(mds_results$C1 - mean1) <= 3*std1 & abs(mds_results$C2 - mean2) <= 3*std2, ]

mds_filtered_3 <- mds_results[rownames(to_keep_3),]


# Créer un fichier avec les individus à supprimer
to_remove <- rownames(mds_results)[!(row.names(mds_results) %in% rownames(to_keep_3))]
mds_to_remove <- mds_results[to_remove,]
ind_to_remove <- mds_to_remove[,1:2]
write.table(ind_to_remove,file='fail_MDS_ind.txt',sep='\t',quote = FALSE,row.names = FALSE)

pdf('PlotMDS_Final.pdf')
# Créer un vecteur pour stocker les couleurs en fonction du statut de chaque individu
colors <- ifelse(row.names(mds_results) %in% rownames(to_keep_3), "black", "red")

# Tracer le nuage de points pour les individus à garder en rond noir
plot(mds_results$C1, mds_results$C2, main = "MDS Plot=", xlab = "Dimension 1", ylab = "Dimension 2", pch = 16, col = colors)


# Afficher le graphique
dev.off()
