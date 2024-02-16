#Script pour identifier les individus à supprimer à cause de la relatdness

related <- read.table('pihat_min0.2_in_founders.genome',header=TRUE)
call_rate <- read.table('plink.imiss',header=TRUE)

# Supprimer les lignes doublons basé sur les premières colonnes, cad les identifiants
verif <- unique(related[, c(2,4)]) #J'ai verifié que y avait pas de doublons (IID1 : A, IID2 :B & IID1 : B, IID2 : A)


# Initialiser une liste pour stocker les groupes d'identifiants
groupes <- list()

# Fonction pour ajouter un nouvel identifiant à un groupe existant ou créer un nouveau groupe
ajouter_a_groupe <- function(groupe, nouvel_iid) {
  # Vérifier si l'identifiant est déjà dans le groupe
  if (!(nouvel_iid %in% groupe)) {
    # Ajouter l'identifiant au groupe
    groupe <- c(groupe, nouvel_iid)
  }
  return(groupe)
}

# Pour chaque ligne unique dans verif, former les groupes d'identifiants
for (i in 1:nrow(verif)) {
  iid1 <- verif$IID1[i]
  iid2 <- verif$IID2[i]
  
  # Vérifier si iid1 et iid2 sont déjà dans un groupe existant
  groupe_trouve <- FALSE
  
  for (j in seq_along(groupes)) {
    groupe <- groupes[[j]]
    
    if (iid1 %in% groupe || iid2 %in% groupe) {
      # Ajouter iid1 et iid2 au groupe existant
      groupe <- ajouter_a_groupe(groupe, iid1)
      groupe <- ajouter_a_groupe(groupe, iid2)
      
      # Mettre à jour le groupe dans la liste
      groupes[[j]] <- groupe
      
      groupe_trouve <- TRUE
      break
    }
  }
  
  # Si iid1 et iid2 ne sont pas dans un groupe existant, créer un nouveau groupe
  if (!groupe_trouve) {
    nouveau_groupe <- c(iid1, iid2)
    groupes <- c(groupes, list(nouveau_groupe))
  }
}

# Afficher les groupes formés
#print(groupes)




# Ensuite on veut supprimer l'individus avec le plus petit call-rate
# Attention certains individus sont présents plusieurs fois
# Vecteur pour stocker les individus à garder
to_keep <- character(0)

# Fonction pour trouver l'individu avec le plus petit F_MISS dans un groupe
trouver_a_supprimer <- function(groupe, call_rate) {
  if (length(groupe) > 1) {
    # Trouver les indices correspondants dans call_rate
    indices <- which(call_rate$IID %in% groupe)
    
    # Récupérer les valeurs de F_MISS
    f_miss_values <- call_rate$F_MISS[indices]
    
    # Trouver l'indice de l'individu avec le plus petit F_MISS
    indice_a_supprimer <- indices[which.min(f_miss_values)]
    
    # Récupérer l'identifiant à supprimer
    iid_a_supprimer <- call_rate$IID[indice_a_supprimer]
    
    return(iid_a_supprimer)
  } else {
    return(NULL)
  }
}

# Pour chaque groupe, trouver l'individu avec le plus petit F_MISS et l'ajouter à to_remove
for (groupe in groupes) {
  iid_a_supprimer <- trouver_a_supprimer(groupe, call_rate)
  
  if (!is.null(iid_a_supprimer)) {
    to_keep <- c(to_keep, iid_a_supprimer)
  }
}

# Soustrait les ind à garder de tous les ID afin d'obtenir ceux à soustraire
ID <- unique(c(verif$IID1,verif$IID2))
to_remove <- ID[!(ID %in% to_keep)]

df_to_remove <- data.frame(FIID = to_remove, IID = to_remove)
# Supprime les lignes avec les mêmes valeurs IID (supprime les ligne avec IID IID du à la concaténation)
df_to_remove <- subset(df_to_remove, !grepl("IID", IID))

write.table(df_to_remove,file='fail_relatdness_ind.txt',sep='\t',row.names = FALSE,quote = FALSE)
