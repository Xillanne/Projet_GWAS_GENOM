#!/bin/bash
# Script qui regroupe les étapes pour supprimer les outliers

path='/home/xillanne/GENOM/QC_functions'
data='/home/xillanne/GENOM/data/Bed_bim_fam/'

output_file='filtering_QC_step_by_step.txt'
to_remove='indtoremove.txt'
map='chr_final.map'
pihat='pihat_min0.2_in_founders.genome'

### On concatène les fichiers ensemble
cd $data
plink --bfile chr5/chr5_8 --bmerge chr2/chr2_8.bed chr2/chr2_8.bim chr2/chr2_8.fam --make-bed --out merged_chr2_chr5
plink --bfile merged_chr2_chr5 --bmerge chr13/chr13_8.bed chr13/chr13_8.bim chr13/chr13_8.fam --make-bed --out merged_chr2_chr5_chr13
plink --bfile merged_chr2_chr5_chr13 --bmerge chr16/chr16_8.bed chr16/chr16_8.bim chr16/chr16_8.fam --make-bed --out merged_dataset
echo -e "merged_chr2_chr5\nmerged_chr2_chr5_chr13\nmerged_dataset" > merge_list.txt
plink --merge-list merge_list.txt --make-bed --out final_merged


# Supprimer le fichier s'il existe déjà
if [ -f "indtoremove.txt" ]; then
    rm "indtoremove.txt"
fi

# Supprimer le fichier s'il existe déjà
if [ -f "chr.map" ]; then
    rm "chr.map"
fi

# Supprimer le fichier s'il existe déjà
if [ -f "pihat_min0.2_in_founders.genome" ]; then
    rm "pihat_min0.2_in_founders.genome"
fi

if [ -f "chr.bim" ]; then
    rm "chr.bim"
fi

if [ -f "chr.bed" ]; then
    rm "chr.bed"
fi


cd $data

   ############### OUTLIER - HÉTÉROZYGOSITÉ #################

# Supprimer les individus avec une non-concordance >0.02.
plink --bfile final_merged --mind 0.02 --make-bed --out final_merged_1
    
    ### Étape 5 ###

    # Générer un graphique de la distribution du taux d'hétérozygotie de vos sujets.
    # Et supprimer les individus dont le taux d'hétérozygotie dévie de plus de 3 écart-types de la moyenne.

    # Les vérifications d'hétérozygotie sont effectuées sur un ensemble de SNPs qui ne sont pas fortement corrélés.
    
    ## ON A PAS DE FICHIER AVEC DES INVERSIONS (je crois)
    # Par conséquent, pour générer une liste de SNPs non (fortement) corrélés, nous excluons les régions d'inversion élevée (inversion.txt [Régions de LD élevée]) et éclaircissons les SNPs à l'aide de la commande --indep-pairwise.
    # Les paramètres  50 5 0.2  correspondent respectivement à : la taille de la fenêtre, le nombre de SNPs pour décaler la fenêtre à chaque étape, et le coefficient de corrélation multiple pour un SNP régressé sur tous les autres SNP simultanément.

plink --bfile final_merged_1 --indep-pairwise 50 5 0.2 --out indepSNP #--exclude inversion.txt --range
    # Remarque, ne supprimez pas le fichier indepSNP.prune.in, nous utiliserons ce fichier dans les étapes ultérieures du tutoriel.

plink --bfile final_merged_1 --extract indepSNP.prune.in --het --out R_check
    # Ce fichier contient votre jeu de données éclairci.

    # Graphique de la distribution du taux d'hétérozygotie
Rscript --no-save "$path/check_heterozygosity_rate.R"

    # Le code suivant génère une liste d'individus dont le taux d'hétérozygotie dévie de plus de 3 écart-types de la moyenne du taux d'hétérozygotie.
    # Pour la manipulation des données, nous recommandons l'utilisation de UNIX. Cependant, lors de calculs statistiques, R peut être plus pratique, d'où l'utilisation de Rscript pour cette étape :
Rscript --no-save "$path/heterozygosity_outliers_list.R"

    # Sortie de la commande ci-dessus : fail-het-qc.txt .

    # Adapter ce fichier pour le rendre compatible avec PLINK en supprimant toutes les guillemets du fichier et en sélectionnant uniquement les deux premières colonnes.
sed 's/"//g' fail-het-qc.txt | awk '{print $1, $2}' >> $to_remove
    
# Supprimer les valeurs aberrantes du taux d'hétérozygotie.
plink --bfile final_merged_1 --remove $to_remove --make-bed --out final_merged_2

# Compter le nombre de SNPs et d'individus après cette étape
num_snps_after=$(wc -l "final_merged_1.bim" | cut -d' ' -f1)
num_indiv_after=$(wc -l "final_merged_1.fam" | cut -d' ' -f1)

# Afficher le nombre de SNPs et d'individus supprimés et écrire dans le fichier
echo -e "##############OUTLIERS ##############\nÉtape 1: Hétérozygotie +- 3std \n Après - SNPs: $num_snps_after, Individus: $num_indiv_after" >> "$output_file"  
    
    
    ### Etape 6 ###
    # On va checker la relatdness : de ce que j'ai vu autant d'ID familiaux que d'ID ind => pas de relations
    # AU cas où on va juste inclure les founders c.a.d les individus sans parents dans le dataset
plink --bfile final_merged_2 --filter-founders --make-bed --out final_merged_3
    
    # Comme on est dans le cas d'une population random, on applique un seuil de 0.2
plink --bfile final_merged_3 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders
    
    # Dans .genomes on retrouve les individus ont un pihat > 0.2. Ce ne sont pas les mêmes pour tous les chromosomes 
    # Y en a même plusieurs très élévés voire à 1
    
    
    # Génére la liste des lowest call rate pour tous les ind
plink --bfile final_merged_3 --missing
    

    
    # Identifie individus à supprimer à cause de relatedness
Rscript --no-save "$path/ind_relatedness.R"

# Supprime les individus avec reladness
plink --bfile final_merged_3 --remove fail_relatdness_ind.txt --make-bed --out final_merged_4
    
# Compter le nombre de SNPs et d'individus après cette étape
num_snps_after=$(wc -l "final_merged_3.bim" | cut -d' ' -f1)
num_indiv_after=$(wc -l "final_merged_3.fam" | cut -d' ' -f1)

    # Afficher le nombre de SNPs et d'individus supprimés et écrire dans le fichier
echo -e "Étape 2: Relatedness pihat>0.2 \n Après - SNPs: $num_snps_after, Individus: $num_indiv_after" >> "$output_file"
    


#On fait MDS
plink --bfile final_merged_4 --extract indepSNP.prune.in --genome --out genome_results
plink --bfile final_merged_4 --read-genome genome_results.genome --cluster --mds-plot 2 --out mds_results
#On visualise
Rscript --no-save "$path/MDS.R"
# Renvoie le fichier fail_MAD_ind.txt qui contient les ind avec une moyenne qui s'écarte de + ou - 2std

# On les supprime
plink --bfile final_merged_4 --remove fail_MDS_ind.txt --make-bed --out chr_clean


# On converti en .ped et .map
plink --bfile chr_clean --recode --out chr_clean

rm *merged*

