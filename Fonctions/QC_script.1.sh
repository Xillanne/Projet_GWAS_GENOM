#!/bin/bash
# Script regroupant les étapes de QC

# Chemin vers les fonctions de QC
path='/home/xillanne/GENOM/QC_functions'
# Chemin vers les données Bed_bim_fam
data='/home/xillanne/GENOM/data/Bed_bim_fam/'

# Créer un fichier pour référencer les SNPs retirés
removed_snps_file='/home/xillanne/GENOM/data/removed_snps.txt'

# Supprimer le fichier s'il existe déjà
if [ -f "$removed_snps_file" ]; then
    rm "$removed_snps_file"
fi

# Supprimer le fichier de sortie s'il existe déjà
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

for i in 2 5 13 16 # Parcourt les chromosomes
do
    cd $data
    file="chr$i"
    mkdir $file
    cd $file
    # Nom du fichier de sortie
    output_file='filtering_QC_step_by_step.txt'
    
    # Supprimer le fichier de sortie s'il existe déjà
    if [ -f "$output_file" ]; then
        rm "$output_file"
    fi

    ### Étape 1 ###

    # Examiner la non-concordance par individu et par SNP et créer des histogrammes.
    plink --bfile "$data$file" --missing
    # sortie : plink.imiss et plink.lmiss, ces fichiers montrent respectivement la proportion de SNPs manquants par individu et la proportion d'individus manquants par SNP.
    
    
    # Générer des graphiques pour visualiser les résultats de la non-concordance.
    Rscript --no-save "$path/hist_miss.R"

    # Supprimer les SNPs et les individus avec des niveaux élevés de non-concordance
    # Les deux commandes de QC suivantes ne supprimeront aucun SNP ou individu.
    # Commencer avec un seuil de 0.2, et après nous pouvons appliquer un seuil plus strict.
    # On peut essayer l'exclusion comme dans l'article de Laurie : on essaie, on calcule la concordance (semble long et compliqué)
    # Comment choisir les seuils ?

    # Supprimer les SNPs avec une non-concordance >0.02.
    plink --bfile "$data$file" --geno 0.02 --make-bed --out "${file}_2"

    

    # Compter le nombre de SNPs et d'individus après cette étape
    num_snps_before=$(wc -l "$data$file.bim" | cut -d' ' -f1)
    num_indiv_before=$(wc -l "$data$file.fam" | cut -d' ' -f1)
    
    num_snps_after=$(wc -l "${file}_3.bim" | cut -d' ' -f1)
    num_indiv_after=$(wc -l "${file}_3.fam" | cut -d' ' -f1)

    # Afficher le nombre de SNPs et d'individus supprimés et écrire dans le fichier
    echo -e "Étape 1: Call-rate 0.02\n Avant - SNPs: $num_snps_before, Individus: $num_indiv_before | Après - SNPs: $num_snps_after, Individus: $num_indiv_after" >> "$output_file"

    ### Étape 2 ###

    # Générer un fichier de type bfile avec uniquement les SNPs autosomiques et supprimer les SNPs avec une fréquence allélique mineure (MAF) faible.

    # Sélectionner uniquement les SNPs autosomiques (c'est-à-dire, des chromosomes 1 à 22).
    awk '{ if ($1 == 2 || $1 == 5 || $1 ==13 || $1 == 16) print $2 }' "${file}_3.bim" > snp.txt
    plink --bfile "${file}_3" --extract snp.txt --make-bed --out "${file}_4"

    # Générer un graphique de la distribution de la fréquence allélique mineure (MAF).
    plink --bfile "${file}_4" --freq --out MAF_check
    Rscript --no-save "$path/MAF_check.R"

    # Supprimer les SNPs avec une fréquence allélique mineure (MAF) faible.
    plink --bfile "${file}_4" --maf 0.05 --make-bed --out "${file}_5"

    # Un seuil MAF conventionnel pour une GWAS régulière est entre 0.01 ou 0.05, selon la taille de l'échantillon.
    # Sur un N=190 857, N large donc un seuil MAF = 0.01, sauf que l'on considère les chromosomes 1 à 1 donc N petit donc 0.05

    # Compter le nombre de SNPs et d'individus après cette étape
    num_snps_after=$(wc -l "${file}_5.bim" | cut -d' ' -f1)
    num_indiv_after=$(wc -l "${file}_5.fam" | cut -d' ' -f1)

    # Afficher le nombre de SNPs et d'individus supprimés et écrire dans le fichier
    echo -e "Étape 2: MAF 0.05\n Après - SNPs: $num_snps_after, Individus: $num_indiv_after" >> "$output_file"

    ### Étape 3 ###

    # Supprimer les SNPs qui ne sont pas en équilibre de Hardy-Weinberg (HWE).
    # Indicateur commun d'erreur de génotypage, peut également indiquer une sélection évolutive
    # Vérifier la distribution des valeurs p de HWE de tous les SNPs.

    plink --bfile "${file}_5" --hardy
    # Sélectionner les SNPs avec une valeur p de HWE inférieure à 0.000001, nécessaire pour l'un des deux graphiques générés par le script R suivant, permet de zoomer sur les SNPs fortement déviants.
    awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
    Rscript --no-save "$path/hwe.R"

    # Par défaut, l'option --hwe de plink ne filtre que pour les témoins.
    # On a des traits quantitatifs (PAS SÛR, tels que le taux de HDL-Cholestérol et le taux d’Apolipoprotéine A1, ApoA1) 
    # Donc, comme recommandé dans l'article, on filtre pour <1e-6
    plink --bfile "${file}_5" --hwe 1e-6 --hwe-all --make-bed --out "${file}_6"
    
    # Compter le nombre de SNPs et d'individus après cette étape
    num_snps_after=$(wc -l "${file}_6.bim" | cut -d' ' -f1)
    num_indiv_after=$(wc -l "${file}_6.fam" | cut -d' ' -f1)

    # Afficher le nombre de SNPs et d'individus supprimés et écrire dans le fichier
    echo -e "Étape 3: HWE 1e-6\n Après - SNPs: $num_snps_after, Individus: $num_indiv_after" >> "$output_file"

    
    # On retire A/T C/G SNPs, c'est à dire les SNP ambigus
    awk '(($5=="C"&&$6=="G")||($5=="G"&&$6=="C")||($5=="A"&&$6=="T")||($5=="T"&&$6=="A")) {print $2} ' "${file}_6.bim" > AT_CG_SNPS.txt
    
    #On les exclue
    plink --bfile "${file}_6" --exclude AT_CG_SNPS.txt --make-bed --out "${file}_7"
    
    #On filtre les SNP non bi-allélique
    plink --bfile "${file}_7" --biallelic-only strict --make-bed --out "${file}_8"

    mkdir Pdf
    mv *.pdf Pdf

done

