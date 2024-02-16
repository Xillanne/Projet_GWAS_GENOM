Projet GENOM 2023-2024 - Master 2 BIM
Données GWAS et variations phénotypiques
Alix BOUTHEROUE-DESMARAIS	Kenza FLIOU 	Yuliya LIM
Encadré par Sophie Garnier


Ce projet a pour but d'étudier les résultats de l’étude GHS (Gutenberg Health Study). Projet B indiqué dans l'énoncé.
Les dossiers .bed. bim, .fam, .ped, .map des chromosomes 2, 5, 13 et 16 de la cohorte de 3300 individus nous ont été fournis. Ainsi que les fichiers phenotype_plink.txt et covar_plink.txt.

PACKAGES REQUIS
PLINK v1.90
shapeit.v2.904.3.10.0
impute_v2.3
snptest_v2.5.6

PIPELINE

Dans le répértoire où se trouve les .bim,.bed et .fam : 

1. QC_script.1.sh : réalise les étapes de la QC liées aux SNP  sur les fichiers .bim, .bed, .fam - Le path des QC_functions doit être modifié

2. Outlier_script.1.sh : réalise les étapes de QC sur les individus - Le path des QC_functions doit être modifié

3. Phasage_script.sh : réalise le phasage

4. Imputation_script.sh : réalise l'imputation, les résultats obtenus doivent être complié pour les fichiers .gen et _info dans des fichiers chrN.impute2 et chrN_info

Dans le répértoire où se trouve les résultats

5. filter_imputation.py : filtre les résultats de l'imputation utilisé chrN.impute2 et chrN_info, les paths doivent être modifié

6. convert_sample.R : doit être modifié et utilisé pour ajouté les colonnes de phenotype_plink.txt et covar_plink.txt au fichier .sample obtenu après l'étape de phasage.

7. snptest_v2.5.6 -data chrN.impute2 chrN.sample -o chrN_QTL.pheno.out -method score -frequentist 1 -pheno pheno -cov_all   - Cmd à lancer dans le terminal, chrN.impute2 est le résultat obtenu après filtrage (5), chrN.sample est obtenu après (6), pheno correspond soit à ApoA1 oit à HDL_CHOL

8. Plot_GWAS.R : utilise le package qqman, sert à plotter les manhattans et QQ plots.




