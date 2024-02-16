#!/bin/bash
#Script qui regroupe les étapes pour l'imputation
#On se place à l'endroit où le résultat de outlier_script.1.sh a été stocké


## PATH 
#Chemin d'accès au génome de référence
#ref_hap='/home/xillanne/GENOM/data/ref_haplotypes/' 
ref_hap='/home/xillanne/GENOM/data/1000GP_Phase3/' 

# Crée un dossier où on va stocker les résultats
mkdir Imputation
cd Imputation

for i in 2 5 13 16
do
	#On  crée un dossier pour chaque chromosome
	mkdir chr$i
	#On extrait les SNP du chromosome d'interet dans des fichiers .ped et .map
	# recode
	plink --bfile ../chr_clean --chr $i  --make-bed --out ./chr$i/unphased_chr$i 
        
        
        
        #Imputation 
        #On a >1000 individus alors on a pas besoin de se réferencer
        shapeit -B  ./chr$i/unphased_chr$i\
        -M $ref_hap'genetic_map_chr'$i'_combined_b37.txt'\
        -O ./chr$i/phased_chr$i\
        -T 8


done


