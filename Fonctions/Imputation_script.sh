#!/bin/bash

#SBATCH --job-name=Imputation_2
#SBATCH --mem=40GB
#SBATCH --account=2319_ens_hts

## Script qui décrit les étapes de l'imputation (après phasage)

## PATH 
#Chemin d'accès au génome de référence
#ref_hap='/home/xillanne/GENOM/data/ref_haplotypes/' 
ref_hap='/home/xillanne/GENOM/data/1000GP_Phase3/' 

cd Imputation

#On parcourt chaque chromosme
for i in 2 5 13 16
do
	#On se place dans le dossier spécifique au chrom
	cd chr$i
	#On récupère la position du dernier SNP 
	len=$(awk '{print $3}' 'phased_chr'$i'.haps' | tail -n 1)
	#ON va parcourir les position
	step=5000000 #5Mb : j'ai enlevé un 0
	k=1
	for start in $(seq 1 $step $len) #de 1 à len max avec des pas de 5Mb
	do     
		end=$((start + step - 1))     
		if [ $end -gt $len ] #Si plus petit que len max
		then         
			end=$len
		fi
		echo $start $end
		impute2 -use_prephased_g \
		 -m $ref_hap'genetic_map_chr'$i'_combined_b37.txt' \
		 -h $ref_hap'1000GP_Phase3_chr'$i'.hap.gz' \
		 -l $ref_hap'1000GP_Phase3_chr'$i'.legend.gz' \
		 -known_haps_g phased_chr$i.haps \
		 -int $start $end \
		 -Ne 20000 \
		 -o 'chr'$i'_chunk_'$k -verbose
		 
		 ((k+=1))
		 done
	 cd ..
done
	 
