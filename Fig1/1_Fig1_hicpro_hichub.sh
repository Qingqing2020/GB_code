#!/bin/bash


## ---- hicpro ----
time /HiC-Pro_3.0.0/bin/HiC-Pro -c ./config_H.txt -i $INPUT/data -o $OUTPUT/Neu_H/

# nohup bash hicpro_H.sh 2>&1 &

## ---- turn hicpro output *.allValidPairs to *.hic file for juicer  ----

for Sample in $(ls $INPUTDIR/data| grep -v "rep")
do  
    
	echo "bash /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh  -i data/$Sample/${Sample}.allValidPairs  -g /HiC-Pro_3.0.0/annotation/chrom_hg19.sizes -j /juicer-1.6/scripts/common/juicer_tools.jar -t tmp -o H_hicpro2juicebox/ "

	echo ""


    bash /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh  -i data/$Sample/${Sample}.allValidPairs  -g /HiC-Pro_3.0.0/annotation/chrom_hg19.sizes -j /juicer-1.6/scripts/common/juicer_tools.jar -t tmp -o H_hicpro2juicebox/

    wait

done


## ---- hichub ----
OUTPUTDIR=./hic_results/hichub_c10q0001
INPUTDIR=./hic_results/H_hicpro2juicebox
## different combination of samples

Samples="Ctrl 4h 2d 4d"

arr=($Samples)

for i in {0..3}
do
	# echo $i
	for j in {0..3}
	do
		# echo $j
		if(( $i<$j ))
		then
			echo ${arr[$i]}
			echo ${arr[$j]}

			# every comparison ,run hichub
			mkdir  $OUTPUTDIR/${arr[$i]}_VS_${arr[$j]}
			cd $OUTPUTDIR/${arr[$i]}_VS_${arr[$j]}
			# step 1: --------hichub convert

			echo "hichub convert -i $INPUTDIR -f HL60_${arr[$i]}.allValidPairs.hic,HL60_${arr[$j]}.allValidPairs.hic -l HL60_${arr[$i]},HL60_${arr[$j]} -r 10000"

			hichub convert -i $INPUTDIR -f HL60_${arr[$i]}.allValidPairs.hic,HL60_${arr[$j]}.allValidPairs.hic -l HL60_${arr[$i]},HL60_${arr[$j]} -r 10000 
			# -----------------------------------------
			# it is strange that the output file is saved where the input file is. There should be an option to direct the output path.!!!!!
			# -----------------------------------------
			# step 2: -------hichub diff

			echo "hichub diff -i $INPUTDIR/Summary_HL60_${arr[$i]}_HL60_${arr[$j]}_Dense_Matrix.txt -l HL60_${arr[$i]},HL60_${arr[$j]} -r 10000 -c 10 -d 1 -p 0.001 -t 8"
			hichub diff -i $INPUTDIR/Summary_HL60_${arr[$i]}_HL60_${arr[$j]}_Dense_Matrix.txt -l HL60_${arr[$i]},HL60_${arr[$j]} -r 10000 -c 10 -d 1 -p 0.001 -t 8 

			
			cd ..

		fi

	done
done
