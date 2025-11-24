#!/bin/bash

# RNAseq pipeline: use hisat2 htseq

GTFFILE=/home/shared_data/Annotation/UCSC/Human_Genome/HG19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2015-07-17-14-32-32/Genes/hg19_2015.gtf

HISAT2Index=/home/qchen/hisat2index/hg19/genome

THREADS=8
# build index
# hisat2-build genome.fa index_name

INPUTDIR=./data/trim_adapter 
OUTPUTDIR=./results


for i in $(ls $INPUTDIR|cat | grep "_1_val_1.fq.gz")
do
	echo ${i}
	while read line
     do
          # echo ${line}
          FileName="$(echo ${line} | cut -d' ' -f1)"
          if [[ ${FileName} = ${i%%_1_val_1.fq.gz} ]]
          then
               SAMPLE="$(echo ${line} | cut -d' ' -f2)"
               echo "${FileName} is ${SAMPLE}"
               echo ""
          fi

     done < ./RNAseq_HL60_SampleInfo.txt

     echo "The running sample is ${SAMPLE}"
     echo ""
 # ----------------- mapping  -----------------
          FASTQFILE1=$i

          FASTQFILE2=${i/_1_val_1.fq.gz/_2_val_2.fq.gz}
        	
	     mkdir $OUTPUTDIR/$SAMPLE
	     cd $OUTPUTDIR/$SAMPLE
          # this is PE data:
	     echo "hisat2 -t -p $THREADS --no-unal --no-mixed --non-deterministic --no-discordant -x $HISAT2Index -1 $INPUTDIR/$FASTQFILE1 -2 $INPUTDIR/$FASTQFILE2 -S ${SAMPLE}.sam"
	     hisat2 -t -p $THREADS --no-unal --no-mixed --non-deterministic --no-discordant -x $HISAT2Index -1 $INPUTDIR/$FASTQFILE1 -2 $INPUTDIR/$FASTQFILE2 -S ${SAMPLE}.sam
	     
	     echo ""
          echo "" 
##----------- sam2bam ---sort ---(to bam to save space) bam2bigiwg----


		samtools view -S ${SAMPLE}.sam -b > ${SAMPLE}.bam

          
          if ! [ -f ${SAMPLE}_sortedn.bam ]
		then
			echo "samtools sort -n -o ${SAMPLE}_sortedn.bam ${SAMPLE}.bam"
			samtools sort -n -o ${SAMPLE}_sortedn.bam ${SAMPLE}.bam

			echo ""
			echo ""
		fi		
		
# ---------------bam2bigiwg-------For visualization--------
          # bam sort
		samtools sort ${SAMPLE}.bam > ${SAMPLE}_sort.bam
		# bam index
		samtools index ${SAMPLE}_sort.bam 
		# bam 2 bigiwg
		bamCoverage -b ${SAMPLE}_sort.bam --binSize 30 --smoothLength 100 --normalizeUsing RPKM -p 10 -o ${SAMPLE}.bw

		rm ./${SAMPLE}.bam
          

## ----------- count -----------------
          # htseq-count -r name -s no -f bam ${SAMPLE}_sortedn.bam $GTFFILE  1>${SAMPLE}.count 2>${SAMPLE}.htseq-count.log
          htseq-count -r name -s no -f bam ${SAMPLE}_sortedn.bam $GTFFILE  1>${SAMPLE}.count 


done

mkdir $OUTPUTDIR/bigwig_files
mkdir $OUTPUTDIR/htseq_count_files


cd $OUTPUTDIR
mv */*.bw bigwig_files
mv */*.count htseq_count_files


