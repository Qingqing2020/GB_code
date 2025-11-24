#!/bin/bash


#################   Start   ###################################################

stime1=`date +"%Y-%m-%d %H:%M:%S"` 

# mkdir ./CEBPA_hicut/results

OUTPUTDIR=./CEBPA_hicut/results
# INPUTDIR=$1
INPUTDIR=./CEBPA_hicut/ANNO_XS01KF2022110123_PM-XS01KF2022110123-307/
# SampleInfo=$2
SampleInfo=./CEBPA_hicut/SampleNameTransform.txt

## !!!  keep in mind INPUTDIR changed after trimming 
# INPUTDIR=./CEBPA_hicut/data/trim_adapter

####################################################################
###          (1)  ------  trim adapter
####################################################################
INPUTDIR=./CEBPA_hicut/ANNO_XS01KF2022110123_PM-XS01KF2022110123-307/
cd $INPUTDIR
## 1.1-----pre_trim fastqc

mkdir -p $OUTPUTDIR/fastqc
fastqc $INPUTDIR/Rawdata/*/*.fq.gz -o $OUTPUTDIR/fastqc
multiqc $OUTPUTDIR/fastqc/*fastqc.zip --ignore *.html -o $OUTPUTDIR/fastqc/

## 1.2-----trim and after_trim fastqc


mkdir $INPUTDIR/trim_adapter
for i in $(ls $INPUTDIR/Rawdata/*|cat | grep "_R1.fq.gz")
do  	
    FASTQFILE1=$i
    FASTQFILE2=${i/_R1.fq.gz/_R2.fq.gz}
        # SAMPLE=${i%%_*}
    echo $FASTQFILE1
    echo $FASTQFILE2
    echo ""
    echo ""

    echo " trim_galore cut adapters started at $(date)"

    echo "trim_galore -q 20 --phred33 --stringency 3 --fastqc --length 30 -e 0.1 --paired $INPUTDIR/$FASTQFILE1 $INPUTDIR/$FASTQFILE2 --gzip -o $INPUTDIR/trim_adapter "

    trim_galore -q 20 --phred33 --stringency 3 --fastqc --length 30 -e 0.1 \
                --paired $INPUTDIR/Rawdata/*/$FASTQFILE1 $INPUTDIR/Rawdata/*/$FASTQFILE2  \
                --gzip -o $INPUTDIR/trim_adapter

    echo " trim_galore cut adapters finished at $(date)"

done

wait

# now still in $INPUTDIR
# mv the fastqc files
mkdir ../results/fastqc_trim
mv trim_adapter/*fastqc* ../results/fastqc_trim
# multiqc 
cd ../results/fastqc_trim
multiqc .


################################################################################### 
# # !!!  INPUTDIR changed !!!
INPUTDIR=./CEBPA_hicut/ANNO_XS01KF2022110123_PM-XS01KF2022110123-307/trim_adapter
# cd $INPUTDIR

NewOUTPUTDIR=./CEBPA_hicut/results_macs2_FE4_q5_nomodel


# ###################################################################################
# ###     (2)  ------ mapping and callpeak separately by individual sample
# ###################################################################################


mkdir $NewOUTPUTDIR/macs2_no_SPMR

for i in $(ls $INPUTDIR | grep "_R1_val_1.fq.gz")
do  

    ## 2.1 ----  rename Sample
	echo ${i}
	while read line
    do 
		# echo ${line}
        Pre_sample="$(echo ${line} | cut -d' ' -f1)"

        if [[ ${Pre_sample} = ${i%_val_1.fq.gz} ]]
        then
			SAMPLE="$(echo ${line} | cut -d' ' -f2)"
            echo "${Pre_sample} is ${SAMPLE}"
        fi

    done < $SampleInfo

    echo "The running sample is ${SAMPLE}"

        mkdir $NewOUTPUTDIR/$SAMPLE

        ##  ------- mapping using bowtie2------
        echo "bowtie2 --mm -p 8 --no-unal --no-mixed --non-deterministic --no-discordant -x $BOWTIEINDEXS -1 $INPUTDIR/$FASTQFILE1 -2 $INPUTDIR/$FASTQFILE2 -S $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam"

        bowtie2 --mm -p 8 --no-unal --no-mixed --non-deterministic --no-discordant -x $BOWTIEINDEXS -1 $INPUTDIR/$FASTQFILE1 -2 $INPUTDIR/$FASTQFILE2 -S $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam

        wait

        echo "samtools view -h -q 10 -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam > $OUTPUTDIR/$SAMPLE/${SAMPLE}.bam"

        samtools view -h -q 10 -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam > $OUTPUTDIR/$SAMPLE/${SAMPLE}.bam
        

        samtools flagstat  $OUTPUTDIR/$SAMPLE/${SAMPLE}.bam > $OUTPUTDIR/bowtie2_stat/${SAMPLE}.bowtie2.stat
        
        # remove sam file because it take too much space -----
        rm $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam

        #  !!!!!	remove blacklist from bam/bed file  !!!!!!
        bedtools intersect -v -a $OUTPUTDIR/$SAMPLE/${SAMPLE}.bam -b $Blacklist_Bed > $OUTPUTDIR/$SAMPLE/${SAMPLE}_rmbl.bam

        # bam--> sorted.bam ----
        echo "samtools sort $OUTPUTDIR/$SAMPLE/${SAMPLE}_rmbl.bam -o $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted.bam"

        samtools sort $OUTPUTDIR/$SAMPLE/${SAMPLE}_rmbl.bam -o $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted.bam
        
        # index  sorted.bam ------
        echo "samtools index $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted.bam"

        samtools index $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted.bam

        wait

        ## -------------- remove duplicates ----------------------
        echo "java -jar $PICARD MarkDuplicates I=$OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted.bam O=$OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam ASSUME_SORTED=TRUE M=$OUTPUTDIR/$SAMPLE/${SAMPLE}_marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true QUIET=true"

        java -jar $PICARD MarkDuplicates I=$OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted.bam O=$OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam ASSUME_SORTED=TRUE M=$OUTPUTDIR/$SAMPLE/${SAMPLE}_marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true QUIET=true


        ReadAfterRmdup=`samtools view -c $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam`
        FragmentAfterRmdup=`expr $ReadAfterRmdup / 2`
        echo " The Fragment count after rmdup is: $FragmentAfterRmdup"
        wait
        
        ## -------------- insert size distribution (after rmdup)----------------------
        echo "java -Xmx8g -jar $PICARD CollectInsertSizeMetrics I=$OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam O=$OUTPUTDIR/picard_InsertSize_out/${SAMPLE}_picard_insert_size_metrics.txt H=$OUTPUTDIR/picard_InsertSize_out/${SAMPLE}_picard_insert_size_histogram.pdf M=0.5"

        java -Xmx8g -jar $PICARD CollectInsertSizeMetrics I=$OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam O=$OUTPUTDIR/picard_InsertSize_out/${SAMPLE}_picard_insert_size_metrics.txt H=$OUTPUTDIR/picard_InsertSize_out/${SAMPLE}_picard_insert_size_histogram.pdf M=0.5


    ## -------------call peak using macs2------------
    
      
    echo "macs2 callpeak -t $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam -f BAMPE -g hs -n ${SAMPLE} -B -q 0.00001 --fe-cutoff 4 --nomodel --extsize 200 --shift 100 --outdir $NewOUTPUTDIR/macs2_no_SPMR"

    macs2 callpeak -t $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam -f BAMPE -g hs -n ${SAMPLE} -B -q 0.00001 --fe-cutoff 4 --nomodel --extsize 200 --shift 100 --outdir $NewOUTPUTDIR/macs2_no_SPMR



    echo "samtools sort -n $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam -o $OUTPUTDIR/$SAMPLE/${SAMPLE}_sortedbyn_rmdup.bam"

    samtools sort -n $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam -o $OUTPUTDIR/$SAMPLE/${SAMPLE}_sortedbyn_rmdup.bam

    # combine pair-end  bam to bed file
    echo ""

    bedtools bamtobed -i $OUTPUTDIR/$SAMPLE/${SAMPLE}_sortedbyn_rmdup.bam -bedpe -mate1 2>/dev/null | awk '{if ($2<$5) $11=$2; else $11=$5; if ($3>$6) $12=$3; else $12=$6; print $1"\t"$11"\t"$12"\t"$9}' > $OUTPUTDIR/$SAMPLE/${SAMPLE}_PE.bed 

wait

done


###################################################################################
###     (3)  ------ pool replicates(bam) and call peak based on condition
###################################################################################


NewOUTPUTDIR=./CEBPA_hicut/results_macs2_FE4_q5_nomodel
mkdir $NewOUTPUTDIR/merged_macs2_no_SPMR


cd $OUTPUTDIR

for i in $(ls |grep -E "_[0-9]" | cut -d '_' -f 1,2 | uniq)
do 
    SAMPLE=$i
    echo $SAMPLE
 ## --------------merge *sorted_rmdup.bam duplicates ----------------------
    echo "ls ${SAMPLE}*/*sorted_rmdup.bam|xargs samtools merge $OUTPUTDIR/merged_macs2_no_SPMR/${SAMPLE}_sorted_rmdup_merged.bam"
    ls ${SAMPLE}*/*sorted_rmdup.bam|xargs samtools merge $OUTPUTDIR/merged_macs2_no_SPMR/${SAMPLE}_sorted_rmdup_merged.bam
    
    wait
    ## -------------call peak using macs2------------


    echo "macs2 callpeak -t $OUTPUTDIR/merged_macs2_no_SPMR/${SAMPLE}_sorted_rmdup_merged.bam -f BAMPE -g hs -n ${SAMPLE}_merged -B -q 0.00001 --fe-cutoff 4 --nomodel --extsize 200 --shift 100  --outdir $NewOUTPUTDIR/merged_macs2_no_SPMR"

    macs2 callpeak -t $OUTPUTDIR/merged_macs2_no_SPMR/${SAMPLE}_sorted_rmdup_merged.bam -f BAMPE -g hs -n ${SAMPLE}_merged -B -q 0.00001 --fe-cutoff 4 --nomodel --extsize 200 --shift 100  --outdir $NewOUTPUTDIR/merged_macs2_no_SPMR

wait

done

###################################################################################
###   (4) ------ bdg2bw_no_SPMR 

###################################################################################

## need to run separately
conda activate local

NewOUTPUTDIR=./CEBPA_hicut/results_macs2_FE4_q5_nomodel

mkdir $NewOUTPUTDIR/bwFiles_no_SPMR

cd $NewOUTPUTDIR


for i in $(ls *macs2_no_SPMR/*pileup.bdg)
do  
    # get the file and sample ---- no need to do this step!!!!
    echo "the file is $i"
    tmp=$(basename $i) ## pay attention use ()
    SAMPLE=${tmp%%_treat_pileup.bdg}
    echo "the sample is $SAMPLE"

    # bdg to bw
    # http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
    echo " bash ~/PhD/LSD1/bdg2bw.sh $i ./hg19.chrom.sizes"
    bash ~/PhD/LSD1/bdg2bw.sh $i ./hg19.chrom.sizes

    wait
done

mv *macs2_no_SPMR/*.bw bwFiles_no_SPMR/

conda deactivate

###################################################################################

### : run  step5  after the above finished
###################################################################################


###  (5) merge pooled peaks and calculate coverage in each sample 
###################################################################################

OUTPUTDIR=./CEBPA_hicut/results
NewOUTPUTDIR=./CEBPA_hicut/results_macs2_FE4_q5_nomodel

# 5.1 --------merge all the **pooled** bed files-------------
# $0 is the whole line of arguments
# OFS: Output Field Separator

cd $NewOUTPUTDIR/merged_macs2_no_SPMR
cat *Peak | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN{OFS="\t"} $0=$0"\t"NR' > CEBPA_merged_pooled_Peaks.bed

# 5.2 -------- calculate the coverage of each sample -------------


Samples="CEBPA_0h_rep2 CEBPA_0h_rep3 CEBPA_0h_rep4 CEBPA_4d_rep2 CEBPA_4d_rep3 CEBPA_4d_rep4 CEBPA_4h_rep2 CEBPA_4h_rep3 CEBPA_4h_rep4"

cd  $NewOUTPUTDIR/macs2_no_SPMR
for i in $Samples; do bedtools coverage -a $NewOUTPUTDIR/merged_macs2_no_SPMR/CEBPA_merged_pooled_Peaks.bed -b $OUTPUTDIR/${i}/${i}_PE.bed > ${i}_coverage.bed; done

# 5.3 ---------------combine the coverage of all the samples together------------
# where does the column # of the results depends also on your input file,
# If you have 3 column in the coverage ( A ), then the f4 ,f6 is the column you want, 
# If you have 4 column in the coverage ( A ), then the f5 ,f7 is the column you want, 


for i in *coverage.bed; do cat ${i} | cut -f5 > ${i}.coverage;done

# just add the peak length to the first column
i=CEBPA_0h_rep2; paste <(cat ${i}_coverage.bed | cut -f7) *coverage > CEBPA_coverage.bed

# 5.4 ---------------add header (column names ) to the file------------


sed -i '1i PeakLength\tCEBPA_0h_rep2\tCEBPA_0h_rep3\tCEBPA_0h_rep4\tCEBPA_4d_rep2\tCEBPA_4d_rep3\tCEBPA_4d_rep4\tCEBPA_4h_rep2\tCEBPA_4h_rep3\tCEBPA_4h_rep4' CEBPA_coverage.bed




# ###########   End   ########################################################################

# stime2=`date +"%Y-%m-%d %H:%M:%S"`
# echo "[$stime1] started"
# echo "[$stime2] done"

##  distiller 
nohup nextflow run distiller-nf -params-file CEBPA.yml > CEBPA.log