#!/bin/bash

#################   Start   ###################################################

OUTPUTDIR=./Neu/ATAC/results
INPUTDIR=./Neu/ATAC/data
SampleInfo=./Neu/ATAC/SampleNameTransform.txt

####################################################################
###          (1)  ------  trim adapter
####################################################################

cd $INPUTDIR
## 1.1-----pre_trim fastqc

mkdir -p $OUTPUTDIR/fastqc
fastqc $INPUTDIR/*.fq.gz -o $OUTPUTDIR/fastqc
multiqc $OUTPUTDIR/fastqc/*fastqc.zip --ignore *.html -o $OUTPUTDIR/fastqc/

## 1.2-----trim and after_trim fastqc

mkdir $INPUTDIR/trim_adapter
for i in $(ls $INPUTDIR|cat | grep "_1.fq.gz")
do  	
    FASTQFILE1=$i
    FASTQFILE2=${i/_1.fq.gz/_2.fq.gz}
        # SAMPLE=${i%%_*}
    echo $FASTQFILE1
    echo $FASTQFILE2
    echo ""
    echo ""

    echo " trim_galore cut adapters started at $(date)"

    echo "trim_galore -q 20 --phred33 --stringency 3 --fastqc --length 30 -e 0.1 --paired $INPUTDIR/$FASTQFILE1 $INPUTDIR/$FASTQFILE2 --gzip -o $INPUTDIR/trim_adapter "

    trim_galore -q 20 --phred33 --stringency 3 --fastqc --length 30 -e 0.1 \
                --paired $INPUTDIR/$FASTQFILE1 $INPUTDIR/$FASTQFILE2  \
                --gzip -o $INPUTDIR/trim_adapter

    echo " trim_galore cut adapters finished at $(date)"

done

wait

mkdir ../results/fastqc_trim
mv trim_adapter/*fastqc* ../results/fastqc_trim
# multiqc 
cd ../results/fastqc_trim
multiqc .

################################################################################### 
# !!!  INPUTDIR changed !!!
INPUTDIR=./Neu/ATAC/data/trim_adapter
cd $INPUTDIR

OUTPUTDIR=./Neu/ATAC/results

PICARD=/home/qchen/Software/picard.jar

##  -- human hg19 BOWTIEINDEXS
BOWTIEINDEXS=/home/shared_data/Annotation/UCSC/Human_Genome/HG19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

Blacklist_Bed=./Blacklist/hg19-blacklist.v2.bed

###################################################################################
###     (2)  ------ mapping and callpeak separately by individual sample
###################################################################################

mkdir $OUTPUTDIR/macs2_no_SPMR
mkdir $OUTPUTDIR/bowtie2_stat
mkdir $OUTPUTDIR/picard_InsertSize_out

for i in $(ls $INPUTDIR | grep "_1_val_1.fq.gz")
do  

    ## 2.1 ----  rename Sample
	echo ${i}
	while read line
    do 
		# echo ${line}
        Pre_sample="$(echo ${line} | cut -d' ' -f1)"
        
        # echo ${Sra}
        # echo ${SAMPLE}
        if [[ ${Pre_sample} = ${i%_1_val_1.fq.gz} ]]
        then
			SAMPLE="$(echo ${line} | cut -d' ' -f2)"
            echo "${Pre_sample} is ${SAMPLE}"
        fi

    done < $SampleInfo

    echo "The running sample is ${SAMPLE}"
	echo ""
	
    FASTQFILE1=$i
    FASTQFILE2=${i/_1_val_1.fq.gz/_2_val_2.fq.gz}
    echo $FASTQFILE1
    echo $FASTQFILE2
    echo ""
    echo ""
    
    mkdir $OUTPUTDIR/$SAMPLE

    ## 2.2 ------- mapping using bowtie2------
    echo "bowtie2 --mm -p 8 --no-unal --no-mixed --non-deterministic --no-discordant -x $BOWTIEINDEXS -1 $INPUTDIR/$FASTQFILE1 -2 $INPUTDIR/$FASTQFILE2 -S $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam"

    bowtie2 --mm -p 8 --no-unal --no-mixed --non-deterministic --no-discordant -x $BOWTIEINDEXS -1 $INPUTDIR/$FASTQFILE1 -2 $INPUTDIR/$FASTQFILE2 -S $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam

    wait

    ## sam2bam ----

    echo "samtools view -h -q 10 -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam > $OUTPUTDIR/$SAMPLE/${SAMPLE}.bam"

    samtools view -h -q 10 -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 $OUTPUTDIR/$SAMPLE/${SAMPLE}.sam > $OUTPUTDIR/$SAMPLE/${SAMPLE}.bam
    
    
    ### get the mapping statistics
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

    wait
    
    ## -------------- insert size distribution (after rmdup)----------------------
    echo "java -Xmx8g -jar $PICARD CollectInsertSizeMetrics I=$OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam O=$OUTPUTDIR/picard_InsertSize_out/${SAMPLE}_picard_insert_size_metrics.txt H=$OUTPUTDIR/picard_InsertSize_out/${SAMPLE}_picard_insert_size_histogram.pdf M=0.5"

    java -Xmx8g -jar $PICARD CollectInsertSizeMetrics I=$OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam O=$OUTPUTDIR/picard_InsertSize_out/${SAMPLE}_picard_insert_size_metrics.txt H=$OUTPUTDIR/picard_InsertSize_out/${SAMPLE}_picard_insert_size_histogram.pdf M=0.5


    ## -------------call peak using macs2------------
    

    echo "macs2 callpeak -t $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam -f BAMPE -g hs -n ${SAMPLE} -B -q 0.05 --outdir $OUTPUTDIR/macs2_no_SPMR"

    macs2 callpeak -t $OUTPUTDIR/$SAMPLE/${SAMPLE}_sorted_rmdup.bam -f BAMPE -g hs -n ${SAMPLE} -B -q 0.05 --outdir $OUTPUTDIR/macs2_no_SPMR


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


OUTPUTDIR=./Neu/ATAC/results
mkdir $OUTPUTDIR/merged_macs2_no_SPMR

cd $OUTPUTDIR

for i in $(ls |grep -E "_[1-9]" | rev | cut -c3- |rev | uniq|grep -v "HL60_Tn5"| sort) 
do 
    SAMPLE=$i
    echo $SAMPLE
 ## --------------merge *sorted_rmdup.bam duplicates ----------------------
    echo "ls ${SAMPLE}*/*sorted_rmdup.bam|xargs samtools merge $OUTPUTDIR/merged_macs2_no_SPMR/${SAMPLE}_sorted_rmdup_merged.bam"
    ls ${SAMPLE}*/*sorted_rmdup.bam|xargs samtools merge $OUTPUTDIR/merged_macs2_no_SPMR/${SAMPLE}_sorted_rmdup_merged.bam
    
    wait

    echo "macs2 callpeak -t $OUTPUTDIR/merged_macs2_no_SPMR/${SAMPLE}_sorted_rmdup_merged.bam -f BAMPE -g hs -n ${SAMPLE}_merged -B -q 0.05 --outdir $OUTPUTDIR/merged_macs2_no_SPMR"

    macs2 callpeak -t $OUTPUTDIR/merged_macs2_no_SPMR/${SAMPLE}_sorted_rmdup_merged.bam -f BAMPE -g hs -n ${SAMPLE}_merged -B -q 0.05 --outdir $OUTPUTDIR/merged_macs2_no_SPMR

wait

done

###################################################################################
###   (4) ------ bdg2bw_no_SPMR (should all be integer in  the bigwig file-- can check use R )
### need to be in conda activated local environment
###################################################################################
# save the bw file ( not normalized ) for future RIP normalization !!

OUTPUTDIR=./Neu/ATAC/results

mkdir $OUTPUTDIR/bwFiles_no_SPMR

cd $OUTPUTDIR

stime1=`date +"%Y-%m-%d %H:%M:%S"`
stime1=`date +"%Y-%m-%d %H:%M:%S"`


for i in $(ls *macs2_no_SPMR/*pileup.bdg)
do  
    # get the file and sample ---- no need to do this step!!!!
    echo "the file is $i"
    tmp=$(basename $i) ## pay attention use ()
    SAMPLE=${tmp%%_treat_pileup.bdg}
    echo "the sample is $SAMPLE"

    # bdg to bw
    # http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
    echo " bash ~/PhD/LSD1/bdg2bw.sh $i ./Neu/hg19.chrom.sizes"
    bash ~/PhD/LSD1/bdg2bw.sh $i ./Neu/hg19.chrom.sizes

    wait
done

mv *macs2_no_SPMR/*.bw bwFiles_no_SPMR/

###################################################################################

mkdir $OUTPUTDIR/merged_macs2_no_SPMR/HL60
mv $OUTPUTDIR/merged_macs2_no_SPMR/HL60_* $OUTPUTDIR/merged_macs2_no_SPMR/HL60
mkdir $OUTPUTDIR/macs2_no_SPMR/HL60
mv $OUTPUTDIR/macs2_no_SPMR/HL60_* $OUTPUTDIR/macs2_no_SPMR/HL60
#######################################################################


###   (5) HL60 ------ merge pooled peaks and calculate coverage in each sample 
###################################################################################


OUTPUTDIR=./Neu/ATAC/results/HL60

# 5.1 --------merge all the **pooled** bed files-------------

cd $OUTPUTDIR/merged_macs2_no_SPMR/HL60
cat *Peak | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN{OFS="\t"} $0=$0"\t"NR' > HL60_merged_pooled_Peaks.bed

# 5.2 -------- calculate the coverage of each sample -------------

Samples="HL60_2d_1 HL60_2d_2 HL60_4d_1 HL60_4d_2 HL60_4h_1 HL60_4h_2 HL60_Ctrl_1 HL60_Ctrl_2"

cd  $OUTPUTDIR/macs2_no_SPMR/HL60
for i in $Samples; do bedtools coverage -a $OUTPUTDIR/merged_macs2_no_SPMR/HL60/HL60_merged_pooled_Peaks.bed -b $OUTPUTDIR/${i}/${i}_PE.bed > ${i}_coverage.bed; done

# 5.3 ---------------combine the coverage of all the samples together------------
# where does the column # of the results depends also on your input file,
# If you have 3 column in the coverage ( A ), then the f4 ,f6 is the column you want, 
# If you have 4 column in the coverage ( A ), then the f5 ,f7 is the column you want, 

# here the  merged_Peaks.bed has 4 columns with the last column as the peak number

for i in *coverage.bed; do cat ${i} | cut -f5 > ${i}.coverage;done

# just add the peak length to the first column
i=HL60_2d_1; paste <(cat ${i}_coverage.bed | cut -f7) *coverage > HL60_coverage.bed

# 5.4 ---------------add header (column names ) to the file------------

## add by myself, then no need to worry about the column names in the R code
sed -i '1i PeakLength\tHL60_2d_1\tHL60_2d_2\tHL60_4d_1\tHL60_4d_2\tHL60_4h_1\tHL60_4h_2\tHL60_Ctrl_1\tHL60_Ctrl_2' HL60_coverage.bed



###########   End   ########################################################################
