#!/bin/bash

#############################################################################
#### profile 
#############################################################################

Peaks_Regions=./ATAC_4h_DEG_table
# diff peaks
# regionsLabel="Stem_specific Term_specific Nochange"
regionsLabel="UpDEG_asso_ATAC DownDEG_asso_ATAC"

region="$Peaks_Regions/Up_DEG_associated_4h_ATAC.bed $Peaks_Regions/Down_DEG_associated_4h_ATAC.bed"


samplesLabel="RA0h RA4h RA2d RA4d"

bws="$PATH/HL60_Ctrl_merged_treat_pileup_scaled.bw $PATH/HL60_4h_merged_treat_pileup_scaled.bw $PATH/HL60_2d_merged_treat_pileup_scaled.bw $PATH/HL60_4d_merged_treat_pileup_scaled.bw"

name="DEG_asso_ATAC4h_profile"


mkdir ./DEG_asso_ATAC4h_profile

CurrDir=./DEG_asso_ATAC4h_profile

mkdir -p $CurrDir/matrix

/opt/miniconda2/bin/computeMatrix reference-point --referencePoint center -b 2000 -a 2000 -R $region -S $bws --binSize 20 -p 8 --numberOfProcessors=max --skipZeros -o $CurrDir/matrix/$name.mat.gz

mkdir -p $CurrDir/png


/opt/miniconda2/bin/plotHeatmap -m $CurrDir/matrix/$name.mat.gz --colorMap Blues  --heatmapHeight 20 --sortUsingSamples 1 --legendLocation best -o $CurrDir/png/$name.png --samplesLabel $samplesLabel --regionsLabel $regionsLabel


## separately
/opt/miniconda2/bin/plotProfile -m $CurrDir/matrix/$name.mat.gz --plotHeight 6 --plotWidth 5 -o $CurrDir/png/${name}_profile.png   --legendLocation best --samplesLabel $samplesLabel --regionsLabel $regionsLabel

## together
/opt/miniconda2/bin/plotProfile -m $CurrDir/matrix/$name.mat.gz --plotHeight 12 --plotWidth 10 -o $CurrDir/png/${name}_size_bigger_profile.png   --legendLocation best --samplesLabel $samplesLabel --regionsLabel $regionsLabel --perGroup
