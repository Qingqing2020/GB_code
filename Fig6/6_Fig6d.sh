#!/bin/bash

#############################################################################
#### profile 
#############################################################################

Peaks_Regions=./CEBPA_hicut/E_P
# diff peaks
# regionsLabel="Stem_specific Term_specific Nochange"
regionsLabel="RA0h_E RA4d_E P"

region="$Peaks_Regions/HL60_Ctrl_Enhancer.bed $Peaks_Regions/HL60_RA4d_Enhancer.bed $Peaks_Regions/HL60_gene_promoters.bed"

samplesLabel="0h 4h 4d"

bw_path=./CEBPA_hicut/results_macs2_FE4_q5_nomodel/bwFiles_no_SPMR/CEBPA_bwNorm_RIP
bws="$bw_path/CEBPA_0h_merged_treat_pileup_scaled.bw $bw_path/CEBPA_4h_merged_treat_pileup_scaled.bw $bw_path/CEBPA_4d_merged_treat_pileup_scaled.bw" 

name="CEBPA_on_ctrl_RA4d_E_P"


mkdir ./CEBPA_hicut/E_P/CEBPA_Profile_Use_18Paper_E_AllP

CurrDir=./CEBPA_hicut/E_P/CEBPA_Profile_Use_18Paper_E_AllP

mkdir -p $CurrDir/matrix

/opt/miniconda2/bin/computeMatrix reference-point --referencePoint center -b 2000 -a 2000 -R $region -S $bws --binSize 20 -p 8 --numberOfProcessors=max --skipZeros -o $CurrDir/matrix/$name.mat.gz

mkdir -p $CurrDir/png


/opt/miniconda2/bin/plotHeatmap -m $CurrDir/matrix/$name.mat.gz --colorMap Blues  --heatmapHeight 20 --sortUsingSamples 1 --legendLocation best -o $CurrDir/png/$name.png --samplesLabel $samplesLabel --regionsLabel $regionsLabel


## separately
/opt/miniconda2/bin/plotProfile -m $CurrDir/matrix/$name.mat.gz --plotHeight 6 --plotWidth 5 -o $CurrDir/png/${name}_profile.png   --legendLocation best --samplesLabel $samplesLabel --regionsLabel $regionsLabel

## together
/opt/miniconda2/bin/plotProfile -m $CurrDir/matrix/$name.mat.gz --plotHeight 12 --plotWidth 10 -o $CurrDir/png/${name}_size_bigger_profile.png   --legendLocation best --samplesLabel $samplesLabel --regionsLabel $regionsLabel --perGroup
