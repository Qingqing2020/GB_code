#!/bin/bash

#############################################################################
#### profile 
#############################################################################


Peaks_Regions=./ATAC/results/DiffPeakClusters2
# diff peaks
regionsLabel="HL60_Ctrl HL60_4h HL60_2d HL60_4d NoChange "

region="$Peaks_Regions/HL60_Ctrl_specific_peaks.bed $Peaks_Regions/HL60_4h_specific_peaks.bed $Peaks_Regions/HL60_2d_specific_peaks.bed $Peaks_Regions/HL60_4d_specific_peaks.bed $Peaks_Regions/NoChangePeak.bed "

samplesLabel="HL60_Ctrl HL60_4h HL60_2d HL60_4d"


bws="./ATAC/results/bwFiles_no_SPMR/HL60_bwNorm_RIP/HL60_Ctrl_merged_treat_pileup_scaled.bw ./ATAC/results/bwFiles_no_SPMR/HL60_bwNorm_RIP/HL60_4h_merged_treat_pileup_scaled.bw ./ATAC/results/bwFiles_no_SPMR/HL60_bwNorm_RIP/HL60_2d_merged_treat_pileup_scaled.bw ./ATAC/results/bwFiles_no_SPMR/HL60_bwNorm_RIP/HL60_4d_merged_treat_pileup_scaled.bw"


name="SpecificPeaks"

CurrDir=./ATAC/results/profile_HL60_NoChange

cd CurrDir
mkdir -p $CurrDir/matrix

computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R $region -S $bws --binSize 20 -p 8 --numberOfProcessors=max --skipZeros -o $CurrDir/matrix/$name.mat.gz

mkdir -p $CurrDir/png


plotHeatmap -m $CurrDir/matrix/$name.mat.gz --colorMap Blues  --heatmapHeight 20 --sortUsingSamples 1 --legendLocation best -o $CurrDir/png/$name.png --samplesLabel $samplesLabel --regionsLabel $regionsLabel
#--zMax 15 


plotProfile -m $CurrDir/matrix/$name.mat.gz --plotHeight 6 --plotWidth 5 -o $CurrDir/png/${name}_profile.png   --legendLocation best --samplesLabel $samplesLabel --regionsLabel $regionsLabel
