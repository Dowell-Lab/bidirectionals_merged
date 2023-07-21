#!/bin/bash

#SBATCH --job-name=test_merge
#SBATCH --output=/scratch/Users/rusi2317/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/rusi2317/e_and_o/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50MB
#SBATCH --partition short
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rutendo.sigauke@colorado.edu

module purge
module load bedtools

##############################
##0: initialize files       ## 
##############################
wd=/scratch/Users/rusi2317/projects/meta_analysis_qc/bidirectionals_merged
bidir=/scratch/Shares/dowell/dbnascent/out/meta_analysis/mumerge/processed_data/all_master_mumerge

#get Paper level muMerged calls
tfit=${bidir}/chr/All_tiers_tfit_hg38_MUMERGE.bed
dreg=${bidir}/All_tiers_dreg_hg38_MUMERGE.bed
out_bidir=${wd}/processed/numSamples_bidirs
mkdir -p ${out_bidir}

#input files from samples_support_bidir.R
tfit_hg38=${out_bidir}/All_tiers_tfit_hg38_MUMERGE_numSamples.bed
dreg_hg38=${out_bidir}/All_tiers_dreg_hg38_MUMERGE_numSamples.bed
#tfit_hg38=${bidir}/All_tiers_tfit_hg38_MUMERGE_numSamples.bed
#dreg_hg38=${bidir}/All_tiers_dreg_hg38_MUMERGE_numSamples.bed
outdir=${wd}/processed/overlaps

mkdir -p ${outdir}

#path to scripts
samples_support_bidir=${wd}/bin/samples_support_bidir.R
dreg_tfit_overlaps=${wd}/bin/dreg_tfit_overlaps.sh

#############################################
#1: Get number of samples supporting a call##
#############################################

Rscript --vanilla ${samples_support_bidir} -b ${tfit} -o ${out_bidir}/
Rscript --vanilla ${samples_support_bidir} -b ${dreg} -o ${out_bidir}/

#############################################
#2: get overlaps                           ##
#############################################

sh ${dreg_tfit_overlaps} -d ${dreg_hg38} -t ${tfit_hg38} -g 'hg38' -o ${outdir}
