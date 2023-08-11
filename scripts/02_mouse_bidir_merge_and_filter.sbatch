#!/bin/bash

#SBATCH --job-name=mouse_merge_and_filter
#SBATCH --output=/scratch/Users/rusi2317/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/rusi2317/e_and_o/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1000MB
#SBATCH --partition short
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rutendo.sigauke@colorado.edu

module purge
module load bedtools/2.28.0

##############################
##0: initialize files       ## 
##############################
wd=/scratch/Shares/dowell/dbnascent/out/meta_analysis/mumerge
bidir=/scratch/Shares/dowell/dbnascent/out/meta_analysis/mumerge/processed_data/all_master_mumerge

#path to scripts
samples_support_bidir=${wd}/bin/samples_support_bidir.R
dreg_tfit_overlaps=${wd}/bin/dreg_tfit_overlaps.sh
merge_tfit_dreg_bidir=${wd}/bin/merge_tfit_dreg_bidir.R
bidir_gene_overlaps=${wd}/bin/bidir_gene_overlaps.sh
gene_transcript_distances=${wd}/bin/gene_transcript_distances.sh
remove_bidirs_overlap_posneg_genes=${wd}/bin/remove_bidirs_overlap_posneg_genes.R
conv_genes=${wd}/bin/conv_genes.R
remove_bidirs_overlap_conv_gene=${wd}/bin/remove_bidirs_overlap_conv_gene.R

#############################################
#1: Get number of samples supporting a call##
#############################################
#get Paper level muMerged calls
tfit=${bidir}/All_tiers_tfit_mm10_MUMERGE.bed
dreg=${bidir}/All_tiers_dreg_mm10_MUMERGE.bed
out_bidir=${wd}/processed_data/numSamples_bidirs
mkdir -p ${out_bidir}

Rscript --vanilla ${samples_support_bidir} -b ${tfit} -o ${out_bidir}/
Rscript --vanilla ${samples_support_bidir} -b ${dreg} -o ${out_bidir}/

#############################################
#2: get overlaps                           ##
#############################################
#input files from samples_support_bidir.R
tfit_num=${out_bidir}/All_tiers_tfit_mm10_MUMERGE_numSamples.bed
dreg_num=${out_bidir}/All_tiers_dreg_mm10_MUMERGE_numSamples.bed
outdir_overlaps=${wd}/processed_data/overlaps

mkdir -p ${outdir_overlaps}

sh ${dreg_tfit_overlaps} -d ${dreg_num} -t ${tfit_num} -g 'mm10' -o ${outdir_overlaps}

#############################################
#3: Filter the mumerge files and save final##
#############################################
#input from dreg_tfit_overlaps.sh
overlaps=${outdir_overlaps}/mm10_overlap.bed
tfit_only=${outdir_overlaps}/mm10_tfit_only.bed
dreg_only=${outdir_overlaps}/mm10_dreg_only.bed
metadata=/scratch/Shares/dowell/dbnascent/out/meta_analysis/mumerge/sample_metadata/230131_db_output.tsv
metadata_qc=/scratch/Shares/dowell/dbnascent/out/meta_analysis/mumerge/sample_metadata/20211026_paper_qc_scores.csv
base_content=/scratch/Shares/dowell/dbnascent/out/meta_analysis/mumerge/processed_data/paper_mumerges/all_calls
genome="mouse"
qc=2
outdir_master=${wd}/processed_data/master_mumerge

mkdir -p ${outdir_master}

Rscript --vanilla ${merge_tfit_dreg_bidir} -b ${overlaps} -t ${tfit_only} -d ${dreg_only} -m ${metadata} -n ${metadata_qc} -c ${base_content} -g ${genome} -q ${qc} -o ${outdir_master}/

############################################
#4: Get the overlaps between bidirectionals# 
#   and genes using bed tools              #
############################################
#input for bidir_gene_overlaps.sh
genes=/scratch/Shares/dowell/genomes/mm10/ncbi/mm10_refseq_transcripts.bed #mm10_refseq_genes.bed
tfit_filt=${outdir_master}/mm10_master_mumerge_qc_gc_len_filtered.bed

sh ${bidir_gene_overlaps} -a ${genes} -b ${tfit_filt} -g 'mm10' -o ${outdir_master}

############################################
#5: Get distances between the X closest    #
############################################
#takes gene transcript annotation bed file
outdir_genes=${wd}/processed_data/gene_distances

mkdir -p ${outdir_genes}

sh ${gene_transcript_distances} -a ${genes} -g 'mm10' -o ${outdir_genes}

############################################
#6: Remove bidirectionals overlapping      #
#   positive and negative genes            #
############################################

bidirs_cleaned=${outdir_master}/mm10_master_mumerge_qc_gc_len_filtered_merged.bed
bidir_gene_overlaps=${outdir_master}/mm10_master_mumerge_qc_gc_len_filtered_merged_geneoverlaps.bed

Rscript --vanilla ${remove_bidirs_overlap_posneg_genes} -a ${bidirs_cleaned} -b ${bidir_gene_overlaps} -g ${genome} -o ${outdir_master}/

###########################################
#7: Get converging genes by 1kb           #
###########################################
all_close_genes=${outdir_genes}/mm10_refseq_transcripts_closest_pos_to_neg.bed 

Rscript --vanilla ${conv_genes} -c ${all_close_genes} -g ${genome} -o ${outdir_master}/

###########################################
#8: Remove bidirectionals overlapping     #
#   converging genes                      #
###########################################
bidirs_no_negpos=${outdir_master}/mm10_master_qc_gc_len_posneg_gene_filt.bed
close_genes=${outdir_master}/mm10_refseq_transcripts_converging_1kb.bed

#get intersect between overlappaing genes. These will be removed from the final bed file with bidirectionals   
bedtools intersect -f 0.5 -wo -a ${bidirs_no_negpos} -b ${close_genes} > ${outdir_master}/mm10_bidir_overlapping_conv_genes.bed

Rscript --vanilla ${remove_bidirs_overlap_conv_gene} -a ${bidirs_no_negpos} -c ${outdir_master}/mm10_bidir_overlapping_conv_genes.bed -g ${genome} -o ${outdir_master}/

###########################################                   
#9: Clean temporary files                 #                   
###########################################                   

rm ${outdir_master}/mm10_bidir_overlapping_conv_genes.bed
rm ${outdir_master}/mm10_master_mumerge_qc_gc_len_filtered.bed
rm ${outdir_master}/mm10_master_mumerge_qc_gc_len_filtered_merged.bed
rm ${outdir_master}/mm10_master_mumerge_qc_gc_len_filtered_merged_geneoverlaps.bed
rm ${outdir_master}/mm10_master_qc_gc_len_posneg_gene_filt.bed
rm ${outdir_master}/mm10_refseq_transcripts_converging_1kb.bed

rm -r ${out_bidir}
rm -r ${outdir_overlaps}
rm -r ${outdir_genes}
