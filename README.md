# bidirectionals_merged

Generates a final merged file from bidirectionals called by Tfit and dREG across multiple samples after paper muMerge.

## Repository stucture

- bin: scripts used in the merging process

- notebooks: jupyter notebooks with analyses used to determine defaults in scripts

## Dependencies

The process was performed in bash and R. Below are the tools and libraries used.

`mumerge version 1.1.0` (https://github.com/Dowell-Lab/mumerge) 

`bedtools version 2.28.0`

`R version 3.6.0`

`data.table version 1.14.2`

## Summary of scripts

This is a multi-step process with starts with all muMerge files from Tfit and dREG. These files are from experiment muMerge (combined by experimental setup), then experiments are muMerged by cell/tissue types.

1. `samples_support_bidir.R` : Get number of samples supporting a bidirectional call.

2. `dreg_tfit_overlaps.sh` : Overlap between dREG and Tfit muMerge calls.

3. `merge_tfit_dreg_bidir.R` : Combine Tfit and dREG calls taking into account the number of sample supporting a bidirectional transcript. This script also filters single calls found in a low quality paper and calls from papers with low GC content. Additionally, regions are filtered by length (a minimum and maximum length filter is applied).

4. `bidir_gene_overlaps.sh` : Sort and merge bidirectional regions. Then overlap the bidirectionals with genes.

5. `gene_transcript_distances.sh` : Get distances of the 2 closest transcript pairs.

6. `conv_genes.R` : Get converging (**--> <--**) gene transcripts from a genome annotations at a minimum distance of 1kb.

7. `remove_bidirs_overlap_posneg_genes.R` : Remove bidirectionals overlapping positive and negative genes.

8. `remove_bidirs_overlap_conv_gene.R` : Remove bidirectionals overlapping converging gene transcripts.

## Example final output

The final output is a bed6 file for bidirectionals with the following columns: `chromosome` ,`start` ,`stop` ,`bidirectional` ,`number of papers a bidirectional was called`, `strand (it is . since bidirectionals are not stranded)`.

```
chr1	3917	4919	dreg	14	.
chr1	6132	6486	dreg	7	.
chr1	7207	7841	dreg	7	.
chr1	13264	13506	tfit	14	.
chr1	16191	16429	tfit	128	.
```
## muMerge Pipeline

![muMerge Pipeline](https://github.com/Dowell-Lab/bidirectionals_merged/blob/main/README_images/muMerge_strategy_version3.png)