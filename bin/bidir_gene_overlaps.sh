while getopts a:b:g:o: flag
do
    case "${flag}" in
        a) genes=${OPTARG};;
        b) bidirectionals=${OPTARG};;
        g) genome=${OPTARG};;
        o) outdir=${OPTARG};;
    esac
done

echo "Gene Transcript Bed File: $genes";
echo "Bidirectional Transcript Bed File: $bidirectionals";
echo "Genome analyzied: $genome"
echo "Output Directory: $outdir";

##using bedtools to sort the output
##and merge overlapping or bookend calls
bedtools sort -i ${bidirectionals} \
    | bedtools merge -d -1 -c 4,5,6 -o distinct,max,distinct -delim "|" -i - \
    > ${outdir}/${genome}_master_mumerge_qc_gc_len_filtered_merged.bed

##get the overlaps with genes and the 
## number of bases overlapping
bedtools intersect -a ${genes} \
    -b ${outdir}/${genome}_master_mumerge_qc_gc_len_filtered_merged.bed \
    -wo > ${outdir}/${genome}_master_mumerge_qc_gc_len_filtered_merged_geneoverlaps.bed
