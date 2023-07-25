while getopts a:g:o: flag
do
    case "${flag}" in
        a) annotations=${OPTARG};;
        g) genome=${OPTARG};;
        o) outdir=${OPTARG};;
    esac
done

echo "Gene Transcript Bed File: $annotations";
echo "Genome ID: $genome"
echo "Output Directory: $outdir";

echo "Get distances bewteen the closest positive and negative trasnscripts."

# get a bed6 file of gene transcripts and remove contig transcripts
# additionally, transcripts are spit into positive and negative strands
cat ${annotations} | grep -v 'chr._\|chr.._' | grep '+' > ${outdir}/${genome}_refseq_transcripts_pos.bed

# get distances of the 5 closest transcript pairs
bedtools closest -S -k 5 -D ref \
    -a ${outdir}/${genome}_refseq_transcripts_pos.bed \
    -b ${annotations} \
    > ${outdir}/${genome}_refseq_transcripts_closest_pos_to_neg.bed
