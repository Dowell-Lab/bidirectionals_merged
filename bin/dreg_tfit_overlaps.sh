while getopts d:t:g:o: flag
do
    case "${flag}" in
        d) dreg=${OPTARG};;
        t) tfit=${OPTARG};;
        g) genome=${OPTARG};;
        o) outdir=${OPTARG};;
    esac
done

echo "dREG File: $dreg";
echo "Tfit File: $tfit";
echo "Genome: $genome";
echo "Output Directory: $outdir"

echo "Overlapping bidirectionals calls with bedtools"

bedtools intersect -wa -wb -a ${tfit} -b ${dreg} > ${outdir}/${genome}_overlap.bed
bedtools intersect -a ${tfit} -b ${dreg} -v > ${outdir}/${genome}_tfit_only.bed
bedtools intersect -a ${dreg} -b ${tfit} -v > ${outdir}/${genome}_dreg_only.bed

echo "DONE!"
