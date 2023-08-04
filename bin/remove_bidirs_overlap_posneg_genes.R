suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(optparse)) ##adding arguments 

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
	    make_option(c("-a", "--bidir"), type="character", default=NULL, 
	          help="Bed file for master bidirectionals", metavar="character"),
		  make_option(c("-b", "--bidirgene"), type="character", default=NULL, 
		        help="Bed file with overlap between bidirectionals and gene transcripts", metavar="character"),  
    make_option(c("-d", "--dist"), type="character", default=1000, 
          help="Distance between pairs", metavar="integer"),
    make_option(c("-f", "--frac"), type="character", default=0.975, 
          help="Fraction of overlap between bidirectionals and genes", metavar="float"),
    make_option(c("-g", "--genome"), type="character", default=NULL, 
          help="Genome type (human or mouse)", metavar="character"),
	  make_option(c("-o", "--out"), type="character", default="./", 
              help="path to output directory [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bidir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-a bedfile of bidirectionals).n", call.=FALSE)
}

#################################################
# 1: Initialize files                          ##
#################################################
max_distance <- as.numeric(opt$dist)
output_folder <- opt$out
bidirs <- data.table::fread(opt$bidir)
overlaps <- data.table::fread(opt$bidirgene)
frac_overlap <- opt$frac

# rename column names for imported data.tables
colnames(bidirs) <- c("chr", "bidir_start", "bidir_stop", 
                    "source", "score", "strand")

colnames(overlaps) <- c("gene_chr", "gene_start", "gene_stop", 
                        "transcript_id", "gene_score", "strand", 
                        "bidir_chr", "bidir_start", "bidir_stop",
                        "source","support","bidir_score", "overlap")

###################################################
# 2: Add column names                            ##
###################################################

# add bidirectional ids based on chr:start-stop-source
overlaps$bidir_id <- paste0(overlaps$bidir_chr,
                            ":", overlaps$bidir_start, 
                            "-", overlaps$bidir_stop, 
                            "-", overlaps$source)

bidirs$bidir_id <- paste0(bidirs$chr, 
                        ":", bidirs$bidir_start, 
                        "-", bidirs$bidir_stop, 
                        "-", bidirs$source)

# get Fraction of gene overlapped
overlaps$gene_length <- overlaps$gene_stop - overlaps$gene_start + 1
overlaps$frac_gene_overlap <- overlaps$overlap/overlaps$gene_length

# get the positive & negative overlaps separated
overlaps_pos <- overlaps[overlaps$strand == "+",]
overlaps_neg <- overlaps[overlaps$strand == "-",]

# get bidirectionas that are found in both negative and positice strands
bidirs_overlap_pos_neg <- intersect(unique(overlaps_neg$bidir_id), unique(overlaps_pos$bidir_id))

overlaps_frac <- unique(overlaps[overlaps$frac_gene_overlap > frac_overlap,]$bidir_id)

# get the new filtered bidirectionals
bidirs_overlap_genes <- bidirs[!bidirs$bidir_id %in% bidirs_overlap_pos_neg,]
bidirs_overlap_genes_frac <- bidirs_overlap_genes[!bidirs_overlap_genes$bidir_id %in% overlaps_frac,]

# save the files
if (opt$genome == "human"){
    data.table::fwrite(bidirs_overlap_genes_frac,
                            col.names = FALSE, 
                            sep='\t',
                            paste0(output_folder,"hg38_master_qc_gc_len_posneg_gene_filt.bed"))

} else {
    data.table::fwrite(bidirs_overlap_genes_frac,
                            col.names = FALSE, 
                            sep='\t',
                            paste0(output_folder,"mm10_master_qc_gc_len_posneg_gene_filt.bed"))
}
