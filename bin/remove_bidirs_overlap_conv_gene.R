suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(optparse)) ##adding arguments 

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
	    make_option(c("-a", "--bidir"), type="character", default=NULL, 
	          help="Bed file for master bidirectionals after filtering (qc, gc, len & pos/neg gene filter)", metavar="character"),
			make_option(c("-c", "--closest"), type="character", default=NULL, 
			      help="Bed file with the closest pos/neg transcripts", metavar="character"),          
    make_option(c("-d", "--dist"), type="character", default=1000, 
          help="Distance between pairs", metavar="integer"),
    make_option(c("-m", "--mindist"), type="character", default=0, 
          help="Distance between pairs", metavar="integer"),
    make_option(c("-g", "--genome"), type="character", default=NULL, 
          help="Genome type (human or mouse)", metavar="character"),
	  make_option(c("-o", "--out"), type="character", default="./", 
              help="path to output directory [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bidir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-a bedfile for bidirectional transcripts).n", call.=FALSE)
}

#################################################
# 1: Initialize files                          ##
#################################################
max_distance <- as.numeric(opt$dist)
min_distance <- as.numeric(opt$mindist)
output_folder <- opt$out
bidirs <- data.table::fread(opt$bidir)
closest <- data.table::fread(opt$closest)

# rename column names for imported data.tables
colnames(bidirs) <- c("chr", "bidir_start", "bidir_stop", 
                    "source", "score", "strand", "bidir_id")

colnames(closest) <- c("chr1", "start1", "stop1", "bidir", "score1", "strand1", "bidir_id", 
                       "chr2", "start2", "stop2", "gene_id","length","dest", "gene_id2", "distance")

###################################################
# 2: Get bidirectionals that overlap run-on txpts #
###################################################
# Get converging regions w/ <=1kb

bidirs_non_conv <- bidirs[!bidirs$bidir_id %in% unique(closest$bidir_id),]

if (opt$genome == "human"){
    
    data.table::fwrite(bidirs_non_conv,
                            col.names = FALSE, 
                            sep='\t',
                            paste0(output_folder,"hg38_master_bidirectionals.bed"))

} else {

    data.table::fwrite(bidirs_non_conv,
                            col.names = FALSE, 
                            sep='\t',
                            paste0(output_folder,"mm10_master_bidirectionals.bed"))
}