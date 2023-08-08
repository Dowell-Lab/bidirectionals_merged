suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(optparse)) ##adding arguments 

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(make_option(c("-c", "--closest"), type="character", default=NULL, 
			      help="Bed file with the closest pos/neg transcripts", metavar="character"),          
    make_option(c("-d", "--dist"), type="integer", default=1000, 
          help="Maximum distance between pairs [default = %default]", metavar="integer"),
    make_option(c("-m", "--mindist"), type="integer", default=0, 
          help="Minimum distance between pairs [default = %default]", metavar="integer"),
    make_option(c("-g", "--genome"), type="character", default=NULL, 
          help="Genome type (human or mouse)", metavar="character"),
	  make_option(c("-o", "--out"), type="character", default="./", 
              help="path to output directory [default = %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$closest)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-c bedfile).n", call.=FALSE)
}

#################################################
# 1: Initialize files                          ##
#################################################
max_distance <- as.numeric(opt$dist)
min_distance <- as.numeric(opt$mindist)
output_folder <- opt$out
closest <- data.table::fread(opt$closest)

# rename column names for imported data.tables
colnames(closest) <- c("chr1", "start1", "stop1", 
                        "transcript_id1", "score1", "strand1",
                       "chr2", "start2", "stop2", 
                       "transcript_id2", "score2", "strand2", 
                       "distance")

#add gene length column
closest$length1 <- closest$stop1 - closest$start1

###################################################
# 2: Get bidirectionals that overlap run-on txpts #
###################################################
# Get converging regions w/ <=1kb
converge <- closest[closest$distance > min_distance & closest$distance <= max_distance,]

# Get converging regions(right most of Pos gene & left most of neg gene)
start <- converge$stop1
end <- converge$start2
distance <- converge$distance
name <- paste0(converge$transcript_id1, "||", converge$transcript_id2)

# now get the center of the distances and add 1kb both directions for region
center <- round((start + end)/2)
start <- center - max_distance
end <- center + max_distance

conv_regions <- data.table::data.table("chr"=converge$chr1, 
                                      "start"=start, 
                                      "stop"=end, 
                                      "name"=name,
                                      "orig_distance"=distance, 
                                      "new_distance"=end-start)

if (opt$genome == "human"){
    
     data.table::fwrite(conv_regions,
                            col.names = FALSE, 
                            sep='\t',
                            paste0(output_folder,"hg38_refseq_transcripts_converging_1kb.bed"))

} else {
    data.table::fwrite(conv_regions,
                            col.names = FALSE, 
                            sep='\t',
                            paste0(output_folder,"mm10_refseq_transcripts_converging_1kb.bed"))
}
