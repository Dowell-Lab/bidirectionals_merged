suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(optparse)) ##adding arguments 

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
	    make_option(c("-b", "--bed"), type="character", default=NULL, 
	          help="mumerge bed files with samples", metavar="character"),
		  make_option(c("-o", "--out"), type="character", default="./", 
              help="path to output directory [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bed)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-b for bed file).n", call.=FALSE)
}

#define paths based in optparse
bed_file <- data.table::fread(opt$bed)
bed_file_name <- sub(pattern = "(.*)\\..*$", 
                    replacement = "\\1",
                    base::basename(opt$bed))
output_folder <- opt$out

#count the number of samples that support a region
bed_file$num_papers <- lapply(strsplit(as.character(bed_file$V4), ','), length)

#save files
final_path <- paste0(output_folder, 
                    bed_file_name, 
                    '_numSamples.bed')
    
data.table::fwrite(bed_file,
		col.names = FALSE,
                    final_path,
                    sep='\t')

#print summary of R session
print("Session Summary")

print(sessionInfo())

print("DONE!")
