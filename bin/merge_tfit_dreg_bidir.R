suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(optparse)) ##adding arguments 

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
	make_option(c("-b", "--bidir_overlap"), type="character", default=NULL,
		help="Bed6 file of overlapping bidirections (first 5 columns are Tfit and followign 5 dREG)", metavar="character"),
	make_option(c("-t", "--tfit"), type="character", default=NULL, 
		help="Tfit only bed file", metavar="character"),
    	make_option(c("-d", "--dreg"), type="character", default=NULL, 
          	help="dREG only bed file", metavar="character"),
    	make_option(c("-m", "--metadata"), type="character", default=NULL, 
          	help="Text file with metadata of samples", metavar="character"),
    	make_option(c("-n", "--metadata_qc"), type="character", default=NULL, 
          	help="Paper QC score (average sample scores)", metavar="character"),
    	make_option(c("-c", "--cg"), type="character", default=NULL, 
          	help="Folder with CG base content files for all papers", metavar="character"),
    	make_option(c("-g", "--genome"), type="character", default=NULL, 
          	help="Genome type (human or mouse)", metavar="character"),
    	make_option(c("-q", "--qc"), type="integer", default=2, 
          	help="Maximum paper QC to filter [default = %default]", metavar="integer"),
	make_option(c("-o", "--out"), type="character", default="./", 
        	help="path to output directory [default = %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bidir_overlap)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-b for bed file).n", call.=FALSE)
}

## Load the data based on optparse
#metadata
meta <- data.table::fread(opt$metadata)
qc_paper <- read.table(opt$metadata_qc, fill = TRUE, sep=',', header=FALSE)

#overlaps
overlap <- data.table::fread(opt$bidir_overlap)
tfit_only <- data.table::fread(opt$tfit)
dreg_only <- data.table::fread(opt$dreg)

#variables and output dir
max_qc <- opt$qc
output_folder <- opt$out

#base composition files
bc_files <- list.files(path=opt$cg, pattern="txt$", full.names=TRUE)
print("Base Composition Files: ")
print(bc_files)

#########################################
## subset the data.tables for overlaps ##
#########################################
# First five columns are Tfit
overlap_tfit <- unique(overlap[,1:5])
overlap_tfit$type <- 'Tfit Overlap'
overlap_tfit$bidirs <- 'Tfit'

# The following 5 are dREG regions
overlap_dreg <- unique(overlap[,6:10])
overlap_dreg$type <- 'dREG Overlap'
overlap_dreg$bidirs <- 'dREG'

## Add type and bidir columns to the *_only data.tables
# Tfit only data
tfit_only$type <- 'Tfit Only'
tfit_only$bidirs <- 'Tfit'
tfit_only_singletons <- subset(tfit_only, V5 == 1) ##get call from 1 sample
tfit_only_singletons$V4 <- as.character(lapply(strsplit(as.character(tfit_only_singletons$V4), '_'), `[`, 1)) #author name without cell type

# dREG only data
dreg_only$type <- 'dREG Only'
dreg_only$bidirs <- 'dREG'
dreg_only_singletons <- subset(dreg_only, V5 == 1) ##get call from 1 sample
dreg_only_singletons$V4 <- as.character(lapply(strsplit(as.character(dreg_only_singletons$V4), '_'), `[`, 1)) #author name without cell type

##########################################
## get metadata by genome type          ##
##########################################

if (opt$genome == "human"){

    meta_genome <- subset(meta, organism=='H. sapiens')

} else {

    meta_genome <- subset(meta, organism=='M. musculus')

}

meta_genome_qc <- unique(meta_genome[,c('paper_id', 'samp_qc_score')])
meta_genome_qc_filtered <- subset(meta_genome_qc, samp_qc_score<=max_qc) 
meta_genome_qc_filtered_papers <- unique(meta_genome_qc_filtered$paper_id)

#############################################
## get GC content data for paper calls     ##
#############################################
##load all base composition files
bc_list <- lapply(bc_files, data.table::fread)
base_composition <- do.call(rbind, bc_list)
base_composition$bidir_caller <- as.character(lapply(lapply((strsplit(as.character(base_composition$id), 
                                                                        '_')),
                                                                         rev), `[`, 1))
base_composition$sample_name <- substr(base_composition$id,1,nchar(base_composition$id)-5)
print(head(base_composition,2))

##Tfit base compositions
base_composition_tfit <- subset(base_composition, bidir_caller=='tfit')

##dREG base composition
base_composition_dreg <- subset(base_composition, bidir_caller=='dreg')

##base composition Tfit and dREG wide
base_composition_tfit_dreg <- merge(base_composition_tfit, base_composition_dreg, by='sample_name')
base_composition_tfit_dreg$author <- as.character(lapply(strsplit(as.character(base_composition_tfit_dreg$sample_name), '_'), `[`, 1))
colnames(base_composition_tfit_dreg) <- c('sample_name','id_tfit','cg_tfit','at_tfit', 'bidir_caller_tfit',
                                         'id_dreg','cg_dreg','at_dreg', 'bidir_caller_dreg','author')

##merge with paper QC data
base_composition_tfit_dreg_qc <- merge(base_composition_tfit_dreg, qc_paper, by.x='author', by.y='V1')

##extract samples with low GC and have a high QC 
##since these will be removed 
tfit_dreg_low_gc <- subset(base_composition_tfit_dreg_qc,
                                  as.numeric(cg_tfit) < 0.49 &
                                  as.numeric(cg_dreg) < 0.49 )

print(head(tfit_dreg_low_gc, 2))
tfit_dreg_low_gc_qc <- subset(tfit_dreg_low_gc, V2 <=2)

#get first author's last name for filtering the samples 
tfit_dreg_low_gc_qc_samples <- unique(as.character(lapply(strsplit(as.character(tfit_dreg_low_gc_qc$author), '2'), `[`, 1)))

#############################################
## subsetting and getting high QC samples  ##
#############################################
#get high QC singletons
#Tfit
tfit_only_singletons_qc <- tfit_only_singletons[tfit_only_singletons$V4 %in% meta_genome_qc_filtered_papers,]

#dREG
dreg_only_singletons_qc <- dreg_only_singletons[dreg_only_singletons$V4 %in% meta_genome_qc_filtered_papers,]


##excluding poor QC singletons...
master_qc1_2 <- rbind(subset(tfit_only, V5 != 1),
                        overlap_tfit,
                        subset(dreg_only, V5 != 1),
                        dreg_only_singletons_qc,
                        tfit_only_singletons_qc,
                        fill=TRUE)

#get the number of low QC and GC papers supporting a call
master_qc1_2$low_gc_num <- lengths(regmatches(master_qc1_2$V4,
                                                   gregexpr(tfit_dreg_low_gc_qc_samples, 
                                                            master_qc1_2$V4)))

#remove those call that are only low QC and GC papers
master_qc1_2_gc_filter <- master_qc1_2[!master_qc1_2$V5==master_qc1_2$low_gc_num,]
print(head(master_qc1_2_gc_filter, 2))

#add a length filter
if (opt$genome == "human"){

    master_qc1_2_gc_filter$width <- as.numeric(master_qc1_2_gc_filter$V3) - as.numeric(master_qc1_2_gc_filter$V2)
    master_qc1_2_gc_len_filter <- subset(master_qc1_2_gc_filter, width > 150 & width < 2500)


} else {

    master_qc1_2_gc_filter$width <- as.numeric(master_qc1_2_gc_filter$V3) - as.numeric(master_qc1_2_gc_filter$V2)
    master_qc1_2_gc_len_filter <- subset(master_qc1_2_gc_filter, width > 150 & width < 2000)

}

master_qc1_2_gc_len_filter$bidir_type <- ifelse(master_qc1_2_gc_len_filter$type == "Tfit Only",
                                                "tfit",
                                                ifelse(master_qc1_2_gc_len_filter$type == "dREG Only","dreg","tfit,dreg"))

master_final <- master_qc1_2_gc_len_filter[, c(1,2,3,10,5)]
master_final$strand <- "."

#There were some calls from EBV. Removing them here
master_final_chr <- subset(master_final, V1 != "chrEBV") #master_final[!master_final$V1 %in% c("chrEBV"),]

if (opt$genome == "human"){

    data.table::fwrite(master_final_chr,
                        col.names = FALSE, 
                        sep='\t',
                        paste0(output_folder,"hg38_master_mumerge_qc_gc_len_filtered.bed"))
                
} else {

    data.table::fwrite(master_final_chr,
                        col.names = FALSE, 
                        sep='\t',
                        paste0(output_folder,"mm10_master_mumerge_qc_gc_len_filtered.bed"))

}
