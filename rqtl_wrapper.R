library(qtl)
library(stringi)

tmp_dir = Sys.getenv("TMPDIR")

args = commandArgs(trailingOnly=TRUE)

# Parsing command line args like this until optparse is installed
geno_file = args[1]
pheno_file = args[2]
cross_file = file.path(tmp_dir, "cross", stri_rand_strings(1, 8)) # Generate randomized filename for cross object

cat('Generating Cross Object\n')
cross_object = geno_to_csvr(geno_file, cross_file)

cat('Calculating genotype probabilities\n')
cross_object = calc.genoprob(cross_object)

trim <- function( x ) { gsub("(^[[:space:]]+|[[:space:]]+$)", "", x) }

getGenoCode <- function(header, name = 'unk'){
    mat = which(unlist(lapply(header,function(x){ length(grep(paste('@',name,sep=''), x)) })) == 1)
    return(trim(strsplit(header[mat],':')[[1]][2]))
}

geno_to_csvr <- function(genotypes, out, phenotype = NULL, sex = NULL, mapping_scale = "Mb", verbose = FALSE){
    header = readLines(genotypes, 40)                                                                                 # Assume a geno header is not longer than 40 lines
    toskip = which(unlist(lapply(header, function(x){ length(grep("Chr\t", x)) })) == 1)-1                            # Major hack to skip the geno headers
    type <- getGenoCode(header, 'type')
    if(type == '4-way'){
    genocodes <- NULL
    } else {
    genocodes <- c(getGenoCode(header, 'mat'), getGenoCode(header, 'het'), getGenoCode(header, 'pat'))             # Get the genotype codes
    }
    genodata <- read.csv(genotypes, sep='\t', skip=toskip, header=TRUE, na.strings=getGenoCode(header,'unk'), colClasses='character', comment.char = '#')
    cat('Genodata:', toskip, " ", dim(genodata), genocodes, '\n')
    if(is.null(phenotype)) phenotype <- runif((ncol(genodata)-4))                                                     # If there isn't a phenotype, generate a random one
    if(is.null(sex)) sex <- rep('m', (ncol(genodata)-4))                                                              # If there isn't a sex phenotype, treat all as males
    outCSVR <- rbind(c('Pheno', '', '', phenotype),                                                                   # Phenotype
                    c('sex', '', '', sex),                                                                           # Sex phenotype for the mice
                    cbind(genodata[,c('Locus','Chr', mapping_scale)], genodata[, 5:ncol(genodata)]))                          # Genotypes
    write.table(outCSVR, file = out, row.names=FALSE, col.names=FALSE,quote=FALSE, sep=',')                           # Save it to a file
    require(qtl)
    if(type == '4-way'){
    cat('Loading in as 4-WAY\n')
    cross = read.cross(file=out, 'csvr', genotypes=NULL, crosstype="4way")                                         # Load the created cross file using R/qtl read.cross
    }else if(type == 'f2'){
    cat('Loading in as F2\n')
    cross = read.cross(file=out, 'csvr', genotypes=genocodes, crosstype="f2")                                       # Load the created cross file using R/qtl read.cross
    }else{
    cat('Loading in as normal\n')
    cross = read.cross(file=out, 'csvr', genotypes=genocodes)                                                       # Load the created cross file using R/qtl read.cross
    }
    if(type == 'riset'){
    cat('Converting to RISELF\n')
    cross <- convert2riself(cross)                                                                # If its a RIL, convert to a RIL in R/qtl
    }
    return(cross)
}