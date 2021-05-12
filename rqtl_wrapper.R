library(optparse)
library(qtl)
library(stringi)

option_list = list(
    make_option(c("-g", "--geno"), type="character", help=".geno file containing a dataset's genotypes"),
    make_option(c("-p", "--pheno"), type="character", help="File containing two columns - sample names and values"),
    make_option(c("-c", "--covar"), type="character", help="File containing covariates - first column sample names, other columns covariate values"),
    make_option(c("--model"), type="character", default="normal", help="Mapping Model - Normal or Non-Parametric"),
    make_option(c("--method"), type="character", default="hk", help="Mapping Method - hk (Haley Knott), ehk (Extended Haley Knott), mr (Marker Regression), em (Expectation-Maximization), imp (Imputation)"),
    make_option(c("-i", "--interval"), action="store_true", default=NULL, help="Use interval mapping"),
    make_option(c("--perm"), type="integer", default=0, help="Number of permutations"),
    make_option(c("-s", "--scale"), type="character", default="mb", help="Mapping scale - Megabases (Mb) or Centimorgans (cM)"),
    make_option(c("--control_marker"), type="character", help="Name of marker (contained in genotype file) to be used as a control"),
    make_option(c("-v", "--verbose"), action="store_true", default=NULL, help="Show extra information")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

verbose_print <- function(...){
    if (!is.null(opt$verbose)) {
        for(item in list(...)){
            cat(item)
        }
        cat("\n")
    }
}

if (is.null(opt$geno) || is.null(opt$pheno)){
    print_help(opt_parser)
    stop("Both a genotype and phenotype file must be provided.", call.=FALSE)
}

tmp_dir = Sys.getenv("TMPDIR")

geno_file = opt$geno
pheno_file = opt$pheno
cross_file = file.path(tmp_dir, "cross", paste(stri_rand_strings(1, 8), ".cross", sep = "")) # Generate randomized filename for cross object

trim <- function( x ) { gsub("(^[[:space:]]+|[[:space:]]+$)", "", x) }

get_geno_code <- function(header, name = 'unk'){
    mat = which(unlist(lapply(header,function(x){ length(grep(paste('@',name,sep=''), x)) })) == 1)
    return(trim(strsplit(header[mat],':')[[1]][2]))
}

geno_to_csvr <- function(genotypes, out, phenotype = NULL, sex = NULL, mapping_scale = "Mb", verbose = FALSE){
    header = readLines(genotypes, 40)                                                                                 # Assume a geno header is not longer than 40 lines
    toskip = which(unlist(lapply(header, function(x){ length(grep("Chr\t", x)) })) == 1)-1                            # Major hack to skip the geno headers
    type <- get_geno_code(header, 'type')
    if(type == '4-way'){
    genocodes <- NULL
    } else {
    genocodes <- c(get_geno_code(header, 'mat'), get_geno_code(header, 'het'), get_geno_code(header, 'pat'))             # Get the genotype codes
    }
    genodata <- read.csv(genotypes, sep='\t', skip=toskip, header=TRUE, na.strings=get_geno_code(header,'unk'), colClasses='character', comment.char = '#')
    verbose_print('Genodata:', toskip, " ", dim(genodata), genocodes, '\n')
    if(is.null(phenotype)) phenotype <- runif((ncol(genodata)-4))                                                     # If there isn't a phenotype, generate a random one
    if(is.null(sex)) sex <- rep('m', (ncol(genodata)-4))                                                              # If there isn't a sex phenotype, treat all as males
    outCSVR <- rbind(c('Pheno', '', '', phenotype),                                                                   # Phenotype
                    c('sex', '', '', sex),                                                                           # Sex phenotype for the mice
                    cbind(genodata[,c('Locus','Chr', mapping_scale)], genodata[, 5:ncol(genodata)]))                          # Genotypes
    write.table(outCSVR, file = out, row.names=FALSE, col.names=FALSE,quote=FALSE, sep=',')                           # Save it to a file
    if(type == '4-way'){
    verbose_print('Loading in as 4-WAY\n')
    cross = read.cross(file=out, 'csvr', genotypes=NULL, crosstype="4way")                                         # Load the created cross file using R/qtl read.cross
    }else if(type == 'f2'){
    verbose_print('Loading in as F2\n')
    cross = read.cross(file=out, 'csvr', genotypes=genocodes, crosstype="f2")                                       # Load the created cross file using R/qtl read.cross
    }else{
    verbose_print('Loading in as normal\n')
    cross = read.cross(file=out, 'csvr', genotypes=genocodes)                                                       # Load the created cross file using R/qtl read.cross
    }
    if(type == 'riset'){
    verbose_print('Converting to RISELF\n')
    cross <- convert2riself(cross)                                                                # If its a RIL, convert to a RIL in R/qtl
    }
    return(cross)
}

gen_pheno_vector_from_file <- function(pheno_file){
    df <- read.table(pheno_file, na.strings = "x", header=TRUE, check.names=FALSE)
    sample_names <- df$Sample
    trait_name <- colnames(df)[2]
    vals <- df[trait_name]

    return(list(trait_name, sample_names, vals))
}

sample_vals = gen_pheno_vector_from_file(pheno_file)
trait_name = sample_vals[1]
samples_vector = unlist(sample_vals[2])
pheno_vector = unlist(sample_vals[3])

verbose_print('Generating cross object\n')
cross_object = geno_to_csvr(geno_file, cross_file)

verbose_print('Calculating genotype probabilities\n')
cross_object = calc.genoprob(cross_object)

verbose_print('Adding phenotype to cross object\n')
cross_object$pheno <- pheno_vector