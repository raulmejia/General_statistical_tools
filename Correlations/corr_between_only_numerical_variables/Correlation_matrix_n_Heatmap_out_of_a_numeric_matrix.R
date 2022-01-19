# This script do some correlation plots over a numeric matrix
#   You can choose between this 

# The structure of your numeric matrix should be tsv
# "name"  "NumericColumn2" "NumericColumn3"
# "ID_A"  0.9 0.2
# "ID_B"  0.6 0.3
# "ID_C"  9.6 0.2

# choose a method of the cor function (stats::cor) or it came from corrplot::cor ?
# my_method = "pearson";  my_method = "kendall" ; my_method = "spearman"

# Example of use:
# Rscript /Path/to/this/nice/script/ .R \
# -m /Path/to/the/numeric/matrix.txt \
# -a "pearson" \
# -l some_label_for_the_results \
# -o /my/output/file.pdf

#################
### installing required packages
#################
if (!require("Hmisc")) {
  install.packages("Hmisc", ask =FALSE)
  library("Hmisc")
}
if (!require("corrplot")) {
  install.packages("corrplot", ask =FALSE)
  library("corrplot")
}
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}

#######################
### defining functions to be used afterwards 
#######################
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
# Column 1 : row names (variable 1 for the correlation test)
# Column 2 : column names (variable 2 for the correlation test)
# Column 3 : the correlation coefficients
# Column 4 : the p-values of the correlations

############################## 
## Data given by the user
##############################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-m", "--matrix", type="character", 
                    help="path to your expression matrix")
parser$add_argument("-a", "--algorithm", type="character", 
                    help="algorithm/method to be used to calculate the correlation")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-o", "--outputfile", type="character", 
                    help="output file where you want to store your correlation plots")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## Reading or preparing the inputs
#############################
mymatrix <-read.table( file=args$matrix, stringsAsFactors = FALSE , check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toys_merged_quantile_norm_by_batch.txt", stringsAsFactors = FALSE, check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Merged/Exp_Mat_GSE115989_MNHK.tsv", stringsAsFactors = FALSE, check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Merged/Exp_Mat_MK_GSE113342LE_GSE115989RJ_MajaL_GSE89880.txt", stringsAsFactors = FALSE, check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swipped_LipCon2-3_for_NorChl2-3.txt", stringsAsFactors = FALSE, check.names = FALSE) 

################
## Variables defined by the user
################
# path_to_my_numerical_matrix <- "/media/rmejia/mountme88/Projects/Maja/Interheart_Project/Data_2/InterheartScoringconsensus_aktuell_1221_header_corrected_nocommas_allvoids_are_NAs_nodouble_points_extracted_cols_of_interest_9n15to27.tsv"
path_to_my_numerical_matrix <-args$matrix
  
# my_choosen_algorithm <- "pearson"
my_choosen_algorithm <- args$algorithm

# my_label <- "Pearson correlation between variables"
my_label <-args$label
  
# path_to_save_the_pdf_graphs <- "/media/rmejia/mountme88/Projects/Maja/Interheart_Project/Data_2/InterheartScoringconsensus_aktuell_1221_header_corrected_nocommas_allvoids_are_NAs_nodouble_points_extracted_cols_of_interest_9n15to27.pdf"
path_to_save_the_pdf_graphs <- args$outputfile

dir.create( dirname( path_to_save_the_pdf_graphs ) , recursive = TRUE ) # creating the output folder if not exists already

#########################
#### The program starts 
#########################
mytable_tsv <- read.table( path_to_my_numerical_matrix , header = TRUE)

#########################
# Quality control of the matrix is it numerical? has it NAs?
#########################
matrix_noNAs <- na.omit(mytable_tsv) # eliminating NAs
matrix_noNAs_makeme_nummat <- as.matrix( matrix_noNAs ) # 
mode(matrix_noNAs_makeme_nummat) <- "numeric"
str(matrix_noNAs_makeme_nummat); dim( matrix_noNAs_makeme_nummat)

######
# Calculating Correlations
######
cor_matrix <- cor(matrix_noNAs_makeme_nummat , method= my_choosen_algorithm  )

# another method to calculate the correlations
rcorred_matrix <- rcorr( matrix_noNAs_makeme_nummat , type= my_choosen_algorithm )
rcorred_matrix$r # correlation coefficients
rcorred_matrix$P # Extract p-values

flattenCorrMatrix( rcorred_matrix$r, rcorred_matrix$P) # handy melted data frame to facilitate plotting 

pdf( file=path_to_save_the_pdf_graphs,  width = 10, height = 7) # let's save the graph in a pdf

corrplot( cor_matrix , type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,
         title=my_label)
#corrplot(rcorred_matrix$r, type="upper", order="hclust", 
#         p.mat = res2$P, sig.level = 0.01, insig = "blank")
col <- colorRampPalette( c("blue", "white", "red") )(20)
heatmap(x = cor_matrix, col = col, symm = TRUE,
        margins = c(15,15) , main = my_label)

dev.off()