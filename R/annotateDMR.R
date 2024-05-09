#' annotateDMR
#'
#' @description
#' This function reads differentially-methylated regions (DMRs) and annotate with hg38 information
#' 
#' @author Tsunghan Hsieh
#'
#' @param input a character string specifying the name and path of the input file (.csv)
#' @param col_names a character vector specifying the column names of choromosome, starting position, ending position
#' Should be in this format: c(chr, start, end, strand)
#' @param output A character string specifying the name and path of output file
#' the default path is as same as the input file
#' the default output file name is the input file name prefixed with "annotated_"
#' @return a dataframe of the output data
#'
#' @export
#'
#' @importFrom here here
#' @importFrom assertthat assert_that
#' @importFrom tidyr separate
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom org.Hs.eg.db
#' @importFrom GenomicRanges
#' @importFrom AnnotationDbi select
#'
#' @examples
#' \dontrun{
#' df <- annotateDMR(
#' input = here::here("./Results/Statistics/DMR/DMR_5hmc_2018_SI_3.csv"),
#' col_names = c("chr", "start", "end", "strand"),
#' output = NULL
#' )
#' }
#' 
## TODO
## Run getDIRWithNoReplicate() to generate DIR files

annotateDMR <- function(input, col_names, output = NULL) {
  # check input ----------------------------------------------------------------
  
  assertthat::assert_that(is.character(input), 
                          msg = "The given argument 'input' is not a string.\n")
  
  assertthat::assert_that(grepl("\\.csv$", input), 
                          msg = "The given file name of 'input' is not a .csv file.\n")
  
  assertthat::assert_that(file.exists(here::here(input)), 
                          msg = "The given file of 'input' does not exist.\n")
  
  assertthat::assert_that(length(col_names) == 4, 
                          msg = "There are not enough column numbers of the given argument 'col_names'.\n")
  
  # check output ---------------------------------------------------------------
  
  ## Generate output file name
  if (is.null(output)) {
    # Use sub() to extract the desired part
    outputfname <- sub(".*/([^.]+)\\.csv", "\\1", input)
    outputfname <- paste0("annotated_", outputfname, ".csv")
    outputdir <- here::here("Results/Statistics/DMR_annotated/")
    output <- here::here(outputdir, outputfname)
    remove(outputfname, outputdir)
  } else {
    assertthat::assert_that(is.character(output), 
                            msg = "The given argument 'output' is not a string.\n")
    
    assertthat::assert_that(grepl("\\.csv$", output), 
                            msg = "The given file name of 'output' is not a .csv file.\n")
    
    assertthat::assert_that(file.exists(here::here(dirname(output))), 
                            msg = "The given file path of 'output' does not exist.\n")
  }
  
  
  # Assign the hg38 objects ----------------------------------------------------
  .GlobalEnv$txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  hgdb <- org.Hs.eg.db
  
  # Extract TSS positions ------------------------------------------------------
  tss_fname <- here::here("Results/Tmp/tss_df.rds") 
  if (file.exists(tss_fname)) {
    tss_df <- readRDS(tss_fname)
  } else {
    tss <- transcripts(txdb, columns=c("tx_id", "gene_id"), use.names=TRUE)
    tss <- resize(tss, width=1, fix='start')
    tss_df <- as.data.frame(tss, row.names = NULL)
    tss_df$transcript <- row.names(tss_df)
  }
  colnames(tss_df) <- c("tx_chr","tx_start","tx_end","tx_width","tx_strand","tx_id","gene_id","transcript")
  
  # Read DMRs ------------------------------------------------------------------
  input_df <- read.csv(input)
  
  ## In case the row.names exists in input file
  input_df <- input_df[, (!names(input_df) %in% c("X"))]
  
  ## Assign unique id to the DMR data frame
  input_df$id <- seq(1:nrow(input_df))
  
  ## Get the position data frame
  ## Re-order the column
  input_df_pos <- input_df[c(col_names, "id")] |>
    dplyr::select(id, col_names[1], col_names[2], col_names[3], col_names[4])
  colnames(input_df_pos) <- c("id", "chr", "start", "end", "strand")
  
  ## Save the methylation data frame
  input_df_meth <- input_df[ , !(names(input_df) %in% col_names)]
  
  # Get the annotated data frame -----------------------------------------------
  
  annotated_df <- .get_annotated_gr(input_df_pos, col_names)
  final_df <- annotated_df[,1:5] |>
    dplyr::left_join(input_df_meth, by = "id")
  final_df <- cbind(final_df, annotated_df[,6:ncol(annotated_df)])
  final_df <- final_df |>
    dplyr::select(-c(id))
  write.csv(final_df, output, row.names = FALSE)
  message("Done! The output is written in: ", output)
  return(final_df)
}


#' .get_annotated_gr
#'
#' @description
#' This function annotate differentially-interacted regions (DIRs) with given genome information
#' 
#' @author Tsunghan Hsieh
#'
#' @param input a character string specifying the input dataframe
#' @return a dataframe of annotated file
#' @examples
#' \dontrun{
#' df <- .get_annotated_gr(
#' df = input_df_pos,
#' col_names = col_names
#' )
#' }
.get_annotated_gr <- function(df, col_names) {
  
  # Check colnames -------------------------------------------------------------
  ## Filter out non-standard chr
  wanted <- paste0("chr", seq(1:23))
  wanted <- c(wanted, "chrX", "chrY")
  df <- df |>
    dplyr::filter(chr %in% wanted)
  
  # Annotate region ------------------------------------------------------------
  ## For creating GRange objects
  ## The id column is used as the identity column
  
  ## Create GRange objects
  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end, strand = df$strand),
    id = df$id
  )
  
  # Finding overlap ------------------------------------------------------------
  features <- transcripts(txdb)
  
  ## Select the id of closet genes 
  ps.idx <- nearest(gr, features, ignore.strand = FALSE)
  
  ## check if any peak not matched to the nearest transcript
  na.idx <- is.na(ps.idx)
  
  ## Remove na 
  if (sum(na.idx) > 0) { ## suggested by Thomas Schwarzl
    ps.idx <- ps.idx[!na.idx]
    gr <- gr[!na.idx]
  }
  
  # Generate the final data frame ----------------------------------------------
  ## Identify the closet gene information
  features_df <- data.frame(txStart = start(features[ps.idx]),
                            txEnd = end(features[ps.idx]),
                            txStrand = strand(features[ps.idx]),
                            tx_id = features[ps.idx]$tx_id,
                            tx_name = features[ps.idx]$tx_name)
  
  
  gr_df <- as.data.frame(gr)
  annotated_df <- cbind(gr_df, features_df)
  annotated_df <- annotated_df |>
    dplyr::mutate(distanceToNearestTss = ifelse(txStrand == "+", txStart - start - width/2, txEnd - start - width/2))
  
  annotated_df <- annotated_df[,6:ncol(annotated_df)]
  
  ## map to the gene ids and symbols
  mapped_genes <- AnnotationDbi::select(txdb, 
                                         keys = annotated_df$tx_name, 
                                         keytype = "TXNAME", 
                                         columns = c("GENEID"))
  
  mapped_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                          keys = mapped_genes$GENEID,
                                          columns = "SYMBOL", 
                                          keytype = "ENTREZID")
  
  annotated_df <- cbind(annotated_df, mapped_genes, mapped_symbols)
  annotated_df <- annotated_df |>
    dplyr::mutate(gene_id = GENEID,
                  symbol_ = SYMBOL) |>
    dplyr::select(-c(TXNAME, GENEID, SYMBOL, ENTREZID))
  
  ## Generate the final data frame
  final_df <- df |>
    dplyr::left_join(annotated_df, by = "id")
  
  return(final_df)
}
