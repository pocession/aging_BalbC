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
#' @importFrom dplyr select filter left_join
#' @importFrom annotatr annotate_regions build_annotations
#' @importFrom GenomicRanges GRanges
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
## Run getDMR2018.Rmd or getDMR2022.Rmd to get DMR data

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
  ## Get annotation
  annotated_df <- .get_annotated_gr(input_df_pos, col_names)
  
  ## Get id and annotation information
  annotated_df_subset <- annotated_df |>
    dplyr::select(id)
  annotated_df_subset <- cbind(annotated_df_subset, annotated_df[,6:ncol(annotated_df)])
  
  final_df <- input_df_pos |>
    dplyr::left_join(input_df_meth, by = "id") |>
    dplyr::left_join(annotated_df_subset, by = "id")
  
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
#' @param col_names a character vector specifying the column names of return dataframe
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
  
  # Create GRange objects ------------------------------------------------------
  ## For creating GRange objects
  ## The id column is used as the identity column
  
  gr <- GenomicRanges::GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end, strand = df$strand),
    id = df$id
  )
  
  # Annotate -------------------------------------------------------------------
  ## Select annotations for intersection with regions
  ## Note inclusion of custom annotation, and use of shortcuts
  annots <- c('hg38_basicgenes')
  
  ## Build the annotations (a single GRanges object)
  annotations <- annotatr::build_annotations(genome = 'hg38', annotations = annots)
  
  dm_annotated <- annotatr::annotate_regions(
    regions = gr,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  dm_annotated <- data.frame(dm_annotated)
  
  ## Only pick up the first annotation
  dm_annotated <- dm_annotated |>
    dplyr::distinct(id, .keep_all = TRUE) |>
    dplyr::select(-c(strand))
  
  ## Add the strand information back
  df_strand <- df |>
    dplyr::select(id, strand)
  
  dm_annotated <- dm_annotated |>
    dplyr::left_join(df_strand, by = "id")
  
  # Generate the final data frame ----------------------------------------------
  final_data <- dm_annotated |>
    dplyr::select(id, seqnames, start, end, strand)
  
  ## Remove the strand column
  dm_annotated <- dm_annotated |>
    dplyr::select(-c(strand))
  
  colnames(final_data) <- c("id", col_names)
  final_data <- final_data |>
    dplyr::left_join(dm_annotated[,5:ncol(dm_annotated)], by = "id")
  
  return(final_data)
}
