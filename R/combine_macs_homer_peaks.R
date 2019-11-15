#' combine peaks from homer and macs2 output
#' @description macs2 callpeak and homer: analyzeChIP-Seq.pl output peaks are combined based on overlaps in genomic regions
#'
#' @param macs_input A character vector, determining path to MACS2 output. Only
#'   \code{ .narrowPeak or .broadPeak  files obtained by macs2 callpeak }
#'   should be used, with standard MACS2 output
#'   (https://github.com/taoliu/MACS). 7th column from macs2 output is referred as fold_change.
#' @param homer_input A character vector, determining path to HOMER output. File
#'   ending with \code{peaks.txt, obtained by analyzeChIP-Seq.pl} should be
#'   used. First 34 lines containing textual information about peaks will be
#'   ignored. 'Clonal Fold Change' column is referred as fold_change
#'
#' @return A tibble, consisting combined targets from MACS2 and HOMER output,
#'   with their classification as unique to MACS/HOMER or common from both
#'   output.Peaks are designated common if they have either partial or complete
#'   overlap.
#'   \code{.bed} file of combined peaks to upload in genome browser and \code{_TotalPeaks.tab} file to input in \code{annotate_peaks()}
#' @export
#' @importFrom readr read_delim
#' @importFrom readr write_delim
#' @import dplyr
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#'
#' @examples
#' \dontrun{
#'
#' homer_file <- "ABC_HOMER_peaks.txt"
#' macs_file <- "ABC_narrow_peaks.narrowPeak"
#' combine_macs_homer_peaks(macs_input=macs_file, homer_input= homer_file)
#' }
#'

combine_macs_homer_peaks <- function(macs_input, homer_input){

          library(dplyr)

          # read homer peaks.txt file
          homer_output <- readr::read_delim(homer_input,delim = "\t" ,skip = 34) %>%
                    dplyr::select(c("chr", "start", "end", `#PeakID`, `Clonal Fold Change`))

          homer_gr <- GenomicRanges::makeGRangesFromDataFrame(homer_output, keep.extra.columns = TRUE)

          names(GenomicRanges::mcols(homer_gr)) <- c("peak_id", "fold_change")
          homer_gr$type="UniqueToHOMER"

          # read macs .narrowpeaks or .broadpeaks file
          macs_output <- readr::read_delim(macs_input,delim = "\t", col_names = FALSE) %>%
                    dplyr::select(c(1,2,3,4,7))

          colnames(macs_output) <- c("chr", "start", "end","peak_id", "fold_change" )

          macs2_gr <- GenomicRanges::makeGRangesFromDataFrame(macs_output, keep.extra.columns = TRUE)
          macs2_gr$type="UniqueToMACS"

          # determine peaks common in macs and homer
          commonInBoth <- GenomicRanges::findOverlaps(homer_gr,macs2_gr, type="any",select = "all", ignore.strand=TRUE)

          CommonInBoth.df = data.frame(macs2_gr[S4Vectors::subjectHits(commonInBoth),], homer_gr[S4Vectors::queryHits(commonInBoth),])

          message("common peaks: ", nrow(CommonInBoth.df))

          # determine peaks unique to macs and homer
          UniqueToMACS <- subset(macs2_gr, !(macs2_gr$peak_id %in% CommonInBoth.df$peak_id))
          UniqueToHOMER <- subset(homer_gr, !(homer_gr$peak_id %in% CommonInBoth.df$peak_id.1))

          common_gr <- CommonInBoth.df %>%
                    dplyr::select(c("seqnames", "start", "end","peak_id", "fold_change" )) %>%
                    dplyr::mutate(type="CommonPeaks") %>%
                    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

          # combine all peaks
          all_peaks_gr <- c(UniqueToMACS,UniqueToHOMER,common_gr)



          all_peaks_gr <- all_peaks_gr %>% tidyr::as_tibble()

         readr::write_delim(x = all_peaks_gr, path = paste(gsub("_peaks.*","", basename(macs_input)), "_TotalPeaks.tab",sep=""),delim = "\t")

         all_peaks_gr %>%
                   dplyr::select(c("seqnames", "start", "end", "peak_id", "fold_change", "strand", "type")) %>%
                   readr::write_delim(path = paste(gsub("_peaks.*","", basename(macs_input)), "_TotalPeaks.bed",sep=""),delim = "\t", col_names = FALSE)

         return(all_peaks_gr)




}




