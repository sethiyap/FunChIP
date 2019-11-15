#' annotate peaks
#' @description assign each peak to target genes based on their location
#'
#' @param feature_file A character vector containing path to the file. Either a
#'   bed file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand} OR
#'   gene feature file with extention \code{.gff or .gff3 }
#' @param peaks_file A character vector containing path to the peaks file.
#'   peaks_file can be either \code {1. macs2 callpeak output (i.e. .narrowPeak
#'   or .broadPeak) 2. _peaks.txt, obtained by Homer, analyzeChIP-Seq.pl 3.
#'   _TotalPeaks.tab, obtained by ChIPFun::combine_macs_homer_peaks()}
#'   For macs, 7th column from macs2 output is referred as fold_change. For Homer output,
#'   first 34 lines containing textual information about peaks will be ignored
#'   and 'Clonal Fold Change' column is referred as fold_change.
#' @param outfile Logical, to write output in excel sheet or not.
#' @param outfileName A character vector to determine name of output file.
#' @return A tibble or output file and pie-chart with target information and classification of
#'   targets with respect to peak location. Targets are classified as \itemize{
#'   \item InPromoter: If the peak lies at 5'-upstream of gene start as shown: (^
#'   ----->), where ^ :peak and -----> :gene in 5'-3' direction.
#'   \item WithInGeneBody: If the peak completely overlaps with the gene-body region
#'   (--^-->)
#'   \item  PartialOverlap_Gene: If the peak partially overlaps at gene
#'   start (-^---->). Further they are classified as
#'   PartialOverlap_Gene_1,2,3... if they overlap partially with more than one
#'   gene, usually the case in compact fungal genomes with many genes
#'   overlapping either on same strand or different strand.(<---^-
#'   ----->)
#'   \item AtGeneEnd: If the peak overlaps with 3'-end of a gene, and there is no gene
#'   in 5'-3' nearby (----^->)
#'   \item NearestTarget: Peaks with nearest gene
#'   promoter at distance > 2000, then assign them to nearest gene irrespective
#'   of strand. (--A--> ^     <--B--           --C-->  ). Here, ie. even-if the
#'   peak is in promoter according to gene-C but in NearestTarget category, peak
#'   will be assigned to gene-A due on its nearest distance }
#' @export
#'
#' @examples
#' \dontrun{
#'
#' bed_file <- "A_nidulans_FGSC_A4_version_s10-m04-r07_cds.bed"
#' peaks_file <- "myc_narrow_TotalPeaks.tab"
#' annotate_peaks(feature_file=bed_file,peaks_file=peaks_file, outfile=TRUE,outfileName="myc_peaks")
#'
#' }
annotate_peaks <- function(feature_file, peaks_file, outfile=FALSE, outfileName="sample"){

          library(dplyr)

          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".gff")){

                    print("Input is genome feature file")
                    gff = GenomicFeatures::makeTxDbFromGFF(feature_file)
                    subject = GenomicFeatures::genes(gff)
                    base::names(subject)=NULL
                    gene_name <- names(GenomicRanges::mcols(subject)[1])
          }
          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".bed")){
                    print("Input is bed file")
                    subject = rtracklayer::import.bed(feature_file)
                    gene_name <- names(GenomicRanges::mcols(subject)[1])
          }

          if(stringr::str_detect(string = as.character(basename(peaks_file)),paste(c(".narrowPeak", ".broadPeak"),collapse = '|'))){

                    print("Peaks file: MACS output")
                    macs_output <- readr::read_delim(peaks_file,delim = "\t", col_names = FALSE) %>%
                              dplyr::select(c(1,2,3,4,7))

                    colnames(macs_output) <- c("chr", "start", "end","peak_id", "fold_change")
                    query <- GenomicRanges::makeGRangesFromDataFrame(macs_output, keep.extra.columns = TRUE)
                    query$type="macs2_input"

          }

          if(stringr::str_detect(string = as.character(basename(peaks_file)),pattern="peaks.txt")){
                    print("Peaks file: HOMER output")

                    homer_output <- readr::read_delim(peaks_file,delim = "\t" ,skip = 34) %>%
                              dplyr::select(c("chr", "start", "end", `#PeakID`, `Clonal Fold Change`))

                    query <- GenomicRanges::makeGRangesFromDataFrame(homer_output, keep.extra.columns = TRUE)
                    names(GenomicRanges::mcols(query)) <- c("peak_id", "fold_change")
                    query$type="homer_input"
          }

          if(stringr::str_detect(string = as.character(basename(peaks_file)),pattern="_TotalPeaks.tab")){
                    print("Peaks file: combine_homer_macs_peak output")

                    peaks_input <- readr::read_delim(peaks_file, col_names = TRUE, delim = "\t")
                    query <- GenomicRanges::makeGRangesFromDataFrame(peaks_input, keep.extra.columns = TRUE)
          }

          ######################################################################################################

          #### preceeding gene to query
          # Step1
          pp = GenomicRanges::precede(query,subject,ignore.strand=F,select="all")
          df.p = data.frame(query[S4Vectors::queryHits(pp),], subject[S4Vectors::subjectHits(pp),])
          nrow(df.p)
          df.p <- tidyr::as_tibble(df.p)

          ####  CDS to query All (i.e. include complete and partial overlap)
          # Step2
          oo.1 = GenomicRanges::findOverlaps(query,subject,ignore.strand=F)
          df.o.1 = data.frame(query[S4Vectors::queryHits(oo.1),], subject[S4Vectors::subjectHits(oo.1),])
          nrow(df.o.1)

          #### overlapping CDS to query within gene body

          # search peak in bed
          # Step3A
          oo.2 = GenomicRanges::findOverlaps(query,subject,ignore.strand=T,select="all",type="within")
          df.o.2 = data.frame(query[S4Vectors::queryHits(oo.2),], subject[S4Vectors::subjectHits(oo.2),])

          # search bed regions in peaks (useful for broad peaks)
          # Step3B
          oo.3 <- GenomicRanges::findOverlaps(subject,query,ignore.strand=T,select="all",type="within")
          df.o.3 <- data.frame(query[S4Vectors::subjectHits(oo.3),], subject[S4Vectors::queryHits(oo.3),])

          df.o.2 <- rbind(df.o.2, df.o.3)

          # remove redundant from Step3A and Step3B
          df.o.2 <- unique(df.o.2)
          message("Number of peaks with complete overlap with the coding region")
          print(nrow(df.o.2))

          ########################
          # subset partial overlap from complete overlap (Step2 and Step3)

          partial.overlap <- subset(df.o.1, !(df.o.1[,6] %in% df.o.2[,6]))
          message("Number of peaks with Partial Overlap")
          print(nrow(partial.overlap))

          ### Filter Peak to get nearest gene
          ### Assign two genes to the partial overlapping peaks
          ### Filter targets having partial overlap at 3'end
          ### Depending upon the distance denote the peaks as PartialOverlap_Gene_1,2,3
          # Step4
          partial.overlap <- tidyr::as_tibble(partial.overlap)
          pp <-  partial.overlap %>%
                    dplyr::mutate(distance =if_else(strand.1=="-", start-end.1, start.1-end)) %>%
                    dplyr::filter(abs(distance)< width) %>%
                    dplyr::group_by(peak_id)  %>%
                    dplyr::mutate(min = order(distance)) %>%
                    dplyr::mutate(Category = paste("PartialOverlap_Gene",min , sep= "_")) %>%
                    dplyr::select(-c(min))


          ### Find the peaks within gene body
          ### Assign gene terminal if the ratio of distance > 0.75 i.e. if the peak is lying at extreme 3' end
          # Step5 (from Step3)

          df.o.2 = as_tibble(df.o.2)
          pGeneBody <- df.o.2 %>%
                    dplyr::mutate(distance =if_else(strand.1=="-",  end.1-end,start-start.1), Ratio = (abs(distance)/width.1), Category=if_else((abs(distance)/width.1)>0.75,"AtGeneEnd","WithInGeneBody")) %>%
                    dplyr::select(-c(Ratio))

          ### Find peaks in promoter region
          ### filter the peaks based on distance

          # Step6 (From Step1)
          pPrecede <- df.p %>%
                    dplyr::mutate(distance =if_else(strand.1=="-", start-end.1, start.1-end),Category="InPromoter")  %>%
                    group_by(peak_id) %>%
                    dplyr::filter(distance==min(distance))

          pPromoter <- pPrecede

          ### Remove peaks from  Promoter if they are already showing partial overlap
          # Step7 (from Step6)
          pPrecede <- dplyr::anti_join(pPrecede, pp, by = "peak_id")

          ### Remove peaks from  Promoter if they are already showing overlap with GeneBody
          # Step8 (from Step5)
          pPrecede <- dplyr::anti_join(pPrecede, pGeneBody, by = "peak_id")

          ### All possible targets for each peak
          # Step6, 7, 8
          df.All <- dplyr::bind_rows(pp,pGeneBody, pPrecede)

          ### Peaks with distance > 2000 in promoter, assign them to nearest gene irrespective of strand (ie --> Peak <-- <-- -->Target(according To Strand))
          # Step9
          df.LongPromoters <- df.All %>%
                    dplyr::mutate(Condition=if_else((distance>2000 & Category=="InPromoter"), "1","0")) %>%
                    dplyr::filter(Condition=="1")

          qq.subsetLongPromoters <-  subset(query, query$peak_id %in% df.LongPromoters$peak_id)

          nn <- GenomicRanges::nearest(qq.subsetLongPromoters,subject,ignore.strand=T,select="all")
          df.near <- data.frame(qq.subsetLongPromoters[S4Vectors::queryHits(nn),], subject[S4Vectors::subjectHits(nn),])
          nrow(df.near)
          df.near <-  tidyr::as_tibble(df.near)
          df.near.1 <- df.near %>%
                    dplyr::mutate(distance =if_else(start > end.1, start-end.1, start.1-end),Category="NearestTarget")

          ### Filter peaks in Promoter if distance > 2000bp from the all targets
          # Step10
          df.All <- dplyr::anti_join(df.All, df.near.1, by="peak_id")

          df.All <- dplyr::bind_rows(df.All, df.near.1)


          ### Filter the peaks AtGeneEnd if they have promoter of next gene nearby
          # Step11
          df.filterGeneTerminal <- dplyr::bind_rows(pGeneBody, pPromoter)  %>%
                    dplyr::filter(Category=="InPromoter"|Category=="AtGeneEnd") %>%
                    dplyr::group_by(peak_id) %>%
                    dplyr::filter(n()>1) %>%
                    dplyr::filter(distance==min(distance))

          ### All unique targets
          #Step12
          pPrecede <- dplyr::anti_join( pPrecede,df.near.1, by="peak_id")
          pGeneBody <- dplyr::anti_join( pGeneBody,df.filterGeneTerminal, by="peak_id")
          df.Unique <- dplyr::bind_rows(df.filterGeneTerminal,pPrecede, pp, pGeneBody, df.near.1) %>% base::unique(.)



          message("Number of Unique peaks")
          print(nrow(df.Unique))

          if(outfile==TRUE){
                    xlsx::write.xlsx(x=as.data.frame(df.Unique), file = paste(outfileName,"annotate_peaks.xlsx", sep="_"), sheetName = "UniqueTargets", append = TRUE,row.names = FALSE)
          }
          else{
                    return(df.Unique)
          }

                           tt = df.Unique %>%
                                     dplyr::mutate(Category=replace(Category,stringr::str_detect(Category, "PartialOverlap_*"), "PartialOverlap")) %>%
                                     dplyr::ungroup() %>%
                                     dplyr::select(Category) %>% dplyr::group_by(Category) %>% dplyr::summarise(count=n())

                           gg = ggplot2::ggplot(tt, ggplot2::aes("",y=count, fill=Category))+
                                     ggplot2::geom_bar(width = 1, stat = "identity")+
                                     ggrepel::geom_text_repel(ggplot2::aes(label = count), position = ggplot2::position_stack(vjust = 0.5),size=5,hjust = 1)+
                                     ggplot2::coord_polar(theta = "y")+
                                     ggplot2::ylab("No. of Peaks")+ggplot2::xlab("")+
                                     ggplot2::ggtitle(outfileName)+
                                     ggplot2::guides(fill = ggplot2::guide_legend(title = ""),color = ggplot2::guide_legend(title = ""))+
                                     ggplot2::theme_void()
                           print(gg)


}
