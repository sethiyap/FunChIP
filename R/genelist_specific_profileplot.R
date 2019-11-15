#' profile plot specific to a genelist
#'
#' @description Profile plot for list of genes to be subsetted from genome feature file.
#'
#' @param feature_file A character vector containing path to the file. Either a
#'   bed file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand} OR
#'   gene feature file with extention \code{.gff or .gff3 }
#' @param genelist A list of genes to be subsetted from the feature file
#' @param bw_files_dir A character vector containing path to files \code{.bw or
#'   .bdg}
#' @param bw_pattern A character vector, to replace a pattern from bw or bdg
#'   file while plotting. Useful when file names are extremely big.
#'   \code{default: NULL}
#' @param max_scale_limit Numeric, maximum limit of the color key. \code{default:60}
#' @param output Logical, to return a output of profile plot or not.
#'   \code{default: TRUE, with output_name,if FALSE; plot will be returned} **
#'   Highly recommended to plot in outfile when .bw files are more than 2.
#' @param promoter Logical, if TRUE plot 1kb up/downstream of start co-ordinate;
#'   if FALSE plot 1kb up/downstream of gene-body \code{default: TRUE}
#' @param output_name A character vector containing name of the output plot.\code{default:Sample}
#'
#'
#' @return Parallel plots of profiles for given gene lists
#' @export
#' @import magrittr
#' @import png
#' @importFrom stringr str_detect
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicRanges mcols
#' @importFrom rtracklayer import.bed
#' @importFrom GenomicFeatures promoters
#' @importFrom broom tidy
#' @importFrom dplyr mutate
#' @importFrom purrr map
#' @import EnrichedHeatmap
#' @import ComplexHeatmap
#' @import grid unit
#'
#'
#' @examples
#' \dontrun{
#'
#' feature_file <- "A_nidulans_FGSC_A4_version_s10-m04-r07_cds.bed"
#' gene_list <- clipr::read_clip()
#' bw_files_dir <- "/Users/Pooja/Documents/Data-Analysis/Others/FangFang/2019/ChIP_Analysis"
#' pattern <- "*_normalized.bw"
#' genelist_specific_profileplot(feature_file, genelist = gene_list,bw_files_dir = bw_files_dir, promoter = FALSE,bw_pattern = pattern)
#' }
#'
genelist_specific_profileplot <- function(feature_file, genelist, bw_files_dir,promoter=TRUE, bw_pattern=NULL, max_scale_limit=60, output=TRUE,output_name="Sample"){


          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".gff")){

                    print("Input is genome feature file")
                    gff = GenomicFeatures::makeTxDbFromGFF(feature_file)
                    subject = GenomicFeatures::genes(gff)
                    base::names(subject)=NULL
                    gene_name <- names(GenomicRanges::mcols(subject)[1])
                    feature_gr <- subset(subject,subject$gene_id %in% genelist)
          }
          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".bed")){
                    print("Input is bed file")
                    subject = rtracklayer::import.bed(feature_file)
                    gene_name <- names(GenomicRanges::mcols(subject)[1])
                    feature_gr <- subset(subject,subject$name %in% genelist)

          }
          if(promoter==TRUE){
                   tss = GenomicFeatures::promoters(feature_gr, upstream = 0, downstream = 1)
                   axis_name = c("-1kb","Start","+1kb")
          }
          else{
                  tss=feature_gr
                  axis_name = c("-1kb","Start","Stop","+1kb")
          }

          ## prepare signal data
          bw_files <- list.files(bw_files_dir, pattern = paste(c("*.bw", "*.bdg"),collapse = '|'), recursive = T, full.names = T)

          if(is.null(bw_pattern)==TRUE){
                    names(bw_files) <- gsub(pattern = paste(c("*.bw", "*.bdg"),collapse = '|'),replacement = "", basename(bw_files))
          }
          else{
                   # "_[[:upper:]]{6,}_.*_normalized.bw*"
                    names(bw_files) <- gsub(pattern = bw_pattern ,replacement = "", basename(bw_files))
          }

          bw_files <- broom::tidy(bw_files)
          print(bw_files)

          ## generate normalised matrix in tidy way


          xx <- bw_files %>%
                    dplyr::mutate(bw = purrr::map(x, function(ii) {
                              rtracklayer::import(ii)
                    })) %>%
                    dplyr::mutate(norm_matrix = purrr::map(bw, function(ii) {
                              nn <- EnrichedHeatmap::normalizeToMatrix(ii, tss, value_column = "score",background = 0,
                                                                       extend = c(1000),w = 50,
                                                                       smooth = TRUE)
                              nn[nn<0]=0
                              return(nn)
                    }))


          get_enrichment_heatmap_list <- function(x, names, titles, ...) {


                    ll <- length(x)

                    ## first heatmap
                    ehml <- EnrichedHeatmap::EnrichedHeatmap(mat = x[[1]], name = names[[1]], column_title = titles[[1]], show_heatmap_legend = T,
                                                             row_order = order(EnrichedHeatmap::enriched_score(x[[1]]), decreasing = TRUE),
                                            use_raster = TRUE, ...)

                    ## several other heatmaps if length of x > 1.
                    if (ll > 1) {
                              for (i in 2:ll) {
                                        print(i)
                                        ehml <- ehml +
                                                  EnrichedHeatmap::EnrichedHeatmap(
                                                            mat = x[[i]],
                                                            name = ifelse(length(names) >= i, names[i], "NA"),
                                                            use_raster = TRUE,
                                                            column_title = ifelse(length(titles) >= i, titles[i], "NA"),
                                                            show_heatmap_legend = ifelse(length(names) >= i, TRUE, FALSE), ...
                                                  ) ## legend will be shown only if the name is given for a heatmap.
                              }
                    }

                    return(ehml)
          }

          split_factor <- (max_scale_limit/6)

          ehm_list <- get_enrichment_heatmap_list(x = xx$norm_matrix,names = xx$names,titles = xx$names,
                                                  cluster_rows = FALSE,
                                                   pos_line = TRUE,
                                                  show_row_names = FALSE,
                                                  axis_name_rot = 90,
                                                  heatmap_legend_param = list(color_bar = "continuous",legend_direction="horizontal"),
                                                  axis_name = axis_name,
                                                  col = circlize::colorRamp2(breaks = seq(0,max_scale_limit, by = split_factor),
                                                                             colors = c("#ffffb2","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026")),
                                                  top_annotation = ComplexHeatmap::HeatmapAnnotation(lines = EnrichedHeatmap::anno_enriched(
                                                            yaxis_facing="inside",yaxis_gp=grid::gpar(fonsize=12))))

          if(output==TRUE){
                    pdf(file=paste(output_name, length(feature_gr), "hm.pdf", sep="_"), width=nrow(bw_files)*2.5, height=11)
                    ComplexHeatmap::draw(ehm_list, heatmap_legend_side = "bottom", gap = grid::unit(1.5, "mm"))
                    dev.off()
          }
          else{
                    return(ComplexHeatmap::draw(ehm_list, heatmap_legend_side = "right", gap = grid::unit(1.5, "mm")))
          }



}
