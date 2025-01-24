#library("usethis")
#usethis::create_project("/scratch/miller.me/CellChat_DB/build-axolotl-DB")

# Source functions for building custom axolotl CellChat DB

# prep_gene_list - build gene lists for downstream database filtering
#
# REQUIRED:
# objects - list of input Seurat objects as RDS files
#
# OPTIONS:
# assay - one of "SCT" or "RNA" (default "RNA")
# contains.spaces - do the gene names in the input gene lists contain spaces? (default TRUE)
# write.outfiles - write gene list to selected directory (default TRUE)
# out.dir - directory to output written files to (default working directory)

prep_gene_list <- function(objects, assay = NULL, contains.spaces = TRUE,
                           write.outfile = TRUE, out.dir = NULL){
  object.list <- objects
  assay <- assay

  writeDir <- out.dir
  if(is.null(writeDir)){
    writeDir <- "./"
  }

  gene.list <- c()
  for(object in 1:length(object.list)){
    input <- object.list[[object]]

    if(is.null(assay)){
      genes.add <- rownames(input@assays[["RNA"]]@features)
      gene.list <- c(gene.list, genes.add)
    } else if(!is.null(assay)){
      if(assay == "RNA"){
        genes.add <- rownames(input@assays[["RNA"]]@features)
        gene.list <- c(gene.list, genes.add)}
      else if(assay == "SCT"){
          genes.add <- input@assays[["SCT"]]@counts@Dimnames[[1]]
          gene.list <- c(gene.list, genes.add)}
      else if(!assay %in% c("SCT", "RNA")){
            stop("The input assay must be one of either RNA or SCT")
    }}

    if(isFALSE(contains.spaces)){
      gene.list <- as.data.frame(gene.list) %>%
        `colnames<-`(c("Genes")) %>%
        mutate(Genes = gsub(" ","",Genes)) %>%
        unique()
      if(isTRUE(write.outfile)){
        write_delim(gene.list,
                    file.path(writeDir, "Genes_for_filter.txt"))
      }
      gene.list <- as.list(gene.list)[[1]]

    } else{
      gene.list <- gene.list %>% unique()

      if(isTRUE(write.outfile)){
        write_delim(as.data.frame(gene.list) %>% `colnames<-`(c("Genes")),
                    file.path(writeDir, "Genes_for_filter.txt"))
      }
    }
    return(gene.list)
  }
}

# prep_gene_map - prepare annotations for mapping axolotl to human for use in downstream mapping
# Alternatively, users can provide their own gene - annotation mapping downstream.
#
#
# REQUIRED:
# gxt.map - tab-delimited file mapping genes to transcripts OR data.frame with columns gene & transcript
# example: AMEX60DD000001 ZFP37[nr]-ZNF568[hs]-AMEX60DD301000001.2
#
# OPTIONS:
# map.header - does the input gene-transcript map have a header line (default FALSE)
# human.priority - prioritize human annotations over non-human (default FALSE)
# drop.nonhuman - ignore any non-human annotations and drop from output (default FALSE)
# contains.spaces - do the genes output by prep_gene_list contain spaces? (default TRUE)
# write.outfiles - write tab-delimited annotation file (default TRUE)
# out.dir - directory to output written files to (default working directory)

# Honestly - might delete all together and require user input of
# CSV of full-transcript / human / other
# e.g.
# ZFP37[nr]|ZNF568[hs]|AMEX60DD301000001.1,ZNF568,ZFP37
prep_gene_map <- function(gxt.map, has.header = FALSE, tx.delim = NULL,
                          human.priority = FALSE, drop.nonhuman = FALSE,
                          contains.spaces = TRUE,
                          write.outfiles = TRUE, out.dir = NULL){

  writeDir <- out.dir
  if(is.null(writeDir)){
    writeDir <- "./"
  }

  if(is.object(gxt.map)){
    gxtMap <- gxt.map %>%
      `colnames<-`(c("geneID", "transcript"))
  } else{
    gxtMap <- read.delim(gxt.map, header = has.header) %>%
      `colnames<-`(c("geneID", "transcript"))
    #if(isTRUE(map.header)){
    #  gxtMap <- read.delim(gxt.map, header = TRUE) %>%
    #    `colnames<-`(c("geneID", "transcript"))
    #} else{
    #  gxtMap <- read.delim(gxt.map, header = FALSE) %>%
    #    `colnames<-`(c("geneID", "transcript"))
    #}
  }

  if(is.null(tx.delim)){
    tx.delim <- "-"
  }

  # Filter out rows with NO annotation - consider just building into below manipulation
  # I don't know that this is adding anything - the counts files DO contain non-annotated genes,
  # but will not contribute to mapping anyway
  # gxtMap <- gxtMap %>%
  #  filter(str_detect(transcript, tx.delim))

  if(isFALSE(contains.spaces)){
    gxtMap <- gxtMap %>% mutate(transcript = gsub(" ", "", transcript))
  }

  # Now need to work in actual manipulation based on preservation / priority of human genes



  return(gxtMap)

}

# Not yet started ------

# prep_synonym_map

prep_synonym_map <- function(){}

# make_complex_input

make_complex_input <- function(from.CC = CellChatDB.human$complex){}

# make_cofactor_input

make_cofactor_input <- function(from.CC = CellChatDB.human$cofactor){}

# make_interaction_input

make_interaction_input <- function(from.CC = CellChatDB.human$interaction){}

# make_geneInfo_input

make_geneInfo_input <- function(from.CC = CellChatDB.human$geneInfo){}

# assemble_database
