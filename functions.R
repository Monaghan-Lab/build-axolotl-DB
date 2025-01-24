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
# contains.spaces - do input feature names contain spaces? (default TRUE)
# write.outfiles - write output to out.dir? (default TRUE)
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

# prep_synonym_map
#
# OPTIONS:
# write.outfile - write output to out.dir? (default TRUE)
# out.dir - directory to output written files to (default working directory)

prep_synonym_map <- function(write.outfile = TRUE, out.dir = NULL){

  geneInfo_input <- CellChatDB.human$geneInfo

  writeDir <- out.dir
  if(is.null(writeDir)){
    writeDir <- "./"
  }

  synonyms <- geneInfo_input %>%
    select(Symbol, Synonym) %>%
    separate_wider_delim(., Synonym,
                         " ", names_sep = "_", too_few = "align_start") %>%
    separate_longer_delim(., Symbol, "; ") # Some are separated by ;, creating new rows for each

  synonyms[] <- lapply(synonyms, function(x) gsub(";[[:space:]]?", "", x))
  synonyms[synonyms == ""] <- NA # ex. RPL9 - "synonyms" ARE ";"

  synonym_count = rowSums(!(is.na(synonyms[,-1])))
  synonyms$count_na <- synonym_count
  synonyms_longform <- synonyms %>%
    filter(count_na != 0) %>% # If a symbol doesn't have any synonyms, it won't help
    select(-c(count_na))  %>%
    as.data.frame() %>%
    pivot_longer(cols = starts_with("Synonym") , values_drop_na = T) %>%
    select(-2) %>%
    unique() %>%
    mutate(Symbol = toupper(Symbol), value = toupper(value))

  if(isTRUE(write.outfile)){
    write.csv(synonyms_longform,
              file.path(writeDir, "Synonym_map.csv"))
  }

  return(synonyms_longform)
}

# make_complex_input
#
# REQUIRED:
# annotation.map - an input data frame of the form full-transcript / human-annotation / other-annotation
# example: ZFP37[nr]-ZNF568[hs]-AMEX60DD301000001.2 / ZNF568 / ZFP37
# filter.genes - the list of genes from prep_gene_list used to filter output
#
# OPTIONS:
# priority - give priority to either human or other annotations, or no priority (default "none")
# filter - how strictly should annotations be filtered using filter.genes? (default "strict")
# With "strict", ALL annotated subunits of a complex must be in filter.genes; "relaxed" requires at least one.
# annotation.cutoff - filter complexes with fewer than this fraction of annotated genes remaining
# i.e. if cutoff = 0.5 (the default), complexes with fewer than 50% of their original # of components annotated are removed
# check.synonyms - attempt to rescue complexes by checking for annotated synonyms? (default FALSE)
  # synonym.map - REQUIRED if check.synonyms = TRUE. Symbol-Synonym map output from prep_synonym_map()
# contains.spaces - do the input feature names contain spaces? (default TRUE)
# write.outfile - write output to out.dir? (default TRUE)
# out.dir - directory to output written files to (default working directory)


make_complex_input <- function(annotation.map = NULL, filter.genes = NULL,
                               priority = c("human", "other", "none"), filter = c("strict", "relaxed", "none"), annotation.cutoff = NULL,
                               check.synonyms = FALSE, synonym.map = NULL,
                               contains.spaces = TRUE,
                               write.outfiles = TRUE, out.dir = NULL){
  complex_input <- CellChatDB.human$complex

  writeDir <- out.dir
  if(is.null(writeDir)){
    writeDir <- "./"
  }



}

# Not yet started ------

# make_cofactor_input
#
# REQUIRED:
# annotation.map - an input data frame of the form full-transcript / human-annotation / other-annotation
# example: ZFP37[nr]-ZNF568[hs]-AMEX60DD301000001.2 / ZNF568 / ZFP37
# filter.genes - the list of genes from prep_gene_list used to filter output
#
# OPTIONS:
# priority - give priority to either human or other annotations, or no priority (default "none")
# filter - how strictly should annotations be filtered using filter.genes? (default "strict")
  # With "strict", ALL annotated components must be in filter.genes; "relaxed" requires at least one.
  # Note: due to limitations with file length handling in R, there is no option for NO filtering here.
# annotation.cutoff - filter cofactors with fewer than this fraction of annotated components remaining
# i.e. if cutoff = 0.5 (the default), cofactors with fewer than 50% of their original # of components annotated are removed
# check.synonyms - attempt to rescue cofactors by checking for annotated synonyms? (default FALSE)
  # synonym.map - REQUIRED if check.synonyms = TRUE. Symbol-Synonym map output from prep_synonym_map()
# contains.spaces - do the input feature names contain spaces? (default TRUE)
# write.outfile - write output to out.dir? (default TRUE)
# out.dir - directory to output written files to (default working directory)

make_cofactor_input <- function(annotation.map = NULL, filter.genes = NULL,
                                priority = c("human", "other", "none"), filter = c("strict", "relaxed"), annotation.cutoff = NULL,
                                check.synonyms = FALSE, synonym.map = NULL,
                                contains.spaces = TRUE,
                                write.outfiles = TRUE, out.dir = NULL){}

# make_interaction_input
# NOTE: BOTH make_complex_input and make_cofactor_input MUST be run before this step
#
# REQUIRED:
# annotation.map - an input data frame of the form full-transcript / human-annotation / other-annotation
# example: ZFP37[nr]-ZNF568[hs]-AMEX60DD301000001.2 / ZNF568 / ZFP37
# filter.genes - the list of genes from prep_gene_list used to filter output
# complexes - output of make_complex_input
# cofactors - output of make_cofactor_input
#
#
# OPTIONS:
# priority - give priority to either human or other annotations, or no priority (default "none")
# filter - how strictly should annotations be filtered using filter.genes? (default "strict")
# With "strict", ALL annotated components must be in filter.genes; "relaxed" requires at least one.
# check.synonyms - attempt to rescue interactions by checking for annotated synonyms? (default FALSE)
# synonym.map - REQUIRED if check.synonyms = TRUE. Symbol-Synonym map output from prep_synonym_map()
# contains.spaces - do the input feature names contain spaces? (default TRUE)
# write.outfile - write output to out.dir? (default TRUE)
# out.dir - directory to output written files to (default working directory)

make_interaction_input <- function(){}

# make_geneInfo_input
#
# REQUIRED:
# annotation.map - an input data frame of the form full-transcript / human-annotation / other-annotation
# example: ZFP37[nr]-ZNF568[hs]-AMEX60DD301000001.2 / ZNF568 / ZFP37
#
# OPTIONS:
# priority - give priority to either human or other annotations, or no priority (default "none")
# filter - should annotated symbols be limited to those present in output from prep_gene_list? (default FALSE)
# filter.genes - REQUIRED if  filter = TRUE. the list of genes from prep_gene_list used to filter output
# check.synonyms - attempt to rescue symbols by checking for annotated synonyms? (default FALSE)
# synonym.map - REQUIRED if check.synonyms = TRUE. Symbol-Synonym map output from prep_synonym_map()
# contains.spaces - do the input feature names contain spaces? (default TRUE)
# write.outfile - write output to out.dir? (default TRUE)
# out.dir - directory to output written files to (default working directory)

make_geneInfo_input <- function(annotation.map = NULL, priority = c("human", "other", "none"),
                                filter = FALSE, filter.genes = NULL,
                                check.synonyms = FALSE, synonym.map = NULL,
                                contains.spaces = TRUE,
                                write.outfile = TRUE, out.dir = NULL){}

# assemble_database


# Prep gene map - probably delete -----
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
