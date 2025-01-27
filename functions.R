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
                               write.outfile = TRUE, out.dir = NULL){

  colnames(annotation.map) <- c("full", "human", "other")
  if(!any(is.na(annotation.map))){
    complex_input[complex_input == ""] <- NA
  }
  emptyLine <- data.frame(full="", human="", other="") # empty line for parsing
  annotation.map <- rbind(annotation.map, emptyLine)

  writeDir <- out.dir
  if(is.null(writeDir)){
    writeDir <- "./"
  }

  complex_input <- CellChatDB.human$complex
  complex_input[complex_input == ""] <- NA
  complex_count = rowSums(!(is.na(complex_input)))
  complex_input$start_na <- complex_count # For later to determine what valid "missingness" is to keep after mapping

  complex_output_1 <- complex_input %>%
    rownames_to_column("rowname") %>%
    mutate(rowname = case_when(rowname %in% c("RA-ALDH1A1", "RA-ALDH1A2", "RA-ALDH1A3") ~ gsub("RA-","", rowname),
                               TRUE ~ rowname)) %>% # To be able to match to interactions
    mutate(subunit_2 = case_when(subunit_2 == "ALDH1A1" ~ NA,
                                 TRUE ~ subunit_2), # ONLY for ALDH1A1 is there two "subunits". ALDH1A2 and ALDH1A3 show only one.
           start_na = case_when(rowname == "ALDH1A1" ~ 1,
                                TRUE ~ start_na)) %>%
    left_join(annotation.map, by=c("subunit_1" = "human"), relationship = "many-to-many") %>%
    rename(subunit_1a = full) %>%
    left_join(annotation.map, by=c("subunit_2" = "human"),relationship = "many-to-many") %>%
    rename(subunit_2a = full)  %>%
    left_join(annotation.map, by=c("subunit_3" = "human"),relationship = "many-to-many") %>%
    rename(subunit_3a = full) %>%
    left_join(annotation.map, by=c("subunit_4" = "human"),relationship = "many-to-many") %>%
    rename(subunit_4a = full) %>%
    left_join(annotation.map, by=c("subunit_5" = "human"),relationship = "many-to-many") %>%
    rename(subunit_5a = full) %>%
    left_join(annotation.map, by=c("subunit_1" = "other"), relationship = "many-to-many") %>%
    rename(subunit_1b = full) %>%
    left_join(annotation.map, by=c("subunit_2" = "other"),relationship = "many-to-many") %>%
    rename(subunit_2b = full)  %>%
    left_join(annotation.map, by=c("subunit_3" = "other"),relationship = "many-to-many") %>%
    rename(subunit_3b = full) %>%
    left_join(annotation.map, by=c("subunit_4" = "other"),relationship = "many-to-many") %>%
    rename(subunit_4b = full) %>%
    left_join(annotation.map, by=c("subunit_5" = "other"),relationship = "many-to-many") %>%
    rename(subunit_5b = full)

  priority <- match.arg(priority)
  if(missing(priority)){
    priority <- "none"
  }

  if(priority == "human"){
    complex_output_1 <- complex_output_1 %>%
    mutate(final_1 = case_when(!is.na(subunit_1a) & is.na(subunit_1b) ~ subunit_1a,
                               !is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1a,
                               is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1b,
                               is.na(subunit_1a) & is.na(subunit_1b) ~ NA),
           final_2 = case_when(!is.na(subunit_2a) & is.na(subunit_2b) ~ subunit_2a,
                               !is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2a,
                               is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2b,
                               is.na(subunit_2) ~ NA,
                               is.na(subunit_2a) & is.na(subunit_2b) ~ NA),
           final_3 = case_when(!is.na(subunit_3a) & is.na(subunit_3b) ~ subunit_3a,
                               !is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3a,
                               is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3b,
                               is.na(subunit_3) ~ NA,
                               is.na(subunit_3a) & is.na(subunit_3b) ~ NA),
           final_4 = case_when(!is.na(subunit_4a) & is.na(subunit_4b) ~ subunit_4a,
                               !is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4a,
                               is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4b,
                               is.na(subunit_4) ~ NA,
                               is.na(subunit_4a) & is.na(subunit_4b) ~ NA),
           final_5 = case_when(!is.na(subunit_5a) & is.na(subunit_5b) ~ subunit_5a,
                               !is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5a,
                               is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5b,
                               is.na(subunit_5) ~ NA,
                               is.na(subunit_5a) & is.na(subunit_5b) ~ NA)) %>%
    select(c("rowname", "start_na", starts_with("final"))) %>%
      unique() %>%
      `colnames<-`(c("rowname", "start_na", "subunit_1", "subunit_2", "subunit_3", "subunit_4", "subunit_5")) } else{
        complex_output_1 <- complex_output_1 %>%
          mutate(final_1 = case_when(!is.na(subunit_1a) & is.na(subunit_1b) ~ subunit_1a,
                                     !is.na(subunit_1a) & !is.na(subunit_1b) ~ paste0(subunit_1a,";",subunit_1b),
                                     is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1b,
                                     is.na(subunit_1a) & is.na(subunit_1b) ~ NA),
                 final_2 = case_when(!is.na(subunit_2a) & is.na(subunit_2b) ~ subunit_2a,
                                     !is.na(subunit_2a) & !is.na(subunit_2b) ~ paste0(subunit_2a,";",subunit_2b),
                                     is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2b,
                                     is.na(subunit_2) ~ NA,
                                     is.na(subunit_2a) & is.na(subunit_2b) ~ NA),
                 final_3 = case_when(!is.na(subunit_3a) & is.na(subunit_3b) ~ subunit_3a,
                                     !is.na(subunit_3a) & !is.na(subunit_3b) ~ paste0(subunit_3a,";",subunit_3b),
                                     is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3b,
                                     is.na(subunit_3) ~ NA,
                                     is.na(subunit_3a) & is.na(subunit_3b) ~ NA),
                 final_4 = case_when(!is.na(subunit_4a) & is.na(subunit_4b) ~ subunit_4a,
                                     !is.na(subunit_4a) & !is.na(subunit_4b) ~ paste0(subunit_4a,";",subunit_4b),
                                     is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4b,
                                     is.na(subunit_4) ~ NA,
                                     is.na(subunit_4a) & is.na(subunit_4b) ~ NA),
                 final_5 = case_when(!is.na(subunit_5a) & is.na(subunit_5b) ~ subunit_5a,
                                     !is.na(subunit_5a) & !is.na(subunit_5b) ~ paste0(subunit_5a,";",subunit_5b),
                                     is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5b,
                                     is.na(subunit_5) ~ NA,
                                     is.na(subunit_5a) & is.na(subunit_5b) ~ NA)) %>%
          select(c("rowname", "start_na", starts_with("final"))) %>%
          separate_longer_delim(final_1,";") %>%
          separate_longer_delim(final_2, ";") %>%
          separate_longer_delim(final_3,";") %>%
          separate_longer_delim(final_4, ";") %>%
          separate_longer_delim(final_5,";") %>%
          unique() %>%
          `colnames<-`(c("rowname", "start_na", "subunit_1", "subunit_2", "subunit_3", "subunit_4", "subunit_5"))
      }

  out_na = rowSums(!(is.na(complex_output_1[,-c(1,2)])))
  complex_output_1$end_na <- out_na

  complex_output_2 <- complex_output_1 %>%
    filter(end_na != 0) %>% # Filtering out rows with no annotations
    select(rowname, start_na, end_na, starts_with("subunit_")) %>%
    mutate(na_diff = (start_na - end_na)/start_na) %>%
    filter(na_diff < annotation.cutoff) %>% # Drop if more than 50% of subunits are not annotated
    select(-c(start_na, end_na, na_diff)) %>%
    unique()

  if(isFALSE(contains.spaces)){
    complex_output_2 <- complex_output_2 %>%
      mutate(subunit_1 = gsub(" \\[", "\\[", subunit_1),
             subunit_2 = case_when(is.na(subunit_2) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", subunit_2)),
             subunit_3 = case_when(is.na(subunit_3) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", subunit_3)),
             subunit_4 = case_when(is.na(subunit_4) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", subunit_4)),
             subunit_5 = case_when(is.na(subunit_5) ~ NA,
                                   TRUE ~ gsub(" \\[","\\[", subunit_5)))
  }

  filter.level <- match.arg(filter)
  if(missing(filter.level)){
    filter.level <- "strict"
  }

  # TO DO: is there any way to make filtering smarter? Less cluttered?
  if(filter.level == "none"){
    complex_output_3 <- complex_output_2
  } else if(filter.level == "relaxed"){
    complex_output_3 <- complex_output_2 %>%
      filter(subunit_1 %in% filter.genes |
               subunit_2 %in% filter.genes |
               subunit_3 %in% filter.genes |
               subunit_4 %in% filter.genes |
               subunit_5 %in% filter.genes)
  } else if(filter.level == "strict"){
    complex_output_3 <- complex_output_2 %>%
      filter( ((!is.na(subunit_1) & is.na(subunit_2) & is.na(subunit_3) & is.na(subunit_4)) & (subunit_1 %in% filter.genes)) |
                ((!is.na(subunit_1) & !is.na(subunit_2) & is.na(subunit_3) & is.na(subunit_4)) & (subunit_1 %in% filter.genes & subunit_2 %in% filter.genes)) |
                ((!is.na(subunit_1) & !is.na(subunit_2) & !is.na(subunit_3) & is.na(subunit_4)) & (subunit_1 %in% filter.genes & subunit_2 %in% filter.genes & subunit_3 %in% filter.genes)) |
                ((!is.na(subunit_1) & !is.na(subunit_2) & !is.na(subunit_3) & !is.na(subunit_4)) & (subunit_1 %in% filter.genes & subunit_2 %in% filter.genes & subunit_3 %in% filter.genes & subunit_4 %in% filter.genes)) |
                ((!is.na(subunit_1) & !is.na(subunit_2) & !is.na(subunit_3) & !is.na(subunit_4) & !is.na(subunit_5)) & (subunit_1 %in% filter.genes & subunit_2 %in% filter.genes & subunit_3 %in% filter.genes & subunit_4 %in% filter.genes & subunit_5 %in% filter.genes)) )
  }

  if(isFALSE(check.synonyms)){
    complex_output <- complex_output_3
    fix_rownames <- complex_output
    fix_rownames$new_rowname <- make.unique(fix_rownames$rowname, sep = ".")
    complex_gene_match <- fix_rownames %>% select(rowname, new_rowname) # For later use with interaction input
    fix_rownames$rowname <- make.unique(fix_rownames$rowname)

    complex_output <- fix_rownames %>%
      column_to_rownames("new_rowname") %>%
      select(-rowname)

    if(isTRUE(write.outfile)){
      write.csv(complex_output, file.path(writeDir, "complex_input_axolotl.csv"))
    }
  }

  # TO DO: Add synonym testing + related complex_gene_match with writing outfile
  if(isTRUE(check.synonyms)){
    if(missing(synonym.map)){
      stop("You must provide a synonym map in order to test for synonyms")
    } else{
      synonym.map <- synonym.map
    }
  }

return(complex_output)
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
                                write.outfile = TRUE, out.dir = NULL){}

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
# OPTIONS:
# priority - give priority to either human or other annotations, or no priority (default "none")
# filter - how strictly should annotations be filtered using filter.genes? (default "strict")
# With "strict", ALL annotated SIMPLE components (ligands and receptors) must be in filter.genes; "relaxed" requires at least one.
# check.synonyms - attempt to rescue interactions by checking for annotated synonyms? (default FALSE)
# synonym.map - REQUIRED if check.synonyms = TRUE. Symbol-Synonym map output from prep_synonym_map()
# contains.spaces - do the input feature names contain spaces? (default TRUE)
# write.outfile - write output to out.dir? (default TRUE)
# out.dir - directory to output written files to (default working directory)

make_interaction_input <- function(annotation.map = NULL, complexes = NULL, cofactors = NULL, priority = c("human", "other", "none"),
                                   filter = c("strict", "relaxed", "none"), filter.genes = NULL,
                                   check.synonyms = FALSE, synonym.map = NULL,
                                   contains.spaces = TRUE,
                                   write.outfile = TRUE, out.dir = NULL){

  interaction_input <- interaction_input %>%
    mutate(ligand = case_when(pathway_name == "RA" ~ gsub("RA-","", ligand),
                              TRUE ~ ligand))
}

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

# Internal functions ------

complex_synonyms_internal <- function(){}

cofactor_synonyms_internal <- function(){}

interaction_synonyms_internal <- function(){}
