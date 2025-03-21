#' Source functions for building custom axolotl CellChat DB
#' TO DO: inherit params for complex, cofactor, interactions; effectively use all the same
#' TO DO: how necessary is contains_spaces after creating filter genes and annotation map?

#' prep_gene_list - build gene lists for downstream database filtering
#'
#' REQUIRED:
#' @param objects - list of input Seurat objects as RDS files
#'
#' OPTIONS:
#' @param assay - one of "SCT" or "RNA" (default "RNA")
#' @param contains.spaces - do the gene names in the input gene lists contain spaces? (default TRUE)
#' @param write.outfile - write gene list to selected directory (default TRUE)
#' @param out.dir - directory to output written files to (default working directory)

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

    if(isTRUE(contains.spaces)){
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

# prep_gene_map - prepare annotations for mapping axolotl to human
#
# REQUIRED:
# gxt.map - tab-delimited file mapping genes to transcripts with "|" separating annotations & transcript number
# e.g.: AMEX60DD000001 ZFP37[nr]|ZNF568[hs]|AMEX60DD301000001.2
# In the event of NO annotations:
# AMEX60DD000100  AMEX60DD201000100.1
#
# OPTIONS:
# apply.filter - restrict genes in map to only those present in gene.list? (default FALSE)
#   REQUIRES gene.list to be supplied - the output of prep_gene_list
# map.header - does the input gene-transcript map have a header line (default FALSE)
# prioritize.human - when both human and other annotations exist, select human as "official" annotation.
#   Otherwise, if no human annotation exists, will use the "other" annotation. (default FALSE - keeps both annotations)
# drop.nonhuman - supercedes 'prioritize.human'; ignore any non-human annotations and drop from output (default FALSE)
# contains.spaces - do the genes output by prep_gene_list contain spaces? (default TRUE)
# write.outfile - write tab-delimited annotation file (default TRUE)
# out.dir - directory to output written files to (default working directory)

prep_gene_map <- function(gxt.map,
                          apply.filter = FALSE, gene.list = NULL,
                          map.header = FALSE,
                          prioritize.human = FALSE, drop.nonhuman = FALSE,
                          contains.spaces = TRUE,
                          write.outfile = TRUE, out.dir = NULL){

  writeDir <- out.dir
  if(is.null(writeDir)){
    writeDir <- "./"
  }

  if(is.character(gxt.map)){
  if(isTRUE(map.header)){
  gxtMap <- read.delim(gxt.map, header = TRUE) %>%
    `colnames<-`(c("geneID", "transcript"))
  } else{
    gxtMap <- read.delim(gxt.map, header = FALSE) %>%
      `colnames<-`(c("geneID", "transcript"))
  }
  } else{
      stop("gxt.map must be a path to a tab-delimited gene-transcript map!")
  }

  if(nrow(gxtMap) == 0){
    stop("Check the format of the input annotation file")
  }

  gxtMap <- gxtMap %>%
    mutate(transcript = case_when(isFALSE(contains.spaces) ~ gsub(" ", "", transcript),
                                  TRUE ~ transcript)) %>%
    filter(str_detect(transcript, "\\|")) %>% # Filtering genes without any annotations
    mutate(full = transcript)

  gxtMap <- gxtMap %>%
    mutate(human = case_when(str_detect(full,"\\[", negate = T) ~ gsub("\\|.*","",full), # ex. "AEN|..."
                                     str_detect(full, "\\[hs\\]", negate = T) & str_detect(full, "\\[nr\\]") ~ "", # ex. LOC102367500[nr]|...
                                     str_detect(full, "\\[nr\\]", negate = T) & str_detect(full, "\\[hs\\]") ~ gsub("\\|.*","",full), # ex. TOP2A[hs]|...
                                     str_detect(full, "\\[nr\\]") & str_detect(full, "\\[hs\\]") ~ gsub(".*\\|(.*)\\|.*", "\\1", full), # ex. ZFP37[nr]|ZNF568[hs]|...
                                     TRUE ~ "TEST"), # To confirm all cases are caught
           other = case_when(str_detect(full,"\\[", negate = T) ~ gsub("\\|.*","", full), # ex. "AEN|..."
                                     str_detect(full, "\\[hs\\]", negate = T) & str_detect(full, "\\[nr\\]") ~ gsub("\\|.*","", full), # ex. LOC102367500[nr]|...
                                     str_detect(full, "\\[nr\\]", negate = T) & str_detect(full, "\\[hs\\]") ~ "", # ex. TOP2A[hs]|...
                                     str_detect(full, "\\[nr\\]") & str_detect(full, "\\[hs\\]") ~ gsub("\\|.*","",full), # ex. ZFP37[nr]|ZNF568[hs]|...
                                     TRUE ~ "TEST")) %>%
    select(full, human, other) %>%
   mutate(full = gsub("\\|", "-", full),
          human = gsub("\\[hs\\]","", human),
          other = gsub("\\[nr\\]","",other))  %>%
    mutate(full = gsub("_", "-", full),
           human = gsub("_","-", human),
           other = gsub("_","-", other)) %>%
    mutate(other = case_when(human == other ~ "",
                             TRUE ~ other)) %>%
    unique()

  if(isTRUE(apply.filter)){
    if(missing(gene.list)){
      stop("Please run prep_gene_list to prepare input for apply.filter.")
    } else{
    filter.genes <- as.data.frame(gene.list) %>% `colnames<-`(c("genes"))
    gxtMap <- gxtMap %>%
      filter(full %in% filter.genes$genes)
    }
  }

  check <- gxtMap[gxtMap$human == "TEST" | gxtMap$other == "TEST",]
  if(nrow(check) > 0){
    stop("Check the format of your input annotation file")
  }

  if(isTRUE(drop.nonhuman) & isFALSE(prioritize.human)){
    gxtMap <- gxtMap %>%
      filter(!human == "") %>%
      mutate(other = "")

    parseline <- c(full = "", human = "", other = "")
    gxtMap <- rbind(gxtMap, parseline)

  } else if(isTRUE(prioritize.human) & isFALSE(drop.nonhuman)){
    gxtMap <- gxtMap %>%
      mutate(human = case_when(human == "" & other != "" ~ other,
                                    TRUE ~ human)) %>%
      mutate(other = "")

    parseline <- c(full = "", human = "", other = "")
    gxtMap <- rbind(gxtMap, parseline)

  } else if(isTRUE(prioritize.human) & isTRUE(drop.nonhuman)){
    message("drop.nonhuman supercedes prioritize.human. All nonhuman annotations will be dropped from output.")
    gxtMap <- gxtMap %>%
      filter(!human == "") %>%
      mutate(other = "")

    parseline <- c(full = "", human = "", other = "")
    gxtMap <- rbind(gxtMap, parseline)

  } else{
    parseline <- c(full = "", human = "", other = "")
    gxtMap <- rbind(gxtMap, parseline)

  }

  if(isTRUE(write.outfile)){
    write_delim(gxtMap, file.path(writeDir, "Transcript_annotation_map.txt"))
  }

  return(gxtMap)
}

# prep_synonym_map - prepare Symbol-synonym map for testing synonymous genes in building inputs
#
# REQUIRED:
# gene.list - the output of prep_gene_list that will be used to filter annotated genes
# gene.map - the output of prep_gene_map that will be used to match human genes to axolotl
#
# OPTIONS:
# write.outfile - write tab-delimited long-form map of Symbol - Synonym (default TRUE)
# out.dir - directory to output written files to (default working directory)

prep_synonym_map <- function(gene.list, gene.map,
                             write.outfile = TRUE,
                             out.dir = NULL){

  writeDir <- out.dir
  if(is.null(writeDir)){
    writeDir <- "./"
  }

  if(missing(gene.list)){
    stop("You must provide a gene list output by prep_gene_list.")
  }

  synonyms <- CellChatDB.human$geneInfo %>%
    select(Symbol, Synonym) %>%
    separate_longer_delim(., Symbol, "; ") %>% # Some are separated by ;, creating new rows for each
    filter(str_count(Symbol, " ") == 0) %>%
    separate_wider_delim(., Synonym, " ", names_sep = "_", too_few = "align_start")
  synonyms[] <- lapply(synonyms, function(x) gsub(";[[:space:]]?", "", x))
  synonyms[synonyms == ""] <- NA # ex. RPL9 - "synonyms" ARE ";"

  synonym_count = rowSums(!(is.na(synonyms[,-1])))
  synonyms$count_na <- synonym_count
  synonyms <- synonyms %>%
    filter(count_na != 0) %>% # If a symbol doesn't have any synonyms, it won't help
    select(-c(count_na)) %>%
    as.data.frame() %>%
    pivot_longer(., cols = starts_with("Synonym"), values_drop_na = T) %>%
    select(-2) %>%
    unique() %>%
    mutate(Symbol = toupper(Symbol),
           value = toupper(value)) %>%
    filter(Symbol != value) %>%
    unique()

  # Input is not clean so trying to deal with a few special cases for completeness
  space_synonyms <- CellChatDB.human$geneInfo %>%
    select(Symbol, Synonym) %>%
    separate_longer_delim(., Symbol, "; ") %>% # Some are separated by ;, creating new rows for each
    filter(str_count(Symbol, " ") != 0) %>%
    filter(!is.na(Synonym))
  space_synonyms_not_symbol <- space_synonyms %>%
    filter(str_count(Synonym, " ") == 0) %>%
    mutate(Symbol = toupper(Symbol),
           Synonym = toupper(Synonym)) %>%
    pivot_longer(., cols = starts_with("Synonym"), values_drop_na = T) %>%
    select(-2) %>%
    unique()

  synonyms <- rbind(synonyms, space_synonyms_not_symbol)

  space_synonyms_symbol <- space_synonyms %>%
    filter(str_count(Synonym, " ") != 0) %>%
    mutate(Symbol = toupper(Symbol), Synonym = toupper(Synonym))

  bind <- space_synonyms_symbol[c(1:3, 5:8, 10),] %>% `colnames<-`(c("Symbol", "value"))
  synonyms <- rbind(synonyms, bind)

  # Only three remain, easy enough to fix manually
  space_synonyms_symbol <- space_synonyms_symbol[c(4,9,11),] %>%
    mutate(Synonym = gsub("-5", "-5;", Synonym)) %>%
    mutate(Synonym = gsub("FUSION", "FUSION;", Synonym)) %>%
    mutate(Synonym = gsub("-ETO", "-ETO;", Synonym)) %>%
    mutate(Synonym = gsub("-MTG8", "-MTG8;", Synonym)) %>%
    separate_longer_delim(., Synonym, "; ") %>%
    mutate(Synonym = gsub(";","", Synonym)) %>%
    `colnames<-`(c("Symbol", "value"))

  synonyms <- rbind(synonyms, space_synonyms_symbol) %>% unique()
  synonyms <- synonyms %>%
    filter(Symbol != value)

  filter.genes <- gene.list
  annotation.map <- gene.map

  synonym_map_new <- synonyms %>%
    left_join(annotation.map, by = c("value" = "human"), relationship = "many-to-many") %>%
    rename(humanAnnot = full) %>%
    left_join(annotation.map, by = c("value" = "other"), relationship = "many-to-many") %>%
    rename(otherAnnot = full) %>%
    left_join(annotation.map, by = c("Symbol" = "human"), relationship = "many-to-many") %>%
    rename(humanAnnot2 = full) %>%
    left_join(annotation.map, by = c("Symbol" = "other"), relationship = "many-to-many") %>%
    rename(otherAnnot2 = full) %>%
    select(-c(human.x, human.y, other.x, other.y)) %>%
    filter(humanAnnot %in% filter.genes | otherAnnot %in% filter.genes | humanAnnot2 %in% filter.genes | otherAnnot2 %in% filter.genes) %>%
    filter(!is.na(humanAnnot) | !is.na(otherAnnot)) %>%
    select(Symbol, value) %>%
    unique()
  synonyms <- synonym_map_new %>% select(Symbol, value) %>% unique()

  if(isTRUE(write.outfile)){
    write_delim(synonyms, file.path(writeDir, "Symbol_synonym_map.txt"))
  }
  return(synonyms)
}

# make_complex_input - prepare complex input of CellChat DB
# REQUIRED:
# gene.list - the output of prep_gene_list that will be used to filter annotated genes
# annotation.map - the output of prep_gene_map that will be used to match human genes to axolotl
#
# OPTIONS:
# priority - in the event that there are both human and other annotations for a gene, give priority to either or none. (default "none")
  # Note: the same priority will be applied if testing for synonyms.
# filter - how to apply filtering genes using gene.list. "Strict" (the default) requires that all annotated subunits appear in gene.list. "Relaxed" requires at least one.
  # Note: the same filter will be applied if testing for synonyms.
# annotation.cutoff - if fewer than this fraction of subunits are matched to axolotl, drop from input.
  # The default 0.5 means that rows will be dropped where < 50% of input subunits are matched to axolotl.
# drop.nonhuman - supercedes 'priority'; ignore any non-human annotations and drop from output (default FALSE)
# check.synonyms - attempt to rescue rows by also testing for synonymous genes? (default FALSE)
# synonyms - REQUIRED if check.synonyms is TRUE. The output of prep_synonym_map that will be used to test for synonyms.
# contains.spaces - do the genes output by prep_gene_list contain spaces? (default TRUE)
# write.outfile - write tab-delimited complex output file (default TRUE)
# out.dir - directory to output written file to (default working directory)

make_complex_input <- function(gene.list, annotation.map,
                               priority = c("human", "other", "none"), filter = c("strict", "relaxed"),
                               annotation.cutoff = NULL,
                               drop.nonhuman = FALSE, check.synonyms = FALSE, synonyms = NULL,
                               contains.spaces = TRUE,
                               write.outfile = TRUE, out.dir = NULL){

  if(missing(filter)){
    filter <- "strict"
  }
  if(!filter %in% c("strict", "relaxed")){
    stop("'filter' must be one of 'strict' or 'relaxed'.")
  }
  if(is.null(annotation.cutoff)){
    annotation.cutoff <- 0.5
  }
  if(!annotation.cutoff >= 0 & !annotation.cutoff <= 1 ){
    stop("'annotation.cutoff' must be a decimal between 0 and 1.")
  } else if(annotation.cutoff == 1){
    warning("An 'annotation.cutoff' of 1 may produce unintended results including an empty output.")
  }
  if(missing(priority)){
    priority <- "none"
  }
  if(!priority %in% c("human", "other", "none")){
    stop("'priority' must be one of 'human', 'other', or 'none'.")
  }

  complex_input <- CellChatDB.human$complex
  complex_input[complex_input == ""] <- NA
  complex_input <- complex_input %>%
    rownames_to_column("rowname") %>%
    mutate(rowname = case_when(rowname %in% c("RA-ALDH1A1", "RA-ALDH1A2", "RA-ALDH1A3") ~ gsub("RA-","", rowname),
                               TRUE ~ rowname),
           subunit_2 = case_when(subunit_2 == "ALDH1A1" ~ NA,
                                 TRUE ~ subunit_2))

  if(isTRUE(check.synonyms)){
    if(is.null(synonyms)){
      stop("You must provide a Symbol-Synonym map output by prep_synonym_map to search for synonymous genes.")
    }
    complex_input_synonyms <- complex_input %>%
      left_join(synonyms, by=c("subunit_1" = "Symbol"), relationship = "many-to-many") %>%
      rename(subunit_1a = value) %>%
      left_join(synonyms, by=c("subunit_2" = "Symbol"),relationship = "many-to-many") %>%
      rename(subunit_2a = value)  %>%
      left_join(synonyms, by=c("subunit_3" = "Symbol"),relationship = "many-to-many") %>%
      rename(subunit_3a = value) %>%
      left_join(synonyms, by=c("subunit_4" = "Symbol"),relationship = "many-to-many") %>%
      rename(subunit_4a = value) %>%
      left_join(synonyms, by=c("subunit_5" = "Symbol"),relationship = "many-to-many") %>%
      rename(subunit_5a = value) %>%
      mutate(subunit_1b = case_when(!is.na(subunit_1) & is.na(subunit_1a) ~ subunit_1,
                                    !is.na(subunit_1) & !is.na(subunit_1a) ~ paste0(subunit_1,";",subunit_1a),
                                    is.na(subunit_1) & is.na(subunit_1a) ~ NA,
                                    TRUE ~ "FIX"),
             subunit_2b = case_when(!is.na(subunit_2) & is.na(subunit_2a) ~ subunit_2,
                                    !is.na(subunit_2) & !is.na(subunit_2a) ~ paste0(subunit_2,";",subunit_2a),
                                    is.na(subunit_2) & is.na(subunit_2a) ~ NA,
                                    TRUE ~ "FIX"),
             subunit_3b = case_when(!is.na(subunit_3) & is.na(subunit_3a) ~ subunit_3,
                                    !is.na(subunit_3) & !is.na(subunit_3a) ~ paste0(subunit_3,";",subunit_3a),
                                    is.na(subunit_3) & is.na(subunit_3a) ~ NA,
                                    TRUE ~ "FIX"),
             subunit_4b = case_when(!is.na(subunit_4) & is.na(subunit_4a) ~ subunit_4,
                                    !is.na(subunit_4) & !is.na(subunit_4a) ~ paste0(subunit_4,";",subunit_4a),
                                    is.na(subunit_4) & is.na(subunit_4a) ~ NA,
                                    TRUE ~ "FIX"),
             subunit_5b = case_when(!is.na(subunit_5) & is.na(subunit_5a) ~ subunit_5,
                                    !is.na(subunit_5) & !is.na(subunit_5a) ~ paste0(subunit_5,";",subunit_5a),
                                    is.na(subunit_5) & is.na(subunit_5a) ~ NA,
                                    TRUE ~ "FIX")) %>%
      select(1, 12:16) %>%
      separate_longer_delim(subunit_1b, ";") %>%
      separate_longer_delim(subunit_2b, ";") %>%
      separate_longer_delim(subunit_3b, ";") %>%
      separate_longer_delim(subunit_4b, ";") %>%
      separate_longer_delim(subunit_5b, ";") %>%
      unique()
    colnames(complex_input_synonyms) <- colnames(complex_input)
    complex_count = rowSums(!(is.na(complex_input_synonyms[,-1])))
    complex_input_synonyms$start_na <- complex_count # For later to determine what valid "missingness" is to keep after mapping
    complex_input <- complex_input_synonyms
  } else{
    complex_count = rowSums(!(is.na(complex_input[,-1])))
    complex_input$start_na <- complex_count # For later to determine what valid "missingness" is to keep after mapping
  }

  if(isTRUE(drop.nonhuman)){
    message("drop.nonhuman supercedes annotation priority. Only human annotations will be searched for.")

    complex_output_1 <- complex_input %>% # Mapping subunits to axolotl genes
      left_join(annotation.map, by=c("subunit_1" = "human"), relationship = "many-to-many") %>%
      rename(final_1 = full) %>%
      left_join(annotation.map, by=c("subunit_2" = "human"),relationship = "many-to-many") %>%
      rename(final_2 = full)  %>%
      left_join(annotation.map, by=c("subunit_3" = "human"),relationship = "many-to-many") %>%
      rename(final_3 = full) %>%
      left_join(annotation.map, by=c("subunit_4" = "human"),relationship = "many-to-many") %>%
      rename(final_4 = full) %>%
      left_join(annotation.map, by=c("subunit_5" = "human"),relationship = "many-to-many") %>%
      rename(final_5 = full)

    complex_output_1 <- complex_output_1 %>%
      select(c("rowname", "start_na", starts_with("final"))) %>%
      unique() %>%
      `colnames<-`(c("rowname", "start_na", "subunit_1", "subunit_2", "subunit_3", "subunit_4", "subunit_5"))

    } else{
    complex_output_1 <- complex_input %>% # Mapping subunits to axolotl genes
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

    if(priority == "human"){
      complex_output_1 <- complex_output_1 %>%
        mutate(final_1 = case_when(!is.na(subunit_1a) & is.na(subunit_1b) ~ subunit_1a,
                                   !is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1a,
                                   is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1b,
                                   is.na(subunit_1a) & is.na(subunit_1b) ~ NA),
               final_2 = case_when(is.na(subunit_2) ~ NA,
                                   !is.na(subunit_2a) & is.na(subunit_2b) ~ subunit_2a,
                                   !is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2a,
                                   is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2b,
                                   is.na(subunit_2a) & is.na(subunit_2b) ~ NA),
               final_3 = case_when(is.na(subunit_3) ~ NA,
                                   !is.na(subunit_3a) & is.na(subunit_3b) ~ subunit_3a,
                                   !is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3a,
                                   is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3b,
                                   is.na(subunit_3a) & is.na(subunit_3b) ~ NA),
               final_4 = case_when(is.na(subunit_4) ~ NA,
                                   !is.na(subunit_4a) & is.na(subunit_4b) ~ subunit_4a,
                                   !is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4a,
                                   is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4b,
                                   is.na(subunit_4a) & is.na(subunit_4b) ~ NA),
               final_5 = case_when(is.na(subunit_5) ~ NA,
                                   !is.na(subunit_5a) & is.na(subunit_5b) ~ subunit_5a,
                                   !is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5a,
                                   is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5b,
                                   is.na(subunit_5a) & is.na(subunit_5b) ~ NA)) %>%
        select(c("rowname", "start_na", starts_with("final"))) %>%
        unique() %>%
        `colnames<-`(c("rowname", "start_na", "subunit_1", "subunit_2", "subunit_3", "subunit_4", "subunit_5"))

    } else if(priority == "other"){
      complex_output_1 <- complex_output_1 %>%
        mutate(final_1 = case_when(!is.na(subunit_1a) & is.na(subunit_1b) ~ subunit_1a,
                                   !is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1b,
                                   is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1b,
                                   is.na(subunit_1a) & is.na(subunit_1b) ~ NA),
               final_2 = case_when(is.na(subunit_2) ~ NA,
                                   !is.na(subunit_2a) & is.na(subunit_2b) ~ subunit_2a,
                                   !is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2b,
                                   is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2b,
                                   is.na(subunit_2a) & is.na(subunit_2b) ~ NA),
               final_3 = case_when(is.na(subunit_3) ~ NA,
                                   !is.na(subunit_3a) & is.na(subunit_3b) ~ subunit_3a,
                                   !is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3b,
                                   is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3b,
                                   is.na(subunit_3a) & is.na(subunit_3b) ~ NA),
               final_4 = case_when(is.na(subunit_4) ~ NA,
                                   !is.na(subunit_4a) & is.na(subunit_4b) ~ subunit_4a,
                                   !is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4b,
                                   is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4b,
                                   is.na(subunit_4a) & is.na(subunit_4b) ~ NA),
               final_5 = case_when(is.na(subunit_5) ~ NA,
                                   !is.na(subunit_5a) & is.na(subunit_5b) ~ subunit_5a,
                                   !is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5b,
                                   is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5b,
                                   is.na(subunit_5a) & is.na(subunit_5b) ~ NA)) %>%
        select(c("rowname", "start_na", starts_with("final"))) %>%
        unique() %>%
        `colnames<-`(c("rowname", "start_na", "subunit_1", "subunit_2", "subunit_3", "subunit_4", "subunit_5"))

    } else if(priority == "none"){
      complex_output_1 <- complex_output_1 %>%
        mutate(final_1 = case_when(!is.na(subunit_1a) & is.na(subunit_1b) ~ subunit_1a,
                                   !is.na(subunit_1a) & !is.na(subunit_1b) ~ paste0(subunit_1a,";",subunit_1b),
                                   is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1b,
                                   is.na(subunit_1a) & is.na(subunit_1b) ~ NA),
               final_2 = case_when(is.na(subunit_2) ~ NA,
                                   !is.na(subunit_2a) & is.na(subunit_2b) ~ subunit_2a,
                                   !is.na(subunit_2a) & !is.na(subunit_2b) ~ paste0(subunit_2a,";",subunit_2b),
                                   is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2b,
                                   is.na(subunit_2a) & is.na(subunit_2b) ~ NA),
               final_3 = case_when(is.na(subunit_3) ~ NA,
                                   !is.na(subunit_3a) & is.na(subunit_3b) ~ subunit_3a,
                                   !is.na(subunit_3a) & !is.na(subunit_3b) ~ paste0(subunit_3a,";",subunit_3b),
                                   is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3b,
                                   is.na(subunit_3a) & is.na(subunit_3b) ~ NA),
               final_4 = case_when(is.na(subunit_4) ~ NA,
                                   !is.na(subunit_4a) & is.na(subunit_4b) ~ subunit_4a,
                                   !is.na(subunit_4a) & !is.na(subunit_4b) ~ paste0(subunit_4a,";",subunit_4b),
                                   is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4b,
                                   is.na(subunit_4a) & is.na(subunit_4b) ~ NA),
               final_5 = case_when(is.na(subunit_5) ~ NA,
                                   !is.na(subunit_5a) & is.na(subunit_5b) ~ subunit_5a,
                                   !is.na(subunit_5a) & !is.na(subunit_5b) ~ paste0(subunit_5a,";",subunit_5b),
                                   is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5b,
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

  }

  out_na = rowSums(!(is.na(complex_output_1[,-c(1,2)])))
  complex_output_1$end_na <- out_na
  complex_output_1 <- complex_output_1 %>%
    select(rowname, start_na, end_na, starts_with("subunit_")) %>%
    mutate(na_diff = end_na/start_na) %>%
    filter(na_diff > annotation.cutoff) %>% # Drop if greater than X fraction of subunits are not annotated
    select(-c(start_na, end_na, na_diff)) %>%
    unique()

  if(isFALSE(contains.spaces)){
    complex_output_1 <- complex_output_1 %>%
      mutate(subunit_1 = gsub(" \\[", "\\[", subunit_1),
             subunit_2 = case_when(is.na(subunit_2) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", subunit_2)),
             subunit_3 = case_when(is.na(subunit_3) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", subunit_3)),
             subunit_4 = case_when(is.na(subunit_4) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", subunit_4)),
             subunit_5 = case_when(is.na(subunit_5) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", subunit_5)))
  }

  if(missing(gene.list)){
      stop("You must provide a gene list for filtering output by prep_gene_list.")
    }
  filter.genes <- as.data.frame(gene.list) %>% `colnames<-`(c("genes"))
  if(filter == "relaxed"){
    complex_output_filter <- complex_output_1 %>%
      filter(subunit_1 %in% filter.genes$genes |
               subunit_2 %in% filter.genes$genes |
               subunit_3 %in% filter.genes$genes |
               subunit_4 %in% filter.genes$genes |
               subunit_5 %in% filter.genes$genes)

  } else if(filter == "strict"){
    complex_output_filter <- complex_output_1 %>%
      filter( ((!is.na(subunit_1) & is.na(subunit_2) & is.na(subunit_3) & is.na(subunit_4)) & is.na(subunit_5) & (subunit_1 %in% filter.genes$genes)) |
                ((!is.na(subunit_1) & !is.na(subunit_2) & is.na(subunit_3) & is.na(subunit_4)) & is.na(subunit_5) & (subunit_1 %in% filter.genes$genes & subunit_2 %in% filter.genes$genes)) |
                ((!is.na(subunit_1) & !is.na(subunit_2) & !is.na(subunit_3) & is.na(subunit_4)) & is.na(subunit_5) & (subunit_1 %in% filter.genes$genes & subunit_2 %in% filter.genes$genes & subunit_3 %in% filter.genes$genes)) |
                ((!is.na(subunit_1) & !is.na(subunit_2) & !is.na(subunit_3) & !is.na(subunit_4)) & is.na(subunit_5) & (subunit_1 %in% filter.genes$genes & subunit_2 %in% filter.genes$genes & subunit_3 %in% filter.genes$genes & subunit_4 %in% filter.genes$genes)) |
                ((!is.na(subunit_1) & !is.na(subunit_2) & !is.na(subunit_3) & !is.na(subunit_4)) & !is.na(subunit_5) & (subunit_1 %in% filter.genes$genes & subunit_2 %in% filter.genes$genes & subunit_3 %in% filter.genes$genes & subunit_4 %in% filter.genes$genes & subunit_5 %in% filter.genes$genes)))

  } else{
    stop("'filter' must be one of 'strict' or 'relaxed'.")
  }

  fix_rownames <- complex_output_filter
  fix_rownames$new_rowname <- make.unique(fix_rownames$rowname,
                                          sep = ".")
 # complex_gene_match <- fix_rownames %>% select(rowname, new_rowname) # For later with interaction input

  fix_rownames$rowname <- make.unique(fix_rownames$rowname)

  complex_output <- fix_rownames %>%
    column_to_rownames("new_rowname") %>%
    select(-rowname)

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write_delim(complex_output, file.path(writeDir, "complex_input_axolotl.tsv"))
  }

 return(complex_output)

}

# make_cofactor_input
# REQUIRED:
# gene.list - the output of prep_gene_list that will be used to filter annotated genes
# annotation.map - the output of prep_gene_map that will be used to match human genes to axolotl
#
# OPTIONS:
# priority - in the event that there are both human and other annotations for a gene, give priority to either or none. (default "none")
# Note: the same priority will be applied if testing for synonyms.
# filter - how to apply filtering genes using gene.list. "Strict" (the default) requires that all annotated cofactors appear in gene.list. "Relaxed" requires at least one.
  # Note 1: the same filter will be applied if testing for synonyms.
  # Note 2: the filter is always "strict" for WNT antagonists owing to their number.
# annotation.cutoff - if fewer than this fraction of cofactors are matched to axolotl, drop from input. (default 0.5)
# drop.nonhuman - supercedes 'priority'; ignore any non-human annotations and drop from output (default FALSE)
# check.synonyms - attempt to rescue rows by also testing for synonymous genes? (default FALSE)
# synonyms - REQUIRED if check.synonyms is TRUE. The output of prep_synonym_map that will be used to test for synonyms.
# contains.spaces - do the genes output by prep_gene_list contain spaces? (default TRUE)
# write.outfile - write tab-delimited complex output file (default TRUE)
# out.dir - directory to output written file to (default working directory)

make_cofactor_input <- function(gene.list, annotation.map,
                                priority = c("human", "other", "none"), filter = c("strict", "relaxed"),
                                annotation.cutoff = NULL,
                                drop.nonhuman = FALSE, check.synonyms = FALSE, synonyms = NULL,
                                contains.spaces = TRUE,
                                write.outfile = TRUE, out.dir = NULL){

  select=dplyr::select
  if(missing(filter)){
    filter <- "strict"
  }
  if(!filter %in% c("strict", "relaxed")){
    stop("'filter' must be one of 'strict' or 'relaxed'.")
  }
  if(is.null(annotation.cutoff)){
    annotation.cutoff <- 0.5
  }
  if(!annotation.cutoff >= 0 & !annotation.cutoff <= 1 ){
    stop("'annotation.cutoff' must be a decimal between 0 and 1.")
  } else if(annotation.cutoff == 1){
    warning("An 'annotation.cutoff' of 1 may produce unintended results including an empty output.")
  }
  if(missing(priority)){
    priority <- "none"
  }
  if(!priority %in% c("human", "other", "none")){
    stop("'priority' must be one of 'human', 'other', or 'none'.")
  }

  cofactor_input <- CellChatDB.human$cofactor
  cofactor_input[cofactor_input == ""] <- NA
  cofactor_input <- cofactor_input %>%
    rownames_to_column("rowname")

  if(missing(gene.list)){
    stop("You must provide a gene list for filtering output by prep_gene_list.")
  }
  filter.genes <- gene.list

  if(isTRUE(check.synonyms)){
    if(is.null(synonyms)){
      stop("You must provide a Symbol-Synonym map output by prep_synonym_map to search for synonymous genes.")
    }
    cofactor_input_synonyms <- cofactor_input %>%
      left_join(synonyms, by=c("cofactor1" = "Symbol"), relationship = "many-to-many") %>%
      rename(cofactor1a = value) %>%
      left_join(synonyms, by=c("cofactor2" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor2a = value)  %>%
      left_join(synonyms, by=c("cofactor3" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor3a = value) %>%
      left_join(synonyms, by=c("cofactor4" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor4a = value) %>%
      left_join(synonyms, by=c("cofactor5" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor5a = value) %>%
      left_join(synonyms, by=c("cofactor6" = "Symbol"), relationship = "many-to-many") %>%
      rename(cofactor6a = value) %>%
      left_join(synonyms, by=c("cofactor7" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor7a = value)  %>%
      left_join(synonyms, by=c("cofactor8" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor8a = value) %>%
      left_join(synonyms, by=c("cofactor9" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor9a = value) %>%
      left_join(synonyms, by=c("cofactor10" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor10a = value) %>%
      left_join(synonyms, by=c("cofactor11" = "Symbol"), relationship = "many-to-many") %>%
      rename(cofactor11a = value) %>%
      left_join(synonyms, by=c("cofactor12" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor12a = value)  %>%
      left_join(synonyms, by=c("cofactor13" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor13a = value) %>%
      left_join(synonyms, by=c("cofactor14" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor14a = value) %>%
      left_join(synonyms, by=c("cofactor15" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor15a = value) %>%
      left_join(synonyms, by=c("cofactor16" = "Symbol"),relationship = "many-to-many") %>%
      rename(cofactor16a = value) %>%
      mutate(cofactor1b = case_when(is.na(cofactor1) ~ NA,
                                    !is.na(cofactor1) & is.na(cofactor1a) ~ cofactor1,
                                    !is.na(cofactor1) & !is.na(cofactor1a) ~ paste0(cofactor1,";",cofactor1a),
                                    TRUE ~ "FIX"),
             cofactor2b = case_when(is.na(cofactor2) ~ NA,
                                    !is.na(cofactor2) & is.na(cofactor2a) ~ cofactor2,
                                    !is.na(cofactor2) & !is.na(cofactor2a) ~ paste0(cofactor2,";",cofactor2a),
                                    TRUE ~ "FIX"),
             cofactor3b = case_when(is.na(cofactor3) ~ NA,
                                    !is.na(cofactor3) & is.na(cofactor3a) ~ cofactor3,
                                    !is.na(cofactor3) & !is.na(cofactor3a) ~ paste0(cofactor3,";",cofactor3a),
                                    TRUE ~ "FIX"),
             cofactor4b = case_when(is.na(cofactor4) ~ NA,
                                    !is.na(cofactor4) & is.na(cofactor4a) ~ cofactor4,
                                    !is.na(cofactor4) & !is.na(cofactor4a) ~ paste0(cofactor4,";",cofactor4a),
                                    TRUE ~ "FIX"),
             cofactor5b = case_when(is.na(cofactor5) ~ NA,
                                    !is.na(cofactor5) & is.na(cofactor5a) ~ cofactor5,
                                    !is.na(cofactor5) & !is.na(cofactor5a) ~ paste0(cofactor5,";",cofactor5a),
                                    TRUE ~ "FIX"),
             cofactor6b = case_when(is.na(cofactor6) ~ NA,
                                    !is.na(cofactor6) & is.na(cofactor6a) ~ cofactor6,
                                    !is.na(cofactor6) & !is.na(cofactor6a) ~ paste0(cofactor6,";",cofactor6a),
                                    TRUE ~ "FIX"),
             cofactor7b = case_when(is.na(cofactor7) ~ NA,
                                    !is.na(cofactor7) & is.na(cofactor7a) ~ cofactor7,
                                    !is.na(cofactor7) & !is.na(cofactor7a) ~ paste0(cofactor7,";",cofactor7a),
                                    TRUE ~ "FIX"),
             cofactor8b = case_when(is.na(cofactor8) ~ NA,
                                    !is.na(cofactor8) & is.na(cofactor8a) ~ cofactor8,
                                    !is.na(cofactor8) & !is.na(cofactor8a) ~ paste0(cofactor8,";",cofactor8a),
                                    TRUE ~ "FIX"),
             cofactor9b = case_when(is.na(cofactor9) ~ NA,
                                    !is.na(cofactor9) & is.na(cofactor9a) ~ cofactor9,
                                    !is.na(cofactor9) & !is.na(cofactor9a) ~ paste0(cofactor9,";",cofactor9a),
                                    TRUE ~ "FIX"),
             cofactor10b = case_when(is.na(cofactor10) ~ NA,
                                     !is.na(cofactor10) & is.na(cofactor10a) ~ cofactor10,
                                     !is.na(cofactor10) & !is.na(cofactor10a) ~ paste0(cofactor10,";",cofactor10a),
                                     TRUE ~ "FIX"),
             cofactor11b = case_when(is.na(cofactor11) ~ NA,
                                     !is.na(cofactor11) & is.na(cofactor11a) ~ cofactor11,
                                     !is.na(cofactor11) & !is.na(cofactor11a) ~ paste0(cofactor11,";",cofactor11a),
                                     TRUE ~ "FIX"),
             cofactor12b = case_when(is.na(cofactor12) ~ NA,
                                     !is.na(cofactor12) & is.na(cofactor12a) ~ cofactor12,
                                     !is.na(cofactor12) & !is.na(cofactor12a) ~ paste0(cofactor12,";",cofactor12a),
                                     TRUE ~ "FIX"),
             cofactor13b = case_when(is.na(cofactor13) ~ NA,
                                     !is.na(cofactor13) & is.na(cofactor13a) ~ cofactor13,
                                     !is.na(cofactor13) & !is.na(cofactor13a) ~ paste0(cofactor13,";",cofactor13a),
                                    TRUE ~ "FIX"),
             cofactor14b = case_when(is.na(cofactor14) ~ NA,
                                     !is.na(cofactor14) & is.na(cofactor14a) ~ cofactor14,
                                     !is.na(cofactor14) & !is.na(cofactor14a) ~ paste0(cofactor14,";",cofactor14a),
                                     TRUE ~ "FIX"),
             cofactor15b = case_when(is.na(cofactor15) ~ NA,
                                     !is.na(cofactor15) & is.na(cofactor15a) ~ cofactor15,
                                     !is.na(cofactor15) & !is.na(cofactor15a) ~ paste0(cofactor15,";",cofactor15a),
                                     TRUE ~ "FIX"),
             cofactor16b = case_when(is.na(cofactor16) ~ NA,
                                     !is.na(cofactor16) & is.na(cofactor16a) ~ cofactor16,
                                     !is.na(cofactor16) & !is.na(cofactor16a) ~ paste0(cofactor16,";",cofactor16a),
                                     TRUE ~ "FIX")) %>%
      select("rowname", ends_with("b")) %>%
      separate_longer_delim(cofactor1b, ";") %>%
      separate_longer_delim(cofactor2b, ";") %>%
      separate_longer_delim(cofactor3b, ";") %>%
      separate_longer_delim(cofactor4b, ";") %>%
      separate_longer_delim(cofactor5b, ";") %>%
      separate_longer_delim(cofactor6b, ";") %>%
      separate_longer_delim(cofactor7b, ";") %>%
      separate_longer_delim(cofactor8b, ";") %>%
      separate_longer_delim(cofactor9b, ";") %>%
      separate_longer_delim(cofactor10b, ";") %>%
      separate_longer_delim(cofactor11b, ";") %>%
      separate_longer_delim(cofactor12b, ";") %>%
      separate_longer_delim(cofactor13b, ";") %>%
      separate_longer_delim(cofactor14b, ";") %>%
      separate_longer_delim(cofactor15b, ";") %>%
      separate_longer_delim(cofactor16b, ";") %>%
      unique()
    colnames(cofactor_input_synonyms) <- colnames(cofactor_input)
    cofactor_count <- rowSums(!(is.na(cofactor_input_synonyms[,-1])))
    cofactor_input_synonyms$start_na <- cofactor_count # For later to determine what valid "missingness" is to keep after mapping
    cofactor_input <- cofactor_input_synonyms
    rm(cofactor_input_synonyms)

  } else{
    cofactor_count <- rowSums(!(is.na(cofactor_input[,-1])))
    cofactor_input$start_na <- cofactor_count
  }

  cofactor_input_wnts <- cofactor_input %>%
    filter(rowname == "WNT antagonist") %>%
    pivot_longer(.,starts_with("cofactor"), names_to = "num") %>%
    left_join(.,annotation_map, by = c("value" = "human")) %>%
    rename(humanAnnot = full) %>%
    left_join(.,annotation_map, by = c("value" = "other")) %>%
    rename(otherAnnot = full) %>%
    filter(humanAnnot %in% filter.genes | otherAnnot %in% filter.genes) %>%
    select(rowname, start_na, num, humanAnnot, otherAnnot) %>%
    unique()

  col.add <- colnames(cofactor_input)[!colnames(cofactor_input) %in% cofactor_input_wnts$num]
  col.add <- col.add[!col.add %in% c("rowname", "start_na")]
  for( i in 1:length(col.add)){
    df <- data.frame(rowname = "WNT antagonist", start_na = 16, num = col.add[[i]],
                     humanAnnot = NA, otherAnnot = NA)
    cofactor_input_wnts <- rbind(cofactor_input_wnts, df)
  }

  if(isTRUE(drop.nonhuman)){
    message("drop.nonhuman supercedes annotation priority. Only human annotations will be searched for.")
    cofactor_output_1 <- cofactor_input %>% # Mapping subunits to axolotl genes
      filter(rowname != "WNT antagonist") %>%
      left_join(annotation.map, by=c("cofactor1" = "human"), relationship = "many-to-many") %>%
      rename(final_1 = full) %>%
      left_join(annotation.map, by=c("cofactor2" = "human"),relationship = "many-to-many") %>%
      rename(final_2 = full)  %>%
      left_join(annotation.map, by=c("cofactor3" = "human"),relationship = "many-to-many") %>%
      rename(final_3 = full) %>%
      left_join(annotation.map, by=c("cofactor4" = "human"),relationship = "many-to-many") %>%
      rename(final_4 = full) %>%
      left_join(annotation.map, by=c("cofactor5" = "human"),relationship = "many-to-many") %>%
      rename(final_5 = full) %>%
      left_join(annotation.map, by=c("cofactor6" = "human"), relationship = "many-to-many") %>%
      rename(final_6 = full) %>%
      left_join(annotation.map, by=c("cofactor7" = "human"),relationship = "many-to-many") %>%
      rename(final_7 = full)  %>%
      left_join(annotation.map, by=c("cofactor8" = "human"),relationship = "many-to-many") %>%
      rename(final_8 = full) %>%
      left_join(annotation.map, by=c("cofactor9" = "human"),relationship = "many-to-many") %>%
      rename(final_9 = full) %>%
      left_join(annotation.map, by=c("cofactor10" = "human"),relationship = "many-to-many") %>%
      rename(final_10 = full) %>%
      left_join(annotation.map, by=c("cofactor11" = "human"), relationship = "many-to-many") %>%
      rename(final_11 = full) %>%
      left_join(annotation.map, by=c("cofactor12" = "human"),relationship = "many-to-many") %>%
      rename(final_12 = full)  %>%
      left_join(annotation.map, by=c("cofactor13" = "human"),relationship = "many-to-many") %>%
      rename(final_13 = full) %>%
      left_join(annotation.map, by=c("cofactor14" = "human"),relationship = "many-to-many") %>%
      rename(final_14 = full) %>%
      left_join(annotation.map, by=c("cofactor15" = "human"),relationship = "many-to-many") %>%
      rename(final_15 = full) %>%
      left_join(annotation.map, by=c("cofactor16" = "human"),relationship = "many-to-many") %>%
      rename(final_16 = full)

    cofactor_output_1 <- cofactor_output_1 %>%
      select(c("rowname", "start_na", starts_with("final_"))) %>%
      unique() %>%
      `colnames<-`(c("rowname", "start_na", "cofactor1", "cofactor2", "cofactor3",
                     "cofactor4", "cofactor5", "cofactor6", "cofactor7", "cofactor8",
                     "cofactor9", "cofactor10", "cofactor11", "cofactor12", "cofactor13",
                     "cofactor14", "cofactor15", "cofactor16"))

    cofactor_input_wnts <- cofactor_input_wnts %>%
      select(-otherAnnot) %>%
      pivot_wider(names_from = num,
                  names_sep = "",
                  values_from = humanAnnot,
                  values_fn = list)

    for( i in 1:ncol(cofactor_input_wnts)){
      col <- colnames(cofactor_input_wnts)[[i]]
      cofactor_input_wnts <- unnest(cofactor_input_wnts, col)
    }

    cofactor_input_wnts <- cofactor_input_wnts %>%
      select(rowname, start_na, colnames(cofactor_input)) %>%
      unique()

    cofactor_output_1 <- rbind(cofactor_output_1, cofactor_input_wnts)

  } else{
    cofactor_output_1 <- cofactor_input %>% # Mapping subunits to axolotl genes
      filter(rowname != "WNT antagonist") %>%
      left_join(annotation.map, by=c("cofactor1" = "human"), relationship = "many-to-many") %>%
      rename(cofactor1a = full) %>%
      left_join(annotation.map, by=c("cofactor2" = "human"),relationship = "many-to-many") %>%
      rename(cofactor2a = full)  %>%
      left_join(annotation.map, by=c("cofactor3" = "human"),relationship = "many-to-many") %>%
      rename(cofactor3a = full) %>%
      left_join(annotation.map, by=c("cofactor4" = "human"),relationship = "many-to-many") %>%
      rename(cofactor4a = full) %>%
      left_join(annotation.map, by=c("cofactor5" = "human"),relationship = "many-to-many") %>%
      rename(cofactor5a = full) %>%
      left_join(annotation.map, by=c("cofactor6" = "human"), relationship = "many-to-many") %>%
      rename(cofactor6a = full) %>%
      left_join(annotation.map, by=c("cofactor7" = "human"),relationship = "many-to-many") %>%
      rename(cofactor7a = full)  %>%
      left_join(annotation.map, by=c("cofactor8" = "human"),relationship = "many-to-many") %>%
      rename(cofactor8a = full) %>%
      left_join(annotation.map, by=c("cofactor9" = "human"),relationship = "many-to-many") %>%
      rename(cofactor9a = full) %>%
      left_join(annotation.map, by=c("cofactor10" = "human"),relationship = "many-to-many") %>%
      rename(cofactor10a = full) %>%
      left_join(annotation.map, by=c("cofactor11" = "human"), relationship = "many-to-many") %>%
      rename(cofactor11a = full) %>%
      left_join(annotation.map, by=c("cofactor12" = "human"),relationship = "many-to-many") %>%
      rename(cofactor12a = full)  %>%
      left_join(annotation.map, by=c("cofactor13" = "human"),relationship = "many-to-many") %>%
      rename(cofactor13a = full) %>%
      left_join(annotation.map, by=c("cofactor14" = "human"),relationship = "many-to-many") %>%
      rename(cofactor14a = full) %>%
      left_join(annotation.map, by=c("cofactor15" = "human"),relationship = "many-to-many") %>%
      rename(cofactor15a = full) %>%
      left_join(annotation.map, by=c("cofactor16" = "human"),relationship = "many-to-many") %>%
      rename(cofactor16a = full) %>%
      left_join(annotation.map, by=c("cofactor1" = "other"), relationship = "many-to-many") %>%
      rename(cofactor1b = full) %>%
      left_join(annotation.map, by=c("cofactor2" = "other"),relationship = "many-to-many") %>%
      rename(cofactor2b = full)  %>%
      left_join(annotation.map, by=c("cofactor3" = "other"),relationship = "many-to-many") %>%
      rename(cofactor3b = full) %>%
      left_join(annotation.map, by=c("cofactor4" = "other"),relationship = "many-to-many") %>%
      rename(cofactor4b = full) %>%
      left_join(annotation.map, by=c("cofactor5" = "other"),relationship = "many-to-many") %>%
      rename(cofactor5b = full) %>%
      left_join(annotation.map, by=c("cofactor6" = "other"), relationship = "many-to-many") %>%
      rename(cofactor6b = full) %>%
      left_join(annotation.map, by=c("cofactor7" = "other"),relationship = "many-to-many") %>%
      rename(cofactor7b = full)  %>%
      left_join(annotation.map, by=c("cofactor8" = "other"),relationship = "many-to-many") %>%
      rename(cofactor8b = full) %>%
      left_join(annotation.map, by=c("cofactor9" = "other"),relationship = "many-to-many") %>%
      rename(cofactor9b = full) %>%
      left_join(annotation.map, by=c("cofactor10" = "other"),relationship = "many-to-many") %>%
      rename(cofactor10b = full) %>%
      left_join(annotation.map, by=c("cofactor11" = "other"), relationship = "many-to-many") %>%
      rename(cofactor11b = full) %>%
      left_join(annotation.map, by=c("cofactor12" = "other"),relationship = "many-to-many") %>%
      rename(cofactor12b = full)  %>%
      left_join(annotation.map, by=c("cofactor13" = "other"),relationship = "many-to-many") %>%
      rename(cofactor13b = full) %>%
      left_join(annotation.map, by=c("cofactor14" = "other"),relationship = "many-to-many") %>%
      rename(cofactor14b = full) %>%
      left_join(annotation.map, by=c("cofactor15" = "other"),relationship = "many-to-many") %>%
      rename(cofactor15b = full) %>%
      left_join(annotation.map, by=c("cofactor16" = "other"),relationship = "many-to-many") %>%
      rename(cofactor16b = full)

    rm(cofactor_input)

    if(priority == "human"){
      cofactor_output_1 <- cofactor_output_1 %>%
        mutate(final_1 = case_when(!is.na(cofactor1a) & is.na(cofactor1b) ~ cofactor1a,
                                   !is.na(cofactor1a) & !is.na(cofactor1b) ~ cofactor1a,
                                   is.na(cofactor1a) & !is.na(cofactor1b) ~ cofactor1b,
                                   is.na(cofactor1a) & is.na(cofactor1b) ~ NA),
               final_2 = case_when(is.na(cofactor2) ~ NA,
                                   !is.na(cofactor2a) & is.na(cofactor2b) ~ cofactor2a,
                                   !is.na(cofactor2a) & !is.na(cofactor2b) ~ cofactor2a,
                                   is.na(cofactor2a) & !is.na(cofactor2b) ~ cofactor2b,
                                   is.na(cofactor2a) & is.na(cofactor2b) ~ NA),
               final_3 = case_when(is.na(cofactor3) ~ NA,
                                   !is.na(cofactor3a) & is.na(cofactor3b) ~ cofactor3a,
                                   !is.na(cofactor3a) & !is.na(cofactor3b) ~ cofactor3a,
                                   is.na(cofactor3a) & !is.na(cofactor3b) ~ cofactor3b,
                                   is.na(cofactor3a) & is.na(cofactor3b) ~ NA),
               final_4 = case_when(is.na(cofactor4) ~ NA,
                                   !is.na(cofactor4a) & is.na(cofactor4b) ~ cofactor4a,
                                   !is.na(cofactor4a) & !is.na(cofactor4b) ~ cofactor4a,
                                   is.na(cofactor4a) & !is.na(cofactor4b) ~ cofactor4b,
                                   is.na(cofactor4a) & is.na(cofactor4b) ~ NA),
               final_5 = case_when(is.na(cofactor5) ~ NA,
                                   !is.na(cofactor5a) & is.na(cofactor5b) ~ cofactor5a,
                                   !is.na(cofactor5a) & !is.na(cofactor5b) ~ cofactor5a,
                                   is.na(cofactor5a) & !is.na(cofactor5b) ~ cofactor5b,
                                   is.na(cofactor5a) & is.na(cofactor5b) ~ NA),
               final_6 = case_when(is.na(cofactor6) ~ NA,
                                   !is.na(cofactor6a) & is.na(cofactor6b) ~ cofactor6a,
                                   !is.na(cofactor6a) & !is.na(cofactor6b) ~ cofactor6a,
                                   is.na(cofactor6a) & !is.na(cofactor6b) ~ cofactor6b,
                                   is.na(cofactor6a) & is.na(cofactor6b) ~ NA),
               final_7 = case_when(is.na(cofactor7) ~ NA,
                                   !is.na(cofactor7a) & is.na(cofactor7b) ~ cofactor7a,
                                   !is.na(cofactor7a) & !is.na(cofactor7b) ~ cofactor7a,
                                   is.na(cofactor7a) & !is.na(cofactor7b) ~ cofactor7b,
                                   is.na(cofactor7a) & is.na(cofactor7b) ~ NA),
               final_8 = case_when(is.na(cofactor8) ~ NA,
                                   !is.na(cofactor8a) & is.na(cofactor8b) ~ cofactor8a,
                                   !is.na(cofactor8a) & !is.na(cofactor8b) ~ cofactor8a,
                                   is.na(cofactor8a) & !is.na(cofactor8b) ~ cofactor8b,
                                   is.na(cofactor8a) & is.na(cofactor8b) ~ NA),
               final_9 = case_when(is.na(cofactor9) ~ NA,
                                   !is.na(cofactor9a) & is.na(cofactor9b) ~ cofactor9a,
                                   !is.na(cofactor9a) & !is.na(cofactor9b) ~ cofactor9a,
                                   is.na(cofactor9a) & !is.na(cofactor9b) ~ cofactor9b,
                                   is.na(cofactor9a) & is.na(cofactor9b) ~ NA),
               final_10 = case_when(is.na(cofactor10) ~ NA,
                                   !is.na(cofactor10a) & is.na(cofactor10b) ~ cofactor10a,
                                   !is.na(cofactor10a) & !is.na(cofactor10b) ~ cofactor10a,
                                   is.na(cofactor10a) & !is.na(cofactor10b) ~ cofactor10b,
                                   is.na(cofactor10a) & is.na(cofactor10b) ~ NA),
               final_11 = case_when(is.na(cofactor11) ~ NA,
                                    !is.na(cofactor11a) & is.na(cofactor11b) ~ cofactor11a,
                                   !is.na(cofactor11a) & !is.na(cofactor11b) ~ cofactor11a,
                                   is.na(cofactor11a) & !is.na(cofactor11b) ~ cofactor11b,
                                   is.na(cofactor11a) & is.na(cofactor11b) ~ NA),
               final_12 = case_when(is.na(cofactor12) ~ NA,
                                   !is.na(cofactor12a) & is.na(cofactor12b) ~ cofactor12a,
                                   !is.na(cofactor12a) & !is.na(cofactor12b) ~ cofactor12a,
                                   is.na(cofactor12a) & !is.na(cofactor12b) ~ cofactor12b,
                                   is.na(cofactor12a) & is.na(cofactor12b) ~ NA),
               final_13 = case_when(is.na(cofactor13) ~ NA,
                                   !is.na(cofactor13a) & is.na(cofactor13b) ~ cofactor13a,
                                   !is.na(cofactor13a) & !is.na(cofactor13b) ~ cofactor13a,
                                   is.na(cofactor13a) & !is.na(cofactor13b) ~ cofactor13b,
                                   is.na(cofactor13a) & is.na(cofactor13b) ~ NA),
               final_14 = case_when(is.na(cofactor14) ~ NA,
                                   !is.na(cofactor14a) & is.na(cofactor14b) ~ cofactor14a,
                                   !is.na(cofactor14a) & !is.na(cofactor14b) ~ cofactor14a,
                                   is.na(cofactor14a) & !is.na(cofactor14b) ~ cofactor14b,
                                   is.na(cofactor14a) & is.na(cofactor14b) ~ NA),
               final_15 = case_when(is.na(cofactor15) ~ NA,
                                   !is.na(cofactor15a) & is.na(cofactor15b) ~ cofactor15a,
                                   !is.na(cofactor15a) & !is.na(cofactor15b) ~ cofactor15a,
                                   is.na(cofactor15a) & !is.na(cofactor15b) ~ cofactor15b,
                                   is.na(cofactor15a) & is.na(cofactor15b) ~ NA),
               final_16 = case_when(is.na(cofactor16) ~ NA,
                                   !is.na(cofactor16a) & is.na(cofactor16b) ~ cofactor16a,
                                   !is.na(cofactor16a) & !is.na(cofactor16b) ~ cofactor16a,
                                   is.na(cofactor16a) & !is.na(cofactor16b) ~ cofactor16b,
                                   is.na(cofactor16a) & is.na(cofactor16b) ~ NA)) %>%
        select(c("rowname", "start_na", starts_with("final"))) %>%
        unique() %>%
        `colnames<-`(c("rowname", "start_na", "cofactor1", "cofactor2",
                       "cofactor3", "cofactor4", "cofactor5", "cofactor6",
                       "cofactor7", "cofactor8", "cofactor9", "cofactor10",
                       "cofactor11", "cofactor12", "cofactor13", "cofactor14",
                       "cofactor15", "cofactor16"))

      cofactor_input_wnts <- cofactor_input_wnts %>%
        mutate(final = case_when(is.na(humanAnnot) & is.na(otherAnnot) ~ NA,
                                 !is.na(humanAnnot) & is.na(otherAnnot) ~ humanAnnot,
                                 is.na(humanAnnot) & !is.na(otherAnnot) ~ otherAnnot,
                                 !is.na(humanAnnot) & !is.na(otherAnnot) ~ humanAnnot)) %>%
        select(rowname, start_na, num, final) %>%
        separate_longer_delim(final, ";") %>%
        unique() %>%
        pivot_wider(names_from = num,
                    names_sep="",
                    values_from = final,
                    values_fn = list)

      for( i in 1:ncol(cofactor_input_wnts)){
        col <- colnames(cofactor_input_wnts)[[i]]
        cofactor_input_wnts <- unnest(cofactor_input_wnts, col)
      }

      cofactor_input_wnts <- cofactor_input_wnts %>%
        select(rowname, start_na, colnames(cofactor_input)) %>%
        unique()

      cofactor_output_1 <- rbind(cofactor_output_1, cofactor_input_wnts)

    } else if(priority == "other"){
      cofactor_output_1 <- cofactor_output_1 %>%
        mutate(final_1 = case_when(!is.na(cofactor1a) & is.na(cofactor1b) ~ cofactor1a,
                                   !is.na(cofactor1a) & !is.na(cofactor1b) ~ cofactor1b,
                                   is.na(cofactor1a) & !is.na(cofactor1b) ~ cofactor1b,
                                   is.na(cofactor1a) & is.na(cofactor1b) ~ NA),
               final_2 = case_when(is.na(cofactor2) ~ NA,
                                   !is.na(cofactor2a) & is.na(cofactor2b) ~ cofactor2a,
                                   !is.na(cofactor2a) & !is.na(cofactor2b) ~ cofactor2b,
                                   is.na(cofactor2a) & !is.na(cofactor2b) ~ cofactor2b,
                                   is.na(cofactor2a) & is.na(cofactor2b) ~ NA),
               final_3 = case_when(is.na(cofactor3) ~ NA,
                                   !is.na(cofactor3a) & is.na(cofactor3b) ~ cofactor3a,
                                   !is.na(cofactor3a) & !is.na(cofactor3b) ~ cofactor3b,
                                   is.na(cofactor3a) & !is.na(cofactor3b) ~ cofactor3b,
                                   is.na(cofactor3a) & is.na(cofactor3b) ~ NA),
               final_4 = case_when(is.na(cofactor4) ~ NA,
                                   !is.na(cofactor4a) & is.na(cofactor4b) ~ cofactor4a,
                                   !is.na(cofactor4a) & !is.na(cofactor4b) ~ cofactor4b,
                                   is.na(cofactor4a) & !is.na(cofactor4b) ~ cofactor4b,
                                   is.na(cofactor4a) & is.na(cofactor4b) ~ NA),
               final_5 = case_when(is.na(cofactor5) ~ NA,
                                   !is.na(cofactor5a) & is.na(cofactor5b) ~ cofactor5a,
                                   !is.na(cofactor5a) & !is.na(cofactor5b) ~ cofactor5b,
                                   is.na(cofactor5a) & !is.na(cofactor5b) ~ cofactor5b,
                                   is.na(cofactor5a) & is.na(cofactor5b) ~ NA),
               final_6 = case_when(is.na(cofactor6) ~ NA,
                                   !is.na(cofactor6a) & is.na(cofactor6b) ~ cofactor6a,
                                   !is.na(cofactor6a) & !is.na(cofactor6b) ~ cofactor6b,
                                   is.na(cofactor6a) & !is.na(cofactor6b) ~ cofactor6b,
                                   is.na(cofactor6a) & is.na(cofactor6b) ~ NA),
               final_7 = case_when(is.na(cofactor7) ~ NA,
                                   !is.na(cofactor7a) & is.na(cofactor7b) ~ cofactor7a,
                                   !is.na(cofactor7a) & !is.na(cofactor7b) ~ cofactor7b,
                                   is.na(cofactor7a) & !is.na(cofactor7b) ~ cofactor7b,
                                   is.na(cofactor7a) & is.na(cofactor7b) ~ NA),
               final_8 = case_when(is.na(cofactor8) ~ NA,
                                   !is.na(cofactor8a) & is.na(cofactor8b) ~ cofactor8a,
                                   !is.na(cofactor8a) & !is.na(cofactor8b) ~ cofactor8b,
                                   is.na(cofactor8a) & !is.na(cofactor8b) ~ cofactor8b,
                                   is.na(cofactor8a) & is.na(cofactor8b) ~ NA),
               final_9 = case_when(is.na(cofactor9) ~ NA,
                                   !is.na(cofactor9a) & is.na(cofactor9b) ~ cofactor9a,
                                   !is.na(cofactor9a) & !is.na(cofactor9b) ~ cofactor9b,
                                   is.na(cofactor9a) & !is.na(cofactor9b) ~ cofactor9b,
                                   is.na(cofactor9a) & is.na(cofactor9b) ~ NA),
               final_10 = case_when(is.na(cofactor10) ~ NA,
                                   !is.na(cofactor10a) & is.na(cofactor10b) ~ cofactor10a,
                                   !is.na(cofactor10a) & !is.na(cofactor10b) ~ cofactor10b,
                                   is.na(cofactor10a) & !is.na(cofactor10b) ~ cofactor10b,
                                   is.na(cofactor10a) & is.na(cofactor10b) ~ NA),
               final_11 = case_when(is.na(cofactor11) ~ NA,
                                   !is.na(cofactor11a) & is.na(cofactor11b) ~ cofactor11a,
                                   !is.na(cofactor11a) & !is.na(cofactor11b) ~ cofactor11b,
                                   is.na(cofactor11a) & !is.na(cofactor11b) ~ cofactor11b,
                                   is.na(cofactor11a) & is.na(cofactor11b) ~ NA),
               final_12 = case_when(is.na(cofactor12) ~ NA,
                                   !is.na(cofactor12a) & is.na(cofactor12b) ~ cofactor12a,
                                   !is.na(cofactor12a) & !is.na(cofactor12b) ~ cofactor12b,
                                   is.na(cofactor12a) & !is.na(cofactor12b) ~ cofactor12b,
                                   is.na(cofactor12a) & is.na(cofactor12b) ~ NA),
               final_13 = case_when(is.na(cofactor13) ~ NA,
                                   !is.na(cofactor13a) & is.na(cofactor13b) ~ cofactor13a,
                                   !is.na(cofactor13a) & !is.na(cofactor13b) ~ cofactor13b,
                                   is.na(cofactor13a) & !is.na(cofactor13b) ~ cofactor13b,
                                   is.na(cofactor13a) & is.na(cofactor13b) ~ NA),
               final_14 = case_when(is.na(cofactor14) ~ NA,
                                   !is.na(cofactor14a) & is.na(cofactor14b) ~ cofactor14a,
                                   !is.na(cofactor14a) & !is.na(cofactor14b) ~ cofactor14b,
                                   is.na(cofactor14a) & !is.na(cofactor14b) ~ cofactor14b,
                                   is.na(cofactor14a) & is.na(cofactor14b) ~ NA),
               final_15 = case_when(is.na(cofactor15) ~ NA,
                                   !is.na(cofactor15a) & is.na(cofactor15b) ~ cofactor15a,
                                   !is.na(cofactor15a) & !is.na(cofactor15b) ~ cofactor15b,
                                   is.na(cofactor15a) & !is.na(cofactor15b) ~ cofactor15b,
                                   is.na(cofactor15a) & is.na(cofactor15b) ~ NA),
               final_16 = case_when(is.na(cofactor16) ~ NA,
                                   !is.na(cofactor16a) & is.na(cofactor16b) ~ cofactor16a,
                                   !is.na(cofactor16a) & !is.na(cofactor16b) ~ cofactor16b,
                                   is.na(cofactor16a) & !is.na(cofactor16b) ~ cofactor16b,
                                   is.na(cofactor16a) & is.na(cofactor16b) ~ NA)) %>%
        select(c("rowname", "start_na", starts_with("final"))) %>%
        unique() %>%
        `colnames<-`(c("rowname", "start_na",
                       "cofactor1", "cofactor2", "cofactor3", "cofactor4", "cofactor5",
                       "cofactor6", "cofactor7", "cofactor8", "cofactor9", "cofactor10",
                       "cofactor11", "cofactor12", "cofactor13", "cofactor14",
                       "cofactor15", "cofactor16"))

      cofactor_input_wnts <- cofactor_input_wnts %>%
        mutate(final = case_when(is.na(humanAnnot) & is.na(otherAnnot) ~ NA,
                                 !is.na(humanAnnot) & is.na(otherAnnot) ~ humanAnnot,
                                 is.na(humanAnnot) & !is.na(otherAnnot) ~ otherAnnot,
                                 !is.na(humanAnnot) & !is.na(otherAnnot) ~ otherAnnot)) %>%
        select(rowname, start_na, num, final) %>%
        separate_longer_delim(final, ";") %>%
        unique() %>%
        pivot_wider(names_from = num,
                    names_sep="",
                    values_from = final,
                    values_fn = list)

      for( i in 1:ncol(cofactor_input_wnts)){
        col <- colnames(cofactor_input_wnts)[[i]]
        cofactor_input_wnts <- unnest(cofactor_input_wnts, col)
      }

      cofactor_input_wnts <- cofactor_input_wnts %>%
        select(rowname, start_na, colnames(cofactor_input)) %>%
        unique()

      cofactor_output_1 <- rbind(cofactor_output_1, cofactor_input_wnts)

    } else if(priority == "none"){
      cofactor_output_1 <- cofactor_output_1 %>%
        mutate(final_1 = case_when(is.na(cofactor1) ~ NA,
                                   !is.na(cofactor1a) & is.na(cofactor1b) ~ cofactor1a,
                                   !is.na(cofactor1a) & !is.na(cofactor1b) ~ paste0(cofactor1a,";",cofactor1b),
                                   is.na(cofactor1a) & !is.na(cofactor1b) ~ cofactor1b,
                                   is.na(cofactor1a) & is.na(cofactor1b) ~ NA),
               final_2 = case_when(is.na(cofactor2) ~ NA,
                                   !is.na(cofactor2a) & is.na(cofactor2b) ~ cofactor2a,
                                   !is.na(cofactor2a) & !is.na(cofactor2b) ~ paste0(cofactor2a,";",cofactor2b),
                                   is.na(cofactor2a) & !is.na(cofactor2b) ~ cofactor2b,
                                   is.na(cofactor2a) & is.na(cofactor2b) ~ NA),
               final_3 = case_when(is.na(cofactor3) ~ NA,
                                   !is.na(cofactor3a) & is.na(cofactor3b) ~ cofactor3a,
                                   !is.na(cofactor3a) & !is.na(cofactor3b) ~ paste0(cofactor3a,";",cofactor3b),
                                   is.na(cofactor3a) & !is.na(cofactor3b) ~ cofactor3b,
                                   is.na(cofactor3a) & is.na(cofactor3b) ~ NA),
               final_4 = case_when(is.na(cofactor4) ~ NA,
                                   !is.na(cofactor4a) & is.na(cofactor4b) ~ cofactor4a,
                                   !is.na(cofactor4a) & !is.na(cofactor4b) ~ paste0(cofactor4a,";",cofactor4b),
                                   is.na(cofactor4a) & !is.na(cofactor4b) ~ cofactor4b,
                                   is.na(cofactor4a) & is.na(cofactor4b) ~ NA),
               final_5 = case_when(is.na(cofactor5) ~ NA,
                                   !is.na(cofactor5a) & is.na(cofactor5b) ~ cofactor5a,
                                   !is.na(cofactor5a) & !is.na(cofactor5b) ~ paste0(cofactor5a,";",cofactor5b),
                                   is.na(cofactor5a) & !is.na(cofactor5b) ~ cofactor5b,
                                   is.na(cofactor5a) & is.na(cofactor5b) ~ NA),
               final_6 = case_when(is.na(cofactor6) ~ NA,
                                   !is.na(cofactor6a) & is.na(cofactor6b) ~ cofactor6a,
                                   !is.na(cofactor6a) & !is.na(cofactor6b) ~ paste0(cofactor6a,";",cofactor6b),
                                   is.na(cofactor6a) & !is.na(cofactor6b) ~ cofactor6b,
                                   is.na(cofactor6a) & is.na(cofactor6b) ~ NA),
               final_7 = case_when(is.na(cofactor7) ~ NA,
                                   !is.na(cofactor7a) & is.na(cofactor7b) ~ cofactor7a,
                                   !is.na(cofactor7a) & !is.na(cofactor7b) ~ paste0(cofactor7a,";",cofactor7b),
                                   is.na(cofactor7a) & !is.na(cofactor7b) ~ cofactor7b,
                                   is.na(cofactor7a) & is.na(cofactor7b) ~ NA),
               final_8 = case_when(is.na(cofactor8) ~ NA,
                                   !is.na(cofactor8a) & is.na(cofactor8b) ~ cofactor8a,
                                   !is.na(cofactor8a) & !is.na(cofactor8b) ~ paste0(cofactor8a,";",cofactor8b),
                                   is.na(cofactor8a) & !is.na(cofactor8b) ~ cofactor8b,
                                   is.na(cofactor8a) & is.na(cofactor8b) ~ NA),
               final_9 = case_when(is.na(cofactor9) ~ NA,
                                   !is.na(cofactor9a) & is.na(cofactor9b) ~ cofactor9a,
                                   !is.na(cofactor9a) & !is.na(cofactor9b) ~ paste0(cofactor9a,";",cofactor9b),
                                   is.na(cofactor9a) & !is.na(cofactor9b) ~ cofactor9b,
                                   is.na(cofactor9a) & is.na(cofactor9b) ~ NA),
               final_10 = case_when(is.na(cofactor10) ~ NA,
                                   !is.na(cofactor10a) & is.na(cofactor10b) ~ cofactor10a,
                                   !is.na(cofactor10a) & !is.na(cofactor10b) ~ paste0(cofactor10a,";",cofactor10b),
                                   is.na(cofactor10a) & !is.na(cofactor10b) ~ cofactor10b,
                                   is.na(cofactor10a) & is.na(cofactor10b) ~ NA),
               final_11 = case_when(is.na(cofactor11) ~ NA,
                                   !is.na(cofactor11a) & is.na(cofactor11b) ~ cofactor11a,
                                   !is.na(cofactor11a) & !is.na(cofactor11b) ~ paste0(cofactor11a,";",cofactor11b),
                                   is.na(cofactor11a) & !is.na(cofactor11b) ~ cofactor11b,
                                   is.na(cofactor11a) & is.na(cofactor11b) ~ NA),
               final_12 = case_when(is.na(cofactor12) ~ NA,
                                   !is.na(cofactor12a) & is.na(cofactor12b) ~ cofactor12a,
                                   !is.na(cofactor12a) & !is.na(cofactor12b) ~ paste0(cofactor12a,";",cofactor12b),
                                   is.na(cofactor12a) & !is.na(cofactor12b) ~ cofactor12b,
                                   is.na(cofactor12a) & is.na(cofactor12b) ~ NA),
               final_13 = case_when(is.na(cofactor13) ~ NA,
                                   !is.na(cofactor13a) & is.na(cofactor13b) ~ cofactor13a,
                                   !is.na(cofactor13a) & !is.na(cofactor13b) ~ paste0(cofactor13a,";",cofactor13b),
                                   is.na(cofactor13a) & !is.na(cofactor13b) ~ cofactor13b,
                                   is.na(cofactor13a) & is.na(cofactor13b) ~ NA),
               final_14 = case_when(is.na(cofactor14) ~ NA,
                                   !is.na(cofactor14a) & is.na(cofactor14b) ~ cofactor14a,
                                   !is.na(cofactor14a) & !is.na(cofactor14b) ~ paste0(cofactor14a,";",cofactor14b),
                                   is.na(cofactor14a) & !is.na(cofactor14b) ~ cofactor14b,
                                   is.na(cofactor14a) & is.na(cofactor14b) ~ NA),
               final_15 = case_when(is.na(cofactor15) ~ NA,
                                   !is.na(cofactor15a) & is.na(cofactor15b) ~ cofactor15a,
                                   !is.na(cofactor15a) & !is.na(cofactor15b) ~ paste0(cofactor15a,";",cofactor15b),
                                   is.na(cofactor15a) & !is.na(cofactor15b) ~ cofactor15b,
                                   is.na(cofactor15a) & is.na(cofactor15b) ~ NA),
               final_16 = case_when(is.na(cofactor16) ~ NA,
                                   !is.na(cofactor16a) & is.na(cofactor16b) ~ cofactor16a,
                                   !is.na(cofactor16a) & !is.na(cofactor16b) ~ paste0(cofactor16a,";",cofactor16b),
                                   is.na(cofactor16a) & !is.na(cofactor16b) ~ cofactor16b,
                                   is.na(cofactor16a) & is.na(cofactor16b) ~ NA)) %>%
        select(c("rowname", "start_na", starts_with("final"))) %>%
        separate_longer_delim(final_1,";") %>%
        separate_longer_delim(final_2, ";") %>%
        separate_longer_delim(final_3,";") %>%
        separate_longer_delim(final_4, ";") %>%
        separate_longer_delim(final_5,";") %>%
        separate_longer_delim(final_6,";") %>%
        separate_longer_delim(final_7, ";") %>%
        separate_longer_delim(final_8,";") %>%
        separate_longer_delim(final_9, ";") %>%
        separate_longer_delim(final_10,";") %>%
        separate_longer_delim(final_11,";") %>%
        separate_longer_delim(final_12, ";") %>%
        separate_longer_delim(final_13,";") %>%
        separate_longer_delim(final_14, ";") %>%
        separate_longer_delim(final_15,";") %>%
        separate_longer_delim(final_16,";") %>%
        unique() %>%
        `colnames<-`(c("rowname", "start_na", "cofactor1", "cofactor2", "cofactor3", "cofactor4", "cofactor5",
                       "cofactor6", "cofactor7", "cofactor8", "cofactor9", "cofactor10",
                       "cofactor11", "cofactor12", "cofactor13", "cofactor14", "cofactor15", "cofactor16"))

      cofactor_input_wnts <- cofactor_input_wnts %>%
        mutate(final = case_when(is.na(humanAnnot) & is.na(otherAnnot) ~ NA,
                                 !is.na(humanAnnot) & is.na(otherAnnot) ~ humanAnnot,
                                 is.na(humanAnnot) & !is.na(otherAnnot) ~ otherAnnot,
                                 !is.na(humanAnnot) & !is.na(otherAnnot) ~ paste0(humanAnnot,";",otherAnnot))) %>%
        select(rowname, start_na, num, final) %>%
        separate_longer_delim(final, ";") %>%
        unique() %>%
        pivot_wider(names_from = num,
                    names_sep="",
                    values_from = final,
                    values_fn = list)

      for( i in 1:ncol(cofactor_input_wnts)){
        col <- colnames(cofactor_input_wnts)[[i]]
        cofactor_input_wnts <- unnest(cofactor_input_wnts, col)
      }

      cofactor_input_wnts <- cofactor_input_wnts %>%
        select(rowname, start_na, colnames(cofactor_input)) %>%
        unique()

      cofactor_output_1 <- rbind(cofactor_output_1, cofactor_input_wnts)
    }

  }

  if(isFALSE(contains.spaces)){
    cofactor_output_1 <- cofactor_output_1 %>%
      mutate(cofactor1 = gsub(" \\[", "\\[", cofactor1),
             cofactor2 = case_when(is.na(cofactor2) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", cofactor2)),
             cofactor3 = case_when(is.na(cofactor3) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", cofactor3)),
             cofactor4 = case_when(is.na(cofactor4) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", cofactor4)),
             cofactor5 = case_when(is.na(cofactor5) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", cofactor5)),
             cofactor6 = case_when(is.na(cofactor6) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", cofactor6)),
             cofactor7 = case_when(is.na(cofactor7) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", cofactor7)),
             cofactor8 = case_when(is.na(cofactor8) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", cofactor8)),
             cofactor9 = case_when(is.na(cofactor9) ~ NA,
                                   TRUE ~ gsub(" \\[", "\\[", cofactor9)),
             cofactor10 = case_when(is.na(cofactor10) ~ NA,
                                    TRUE ~ gsub(" \\[", "\\[", cofactor10)),
             cofactor11 = case_when(is.na(cofactor11) ~ NA,
                                    TRUE ~ gsub(" \\[", "\\[", cofactor11)),
             cofactor12 = case_when(is.na(cofactor12) ~ NA,
                                    TRUE ~ gsub(" \\[", "\\[", cofactor12)),
             cofactor13 = case_when(is.na(cofactor13) ~ NA,
                                    TRUE ~ gsub(" \\[", "\\[", cofactor13)),
             cofactor14 = case_when(is.na(cofactor14) ~ NA,
                                    TRUE ~ gsub(" \\[", "\\[", cofactor14)),
             cofactor15 = case_when(is.na(cofactor15) ~ NA,
                                    TRUE ~ gsub(" \\[", "\\[", cofactor15)),
             cofactor16 = case_when(is.na(cofactor16) ~ NA,
                                    TRUE ~ gsub(" \\[", "\\[", cofactor16)))
  }

  # Now filtering
  if(filter == "strict"){
    cofactor_output_filter <- cofactor_output_1 %>%
      filter(!is.na(cofactor1) & cofactor1 %in% filter.genes) %>%
      filter((is.na(cofactor2)) |
               (is.na(cofactor3) & cofactor2 %in% filter.genes) |
               (is.na(cofactor4) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes) |
               (is.na(cofactor5) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes) |
               (is.na(cofactor6) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes) |
               (is.na(cofactor7) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes) |
               (is.na(cofactor8) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes) |
               (is.na(cofactor9) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes & cofactor8 %in% filter.genes) |
               (is.na(cofactor10) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes & cofactor8 %in% filter.genes &
                  cofactor9 %in% filter.genes) |
               (is.na(cofactor11) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes & cofactor8 %in% filter.genes &
                  cofactor9 %in% filter.genes & cofactor10 %in% filter.genes) |
               (is.na(cofactor12) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes & cofactor8 %in% filter.genes &
                  cofactor9 %in% filter.genes & cofactor10 %in% filter.genes & cofactor11 %in% filter.genes) |
               (is.na(cofactor13) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes & cofactor8 %in% filter.genes &
                 cofactor9 %in% filter.genes & cofactor10 %in% filter.genes & cofactor11 %in% filter.genes & cofactor12 %in% filter.genes) |
             (is.na(cofactor14) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes & cofactor8 %in% filter.genes &
                  cofactor9 %in% filter.genes & cofactor10 %in% filter.genes & cofactor11 %in% filter.genes & cofactor12 %in% filter.genes &
                  cofactor13 %in% filter.genes) |
               (!is.na(cofactor15) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                  cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes & cofactor8 %in% filter.genes &
                  cofactor9 %in% filter.genes & cofactor10 %in% filter.genes & cofactor11 %in% filter.genes & cofactor12 %in% filter.genes &
                  cofactor13 %in% filter.genes & cofactor14 %in% filter.genes) |
             (!is.na(cofactor16) & cofactor2 %in% filter.genes & cofactor3 %in% filter.genes & cofactor4 %in% filter.genes &
                cofactor5 %in% filter.genes & cofactor6 %in% filter.genes & cofactor7 %in% filter.genes & cofactor8 %in% filter.genes &
                cofactor9 %in% filter.genes & cofactor10 %in% filter.genes & cofactor11 %in% filter.genes & cofactor12 %in% filter.genes &
                cofactor13 %in% filter.genes & cofactor14 %in% filter.genes & cofactor15 %in% filter.genes))

  } else if(filter == "relaxed"){
    cofactor_output_filter <- cofactor_output_1 %>%
      filter(cofactor1 %in% filter.genes |
               cofactor2 %in% filter.genes |
               cofactor3 %in% filter.genes |
               cofactor4 %in% filter.genes |
               cofactor5 %in% filter.genes |
               cofactor6 %in% filter.genes |
               cofactor7 %in% filter.genes |
               cofactor8 %in% filter.genes |
               cofactor9 %in% filter.genes |
               cofactor10 %in% filter.genes |
               cofactor11 %in% filter.genes |
               cofactor12 %in% filter.genes |
               cofactor13 %in% filter.genes |
               cofactor14 %in% filter.genes |
               cofactor15 %in% filter.genes |
               cofactor16 %in% filter.genes)
  }

  cofactor_output_2 <- cofactor_output_filter %>% unique()

  out_na = rowSums(!(is.na(cofactor_output_2[,-c(1,2)])))
  cofactor_output_2$end_na <- out_na
  cofactor_output_2 <- cofactor_output_2 %>%
    select(rowname, start_na, end_na, starts_with("cofactor")) %>%
    mutate(na_diff = end_na/start_na) %>%
    filter(na_diff >= annotation.cutoff) %>% # Drop if greater than X fraction of subunits are not annotated
    select(-c(start_na, end_na, na_diff)) %>%
    unique()

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write_delim(cofactor_output_2, file.path(writeDir, "cofactor_input_axolotl.tsv"))
  }
  return(cofactor_output_2)

}


# make_interaction_input
# REQUIRED:
# gene.list - the output of prep_gene_list that will be used to filter annotated genes
# annotation.map - the output of prep_gene_map that will be used to match human genes to axolotl
# complex - the output of make_complex_input
# cofactor - the output of make_cofactor_input
#
# OPTIONS:
# priority - in the event that there are both human and other annotations for a gene, give priority to either or none. (default "none")
# Note: the same priority will be applied if testing for synonyms.
# filter - how to apply filtering genes using gene.list. "Strict" (the default) requires that all annotated cofactors appear in gene.list. "Relaxed" requires at least one.
# Note 1: the same filter will be applied if testing for synonyms.
# Note 2: the filter is always "strict" for WNT antagonists owing to their number.
# annotation.cutoff - if fewer than this fraction of cofactors are matched to axolotl, drop from input. (default 0.5)
# drop.nonhuman - supercedes 'priority'; ignore any non-human annotations and drop from output (default FALSE)
# check.synonyms - attempt to rescue rows by also testing for synonymous genes? (default FALSE)
# synonyms - REQUIRED if check.synonyms is TRUE. The output of prep_synonym_map that will be used to test for synonyms.
# contains.spaces - do the genes output by prep_gene_list contain spaces? (default TRUE)
# write.outfile - write tab-delimited complex output file (default TRUE)
# out.dir - directory to output written file to (default working directory)

make_interaction_input <- function(gene.list, annotation.map,
                                   complex, cofactor,
                                   priority = c("human", "other", "none"), filter = c("strict", "relaxed"),
                                   annotation.cutoff = NULL,
                                   drop.nonhuman = FALSE, check.synonyms = FALSE, synonyms = NULL,
                                   contains.spaces = TRUE,
                                   write.outfile = TRUE, out.dir = NULL){

  # To do ---
  # inherit params from complex, cofactor
  # Complex interactions

  if(missing(filter)){
    filter <- "strict"
  }
  if(!filter %in% c("strict", "relaxed")){
    stop("'filter' must be one of 'strict' or 'relaxed'.")
  }
  #if(is.null(annotation.cutoff)){
  #  annotation.cutoff <- 0.5
  #}
  #if(!annotation.cutoff >= 0 & !annotation.cutoff <= 1 ){
  #  stop("'annotation.cutoff' must be a decimal between 0 and 1.")
  #} else if(annotation.cutoff == 1){
  #  warning("An 'annotation.cutoff' of 1 may produce unintended results including an empty output.")
  #}
  if(missing(priority)){
    priority <- "none"
  }
  if(!priority %in% c("human", "other", "none")){
    stop("'priority' must be one of 'human', 'other', or 'none'.")
  }
  if(missing(cofactor)){
    stop("You must provide the output of make_cofactor_input.")
  }
  if(missing(complex)){
    stop("You must provide the output of make_complex_input.")
  }

  interaction_input <- CellChatDB.human$interaction
  interaction_input[interaction_input == ""] <- NA
  interaction_input <- interaction_input %>%
    rownames_to_column("rowname")

  filter.genes <- gene.list

  if(isTRUE(check.synonyms)){
    if(is.null(synonyms)){
      stop("You must provide the Symbol-Synonym map output by prep_synonym_map in order to check for synonyms.")
    }
    # Can only check for synonyms with non-complex ligands/receptors
    interaction_input_synonyms <- interaction_input %>%
      left_join(synonyms, by=c("ligand" = "Symbol"), relationship = "many-to-many") %>%
      rename(ligand1a = value) %>%
      left_join(synonyms, by=c("receptor" = "Symbol"),relationship = "many-to-many") %>%
      rename(receptor1a = value) %>%
      mutate(ligand1b = case_when(!is.na(ligand1a) ~ paste0(ligand, ";", ligand1a),
                                  TRUE ~ ligand),
             receptor1b = case_when(!is.na(receptor1a) ~ paste0(receptor,";",receptor1a),
                                    TRUE ~ receptor)) %>%
      mutate(ligand = ligand1b, receptor = receptor1b) %>%
      select(colnames(interaction_input)) %>%
      separate_longer_delim(ligand,";") %>%
      separate_longer_delim(receptor,";") %>%
      unique()

    interaction_input <- interaction_input_synonyms
    rm(interaction_input_synonyms)
  }

  if(isTRUE(drop.nonhuman)){
    message("drop.nonhuman supercedes annotation priority. Only human annotations will be searched for.")
    # Simple interactions first - one ligand, one receptor
    interaction_input_simple <- interaction_input %>%
      dplyr::filter(!(ligand %in% rownames(complex) & receptor %in% rownames(complex))) %>%
      left_join(annotation.map, by=c("ligand" = "human"), relationship = "many-to-many") %>%
      rename(ligand_1a = full) %>%
      left_join(annotation.map, by=c("receptor" = "human"),relationship = "many-to-many") %>%
      rename(receptor_1a = full) %>%
      mutate(ligand = ligand_1a, receptor = receptor_1a) %>%
      select(c(colnames(interaction_input))) %>%
      filter(!is.na(ligand) & !is.na(receptor)) %>%
      unique()

    interaction_input_complex <- interaction_input %>%
      dplyr::filter(ligand %in% rownames(complex) | receptor %in% rownames(complex)) %>%
      left_join(annotation.map, by=c("ligand" = "human"), relationship = "many-to-many") %>%
      rename(ligand_1a = full) %>%
      left_join(annotation.map, by=c("receptor" = "human"),relationship = "many-to-many") %>%
      rename(receptor_1a = full) %>%
      mutate(final_ligand = case_when(!is.na(ligand_1a) ~ ligand_1a,
                                      is.na(ligand_1a) & ligand %in% rownames(complex) ~ ligand,
                                      is.na(ligand_1a) & !ligand %in% rownames(complex) ~ NA,
                                      TRUE ~ "FIX"),
             final_receptor = case_when(!is.na(receptor_1a) ~ receptor_1a,
                                        is.na(receptor_1a) & receptor %in% rownames(complex) ~ receptor,
                                        is.na(receptor_1a) & !receptor %in% rownames(complex) ~ NA,
                                        TRUE ~ "FIX")) %>%
      mutate(ligand = final_ligand, receptor = final_receptor) %>%
      dplyr::filter(!is.na(ligand)) %>%
      dplyr::filter(!is.na(receptor)) %>%
      select(colnames(interaction_input)) %>%
      unique()

  } else{
    interaction_input_simple <- interaction_input %>%
      dplyr::filter(!(ligand %in% rownames(complex) & receptor %in% rownames(complex))) %>%
      left_join(annotation.map, by=c("ligand" = "human"), relationship = "many-to-many") %>%
      rename(ligand_1a = full) %>%
      left_join(annotation.map, by=c("receptor" = "human"),relationship = "many-to-many") %>%
      rename(receptor_1a = full)  %>%
      left_join(annotation.map, by=c("ligand" = "other"), relationship = "many-to-many") %>%
      rename(ligand_1b = full) %>%
      left_join(annotation.map, by=c("receptor" = "other"),relationship = "many-to-many") %>%
      rename(receptor_1b = full)

    if(priority == "human"){
      interaction_input_simple <- interaction_input_simple %>%
        mutate(final_ligand = case_when(!is.na(ligand_1a) & is.na(ligand_1b) ~ ligand_1a,
                                        !is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1a,
                                        is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                        is.na(ligand_1a) & is.na(ligand_1b) ~ NA),
               final_receptor = case_when(!is.na(receptor_1a) & is.na(receptor_1b) ~ receptor_1a,
                                          !is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1a,
                                          is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                          is.na(receptor_1a) & is.na(receptor_1b) ~ NA)) %>%
        mutate(ligand = final_ligand, receptor = final_receptor) %>%
        select(c(colnames(interaction_input))) %>%
        filter(!is.na(ligand) & !is.na(receptor)) %>%
        unique()

      interaction_input_complex <- interaction_input %>%
        dplyr::filter(ligand %in% rownames(complex) | receptor %in% rownames(complex)) %>%
        left_join(annotation.map, by=c("ligand" = "human"), relationship = "many-to-many") %>%
        rename(ligand_1a = full) %>%
        left_join(annotation.map, by=c("receptor" = "human"),relationship = "many-to-many") %>%
        rename(receptor_1a = full) %>%
        left_join(annotation.map, by=c("ligand" = "other"), relationship = "many-to-many") %>%
        rename(ligand_1b = full) %>%
        left_join(annotation.map, by=c("receptor" = "other"),relationship = "many-to-many") %>%
        rename(receptor_1b = full) %>%
        mutate(final_ligand = case_when(!is.na(ligand_1a) & is.na(ligand_1b) ~ ligand_1a,
                                        !is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1a,
                                        is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                        is.na(ligand_1a) & is.na(ligand_1b) ~ NA,
                                        TRUE ~ "FIX"),
               final_receptor = case_when(!is.na(receptor_1a) & is.na(receptor_1b) ~ receptor_1a,
                                          !is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1a,
                                          is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                          is.na(receptor_1a) & is.na(receptor_1b) ~ NA,
                                          TRUE ~ "FIX")) %>%
        mutate(ligand = ifelse(!(ligand %in% rownames(complex)), final_ligand, ligand)) %>%
        mutate(receptor = ifelse(!(receptor %in% rownames(complex)), final_receptor, receptor)) %>%
        unique() %>%
        dplyr::filter(!is.na(ligand)) %>%
        dplyr::filter(!is.na(receptor)) %>%
        select(colnames(interaction_input)) %>%
        unique()

    } else if(priority == "other"){
      interaction_input_simple <- interaction_input_simple %>%
        mutate(final_ligand = case_when(!is.na(ligand_1a) & is.na(ligand_1b) ~ ligand_1a,
                                        !is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                        is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                        is.na(ligand_1a) & is.na(ligand_1b) ~ NA),
               final_receptor = case_when(!is.na(receptor_1a) & is.na(receptor_1b) ~ receptor_1a,
                                          !is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                          is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                          is.na(receptor_1a) & is.na(receptor_1b) ~ NA)) %>%
        mutate(ligand = final_ligand, receptor = final_receptor) %>%
        select(c(colnames(interaction_input))) %>%
        filter(!is.na(ligand) & !is.na(receptor)) %>%
        unique()

      # Complex interactions next - ligand or receptor or both are a complex
      interaction_input_complex <- interaction_input %>%
        dplyr::filter(ligand %in% rownames(complex) | receptor %in% rownames(complex)) %>%
        left_join(annotation.map, by=c("ligand" = "human"), relationship = "many-to-many") %>%
        rename(ligand_1a = full) %>%
        left_join(annotation.map, by=c("receptor" = "human"),relationship = "many-to-many") %>%
        rename(receptor_1a = full) %>%
        left_join(annotation.map, by=c("ligand" = "other"), relationship = "many-to-many") %>%
        rename(ligand_1b = full) %>%
        left_join(annotation.map, by=c("receptor" = "other"),relationship = "many-to-many") %>%
        rename(receptor_1b = full) %>%
        mutate(final_ligand = case_when(!is.na(ligand_1a) & is.na(ligand_1b) ~ ligand_1a,
                                        !is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                        is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                        is.na(ligand_1a) & is.na(ligand_1b) ~ NA,
                                        TRUE ~ "FIX"),
               final_receptor = case_when(!is.na(receptor_1a) & is.na(receptor_1b) ~ receptor_1a,
                                          !is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                          is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                          is.na(receptor_1a) & is.na(receptor_1b) ~ NA,
                                          TRUE ~ "FIX")) %>%
        mutate(ligand = ifelse(!(ligand %in% rownames(complex)), final_ligand, ligand)) %>%
        mutate(receptor = ifelse(!(receptor %in% rownames(complex)), final_receptor, receptor)) %>%
        unique() %>%
        dplyr::filter(!is.na(ligand)) %>%
        dplyr::filter(!is.na(receptor)) %>%
        select(colnames(interaction_input)) %>%
        unique()

    } else if(priority == "none"){
      interaction_input_simple <- interaction_input_simple %>%
        mutate(final_ligand = case_when(!is.na(ligand_1a) & is.na(ligand_1b) ~ ligand_1a,
                                        !is.na(ligand_1a) & !is.na(ligand_1b) ~ paste0(ligand_1a,";",ligand_1b),
                                        is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                        is.na(ligand_1a) & is.na(ligand_1b) ~ NA),
               final_receptor = case_when(!is.na(receptor_1a) & is.na(receptor_1b) ~ receptor_1a,
                                          !is.na(receptor_1a) & !is.na(receptor_1b) ~ paste0(receptor_1a,";",receptor_1b),
                                          is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                          is.na(receptor_1a) & is.na(receptor_1b) ~ NA)) %>%
        mutate(ligand = final_ligand, receptor = final_receptor) %>%
        select(c(colnames(interaction_input))) %>%
        separate_longer_delim(ligand,";") %>%
        separate_longer_delim(receptor, ";") %>%
        filter(!is.na(ligand) & !is.na(receptor)) %>%
        unique()

      # Complex interactions next - ligand or receptor or both are a complex
      interaction_input_complex <- interaction_input %>%
        dplyr::filter(ligand %in% rownames(complex) | receptor %in% rownames(complex)) %>%
        left_join(annotation.map, by=c("ligand" = "human"), relationship = "many-to-many") %>%
        rename(ligand_1a = full) %>%
        left_join(annotation.map, by=c("receptor" = "human"),relationship = "many-to-many") %>%
        rename(receptor_1a = full) %>%
        left_join(annotation.map, by=c("ligand" = "other"), relationship = "many-to-many") %>%
        rename(ligand_1b = full) %>%
        left_join(annotation.map, by=c("receptor" = "other"),relationship = "many-to-many") %>%
        rename(receptor_1b = full) %>%
        mutate(final_ligand = case_when(!is.na(ligand_1a) & is.na(ligand_1b) ~ ligand_1a,
                                        !is.na(ligand_1a) & !is.na(ligand_1b) ~ paste0(ligand_1a,";",ligand_1b),
                                        is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                        is.na(ligand_1a) & is.na(ligand_1b) ~ NA,
                                        TRUE ~ "FIX"),
               final_receptor = case_when(!is.na(receptor_1a) & is.na(receptor_1b) ~ receptor_1a,
                                          !is.na(receptor_1a) & !is.na(receptor_1b) ~ paste0(receptor_1a,";",receptor_1b),
                                          is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                          is.na(receptor_1a) & is.na(receptor_1b) ~ NA,
                                          TRUE ~ "FIX")) %>%
        mutate(ligand = ifelse(!(ligand %in% rownames(complex)), final_ligand, ligand)) %>%
        mutate(receptor = ifelse(!(receptor %in% rownames(complex)), final_receptor, receptor)) %>%
        separate_longer_delim(ligand,";") %>%
        separate_longer_delim(receptor, ";") %>%
        unique() %>%
        dplyr::filter(!is.na(ligand)) %>%
        dplyr::filter(!is.na(receptor)) %>%
        select(colnames(interaction_input)) %>%
        unique()

    }
  }

  # Initial filter on simple input only
  if(filter == "strict"){
    interaction_input_simple <- interaction_input_simple %>%
      filter(ligand %in% filter.genes & receptor %in% filter.genes)
  } else if(filter == "relaxed"){
    interaction_input_simple <- interaction_input_simple %>%
      filter(ligand %in% filter.genes | receptor %in% filter.genes)
  }

  # Now mapping ligands and receptors to rownames of complex input
  interaction_input_complex_rownameMerge <- interaction_input_complex %>%
    left_join(complex_gene_match, by=c("ligand" = "rowname"), relationship = "many-to-many") %>%
    rename(ligand_match = new_rowname) %>%
    left_join(complex_gene_match, by=c("receptor" = "rowname"), relationship = "many-to-many") %>%
    rename(receptor_match = new_rowname) %>%
    mutate(final_ligand = case_when((!ligand %in% complex_gene_match$rowname) & is.na(ligand_match) ~ ligand,
                                    (ligand %in% complex_gene_match$rowname) & !is.na(ligand_match) ~ ligand_match,
                                    TRUE ~ "FIX"),
           final_receptor = case_when((!receptor %in% complex_gene_match$rowname) & is.na(receptor_match) ~ receptor,
                                      (receptor %in% complex_gene_match$rowname) & !is.na(receptor_match) ~ receptor_match,
                                      TRUE ~ "FIX")) %>%
    mutate(ligand = final_ligand, receptor = final_receptor)

  interaction_input_complex_rownameMerge <- interaction_input_complex_rownameMerge %>%
    select(colnames(interaction_input_complex)) %>%
    unique()

  interaction_input_bind <- rbind(interaction_input_simple, interaction_input_complex_rownameMerge)
  interaction_input_bind <- interaction_input_bind %>%
    mutate(ligand = gsub(" \\[", "\\[", ligand),
           receptor = gsub(" \\[", "\\[", receptor)) %>%
    unique()

  # Another filter - for both strict and relaxed
  interaction_input_bind_test <- interaction_input_bind %>%
    filter((ligand %in% rna.test.genes | receptor %in% rna.test.genes) |
             (ligand %in% rownames(complex_output) & receptor %in% rna.test.genes) |
             (ligand %in% rna.test.genes & receptor %in% rownames(complex_output)) |
             (ligand %in% rownames(complex_output) & receptor %in% rownames(complex_output)))

  interaction_input_bind <- interaction_input_bind_test

  # Map interaction input ant/agonists, co-receptors to cofactor input -----
  interaction_input_cofactor_rownameMerge <- interaction_input_bind %>%
    mutate(agonist = case_when(agonist == "" ~ NA,
                               TRUE ~ agonist),
           antagonist = case_when(antagonist == "" ~ NA,
                                  TRUE ~ antagonist),
           co_A_receptor = case_when(co_A_receptor == "" ~ NA,
                                     TRUE ~ co_A_receptor),
           co_I_receptor = case_when(co_I_receptor == "" ~ NA,
                                     TRUE ~ co_I_receptor))

  # Nothing to map in these - will just tack on after
  interaction_input_cofactor_rownameMerge_noCofactors <- interaction_input_cofactor_rownameMerge %>%
    filter(is.na(agonist) & is.na(antagonist) & is.na(co_A_receptor) & is.na(co_I_receptor))

  # Need to map to cofactor input
  interaction_input_cofactor_rownameMerge_withCofactors <- interaction_input_cofactor_rownameMerge %>%
    filter(!is.na(agonist) | !is.na(antagonist) | !is.na(co_A_receptor) | !is.na(co_I_receptor))

  interaction_input_cofactor_rownameMerge_Full <- interaction_input_cofactor_rownameMerge_withCofactors %>%
    left_join(cofactor_gene_match, by=c("agonist" = "interaction"), relationship = "many-to-many") %>%
    rename(agonist_match = rowname) %>%
    left_join(cofactor_gene_match, by=c("antagonist" = "interaction"), relationship = "many-to-many") %>%
    rename(antagonist_match = rowname) %>%
    left_join(cofactor_gene_match, by=c("co_A_receptor" = "interaction"), relationship = "many-to-many") %>%
    rename(co_A_match = rowname) %>%
    left_join(cofactor_gene_match, by=c("co_I_receptor" = "interaction"), relationship = "many-to-many") %>%
    rename(co_I_match = rowname) %>%
    mutate(agonist = agonist_match,
           antagonist = antagonist_match,
           co_A_receptor = co_A_match,
           co_I_receptor = co_I_match) %>%
    select(colnames(interaction_input_cofactor_rownameMerge_withCofactors)) %>%
    unique()

  interaction_input_cofactors_matched <- rbind(interaction_input_cofactor_rownameMerge_Full,
                                               interaction_input_cofactor_rownameMerge_noCofactors) %>%
    unique()

  interaction_input_cofactors_matched$interaction_name <- make.unique(interaction_input_cofactors_matched$interaction_name,
                                                                      sep = ".")
  rownames(interaction_input_cofactors_matched) <- NULL

  interaction_output <- interaction_input_cofactors_matched %>%
    column_to_rownames("interaction_name") %>%
    unique()

  print( paste0("final number of rows: ", nrow(interaction_output)))

  gc()

  write.csv(interaction_output, file = "./for_stats/Interaction.RNAFilter.HumanPriority.Relaxed.NoSynonyms.csv")

  # Finals
  write.csv(complex_output, "./complex_input_12NOV24_SCTFiltered.csv")
  write.csv(cofactor_output, "./cofactor_input_12NOV24_SCTFiltered.csv")
  write.csv(interaction_output, "./interaction_input_12NOV24_SCTFiltered.csv")

}

# Not yet started -----

# make_geneInfo_input

make_geneInfo_input <- function(from.CC = CellChatDB.human$geneInfo)

# assemble_database

# build_DB

build_DB <- function(use.defaults = TRUE, ...){}
