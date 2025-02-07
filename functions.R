# Source functions for building custom axolotl CellChat DB

# prep_gene_list - build gene lists for downstream database filtering
#
# REQUIRED:
# objects - list of input Seurat objects as RDS files
#
# OPTIONS:
# assay - one of "SCT" or "RNA" (default "RNA")
# contains.spaces - do the gene names in the input gene lists contain spaces? (default TRUE)
# write.outfile - write gene list to selected directory (default TRUE)
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

# prep_synonym_map
#
# OPTIONS:
# write.outfile - write tab-delimited long-form map of Symbol - Synonym (default TRUE)
# out.dir - directory to output written files to (default working directory)

prep_synonym_map <- function(write.outfile = TRUE, out.dir = NULL){

  writeDir <- out.dir
  if(is.null(writeDir)){
    writeDir <- "./"
  }

  synonyms <- CellChatDB.human$geneInfo %>%
    select(Symbol, Synonym) %>%
    separate_wider_delim(., Synonym, " ", names_sep = "_", too_few = "align_start") %>% # Fill empties with NA
    separate_longer_delim(., Symbol, "; ") # Some are separated by ;, creating new rows for each
  synonyms[] <- lapply(synonyms, function(x) gsub(";[[:space:]]?", "", x))
  synonyms[synonyms == ""] <- NA # ex. RPL9 - "synonyms" ARE ";"

  synonym_count = rowSums(!(is.na(synonyms[,-1])))
  synonyms$count_na <- synonym_count
  synonyms <- synonyms %>%
    filter(count_na != 0) %>% # If a symbol doesn't have any synonyms, it won't help
    select(-c(count_na)) %>%
    as.data.frame()

  synonyms_long <- pivot_longer(synonyms, cols = starts_with("Synonym"), values_drop_na = T) %>%
    select(-2) %>%
    unique() %>%
    mutate(Symbol = toupper(Symbol),
           value = toupper(value))

  if(isTRUE(write.outfile)){
    write_delim(synonyms_long, file.path(writeDir, "Symbol_synonym_map.txt"))
  }
  return(synonyms_long)
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
    stop("'annotation.cutoff' must be a decimal between 0 and 1 (inclusive).")
  } else if(annotation.cutoff == 0){
    warning("An 'annotation.cutoff' of 0 may produce unintended consequences including an empty output.")
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
                                    is.na(subunit_1) & !is.na(subunit_1a) ~ subunit_1a,
                                    is.na(subunit_1) & is.na(subunit_1a) ~ NA,
                                    TRUE ~ "FIX"),
             subunit_2b = case_when(!is.na(subunit_2) & is.na(subunit_2a) ~ subunit_2,
                                    !is.na(subunit_2) & !is.na(subunit_2a) ~ paste0(subunit_2,";",subunit_2a),
                                    is.na(subunit_2) & !is.na(subunit_2a) ~ subunit_2a,
                                    is.na(subunit_2) & is.na(subunit_2a) ~ NA,
                                    TRUE ~ "FIX"),
             subunit_3b = case_when(!is.na(subunit_3) & is.na(subunit_3a) ~ subunit_3,
                                    !is.na(subunit_3) & !is.na(subunit_3a) ~ paste0(subunit_3,";",subunit_3a),
                                    is.na(subunit_3) & !is.na(subunit_3a) ~ subunit_3a,
                                    is.na(subunit_3) & is.na(subunit_3a) ~ NA,
                                    TRUE ~ "FIX"),
             subunit_4b = case_when(!is.na(subunit_4) & is.na(subunit_4a) ~ subunit_4,
                                    !is.na(subunit_4) & !is.na(subunit_4a) ~ paste0(subunit_4,";",subunit_4a),
                                    is.na(subunit_4) & !is.na(subunit_4a) ~ subunit_4a,
                                    is.na(subunit_4) & is.na(subunit_4a) ~ NA,
                                    TRUE ~ "FIX"),
             subunit_5b = case_when(!is.na(subunit_5) & is.na(subunit_5a) ~ subunit_5,
                                    !is.na(subunit_5) & !is.na(subunit_5a) ~ paste0(subunit_5,";",subunit_5a),
                                    is.na(subunit_5) & !is.na(subunit_5a) ~ subunit_5a,
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
    filter(na_diff >= annotation.cutoff) %>% # Drop if greater than X fraction of subunits are not annotated
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
  complex_gene_match <- fix_rownames %>% select(rowname, new_rowname) # For later with interaction input

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

# check_synonyms_internal
# Might not need many of these parameters & just take from input functions
# probably retain: input, input_type
# Might be obsolete if just generating inputs pre-mapped with synonyms
check_synonyms <- function(input, synonym_map,
                                    input_type = c("complex", "interaction", "cofactor", "geneInfo"),
                                    priority = c("human", "other", "none"), filter = c("strict", "relaxed"),
                                    annotation.cutoff = NULL, drop.nonhuman = FALSE,
                                    contains.spaces = TRUE,
                                    write.outfile = TRUE, out.dir = NULL){

  if(input_type %in% c("complex")){

  } else if(input_type %in% c("cofactor")){

  } else if(input_type %in% c("interaction")){

  } else if(input_type %in% c("geneInfo")){

  }

}

# Not yet started -----
# make_cofactor_input

make_cofactor_input <- function(from.CC = CellChatDB.human$cofactor){}

# make_interaction_input

make_interaction_input <- function(from.CC = CellChatDB.human$interaction){}

# make_geneInfo_input

make_geneInfo_input <- function(from.CC = CellChatDB.human$geneInfo)

# assemble_database

# build_DB

build_DB <- function(use.defaults = TRUE, ...){}
