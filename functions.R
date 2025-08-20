# prep_gene_list - build gene lists for downstream database filtering
#
# REQUIRED:
# objects - list of input Seurat objects
#
# OPTIONS:
# assay - the assay to pull feature names from (default "RNA")
# write.outfile - write gene list to selected directory (default TRUE)
# out.dir - directory to output written files to (defaults to working directory)

prep_gene_list <- function(objects,
                           assay = NULL,
                           write.outfile = TRUE,
                           out.dir = NULL){
  if(missing(objects)){
    stop("You must provide at least one Seurat object.")
  }

  object.list <- objects

  gene.list <- c()
  for(object in 1:length(object.list)){
    input <- object.list[[object]]

    if(missing(assay)){
      DefaultAssay(input) <- "RNA"
      genes.add <- dimnames(input)[1]
      gene.list <- c(gene.list, genes.add)
    } else{
      DefaultAssay(input) <- assay
      genes.add <- dimnames(input)[1]
      gene.list <- c(gene.list, genes.add)
    }
  }

  gene.list <- gene.list %>% unlist() %>% unique()

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write_delim(as.data.frame(gene.list) %>% `colnames<-`(c("Genes")),
                file.path(writeDir, "Genes_for_filtering.txt"))
  }
  return(gene.list)
}

# prep_gene_map - prepare annotations for mapping axolotl to human
#
# REQUIRED:
# gxt.map - tab-delimited file mapping transcripts to gene annotations; provided in repository as AmexV6_TxAnnot.tsv
# gene.list - the output of prep_gene_list
#
# OPTIONS:
# contains.spaces - do the genes output by prep_gene_list contain spaces? (default FALSE)
#   - For example, "ZFP37 [nr]|ZNF568 [hs]|AMEX60DD301000001.1" vs. "ZFP37[nr]|ZNF568[hs]|AMEX60DD301000001.1"
# drop.nonhuman - will not check for any non-human annotations (default FALSE)
# write.outfile - write tab-delimited annotation file (default TRUE)
# out.dir - directory to output written files to (default working directory)

prep_gene_map <- function(gxt.map,
                          gene.list = NULL,
                          contains.spaces = FALSE,
                          drop.nonhuman = FALSE,
                          write.outfile = TRUE, out.dir = NULL){

  if(missing(gxt.map)){
    stop("You must provide the path to the 'AmexV6_TxAnnot.tsv' file provided in this repository.")
  }

  gxtMap <- read.delim(gxt.map, header = TRUE, sep = " ") %>%
    `colnames<-`(c("full", "human", "other")) %>%
    mutate(full = case_when(isTRUE(contains.spaces) ~ gsub("\\[", " \\[", full),
                            TRUE ~ full))

  if(is.null(gene.list)){
    stop("Please provide the output of prep_gene_list.")
  } else{
    filter.genes <- as.data.frame(gene.list) %>% `colnames<-`(c("genes"))
    gxtMap <- gxtMap %>%
      filter(full %in% filter.genes$genes)
  }

  if(isTRUE(drop.nonhuman)){
    message("drop.nonhuman supercedes prioritize.human. All nonhuman annotations will be dropped from output.")
    gxtMap <- gxtMap %>%
      filter(!human == "") %>%
      mutate(other = "")
  }

  parseline <- c(full = "", human = "", other = "")
  gxtMap <- rbind(gxtMap, parseline)

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write_delim(gxtMap, file.path(writeDir, "geneAnnotations.txt"))
  }

  return(gxtMap)
}

# prep_synonym_map - prepare Symbol-synonym map for testing synonymous genes in building inputs
#
# REQUIRED:
# gene.map - the output of prep_gene_map that will be used to match annotations to axolotl genes
#
# OPTIONS:
# write.outfile - write tab-delimited long-form map of Symbol - Synonym (default TRUE)
# out.dir - directory to output written files to (default working directory)

prep_synonym_map <- function(gene.map = NULL,
                             write.outfile = TRUE,
                             out.dir = NULL){

  if(is.null(gene.map)){
    stop("You must provide the annotation map output by prep_gene_map.")
  }

  geneInfo_input <- CellChatDB.human$geneInfo %>%
    mutate(Symbol = toupper(Symbol),
           Synonym = toupper(Synonym)) %>%
    separate_longer_delim(Symbol,
                          "; ") %>% # ex. Q8IVW1
    separate_longer_delim(Symbol, ";") %>% # ex. M1V481
    filter(!is.na(EntryID.uniprot))

  geneInfo_toSub <- geneInfo_input[,2:3] # only want to edit these columns
  geneInfo_toSub[] <- lapply(geneInfo_toSub, function(x) gsub(";", "", x))
  geneInfo_toSub[geneInfo_toSub == ""] <- NA # ex. RPL9 - "synonyms" ARE ";"
  geneInfo_toSub <- geneInfo_toSub %>%
    mutate(Synonym = case_when(str_detect(Synonym, "^\\s+$") ~ NA,
                               TRUE ~ Synonym)) %>%
    mutate(Synonym = gsub("\\s+$", "", Synonym))

  geneInfo_noSub <- geneInfo_input[,-c(2:3)] # need to replace these columns with toSub values
  geneInfo_input <- cbind(geneInfo_toSub, geneInfo_noSub) %>% select(colnames(geneInfo_input))

  # Collapsing duplicated UNIPROT IDs into one gene with XT symbols as synonyms
  duplicates <- rbind(geneInfo_input[duplicated(geneInfo_input$EntryID.uniprot),],
                      geneInfo_input[duplicated(geneInfo_input$EntryID.uniprot, fromLast = TRUE),]) %>%
    unique() %>%
    arrange(Symbol) %>%
    group_by(EntryID.uniprot, Synonym, Reviewed, EntryName, ProteinName, Keywords, Location) %>%
    summarise(Symbol = paste(Symbol, collapse = " ")) %>%
    select(colnames(geneInfo_input)) %>%
    mutate(Synonym = case_when(str_detect(Symbol," ") & !is.na(Synonym) ~ paste(Synonym, gsub("^.*? ", "", Symbol)),
                               str_detect(Symbol," ") & is.na(Synonym) ~ gsub("^.*? ", "", Symbol),
                               str_detect(Symbol," ",negate = TRUE) & is.na(Synonym) ~ NA,
                               TRUE ~ Synonym)) %>%
    mutate(Synonym = gsub("^ ?","",Synonym)) %>%
    mutate(Synonym = gsub("\\s+"," ",Synonym)) %>%
    mutate(Symbol = gsub(" .*$","", Symbol))

  geneInfo_input <- rbind(geneInfo_input, duplicates) %>% unique() %>% arrange(Symbol)
  rm(duplicates)

  synonyms <- geneInfo_input %>%
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

  # Input is not clean, convoluted but trying to deal with a few special cases for completeness
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
  rm(space_synonyms_not_symbol)

  space_synonyms_symbol <- space_synonyms %>%
    filter(str_count(Synonym, " ") != 0) %>%
    mutate(Symbol = toupper(Symbol), Synonym = toupper(Synonym))
  rm(space_synonyms)

  bind <- space_synonyms_symbol[c(1:3, 5:8, 10),] %>% `colnames<-`(c("Symbol", "value"))
  synonyms <- rbind(synonyms, bind)
  rm(bind)

  # Only three remain, easy enough to fix manually
  space_synonyms_symbol <- space_synonyms_symbol[c(4,9,11),] %>%
    mutate(Synonym = gsub("-5", "-5;", Synonym)) %>%
    mutate(Synonym = gsub("FUSION", "FUSION;", Synonym)) %>%
    mutate(Synonym = gsub("-ETO", "-ETO;", Synonym)) %>%
    mutate(Synonym = gsub("-MTG8", "-MTG8;", Synonym)) %>%
    separate_longer_delim(., Synonym, "; ") %>%
    mutate(Synonym = gsub(";","", Synonym)) %>%
    `colnames<-`(c("Symbol", "value"))

  synonyms <- rbind(synonyms, space_synonyms_symbol) %>%
    unique() %>%
    filter(Symbol != value)
  rm(space_synonyms_symbol)

  annotation.map <- gene.map

  synonyms <- synonyms %>%
    left_join(annotation.map, by = c("value" = "human"), relationship = "many-to-many") %>%
    rename(humanAnnot = full) %>%
    left_join(annotation.map, by = c("value" = "other"), relationship = "many-to-many") %>%
    rename(otherAnnot = full) %>%
    left_join(annotation.map, by = c("Symbol" = "human"), relationship = "many-to-many") %>%
    rename(humanAnnot2 = full) %>%
    left_join(annotation.map, by = c("Symbol" = "other"), relationship = "many-to-many") %>%
    rename(otherAnnot2 = full) %>%
    select(-c(human.x, human.y, other.x, other.y)) %>%
    filter(!is.na(humanAnnot) | !is.na(otherAnnot)) %>%
    select(Symbol, value) %>%
    unique()

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write_delim(synonyms, file.path(writeDir, "Symbol_synonym_map.txt"))
  }
  return(synonyms)
}


# make_complex_input - prepare complex input of database
# REQUIRED:
# annotation.map - the output of prep_gene_map that will be used to match human genes to axolotl
#
# OPTIONS:
# prioritize.human - in the event that there are both human and other annotations for a gene, select human. (default FALSE)
# annotation.cutoff - if fewer than this fraction of subunits are matched to axolotl, drop from input.
#   - The default 0.5 means that rows will be dropped where < 50% of input subunits are matched to axolotl.
# check.synonyms - attempt to rescue rows by also testing for synonymous genes? (default FALSE)
#   - Requires synonyms if TRUE; the output of prep_synonym_map
# write.outfile - write tab-delimited complex output file (default TRUE)
# out.dir - directory to output written file to (default working directory)

make_complex_input <- function(annotation.map = NULL,
                               prioritize.human = FALSE,
                               annotation.cutoff = NULL,
                               check.synonyms = FALSE, synonyms = NULL,
                               write.outfile = TRUE, out.dir = NULL){

  if(is.null(annotation.map)){
    stop("You must provide the annotation map output by prep_gene_map.")
  }
  if(is.null(annotation.cutoff)){
    annotation.cutoff <- 0.5
  }
  if(!annotation.cutoff >= 0 & !annotation.cutoff <= 1 ){
    stop("'annotation.cutoff' must be a decimal between 0 and 1.")
  }

  complex_input <- CellChatDB.human$complex
  complex_input[complex_input == ""] <- NA
  complex_input <- complex_input %>%
    rownames_to_column("rowname") %>%
    mutate(rowname = case_when(rowname %in% c("RA-ALDH1A1", "RA-ALDH1A2", "RA-ALDH1A3") ~ gsub("RA-","", rowname), # This is the only special case where this is needed
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
                                    is.na(subunit_1) & is.na(subunit_1a) ~ NA),
             subunit_2b = case_when(!is.na(subunit_2) & is.na(subunit_2a) ~ subunit_2,
                                    !is.na(subunit_2) & !is.na(subunit_2a) ~ paste0(subunit_2,";",subunit_2a),
                                    is.na(subunit_2) & is.na(subunit_2a) ~ NA),
             subunit_3b = case_when(!is.na(subunit_3) & is.na(subunit_3a) ~ subunit_3,
                                    !is.na(subunit_3) & !is.na(subunit_3a) ~ paste0(subunit_3,";",subunit_3a),
                                    is.na(subunit_3) & is.na(subunit_3a) ~ NA),
             subunit_4b = case_when(!is.na(subunit_4) & is.na(subunit_4a) ~ subunit_4,
                                    !is.na(subunit_4) & !is.na(subunit_4a) ~ paste0(subunit_4,";",subunit_4a),
                                    is.na(subunit_4) & is.na(subunit_4a) ~ NA),
             subunit_5b = case_when(!is.na(subunit_5) & is.na(subunit_5a) ~ subunit_5,
                                    !is.na(subunit_5) & !is.na(subunit_5a) ~ paste0(subunit_5,";",subunit_5a),
                                    is.na(subunit_5) & is.na(subunit_5a) ~ NA)) %>%
      select(1, 12:16) %>%
      separate_longer_delim(subunit_1b,";") %>%
      separate_longer_delim(subunit_2b,";") %>%
      separate_longer_delim(subunit_3b,";") %>%
      separate_longer_delim(subunit_4b,";") %>%
      separate_longer_delim(subunit_5b,";") %>%
      unique()
    colnames(complex_input_synonyms) <- colnames(complex_input)
    complex_count = rowSums(!(is.na(complex_input_synonyms[,-1])))
    complex_input_synonyms$start_na <- complex_count # For later to determine what valid "missingness" is to keep after mapping
    complex_input <- complex_input_synonyms
    rm(complex_input_synonyms)
  } else{
    complex_count = rowSums(!(is.na(complex_input[,-1])))
    complex_input$start_na <- complex_count # For later to determine what valid "missingness" is to keep after mapping
  }

  complex_output_1 <- complex_input %>% # Mapping subunits to axolotl genes
    left_join(annotation.map, by=c("subunit_1" = "human"), relationship = "many-to-many") %>%
    rename(subunit_1a = full) %>%
    left_join(annotation.map, by=c("subunit_2" = "human"), relationship = "many-to-many") %>%
    rename(subunit_2a = full)  %>%
    left_join(annotation.map, by=c("subunit_3" = "human"), relationship = "many-to-many") %>%
    rename(subunit_3a = full) %>%
    left_join(annotation.map, by=c("subunit_4" = "human"), relationship = "many-to-many") %>%
    rename(subunit_4a = full) %>%
    left_join(annotation.map, by=c("subunit_5" = "human"), relationship = "many-to-many") %>%
    rename(subunit_5a = full) %>%
    left_join(annotation.map, by=c("subunit_1" = "other"), relationship = "many-to-many") %>%
    rename(subunit_1b = full) %>%
    left_join(annotation.map, by=c("subunit_2" = "other"), relationship = "many-to-many") %>%
    rename(subunit_2b = full)  %>%
    left_join(annotation.map, by=c("subunit_3" = "other"), relationship = "many-to-many") %>%
    rename(subunit_3b = full) %>%
    left_join(annotation.map, by=c("subunit_4" = "other"), relationship = "many-to-many") %>%
    rename(subunit_4b = full) %>%
    left_join(annotation.map, by=c("subunit_5" = "other"), relationship = "many-to-many") %>%
    rename(subunit_5b = full) %>%
    mutate(final_1 = case_when(!is.na(subunit_1a) & is.na(subunit_1b) ~ subunit_1a,
                               !is.na(subunit_1a) & !is.na(subunit_1b) & isTRUE(prioritize.human) ~ subunit_1a,
                               !is.na(subunit_1a) & !is.na(subunit_1b) & isFALSE(prioritize.human) ~ paste0(subunit_1a,";",subunit_1b),
                               is.na(subunit_1a) & !is.na(subunit_1b) ~ subunit_1b,
                               is.na(subunit_1a) & is.na(subunit_1b) ~ NA),
           final_2 = case_when(is.na(subunit_2) ~ NA,
                               !is.na(subunit_2a) & is.na(subunit_2b) ~ subunit_2a,
                               !is.na(subunit_2a) & !is.na(subunit_2b) & isTRUE(prioritize.human) ~ subunit_2a,
                               !is.na(subunit_2a) & !is.na(subunit_2b) & isFALSE(prioritize.human) ~ paste0(subunit_2a,";",subunit_2b),
                               is.na(subunit_2a) & !is.na(subunit_2b) ~ subunit_2b,
                               is.na(subunit_2a) & is.na(subunit_2b) ~ NA),
           final_3 = case_when(is.na(subunit_3) ~ NA,
                               !is.na(subunit_3a) & is.na(subunit_3b) ~ subunit_3a,
                               !is.na(subunit_3a) & !is.na(subunit_3b) & isTRUE(prioritize.human) ~ subunit_3a,
                               !is.na(subunit_3a) & !is.na(subunit_3b) & isFALSE(prioritize.human) ~ paste0(subunit_3a,";",subunit_3b),
                               is.na(subunit_3a) & !is.na(subunit_3b) ~ subunit_3b,
                               is.na(subunit_3a) & is.na(subunit_3b) ~ NA),
           final_4 = case_when(is.na(subunit_4) ~ NA,
                               !is.na(subunit_4a) & is.na(subunit_4b) ~ subunit_4a,
                               !is.na(subunit_4a) & !is.na(subunit_4b) & isTRUE(prioritize.human) ~ subunit_4a,
                               !is.na(subunit_4a) & !is.na(subunit_4b) & isFALSE(prioritize.human) ~ paste0(subunit_4a,";",subunit_4b),
                               is.na(subunit_4a) & !is.na(subunit_4b) ~ subunit_4b,
                               is.na(subunit_4a) & is.na(subunit_4b) ~ NA),
           final_5 = case_when(is.na(subunit_5) ~ NA,
                               !is.na(subunit_5a) & is.na(subunit_5b) ~ subunit_5a,
                               !is.na(subunit_5a) & !is.na(subunit_5b) & isTRUE(prioritize.human) ~ subunit_5a,
                               !is.na(subunit_5a) & !is.na(subunit_5b) & isFALSE(prioritize.human) ~ paste0(subunit_5a,";",subunit_5b),
                               is.na(subunit_5a) & !is.na(subunit_5b) ~ subunit_5b,
                               is.na(subunit_5a) & is.na(subunit_5b) ~ NA)) %>%
    select(c("rowname", "start_na", starts_with("final"))) %>%
    separate_longer_delim(final_1, ";") %>%
    separate_longer_delim(final_2, ";") %>%
    separate_longer_delim(final_3, ";") %>%
    separate_longer_delim(final_4, ";") %>%
    separate_longer_delim(final_5, ";") %>%
    unique() %>%
    `colnames<-`(c("rowname", "start_na", colnames(CellChatDB.human$complex)))

  out_na = rowSums(!(is.na(complex_output_1[,-c(1,2)])))
  complex_output_1$end_na <- out_na
  complex_output_1 <- complex_output_1 %>%
    select(rowname, start_na, end_na, starts_with("subunit_")) %>%
    mutate(na_diff = end_na/start_na) %>%
    filter(na_diff >= annotation.cutoff) %>% # Drop if greater than X fraction of subunits are not annotated
    select(-c(start_na, end_na, na_diff)) %>%
    unique()

  # left-justifying column values
  complex_output_1 <- do.call(rbind,
                              lapply(1:nrow(complex_output_1),
                                     function(x) {z = complex_output_1[x,][order(is.na(complex_output_1[x,]))]; colnames(z)<-colnames(complex_output_1); return(z) } ))

  complex_output_1 <- complex_output_1
  complex_output_1$new_rowname <- make.unique(complex_output_1$rowname,
                                              sep = ".")

  complex_output_1$rowname <- make.unique(complex_output_1$rowname)

  complex_output_1 <- complex_output_1 %>%
    column_to_rownames("new_rowname") %>%
    select(-rowname)

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write.table(complex_output_1, file = file.path(writeDir, "complex_input_axolotl.tsv"), quote = FALSE)
  }

  return(complex_output_1)

}

# make_cofactor_input
# REQUIRED:
# annotation.map - the output of prep_gene_map that will be used to match human genes to axolotl
#
# OPTIONS:
# prioritize.human - in the event that there are both human and other annotations for a gene, select human. (default FALSE)
# annotation.cutoff - if fewer than this fraction of cofactors are matched to axolotl, drop from input. (default 0.5)
# check.synonyms - attempt to rescue rows by also testing for synonymous genes? (default FALSE)
#   - Requires synonyms if TRUE; the output of prep_synonym_map
# write.outfile - write tab-delimited complex output file (default TRUE)
# out.dir - directory to output written file to (default working directory)

make_cofactor_input <- function(annotation.map = NULL,
                                prioritize.human = FALSE,
                                annotation.cutoff = NULL,
                                check.synonyms = FALSE, synonyms = NULL,
                                write.outfile = TRUE, out.dir = NULL){

  select = dplyr::select

  if(is.null(annotation.map)){
    stop("You must provide the annotation map output by prep_gene_map.")
  }
  if(is.null(annotation.cutoff)){
    annotation.cutoff <- 0.5
  }
  if(!annotation.cutoff >= 0 & !annotation.cutoff <= 1 ){
    stop("'annotation.cutoff' must be a decimal between 0 and 1.")
  }

  cofactor_input <- CellChatDB.human$cofactor
  cofactor_input[cofactor_input == ""] <- NA
  cofactor_input <- cofactor_input %>%
    rownames_to_column("rowname")

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
                                    !is.na(cofactor1) & !is.na(cofactor1a) ~ paste0(cofactor1,";",cofactor1a)),
             cofactor2b = case_when(is.na(cofactor2) ~ NA,
                                    !is.na(cofactor2) & is.na(cofactor2a) ~ cofactor2,
                                    !is.na(cofactor2) & !is.na(cofactor2a) ~ paste0(cofactor2,";",cofactor2a)),
             cofactor3b = case_when(is.na(cofactor3) ~ NA,
                                    !is.na(cofactor3) & is.na(cofactor3a) ~ cofactor3,
                                    !is.na(cofactor3) & !is.na(cofactor3a) ~ paste0(cofactor3,";",cofactor3a)),
             cofactor4b = case_when(is.na(cofactor4) ~ NA,
                                    !is.na(cofactor4) & is.na(cofactor4a) ~ cofactor4,
                                    !is.na(cofactor4) & !is.na(cofactor4a) ~ paste0(cofactor4,";",cofactor4a)),
             cofactor5b = case_when(is.na(cofactor5) ~ NA,
                                    !is.na(cofactor5) & is.na(cofactor5a) ~ cofactor5,
                                    !is.na(cofactor5) & !is.na(cofactor5a) ~ paste0(cofactor5,";",cofactor5a)),
             cofactor6b = case_when(is.na(cofactor6) ~ NA,
                                    !is.na(cofactor6) & is.na(cofactor6a) ~ cofactor6,
                                    !is.na(cofactor6) & !is.na(cofactor6a) ~ paste0(cofactor6,";",cofactor6a)),
             cofactor7b = case_when(is.na(cofactor7) ~ NA,
                                    !is.na(cofactor7) & is.na(cofactor7a) ~ cofactor7,
                                    !is.na(cofactor7) & !is.na(cofactor7a) ~ paste0(cofactor7,";",cofactor7a)),
             cofactor8b = case_when(is.na(cofactor8) ~ NA,
                                    !is.na(cofactor8) & is.na(cofactor8a) ~ cofactor8,
                                    !is.na(cofactor8) & !is.na(cofactor8a) ~ paste0(cofactor8,";",cofactor8a)),
             cofactor9b = case_when(is.na(cofactor9) ~ NA,
                                    !is.na(cofactor9) & is.na(cofactor9a) ~ cofactor9,
                                    !is.na(cofactor9) & !is.na(cofactor9a) ~ paste0(cofactor9,";",cofactor9a)),
             cofactor10b = case_when(is.na(cofactor10) ~ NA,
                                     !is.na(cofactor10) & is.na(cofactor10a) ~ cofactor10,
                                     !is.na(cofactor10) & !is.na(cofactor10a) ~ paste0(cofactor10,";",cofactor10a)),
             cofactor11b = case_when(is.na(cofactor11) ~ NA,
                                     !is.na(cofactor11) & is.na(cofactor11a) ~ cofactor11,
                                     !is.na(cofactor11) & !is.na(cofactor11a) ~ paste0(cofactor11,";",cofactor11a)),
             cofactor12b = case_when(is.na(cofactor12) ~ NA,
                                     !is.na(cofactor12) & is.na(cofactor12a) ~ cofactor12,
                                     !is.na(cofactor12) & !is.na(cofactor12a) ~ paste0(cofactor12,";",cofactor12a)),
             cofactor13b = case_when(is.na(cofactor13) ~ NA,
                                     !is.na(cofactor13) & is.na(cofactor13a) ~ cofactor13,
                                     !is.na(cofactor13) & !is.na(cofactor13a) ~ paste0(cofactor13,";",cofactor13a)),
             cofactor14b = case_when(is.na(cofactor14) ~ NA,
                                     !is.na(cofactor14) & is.na(cofactor14a) ~ cofactor14,
                                     !is.na(cofactor14) & !is.na(cofactor14a) ~ paste0(cofactor14,";",cofactor14a)),
             cofactor15b = case_when(is.na(cofactor15) ~ NA,
                                     !is.na(cofactor15) & is.na(cofactor15a) ~ cofactor15,
                                     !is.na(cofactor15) & !is.na(cofactor15a) ~ paste0(cofactor15,";",cofactor15a)),
             cofactor16b = case_when(is.na(cofactor16) ~ NA,
                                     !is.na(cofactor16) & is.na(cofactor16a) ~ cofactor16,
                                     !is.na(cofactor16) & !is.na(cofactor16a) ~ paste0(cofactor16,";",cofactor16a))) %>%
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

  # Need to take care of WNTs specially because there are so many combinations
  cofactor_input_wnts <- cofactor_input %>%
    filter(rowname == "WNT antagonist") %>%
    pivot_longer(.,starts_with("cofactor"), names_to = "num") %>%
    left_join(.,annotation.map, by = c("value" = "human")) %>%
    rename(humanAnnot = full) %>%
    left_join(.,annotation.map, by = c("value" = "other")) %>%
    rename(otherAnnot = full) %>%
    select(rowname, start_na, num, humanAnnot, otherAnnot) %>%
    unique()

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
    rename(cofactor16b = full) %>%
    mutate(final_1 = case_when(!is.na(cofactor1a) & is.na(cofactor1b) ~ cofactor1a,
                               !is.na(cofactor1a) & !is.na(cofactor1b) & isTRUE(prioritize.human) ~ cofactor1a,
                               !is.na(cofactor1a) & !is.na(cofactor1b) & isFALSE(prioritize.human) ~ paste0(cofactor1a,";",cofactor1b),
                               is.na(cofactor1a) & !is.na(cofactor1b) ~ cofactor1b,
                               is.na(cofactor1a) & is.na(cofactor1b) ~ NA),
           final_2 = case_when(is.na(cofactor2) ~ NA,
                               !is.na(cofactor2a) & is.na(cofactor2b) ~ cofactor2a,
                               !is.na(cofactor2a) & !is.na(cofactor2b) & isTRUE(prioritize.human) ~ cofactor2a,
                               !is.na(cofactor2a) & !is.na(cofactor2b) & isFALSE(prioritize.human) ~ paste0(cofactor2a,";",cofactor2b),
                               is.na(cofactor2a) & !is.na(cofactor2b) ~ cofactor2b,
                               is.na(cofactor2a) & is.na(cofactor2b) ~ NA),
           final_3 = case_when(is.na(cofactor3) ~ NA,
                               !is.na(cofactor3a) & is.na(cofactor3b) ~ cofactor3a,
                               !is.na(cofactor3a) & !is.na(cofactor3b) & isTRUE(prioritize.human) ~ cofactor3a,
                               !is.na(cofactor3a) & !is.na(cofactor3b) &isFALSE(prioritize.human) ~ paste0(cofactor3a,";",cofactor3b),
                               is.na(cofactor3a) & !is.na(cofactor3b) ~ cofactor3b,
                               is.na(cofactor3a) & is.na(cofactor3b) ~ NA),
           final_4 = case_when(is.na(cofactor4) ~ NA,
                               !is.na(cofactor4a) & is.na(cofactor4b) ~ cofactor4a,
                               !is.na(cofactor4a) & !is.na(cofactor4b) & isTRUE(prioritize.human) ~ cofactor4a,
                               !is.na(cofactor4a) & !is.na(cofactor4b) & isFALSE(prioritize.human) ~ paste0(cofactor4a,";",cofactor4b),
                               is.na(cofactor4a) & !is.na(cofactor4b) ~ cofactor4b,
                               is.na(cofactor4a) & is.na(cofactor4b) ~ NA),
           final_5 = case_when(is.na(cofactor5) ~ NA,
                               !is.na(cofactor5a) & is.na(cofactor5b) ~ cofactor5a,
                               !is.na(cofactor5a) & !is.na(cofactor5b) & isTRUE(prioritize.human) ~ cofactor5a,
                               !is.na(cofactor5a) & !is.na(cofactor5b) & isFALSE(prioritize.human) ~ paste0(cofactor5a,";",cofactor5b),
                               is.na(cofactor5a) & !is.na(cofactor5b) ~ cofactor5b,
                               is.na(cofactor5a) & is.na(cofactor5b) ~ NA),
           final_6 = case_when(is.na(cofactor6) ~ NA,
                               !is.na(cofactor6a) & is.na(cofactor6b) ~ cofactor6a,
                               !is.na(cofactor6a) & !is.na(cofactor6b) & isTRUE(prioritize.human) ~ cofactor6a,
                               !is.na(cofactor6a) & !is.na(cofactor6b) & isFALSE(prioritize.human) ~ paste0(cofactor6a,";",cofactor6b),
                               is.na(cofactor6a) & !is.na(cofactor6b) ~ cofactor6b,
                               is.na(cofactor6a) & is.na(cofactor6b) ~ NA),
           final_7 = case_when(is.na(cofactor7) ~ NA,
                               !is.na(cofactor7a) & is.na(cofactor7b) ~ cofactor7a,
                               !is.na(cofactor7a) & !is.na(cofactor7b) & isTRUE(prioritize.human) ~ cofactor7a,
                               !is.na(cofactor7a) & !is.na(cofactor7b) & isFALSE(prioritize.human) ~ paste0(cofactor7a,";",cofactor7b),
                               is.na(cofactor7a) & !is.na(cofactor7b) ~ cofactor7b,
                               is.na(cofactor7a) & is.na(cofactor7b) ~ NA),
           final_8 = case_when(is.na(cofactor8) ~ NA,
                               !is.na(cofactor8a) & is.na(cofactor8b) ~ cofactor8a,
                               !is.na(cofactor8a) & !is.na(cofactor8b) & isTRUE(prioritize.human) ~ cofactor8a,
                               !is.na(cofactor8a) & !is.na(cofactor8b) & isFALSE(prioritize.human) ~ paste0(cofactor8a,";",cofactor8b),
                               is.na(cofactor8a) & !is.na(cofactor8b) ~ cofactor8b,
                               is.na(cofactor8a) & is.na(cofactor8b) ~ NA),
           final_9 = case_when(is.na(cofactor9) ~ NA,
                               !is.na(cofactor9a) & is.na(cofactor9b) ~ cofactor9a,
                               !is.na(cofactor9a) & !is.na(cofactor9b) & isTRUE(prioritize.human) ~ cofactor9a,
                               !is.na(cofactor9a) & !is.na(cofactor9b) & isFALSE(prioritize.human) ~ paste0(cofactor9a,";",cofactor9b),
                               is.na(cofactor9a) & !is.na(cofactor9b) ~ cofactor9b,
                               is.na(cofactor9a) & is.na(cofactor9b) ~ NA),
           final_10 = case_when(is.na(cofactor10) ~ NA,
                                !is.na(cofactor10a) & is.na(cofactor10b) ~ cofactor10a,
                                !is.na(cofactor10a) & !is.na(cofactor10b) & isTRUE(prioritize.human) ~ cofactor10a,
                                !is.na(cofactor10a) & !is.na(cofactor10b) & isFALSE(prioritize.human) ~ paste0(cofactor10a,";",cofactor10b),
                                is.na(cofactor10a) & !is.na(cofactor10b) ~ cofactor10b,
                                is.na(cofactor10a) & is.na(cofactor10b) ~ NA),
           final_11 = case_when(is.na(cofactor11) ~ NA,
                                !is.na(cofactor11a) & is.na(cofactor11b) ~ cofactor11a,
                                !is.na(cofactor11a) & !is.na(cofactor11b) & isTRUE(prioritize.human) ~ cofactor11a,
                                !is.na(cofactor11a) & !is.na(cofactor11b) & isFALSE(prioritize.human) ~ paste0(cofactor11a,";",cofactor11b),
                                is.na(cofactor11a) & !is.na(cofactor11b) ~ cofactor11b,
                                is.na(cofactor11a) & is.na(cofactor11b) ~ NA),
           final_12 = case_when(is.na(cofactor12) ~ NA,
                                !is.na(cofactor12a) & is.na(cofactor12b) ~ cofactor12a,
                                !is.na(cofactor12a) & !is.na(cofactor12b) & isTRUE(prioritize.human) ~ cofactor12a,
                                !is.na(cofactor12a) & !is.na(cofactor12b) & isFALSE(prioritize.human) ~ paste0(cofactor12a,";",cofactor12b),
                                is.na(cofactor12a) & !is.na(cofactor12b) ~ cofactor12b,
                                is.na(cofactor12a) & is.na(cofactor12b) ~ NA),
           final_13 = case_when(is.na(cofactor13) ~ NA,
                                !is.na(cofactor13a) & is.na(cofactor13b) ~ cofactor13a,
                                !is.na(cofactor13a) & !is.na(cofactor13b) & isTRUE(prioritize.human) ~ cofactor13a,
                                !is.na(cofactor13a) & !is.na(cofactor13b) & isFALSE(prioritize.human) ~ paste0(cofactor13a,";",cofactor13b),
                                is.na(cofactor13a) & !is.na(cofactor13b) ~ cofactor13b,
                                is.na(cofactor13a) & is.na(cofactor13b) ~ NA),
           final_14 = case_when(is.na(cofactor14) ~ NA,
                                !is.na(cofactor14a) & is.na(cofactor14b) ~ cofactor14a,
                                !is.na(cofactor14a) & !is.na(cofactor14b) & isTRUE(prioritize.human) ~ cofactor14a,
                                !is.na(cofactor14a) & !is.na(cofactor14b) & isFALSE(prioritize.human) ~ paste0(cofactor14a,";",cofactor14b),
                                is.na(cofactor14a) & !is.na(cofactor14b) ~ cofactor14b,
                                is.na(cofactor14a) & is.na(cofactor14b) ~ NA),
           final_15 = case_when(is.na(cofactor15) ~ NA,
                                !is.na(cofactor15a) & is.na(cofactor15b) ~ cofactor15a,
                                !is.na(cofactor15a) & !is.na(cofactor15b) & isTRUE(prioritize.human) ~ cofactor15a,
                                !is.na(cofactor15a) & !is.na(cofactor15b) & isFALSE(prioritize.human) ~ paste0(cofactor15a,";",cofactor15b),
                                is.na(cofactor15a) & !is.na(cofactor15b) ~ cofactor15b,
                                is.na(cofactor15a) & is.na(cofactor15b) ~ NA),
           final_16 = case_when(is.na(cofactor16) ~ NA,
                                !is.na(cofactor16a) & is.na(cofactor16b) ~ cofactor16a,
                                !is.na(cofactor16a) & !is.na(cofactor16b) & isTRUE(prioritize.human) ~ cofactor16a,
                                !is.na(cofactor16a) & !is.na(cofactor16b) & isFALSE(prioritize.human) ~ paste0(cofactor16a,";",cofactor16b),
                                is.na(cofactor16a) & !is.na(cofactor16b) ~ cofactor16b,
                                is.na(cofactor16a) & is.na(cofactor16b) ~ NA)) %>%
    select(c("rowname", "start_na", starts_with("final"))) %>%
    separate_longer_delim(final_1, ";") %>%
    separate_longer_delim(final_2, ";") %>%
    separate_longer_delim(final_3, ";") %>%
    separate_longer_delim(final_4, ";") %>%
    separate_longer_delim(final_5, ";") %>%
    separate_longer_delim(final_6, ";") %>%
    separate_longer_delim(final_7, ";") %>%
    separate_longer_delim(final_8, ";") %>%
    separate_longer_delim(final_9, ";") %>%
    separate_longer_delim(final_10, ";") %>%
    separate_longer_delim(final_11, ";") %>%
    separate_longer_delim(final_12, ";") %>%
    separate_longer_delim(final_13, ";") %>%
    separate_longer_delim(final_14, ";") %>%
    separate_longer_delim(final_15, ";") %>%
    separate_longer_delim(final_16, ";") %>%
    unique() %>%
    `colnames<-`(c("rowname", "start_na", colnames(CellChatDB.human$cofactor)))

  cofactor_input_wnts <- cofactor_input_wnts %>%
    mutate(final = case_when(is.na(humanAnnot) & is.na(otherAnnot) ~ NA,
                             !is.na(humanAnnot) & is.na(otherAnnot) ~ humanAnnot,
                             !is.na(humanAnnot) & !is.na(otherAnnot) & isTRUE(prioritize.human) ~ humanAnnot,
                             !is.na(humanAnnot) & !is.na(otherAnnot) & isFALSE(prioritize.human) ~ paste0(humanAnnot,";",otherAnnot),
                             is.na(humanAnnot) & !is.na(otherAnnot) ~ otherAnnot)) %>%
    select(rowname, start_na, num, final) %>%
    separate_longer_delim(final, ";") %>%
    unique() %>%
    pivot_wider(names_from = num,
                names_sep="",
                values_from = final,
                values_fn = list)

  for( i in 1:ncol(cofactor_input_wnts)){
    col <- colnames(cofactor_input_wnts)[[i]]
    cofactor_input_wnts <- unnest(cofactor_input_wnts, all_of(col))
  }

  cofactor_input_wnts <- cofactor_input_wnts %>%
    select(rowname, start_na, colnames(cofactor_input)) %>%
    unique()

  cofactor_output_1 <- rbind(cofactor_output_1, cofactor_input_wnts)

  rm(cofactor_input, cofactor_input_wnts)

  # left-justifying column values
  cofactor_output_1 <- do.call(rbind,
                               lapply(1:nrow(cofactor_output_1),
                                      function(x) {z = cofactor_output_1[x,][order(is.na(cofactor_output_1[x,]))]; colnames(z)<-colnames(cofactor_output_1); return(z) } ))

  out_na = rowSums(!(is.na(cofactor_output_1[,-c(1,2)])))
  cofactor_output_1$end_na <- out_na
  cofactor_output_1 <- cofactor_output_1 %>%
    select(rowname, start_na, end_na, starts_with("cofactor")) %>%
    mutate(na_diff = end_na/start_na) %>%
    filter(na_diff >= annotation.cutoff) %>% # Drop if greater than X fraction of subunits are not annotated
    select(-c(start_na, end_na, na_diff)) %>%
    unique()

  rownames(cofactor_output_1) <- NULL
  cofactor_output_1$interaction <- make.unique(cofactor_output_1$rowname, sep = ".")
  cofactor_output_1 <- cofactor_output_1 %>%
    column_to_rownames("interaction") %>%
    select(-rowname)

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write.table(cofactor_output_1, file = file.path(writeDir, "cofactor_input_axolotl.tsv"), quote = FALSE)
  }
  return(cofactor_output_1)

}

# make_interaction_input
# REQUIRED:
# annotation.map - the output of prep_gene_map that will be used to match human genes to axolotl
# complex - the output of make_complex_input
# cofactor - the output of make_cofactor_input
#
# OPTIONS:
# prioritize.human - in the event that there are both human and other annotations for a gene, select human. (default FALSE)
# check.synonyms - attempt to rescue rows by also testing for synonymous genes? (default FALSE)
#   - Requires synonyms if TRUE; the output of prep_synonym_map
# write.outfile - write tab-delimited complex output file (default TRUE)
# out.dir - directory to output written file to (default working directory)

make_interaction_input <- function(annotation.map = NULL,
                                   complex = NULL, cofactor = NULL,
                                   prioritize.human = FALSE,
                                   check.synonyms = FALSE, synonyms = NULL,
                                   write.outfile = TRUE, out.dir = NULL){

  if(is.null(annotation.map) | is.null(complex) | is.null(cofactor)){
    stop("Please provide all required inputs.")
  }

  interaction_input <- CellChatDB.human$interaction
  interaction_input[interaction_input == ""] <- NA
  interaction_input <- interaction_input %>%
    rownames_to_column("rowname")

  if(isTRUE(check.synonyms)){
    if(is.null(synonyms)){
      stop("You must provide the Symbol-Synonym map output by prep_synonym_map in order to check for synonyms.")
    }
    # Can only check for synonyms with non-complex ligands/receptors
    interaction_input_synonyms <- interaction_input %>%
      left_join(synonyms, by=c("ligand" = "Symbol"), relationship = "many-to-many") %>%
      rename(ligand1a = value) %>%
      left_join(synonyms, by=c("receptor" = "Symbol"), relationship = "many-to-many") %>%
      rename(receptor1a = value) %>%
      mutate(ligand1b = case_when(!is.na(ligand1a) ~ paste0(ligand, ";", ligand1a),
                                  TRUE ~ ligand),
             receptor1b = case_when(!is.na(receptor1a) ~ paste0(receptor,";",receptor1a),
                                    TRUE ~ receptor)) %>%
      separate_longer_delim(ligand1b, ";") %>%
      separate_longer_delim(receptor1b, ";") %>%
      mutate(ligand = ligand1b, receptor = receptor1b) %>%
      select(colnames(interaction_input)) %>%
      unique()

    interaction_input <- interaction_input_synonyms
    rm(interaction_input_synonyms)

  }

  interaction_input_simple <- interaction_input %>%
    filter(!ligand %in% rownames(complex) & !receptor %in% rownames(complex)) %>%
    left_join(annotation.map, by=c("ligand" = "human"), relationship = "many-to-many") %>%
    rename(ligand_1a = full) %>%
    left_join(annotation.map, by=c("receptor" = "human"),relationship = "many-to-many") %>%
    rename(receptor_1a = full)  %>%
    left_join(annotation.map, by=c("ligand" = "other"), relationship = "many-to-many") %>%
    rename(ligand_1b = full) %>%
    left_join(annotation.map, by=c("receptor" = "other"),relationship = "many-to-many") %>%
    rename(receptor_1b = full) %>%
    mutate(final_ligand = case_when(!is.na(ligand_1a) & is.na(ligand_1b) ~ ligand_1a,
                                    !is.na(ligand_1a) & !is.na(ligand_1b) & isTRUE(prioritize.human) ~ ligand_1a,
                                    !is.na(ligand_1a) & !is.na(ligand_1b) & isFALSE(prioritize.human) ~ paste0(ligand_1a,";",ligand_1b),
                                    is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                    is.na(ligand_1a) & is.na(ligand_1b) ~ NA),
           final_receptor = case_when(!is.na(receptor_1a) & is.na(receptor_1b) ~ receptor_1a,
                                      !is.na(receptor_1a) & !is.na(receptor_1b) & isTRUE(prioritize.human) ~ receptor_1a,
                                      !is.na(receptor_1a) & !is.na(receptor_1b) & isFALSE(prioritize.human) ~ paste0(receptor_1a,";",receptor_1b),
                                      is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                      is.na(receptor_1a) & is.na(receptor_1b) ~ NA)) %>%
    separate_longer_delim(final_ligand, ";") %>%
    separate_longer_delim(final_receptor, ";") %>%
    mutate(ligand = final_ligand, receptor = final_receptor) %>%
    filter(!is.na(ligand) & !is.na(receptor)) %>%
    select(c(colnames(interaction_input))) %>%
    unique()

  interaction_input_complex <- interaction_input %>%
    filter(ligand %in% rownames(complex) | receptor %in% rownames(complex)) %>%
    left_join(annotation.map, by=c("ligand" = "human"), relationship = "many-to-many") %>%
    rename(ligand_1a = full) %>%
    left_join(annotation.map, by=c("receptor" = "human"), relationship = "many-to-many") %>%
    rename(receptor_1a = full) %>%
    left_join(annotation.map, by=c("ligand" = "other"), relationship = "many-to-many") %>%
    rename(ligand_1b = full) %>%
    left_join(annotation.map, by=c("receptor" = "other"), relationship = "many-to-many") %>%
    rename(receptor_1b = full) %>%
    mutate(final_ligand = case_when(!is.na(ligand_1a) & is.na(ligand_1b) ~ ligand_1a,
                                    !is.na(ligand_1a) & !is.na(ligand_1b) & isTRUE(prioritize.human) ~ ligand_1a,
                                    !is.na(ligand_1a) & !is.na(ligand_1b) & isFALSE(prioritize.human) ~ paste0(ligand_1a,";",ligand_1b),
                                    is.na(ligand_1a) & !is.na(ligand_1b) ~ ligand_1b,
                                    is.na(ligand_1a) & is.na(ligand_1b) ~ NA),
           final_receptor = case_when(!is.na(receptor_1a) & is.na(receptor_1b) ~ receptor_1a,
                                      !is.na(receptor_1a) & !is.na(receptor_1b) & isTRUE(prioritize.human) ~ receptor_1a,
                                      !is.na(receptor_1a) & !is.na(receptor_1b) & isFALSE(prioritize.human) ~ paste0(receptor_1a,";",receptor_1b),
                                      is.na(receptor_1a) & !is.na(receptor_1b) ~ receptor_1b,
                                      is.na(receptor_1a) & is.na(receptor_1b) ~ NA)) %>%
    separate_longer_delim(final_ligand, ";") %>%
    separate_longer_delim(final_receptor, ";") %>%
    mutate(ligand = ifelse(!(ligand %in% rownames(complex)), final_ligand, ligand)) %>%
    mutate(receptor = ifelse(!(receptor %in% rownames(complex)), final_receptor, receptor)) %>%
    filter(!is.na(ligand) & !is.na(receptor)) %>%
    select(colnames(interaction_input)) %>%
    unique()

  rowname_match <- complex %>%
    rownames_to_column("rowname") %>%
    mutate(rowname = gsub("\\..$", "", rowname)) # "rowname" is now actually interaction name
  rowname_match$new_rowname <- make.unique(rowname_match$rowname, sep = ".") # "new_rowname" is interaction name made unique
  complex_gene_match <- rowname_match %>% select(rowname, new_rowname)

  # Now mapping ligands and receptors to rownames of complex input
  interaction_input_complex_rownameMerge <- interaction_input_complex %>%
    left_join(complex_gene_match, by=c("ligand" = "rowname"), relationship = "many-to-many") %>%
    rename(ligand_match = new_rowname) %>%
    left_join(complex_gene_match, by=c("receptor" = "rowname"), relationship = "many-to-many") %>%
    rename(receptor_match = new_rowname) %>%
    mutate(final_ligand = case_when((!ligand %in% complex_gene_match$rowname) & is.na(ligand_match) ~ ligand,
                                    (ligand %in% complex_gene_match$rowname) & !is.na(ligand_match) ~ ligand_match),
           final_receptor = case_when((!receptor %in% complex_gene_match$rowname) & is.na(receptor_match) ~ receptor,
                                      (receptor %in% complex_gene_match$rowname) & !is.na(receptor_match) ~ receptor_match)) %>%
    mutate(ligand = final_ligand, receptor = final_receptor) %>%
    filter(!is.na(ligand) & !is.na(receptor)) %>%
    select(colnames(interaction_input_complex)) %>%
    unique()

  interaction_input_bind <- rbind(interaction_input_simple, interaction_input_complex_rownameMerge)

  # Map interaction input ant/agonists, co-receptors to cofactor input
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
    filter(if_all(c(agonist, antagonist, co_A_receptor, co_I_receptor), is.na))

  # Need to map to cofactor input
  interaction_input_cofactor_rownameMerge_withCofactors <- interaction_input_cofactor_rownameMerge %>%
    filter(!if_all(c(agonist, antagonist, co_A_receptor, co_I_receptor), is.na))

  cofactor_gene_match <- as.data.frame(rownames(cofactor))
  colnames(cofactor_gene_match) <- c("row")
  cofactor_gene_match$interaction <- gsub("\\..$", "", cofactor_gene_match$row)

  interaction_input_cofactor_rownameMerge_Full <- interaction_input_cofactor_rownameMerge_withCofactors %>%
    left_join(cofactor_gene_match, by=c("agonist" = "interaction"), relationship = "many-to-many") %>%
    rename(agonist_match = row) %>%
    left_join(cofactor_gene_match, by=c("antagonist" = "interaction"), relationship = "many-to-many") %>%
    rename(antagonist_match = row) %>%
    left_join(cofactor_gene_match, by=c("co_A_receptor" = "interaction"), relationship = "many-to-many") %>%
    rename(co_A_match = row) %>%
    left_join(cofactor_gene_match, by=c("co_I_receptor" = "interaction"), relationship = "many-to-many") %>%
    rename(co_I_match = row) %>%
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
    rename(interaction_name = rowname) %>%
    unique()

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write.table(interaction_output, file = file.path(writeDir, "interaction_input_axolotl.tsv"), quote = FALSE)
  }
  return(interaction_output)
}

# Gene info input
# REQUIRED:
# annotation.map - the output of prep_gene_map
# synonyms - the output of prep_synonym_map
# complex - the output of make_complex_input
# cofactor - the output of make_cofactor_input
# interaction - the output of make_interaction_input
#
# OPTIONAL:
# write.outfile - write tab-delimited complex output file (default TRUE)
# out.dir - directory to output written file to (default working directory)

make_geneInfo_input <- function(annotation.map = NULL, synonyms = NULL,
                                complex = NULL, cofactor = NULL, interaction = NULL,
                                write.outfile = TRUE, out.dir = NULL){

  if(is.null(annotation.map) | is.null(synonyms) | is.null(complex) | is.null(cofactor) | is.null(interaction)){
    stop("Please provide all required inputs.")
  }

  geneInfo_input <- CellChatDB.human$geneInfo %>%
    mutate(Symbol = toupper(Symbol),
           Synonym = toupper(Synonym)) %>%
    separate_longer_delim(Symbol,
                          "; ") %>% # ex. Q8IVW1
    separate_longer_delim(Symbol, ";") %>% # ex. M1V481
    filter(!is.na(EntryID.uniprot))

  geneInfo_toSub <- geneInfo_input[,2:3] # only want to edit these columns
  geneInfo_toSub[] <- lapply(geneInfo_toSub, function(x) gsub(";", "", x))
  geneInfo_toSub[geneInfo_toSub == ""] <- NA # ex. RPL9 - "synonyms" ARE ";"
  geneInfo_toSub <- geneInfo_toSub %>%
    mutate(Synonym = case_when(str_detect(Synonym, "^\\s+$") ~ NA,
                               TRUE ~ Synonym)) %>%
    mutate(Synonym = gsub("\\s+$", "", Synonym))

  geneInfo_noSub <- geneInfo_input[,-c(2:3)] # need to replace these columns with toSub values
  geneInfo_input <- cbind(geneInfo_toSub, geneInfo_noSub) %>% select(colnames(geneInfo_input))

  # Collapsing duplicated UNIPROT IDs into one gene with XT symbols as synonyms
  duplicates <- rbind(geneInfo_input[duplicated(geneInfo_input$EntryID.uniprot),],
                      geneInfo_input[duplicated(geneInfo_input$EntryID.uniprot, fromLast = TRUE),]) %>%
    unique() %>%
    arrange(Symbol) %>%
    group_by(EntryID.uniprot, Synonym, Reviewed, EntryName, ProteinName, Keywords, Location) %>%
    summarise(Symbol = paste(Symbol, collapse = " ")) %>%
    select(colnames(geneInfo_input)) %>%
    mutate(Synonym = case_when(str_detect(Symbol," ") & !is.na(Synonym) ~ paste(Synonym, gsub("^.*? ", "", Symbol)),
                               str_detect(Symbol," ") & is.na(Synonym) ~ gsub("^.*? ", "", Symbol),
                               str_detect(Symbol," ",negate = TRUE) & is.na(Synonym) ~ NA,
                               TRUE ~ Synonym)) %>%
    mutate(Synonym = gsub("^ ?","",Synonym)) %>%
    mutate(Synonym = gsub("\\s+"," ",Synonym)) %>%
    mutate(Symbol = gsub(" .*$","", Symbol))

  geneInfo_input <- rbind(geneInfo_input, duplicates) %>% unique() %>% arrange(Symbol)
  rm(duplicates)

  geneInfo_output <- geneInfo_input %>%
    left_join(annotation.map, by=c("Symbol" = "human"), relationship = "many-to-many") %>%
    rename(human_symbol = Symbol) %>%
    rename(Symbol = full) %>% # Taking human annotations at first pass
    select(colnames(geneInfo_input), human_symbol) %>%
    unique() %>%
    left_join(annotation.map, by=c("human_symbol" = "other"), relationship = "many-to-many") %>%
    rename(NR_Symbol = full) %>% # NR annotations at second pass
    unique() %>%
    filter(!if_all(c(Symbol, NR_Symbol, Synonym), is.na)) %>% # Drop if there's no human or NR annotations and no synonyms
    mutate(Symbol = case_when(is.na(Symbol) & is.na(NR_Symbol) ~ NA,
                              !is.na(Symbol) & is.na(NR_Symbol) ~ Symbol,
                              is.na(Symbol) & !is.na(NR_Symbol) ~ NR_Symbol,
                              !is.na(Symbol) & !is.na(NR_Symbol) ~ paste0(Symbol, ";", NR_Symbol))) %>%
    select(colnames(geneInfo_input), human_symbol) %>%
    separate_longer_delim(Symbol, ";") %>%
    unique() %>%
    filter(!is.na(Symbol)) %>%
    select(colnames(geneInfo_input))

  cofactor.genes <- cofactor %>%
    pivot_longer(., 1:16) %>%
    select(value) %>%
    unique() %>%
    `colnames<-`(c("gene"))
  complex.genes <- complex %>%
    pivot_longer(., 1:5) %>%
    select(value) %>%
    unique() %>%
    `colnames<-`(c("gene"))
  ligands <- interaction %>%
    select(ligand) %>%
    unique() %>%
    filter(!ligand %in% rownames(complex))%>%
    `colnames<-`(c("gene"))
  receptors <- interaction %>%
    select(receptor) %>%
    unique() %>%
    filter(!receptor %in% rownames(complex))%>%
    `colnames<-`(c("gene"))
  all.genes <- rbind(cofactor.genes, complex.genes, ligands, receptors) %>%
    unique() %>%
    drop_na()

  missing <- all.genes[!all.genes$gene %in% geneInfo_output$Symbol,] %>%
    left_join(annotation.map, by = c("gene" = "full")) %>%
    left_join(synonyms, by = c("human" = "value")) %>%
    left_join(synonyms, by = c("other" = "value"))

  missing[missing == ""] <- NA

  missing <- missing %>%
    mutate(final_synonym = case_when(!is.na(human) & is.na(other) ~ human,
                                     is.na(human) & !is.na(other) ~ other,
                                     !is.na(human) & !is.na(other) ~ paste0(human,";",other)),
           final_symbol = case_when(!is.na(Symbol.x) & is.na(Symbol.y) ~ Symbol.x,
                                    is.na(Symbol.x) & !is.na(Symbol.y) ~ Symbol.y,
                                    !is.na(Symbol.x) & !is.na(Symbol.y) ~ paste0(Symbol.x,";",Symbol.y))) %>%
    separate_longer_delim(final_synonym, ";") %>%
    separate_longer_delim(final_symbol, ";") %>%
    select(gene, final_synonym, final_symbol) %>%
    unique() %>%
    filter(final_synonym %in% synonyms$value) %>%
    left_join(., geneInfo_input, by = c("final_symbol" = "Symbol"), relationship = "many-to-many") %>%
    filter(!is.na(Synonym)) %>%
    mutate(Symbol = case_when(str_detect(Synonym, " ", negate = TRUE) ~ gene,
                              TRUE ~ final_symbol),
           Synonym = case_when(str_detect(Synonym," ", negate = TRUE) ~ final_symbol,
                               TRUE ~ Synonym)) %>%
    mutate(Synonym = case_when(str_detect(Synonym, " ") ~ paste0(str_replace_all(Synonym, final_synonym, ""), Symbol),
                               TRUE ~ Synonym),
           Symbol = gene) %>%
    select(colnames(geneInfo_input))

  geneInfo_output <- rbind(geneInfo_output, missing) %>% unique()

  if(isTRUE(write.outfile)){
    writeDir <- out.dir
    if(is.null(writeDir)){
      writeDir <- "./"
    }
    write.table(geneInfo_output, file = file.path(writeDir, "geneInfo_input_axolotl.tsv"), quote = FALSE)
  }
  return(geneInfo_output)
}

# Assemble prepared inputs into final database
# REQUIRED:
# complex - the output of make_complex_input
# interaction - the output of make_interaction_input
# cofactor - the output of make_cofactor_input
# geneInfo - the output of make_geneInfo_input

assemble_database <- function(complex = NULL,
                              interaction = NULL,
                              cofactor = NULL,
                              geneInfo = NULL){

  if(is.null(complex) | is.null(interaction) | is.null(cofactor) | is.null(geneInfo)){
    stop("You must provide all four inputs to assemble the database.")
  }

  options(stringsAsFactors = FALSE)

  complex[is.na(complex)] <- ""

  interaction[is.na(interaction)] <- ""
  interaction$interaction_name <- rownames(interaction)

  cofactor[is.na(cofactor)] <- ""

  database <- list(complex = complex,
                   interaction = interaction,
                   cofactor = cofactor,
                   geneInfo = geneInfo)

  return(database)

}
