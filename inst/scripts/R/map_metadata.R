#!usr/bin/env R

# Author: Sean Maden
#
# Functions for metadata mapping with regular expressions
#
# Notes:
# * For application examples, see following scripts:
#     -- recount3_dataset_properties.R 

# 

map_libprep <- function(cdi, 
                        libprep.cnv = c("sra.library_construction_protocol", 
                                        "sra.study_abstract"), 
                        grepl.vectors = c(grepl.pattv, 
                                          grepl.pattv.polya, 
                                          grepl.pattv.ribozero), 
                        cnv = c("libprep.available", "pa.true", "rz.true")){
  # map_libprep
  # 
  # Gets the libprep condition outcome from srp coldata.
  # 
  # arguments:
  # cdi : coldata from rse
  # libprep.cnv : vector of variables, colnames in cdi to inspect
  # grepl.vectors : grepl pattern vectors corresponding to cnv items
  # cnv : new variable names corresponding to grepl.vectors items
  # 
  # returns:
  # dfi, data.frame of condition outcomes
  # 
  dfi <- do.call(cbind, lapply(grepl.vectors, function(gi){
    apply(do.call(cbind, lapply(libprep.cnv, function(ci){
      grepl(gi, cdi[,ci])
    })), 1, function(ri){
      ifelse(length(which(ri))>0, TRUE, FALSE)
    })
  }))
  colnames(dfi) <- cnv
  return(dfi)
}


get_grepl <- function(strv){
  paste0(unlist(lapply(strv, function(ii){
    paste0(".*", ii, ".*")})), collapse = "|")
}
grepl.strv.polya <- c("poly\\(A\\)", "polyA", "poly-A", "poly A")
grepl.strv.ribozero <- c("ribo z", "riboz", "riboZ", "RiboZ")
grepl.strv <- c(grepl.strv.polya, grepl.strv.ribozero)
grepl.pattv <- get_grepl(c(grepl.strv.polya, grepl.strv.ribozero))
grepl.pattv.polya <- get_grepl(grepl.strv.polya)
grepl.pattv.ribozero <- get_grepl(grepl.strv.ribozero)


do_mappings <- function(){
  lvar <- list(cell_info = c("cell"),
               disease_info = c("disease","pathology","case"),
               sex = c("sex","gender"), tissue = c("tissue","region"),
               demographics = c("race","ethnicity"), age = c("age"),
               identifier = c("id", "ID"), death_details = c("death"),
               library = c("library"))
  
  # get regex patterns
  start.str <- ".*"; end.str <- ".*"
  lgrep <- lapply(lvar, function(li){
    paste0(unlist(lapply(li, function(ii){
      charv <- unlist(strsplit(ii, "")); letter <- charv[1] # get first letter
      which.letter <- paste0("(", letter, "|", toupper(letter), ")") # 1st chr cap
      notfirst.str <- paste0(charv[2:length(charv)], collapse = "") # chr remain
      char.str <- paste0(which.letter, notfirst.str) # new str
      start.str <- ifelse(ii == "age", paste0(start.str, "(^|\\_|!s)"), start.str)
      paste0(start.str, char.str, end.str, collapse = "")
    })), collapse = "|")
  })
  
  # get mappings
  # get colname mappings by srp
  lmap <- lapply(lsr, function(sri){
    df.sri <- sri[,grepl("^sra_attribute.*",colnames(sri))]
    cnv <- colnames(df.sri)
    lmapi <- lapply(seq(length(lgrep)), function(ii){
      message(ii)
      varname <- names(lgrep)[ii]; patti <- lgrep[[ii]]
      valuev <- rep("NA", nrow(df.sri)) # default variable values
      cnv[grepl(patti, cnv)]
    })
    names(lmapi) <- names(lgrep)
    lmapi
  })
  names(lmap) <- names(lsr)
  # get mapped values by matched colnames
  dfmap <- do.call(rbind, lapply(seq(length(lsr)), function(ii){
    srpi <- names(lsr)[ii]; sri <- lsr[[ii]]
    df.sri <- sri[,grepl("^sra_attribute.*",colnames(sri))]
    cnv <- colnames(df.sri)
    dfmapi <- do.call(cbind, lapply(seq(length(lgrep)), function(ii){
      varname <- names(lgrep)[ii]; patti <- lgrep[[ii]]
      valuev <- rep("NA", nrow(df.sri)) # default variable values
      which.cn <- grepl(patti, cnv)
      # replace valuev if any matching colnames
      if(length(which.cn) > 0){ 
        # uniform formatting of matching values
        valuev <- as.character(apply(df.sri[,which.cn,drop=F], 1, function(ri){
          paste0(gsub("'", "", 
                      gsub(" ", "_", tolower(ri))),
                 collapse = ";")
        }))
      }
      dfi <- data.frame(var=valuev); colnames(dfi) <- varname; dfi
    }))
    dfmapi$srp <- srpi; dfmapi$run_id <- sri$external_id; dfmapi
  }))
  
  # save mappings
  lmd <- list(varname.mappings = lmap, dfmap = dfmap)
  lmd.fname <- "list-metadata_dorsolateral-human_recount3.rda"
  save(lmd, file = file.path(dir.name, lmd.fname))
}