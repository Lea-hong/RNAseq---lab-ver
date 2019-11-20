#tximport

tximport <- function(files,
                     type=c("none","salmon","sailfish","kallisto","rsem","stringtie"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=c("no","scaledTPM","lengthScaledTPM"),
                     tx2gene=NULL,
                     varReduce=FALSE,
                     dropInfReps=FALSE,
                     ignoreTxVersion=FALSE,
                     ignoreAfterBar=FALSE,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     importer=NULL,
                     existenceOptional=FALSE,
                     readLength=75) {
  
  # inferential replicate importer
  infRepImporter <- NULL
  
  type <- match.arg(type, c("none","salmon","sailfish","kallisto","rsem","stringtie"))
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  
  if (!existenceOptional) stopifnot(all(file.exists(files)))
  if (!txIn & txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
  
  stopifnot(length(files) > 0)
  kallisto.h5 <- basename(files[1]) == "abundance.h5"
  if (type == "kallisto" & !kallisto.h5) {
    message("Note: importing `abundance.h5` is typically faster than `abundance.tsv`")
  }
  
  if (type=="rsem" & txIn & grepl("genes", files[1])) {
    message("It looks like you are importing RSEM genes.results files, setting txIn=FALSE")
    txIn <- FALSE
  }
  
  readrStatus <- FALSE
  if (is.null(importer) & !kallisto.h5) {
    if (!requireNamespace("readr", quietly=TRUE)) {
      message("reading in files with read.delim (install 'readr' package for speed up)")
      importer <- read.delim
    } else {
      message("reading in files with read_tsv")
      readrStatus <- TRUE
    }
  }
  
  # salmon/sailfish presets
  if (type %in% c("salmon","sailfish")) {
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    if (readrStatus & is.null(importer)) {
      col.types <- readr::cols(
        readr::col_character(),readr::col_integer(),readr::col_double(),readr::col_double(),readr::col_double()
      )
      importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
    }
    infRepImporter <- if (dropInfReps) { NULL } else { function(x) readInfRepFish(x, type) }
  }
  
  # kallisto presets
  if (type == "kallisto") {
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    if (kallisto.h5) {
      importer <- read_kallisto_h5
    } else if (readrStatus & is.null(importer)) {
      col.types <- readr::cols(
        readr::col_character(),readr::col_integer(),readr::col_double(),readr::col_double(),readr::col_double()
      )
      importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
    }
    infRepImporter <- if (dropInfReps) { NULL } else { readInfRepKallisto }
  }
  
  # rsem presets
  if (type == "rsem") {
    if (txIn) {
      txIdCol <- "transcript_id"
      abundanceCol <- "TPM"
      countsCol <- "expected_count"
      lengthCol <- "effective_length"
      if (readrStatus & is.null(importer)) {
        col.types <- readr::cols(
          readr::col_character(),readr::col_character(),readr::col_integer(),readr::col_double(),
          readr::col_double(),readr::col_double(),readr::col_double(),readr::col_double()
        )
        importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
      }
    } else {
      geneIdCol <- "gene_id"
      abundanceCol <- "TPM"
      countsCol <- "expected_count"
      lengthCol <- "effective_length"
      if (readrStatus & is.null(importer)) {
        col.types <- readr::cols(
          readr::col_character(),readr::col_character(),readr::col_double(),readr::col_double(),
          readr::col_double(),readr::col_double(),readr::col_double()
        )
        importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
      }
    }
  }
  
  if (type == c("stringtie")) {
    txIdCol <- "t_name"
    geneIdCol <- "gene_name"
    abundanceCol <- "FPKM"
    countsCol <- "cov"
    lengthCol <- "length"
    if (readrStatus & is.null(importer)) {
      col.types <- readr::cols(
        readr::col_character(),readr::col_character(),readr::col_character(),readr::col_integer(),readr::col_integer(),readr::col_character(),readr::col_integer(),readr::col_integer(),readr::col_character(),readr::col_character(),readr::col_double(),readr::col_double()
      )
      importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
    }
  }
  
  infRepType <- "none"
  if (type %in% c("salmon", "sailfish", "kallisto") & !dropInfReps) {
    infRepType <- if (varReduce) { "var" } else { "full" }
  }
  
  # if input is tx-level (this is every case but RSEM gene.results files)
  if (txIn) {
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      
      raw <- as.data.frame(importer(files[i]))
      
      # if we expect inferential replicate info
      repInfo <- NULL
      if (infRepType != "none") {
        repInfo <- infRepImporter(dirname(files[i]))
        # if we didn't find inferential replicate info
        if (is.null(repInfo)) {
          infRepType <- "none"
        }
      }
      
      # if external tx2gene table not provided, send user to vignette
      if (is.null(tx2gene) & !txOut) {
        summarizeFail() # ...long message in helper.R
      } else {
        # e.g. Salmon and kallisto do not include the gene ID, need an external table
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          txId <- raw[[txIdCol]]
        } else {
          stopifnot(all(txId == raw[[txIdCol]]))
        }
      }
      
      # create empty matrices
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[txIdCol]]
        colnames(mat) <- names(files)
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
        if (infRepType == "var") {
          varMatTx <- mat
        } else if (infRepType == "full") {
          infRepMatTx <- list()
        }
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
      if (infRepType == "var") {
        varMatTx[,i] <- repInfo$vars
      } else if (infRepType == "full") {
        infRepMatTx[[i]] <- repInfo$reps
      }
    }
    
    # propagate names to inferential replicate list
    if (infRepType == "full") {
      names(infRepMatTx) <- names(files)
    }
    
    message("")
    
    # if there is no information about inferential replicates
    if (infRepType == "none") {
      txi <- list(abundance=abundanceMatTx,
                  counts=countsMatTx,
                  length=lengthMatTx,
                  countsFromAbundance=countsFromAbundance)
    } else if (infRepType == "var") {
      # if we're keeping only the variance from inferential replicates
      txi <- list(abundance=abundanceMatTx,
                  counts=countsMatTx,
                  variance=varMatTx,
                  length=lengthMatTx,
                  countsFromAbundance=countsFromAbundance)
    } else if (infRepType == "full") {
      # if we're keeping the full samples from inferential replicates
      txi <- list(abundance=abundanceMatTx,
                  counts=countsMatTx,
                  infReps=infRepMatTx,
                  length=lengthMatTx,
                  countsFromAbundance=countsFromAbundance)
    }
    
    # stringtie outputs coverage, here we turn into counts
    if (type == "stringtie") {
      # here "counts" is still just coverage, this formula gives back original counts
      txi$counts <- txi$counts * txi$length / readLength
    }
    
    if (type == "rsem") {
      # protect against 0 bp length transcripts
      txi$length[txi$length < 1] <- 1
    }
    
    # if the user requested just the transcript-level data, return it now
    if (txOut) {
      if (countsFromAbundance != "no") {
        txi$counts <- makeCountsFromAbundance(txi$counts, txi$abundance, txi$length, countsFromAbundance)
      }
      return(txi)
    }
    
    
    # otherwise, summarize to the gene-level
    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion, ignoreAfterBar, countsFromAbundance)
    return(txiGene)  
    
    
    # else, not txIn...
  } else {
    # RSEM already has gene-level summaries
    # so we just combine the gene-level summaries across files
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      out <- capture.output({
        raw <- as.data.frame(importer(files[i]))
      }, type="message")
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% names(raw)))
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[geneIdCol]]
        colnames(mat) <- names(files)
        abundanceMat <- mat
        countsMat <- mat
        lengthMat <- mat
      }
      abundanceMat[,i] <- raw[[abundanceCol]]
      countsMat[,i] <- raw[[countsCol]]
      lengthMat[,i] <- raw[[lengthCol]]
    }
  } 
  message("")
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance="no"))
}