
library(stringr)
library(Seurat)
library(R6)

bulkObject <- R6Class("BulkObject", list(
   colname = NULL, 
   cell_type = NULL,
   sample = NULL,
   debug = FALSE,
   name = NULL,
   rs = NULL,
   initialize = function(sample, cell_type, colname = "MinorSep",debug = FALSE) {
    stopifnot(is.character(cell_type), length(cell_type) == 1)
    stopifnot(is.character(sample), length(sample) == 1)
    
    self$sample <- sample
    self$cell_type <- cell_type
    self$colname <- colname
    self$debug <- debug
  },
  print = function(...) {
    cat("BulkObject: \n")
    cat("  Sample:    ", self$sample, "\n", sep = "")
    cat("  Cell Type: ", self$cell_type, "\n", sep = "")
    cat("  Column:    ", self$colname, "\n", sep = "")
    cat("  Name:      ", self$name, "\n", sep = "")
    invisible(self)
  },
  bulk = function(data, cells){
    c <- rownames(cells)[cells[,self$colname] %in% self$cell_type]
    if ( self$debug ){
      print(paste("c:",length(c)))
    }
    data2 <- data[,grepl(self$sample,colnames(data))]
    if ( self$debug ){
      print(paste("data2:",dim(data2)))
    }
    data2a <- data2[,colnames(data2) %in% c]
    if ( is.vector(data2a)) {
      self$rs <- data2a
    } else if ( dim(data2a)[2] == 0){
      self$rs <- rep.int(0, dim(data2a)[1])
    } else {
      self$rs <- rowSums(data2a)
    }
    self$name <- paste(self$cell_type, self$sample, sep="_")
    #ob <- list("rs" = rs2, "name" = name)
    #return(ob)
  }
 )
)   

SeuratPseudoBulk <- R6Class("SeuratPseudoBulk", list(
   inS = NULL,
   colname = NULL,
   outS = NULL,
   data = NULL,
   cells = NULL,
   rs = NULL,
   debug = FALSE,
   initialize = function(inS, colname = "MinorSep", debug = FALSE) {
    # stopifnot(is(self$inS,"Seurat"))
    
    self$inS <- inS
    self$colname <- colname
    self$debug <- debug
  },
  pseudobulk = function(){
    samples <- unique(str_split(colnames(self$data),"_",simplify=TRUE)[,1])
    cell_types <- unique(self$cells[,self$colname])
  
    message(samples[1])
    ob <- bulkObject$new(samples[1], cell_types[1], self$colname, self$debug)
  
    ob$bulk(self$data, self$cells)
    # rs <- unlist(ob["rs"])
    # names <- ob["name"]
    self$rs <- unlist(ob$rs)
    if ( self$debug ){
      print(paste("rs:",dim(self$rs)))
    }
    names <- ob$name
    print(ob)
  
    for ( y in seq(2,length(cell_types))){
      if ( self$debug ){ message(paste(y, cell_types[y], sep=": ")) }
      ob <- bulkObject$new(samples[1], cell_types[y], self$colname, self$debug)
      ob$bulk(self$data, self$cells)
      # ob <- .bulk(data, cells, samples[1], cell_types[y], colname)
      # names <- append(names,ob["name"])
      # rs <- cbind(rs,unlist(ob["rs"]))
      names <- append(names,ob$name)
      self$rs <- cbind(self$rs,unlist(ob$rs))
      if ( self$debug){ print(paste("rs:",dim(self$rs)))}
    }
  
    for ( i in seq(2,length(samples)) ){
      message(samples[i])
      for ( y in seq(1,length(cell_types))){
        ob <- bulkObject$new(samples[i], cell_types[y], self$colname, self$debug)
        ob$bulk(self$data, self$cells)
        names <- append(names,ob$name)
        self$rs <- cbind(self$rs,unlist(ob$rs))
        if ( self$debug){ print(paste("rs:",dim(self$rs)))}
      #  ob <- .bulk(data, cells, samples[i], cell_types[y], colname)
      #  names <- append(names,ob["name"])
      #  rs <- cbind(rs,unlist(ob["rs"]))
      }
    }
    colnames(self$rs) <- names
    rownames(self$rs) <- gsub("^rs.","",rownames(self$rs))
    self$rs <- as.matrix(self$rs)
    #return(rs)
  },
  seurat_pseudobulk = function(){
  
    if ( !is(self$inS,"Seurat")){
      stop("input is not a Seurat object!")
    }
  
    # Get meta data
    self$cells <- self$inS@meta.data
    # rownames(cells) <- cells$X
  
    # Get count data
    self$data <- GetAssayData(self$inS, slot="counts")
  
    # pseudo bulk counts
    self$pseudobulk()
  
    # create new metadata
    names <- colnames(self$rs)
    meta <- cbind(str_split(names,"_",simplify=TRUE)[,2], str_split(names,"_",simplify=TRUE)[,1])
    meta <- as.data.frame(meta)
    rownames(meta) <- names
    colnames(meta) <- c("Sample_ID", "Cell_type")
  
    # Create new Seurat object
    self$outS <- CreateSeuratObject(self$rs, meta.data = meta)
  
    # return(outS)
  }

 )
)

frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
wd <- dirname(frame_files[[length(frame_files)]])

seuratClass <- readRDS(paste0(wd,"/scRNASeq.Rda"))

pseudoClass <- SeuratPseudoBulk$new(seuratClass,"SingleR.labels")
pseudoClass$seurat_pseudobulk()

saveRDS(pseudoClass$outS,paste0(wd,"/scRNASeq_Pseudo_test2.Rds"))


