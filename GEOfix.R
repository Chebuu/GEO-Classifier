library(GEOquery)
getDirListing <- function(url) {
  # https://github.com/seandavi/GEOquery/issues/38
  
  # alexvpickering commented on Jul 31, 2017 • 
  # EDIT: now compatible with getGEOSuppFiles (fixing same error)
  
  # I think the issue is with getDirListing. When processing HTML content, 
  # there can be links other than the series_matrix.txt.gz that seem to be 
  # related to determining the file size. In ronammar's case, the implicated link is 
  # '/geo/series/GSE94nnn/GSE94802/'. As a result, a download of 
  # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94802/matrix//geo/series/GSE94nnn/GSE94802/ 
  # is attempted and fails because it doesn't exist. To fix, just make sure the only links returned 
  # don't end with a forward slash (this approach also fixes the same download errors for getGEOSuppFiles). 
  # The following fixed it for me (added line ending with # !!!):
  
  message(url)
  # Takes a URL and returns a character vector of filenames
  a <- RCurl::getURL(url)
  if( grepl("<HTML", a, ignore.case=T) ){ # process HTML content
    doc <- XML::htmlParse(a)
    links <- XML::xpathSApply(doc, "//a/@href")
    links <- links[!grepl('/$', links)]   # !!!
    XML::free(doc)
    b <- as.matrix(links)
    message('OK')
  } else { # standard processing of txt content
    tmpcon <- textConnection(a, "r")
    b <- read.table(tmpcon)
    close(tmpcon)
  }
  b <- as.character(b[,ncol(b)])
  return(b)
}

getGEO.simple <- function(GEO = NULL, destdir = tempdir(), GSElimits = NULL, AnnotGPL = FALSE, getGPL = TRUE) {
  # TODO: Args should be passed to getGEOfile. Some parameters from the original functions are not used.
  # This function replaces getGEO after setting getDirListing above.
  return(
    parseGEO(
      getGEOfile(GEO), 
      GSElimits, 
      destdir, 
      AnnotGPL, 
      getGPL
    )
  )
}

# EXAMPLE
# gfile <- getGEO.simple('GSE72078')
# pData(gfile)
# GDS2eSet(gfile)
# gfile@gsms
# names(Meta(gfile))

