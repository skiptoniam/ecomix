#' @title data.prep objects
#' @rdname data.prep
#' @name table2pam
#' @description prepares table of occurrence records or precence-absence matrix for use in \code{bbgdm}.
#' @param x data to be imported into R. A table of species occurrences and site records.
#' @param site.id Name of column that contains site names /IDs.
#' @param sp.id Name of column that contains species information/ID.
#' @param abund LOGICAL if TRUE calculate abundances.
#' @param abund.col Name of column that contains abundance/counts data.
#' @param siteXsp LOGICAL if TRUE Returns sites as rows, sites as columns, if FALSE returns species as rows and sites as cols.
#' @return A species by sites matrix for use in dissimilarity calculation
#' @export


table2pam <- function (x, site.id = "site.id", sp.id = "sp.id", abund = FALSE, abund.col = " No.of.specimens",siteXsp=TRUE)
{
  a <- site.id
  nr <- length(levels(as.factor(x[, a])))
  rn <- levels(as.factor(x[, a]))
  z <- sp.id
  cn <- levels(as.factor(x[, z]))
  nc <- length(cn)
  nm <- matrix(0, nr, nc, dimnames = list(rn, cn))
  for (i in 1:length(x[, 1])) {
    m <- as.character(x[i, a])
    n <- as.character(x[i, z])
    if (is.na(m) == TRUE | is.null(m) == TRUE | is.na(n) ==
        TRUE | is.null(n) == TRUE)
      (next)(i)
    if (m == "" | m == " " | n == "" | n == " ")
      (next)(i)
    if (abund == TRUE)
      nm[m, n] <- nm[m, n] + x[i, abund.col]
    else nm[m, n] <- 1
  }
  fm <- nm[rowSums(nm) > 0, ]
  if(siteXsp){ return(as.matrix(fm))
  } else {
    return(as.matrix(t(fm)))
  }
}

