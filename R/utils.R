
# Copy-paste of the code from {wbData}, so we don't need the package
i2s <- function (gene_id, geneIDs, warn_missing = TRUE) {

  res <- geneIDs$name[match(gene_id, geneIDs$gene_id, incomparables = NA)]
  if (warn_missing && any(is.na(res))) {
    warning("i2s: ", sum(is.na(res)), " gene IDs could not be converted. NA are returned.")
  }
  res
}



replace_gene_ids <- function(x, gids) {
  key <- tolower(x)
  idx <- match(key, tolower(gids$name))                     # "hbl-1", "HBL-1"
  idx[is.na(idx)] <- match(key[is.na(idx)], tolower(gids$sequence))  # "F13D11.2"
  idx[is.na(idx)] <- match(key[is.na(idx)], tolower(gids$gene_id))   # "WBGene00001851"
  gids$name[idx]                                            # NA if nothing matched
}
