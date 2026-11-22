
# Copy-paste of the code from {wbData}, so we don't need the package
i2s <- function (gene_id, geneIDs, warn_missing = TRUE) {

  res <- geneIDs$name[match(gene_id, geneIDs$gene_id, incomparables = NA)]
  if (warn_missing && any(is.na(res))) {
    warning("i2s: ", sum(is.na(res)), " gene IDs could not be converted. NA are returned.")
  }
  res
}



replace_gene_ids <- function(gene_list, gids){

  is_id <- !(gene_list %in% gids$name) & (gene_list %in% gids$gene_id)

  gene_list[is_id] <- i2s( gene_list[is_id], gids )

  gene_list

}
