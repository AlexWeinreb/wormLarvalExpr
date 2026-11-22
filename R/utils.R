
replace_gene_ids <- function(x, gids) {

  key <- tolower(x)
  idx <- match(key, tolower(gids$symbol))
  idx[is.na(idx)] <- match(key[is.na(idx)], tolower(gids$sequence))
  idx[is.na(idx)] <- match(key[is.na(idx)], tolower(gids$gene_id))

  gids$name[idx] # NA if nothing matched
}
