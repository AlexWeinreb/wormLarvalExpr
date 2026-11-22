# Prepare dataset
# download it if needed, normalize, save


# Inits ----


library(tidyverse)


# for lightweight worm ID conversions
# pak::pak("AlexWeinreb/wbData")



#~ Obtain data ----


# dir.create("data-raw")
# download.file(
#   "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130811/suppl/GSE130811_expr_mRNA_CE10_coding.tab.gz",
#   destfile = "data-raw/GSE130811_expr_mRNA_CE10_coding.tab.gz"
# )

gids_298 <- wbData::wb_load_gene_ids(298)



#~ load ----


dat <- read_tsv("data-raw/GSE130811_expr_mRNA_CE10_coding.tab.gz",
                col_names = c("gene_id",
                              read_lines("data-raw/GSE130811_expr_mRNA_CE10_coding.tab.gz", n_max = 1L) |> str_split_1("\\t")),
                skip = 1) |>
  select(-ends_with(".2")) |>
  pivot_longer(cols = -c(1,2),
               names_to = "time",
               names_pattern = c("^([0-9]{1,2})hr$")) |>
  mutate(time = time |> str_remove("hr$") |> as.numeric(),
         gene_name = wbData::i2s(gene_id, gids_298))









# Normalization ----

# TPM

df2mat <- function(df, xcol, ycol, valcol){
  df |>
    select({{xcol}}, {{ycol}}, {{valcol}}) |>
    pivot_wider(names_from = {{ycol}}, values_from = {{valcol}}) |>
    column_to_rownames(xcol) |>
    as.matrix()
}


norm_factors <- tibble(
  time = unique(dat$time),
  tmm_factor = dat |>
    df2mat("gene_id", "time", "value") |>
    edgeR::DGEList() |>
    edgeR::calcNormFactors(method = "TMM") |>
    chuck("samples", "norm.factors")
)

dat_norm <- dat |>
  left_join(norm_factors,
            by = "time") |>
  mutate(nf = sum(value / width), .by = time) |>
  mutate(
    TPM = 1e6 * value / (width * nf),
    TPM_TMM = TPM / tmm_factor
  ) |>
  select(gene_name, gene_id, time, TPM_TMM)


saveRDS(dat_norm,
        "data/dat_norm.rds")
saveRDS(gids_298,
        "data/gids_298.rds")


