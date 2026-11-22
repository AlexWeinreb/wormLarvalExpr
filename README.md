

## About wormLarvalExpr

This app plots bulk RNA-seq expression of *C. elegans* genes across larval
development (L1 to adult). Enter one or more genes, separated by commas, as
gene names (e.g. `mam-2`), sequence names (e.g. `C26B9.7`) or WormBase IDs.

**Note:** many genes show a dip at 40 h. It appears in both replicates
and probably reflects a technical batch effect rather than biology

### Data

All expression data come from:

> Meeuse MWM, Hauser YP, Morales Moya LJ, et al. (2020). Developmental function
> and state transitions of a gene expression oscillator in *Caenorhabditis
> elegans*. *Molecular Systems Biology* 16: e9498.
> [doi:10.15252/msb.20209498](https://doi.org/10.15252/msb.20209498)

The data belong to the original authors. Please cite their paper if you use
them.

Larval stage boundaries are approximate.

### Preprocessing

Expression values were taken from GEO series
[GSE130811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130811),
file `GSE130811_expr_mRNA_CE10_coding.tab.gz` (protein-coding genes, ce10
assembly). The last seven time points were sequenced in duplicate; only the first
replicate is kept here.

For each sample, counts were divided by gene length and rescaled to sum to one
million (TPM). TMM normalization factors were then computed with edgeR
(`calcNormFactors()`), and each sample's TPM was divided by its factor.
WormBase gene IDs were converted to gene names with WormBase release WS298
(through [wbData](https://github.com/AlexWeinreb/wbData)). The plot shows
log10(1 + TPM) against developmental time in hours.


### Source code

Full source code is available at:
<https://github.com/AlexWeinreb/wormLarvalExpr>

See `data-raw/` for the preprocessing steps.
