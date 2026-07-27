# =============================================================================
# GENCODE gene models: spans (one row per gene) and exons (merged per gene).
#
# Shared by tools/make_gene_models.R (regenerate the committed reference tables)
# and Fig5/00_setup.R (resolve genes for a fine-mapping window at build time).
#
# Both callers want the same 5-column shape, which is what Fig5's gene-model
# track reads:  gene, chrom, start, end, strand   (1-based inclusive, as GENCODE
# reports it; headerless on disk).
#
# Exons are reduced per gene: GENCODE lists them once PER TRANSCRIPT, so without
# merging a gene's shared exons are drawn dozens of times on top of each other.
#
# rtracklayer/GenomicRanges are used namespaced and loaded lazily -- attaching
# rtracklayer lets S4Vectors/BiocGenerics mask dplyr verbs (rename/slice/desc),
# which is why Fig5/00_setup.R keeps it unattached too.
# =============================================================================

# A full GENCODE GTF takes ~1-2 min to parse, so hold the parsed records for the
# session, keyed by path. Callers asking for several windows pay that once.
.gtf_cache <- new.env(parent = emptyenv())

.gtf_records <- function(gtf) {
  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("package 'rtracklayer' is required to read a GTF")
  key <- normalizePath(gtf, mustWork = TRUE)
  if (is.null(.gtf_cache[[key]])) {
    message("  parsing GTF (once per session): ", basename(key))
    .gtf_cache[[key]] <- rtracklayer::import(gtf, feature.type = c("gene", "exon"))
  }
  .gtf_cache[[key]]
}

.gr_to_dt <- function(g) data.table::data.table(
  gene   = g$gene_name,
  chrom  = as.character(GenomicRanges::seqnames(g)),
  start  = GenomicRanges::start(g),
  end    = GenomicRanges::end(g),
  strand = as.character(GenomicRanges::strand(g)))

# Extract gene models from `gtf`, selected EITHER by name (`genes`) or by window
# (`chrom` + `start`/`end`, in bp). Returns list(spans =, exons =), each a
# data.table sorted by chrom then start.
#
# `gene_types` filters on GENCODE's gene_type. The default keeps protein-coding
# genes only: a window drawn from a real GTF otherwise fills with lncRNAs,
# antisense transcripts and pseudogenes, which turns the gene track into an
# unreadable stack. Pass NULL to keep every biotype.
extract_gene_models <- function(gtf, genes = NULL, chrom = NULL,
                                start = NULL, end = NULL,
                                gene_types = "protein_coding") {
  if (is.null(genes) && is.null(chrom))
    stop("give either `genes` or a `chrom`/`start`/`end` window")

  gr <- .gtf_records(gtf)
  if (!is.null(gene_types) && !is.null(gr$gene_type))
    gr <- gr[gr$gene_type %in% gene_types | gr$gene_name %in% genes]
  gr <- gr[!is.na(gr$gene_name)]

  if (!is.null(genes)) {
    gr <- gr[gr$gene_name %in% genes]
    if (!length(gr)) stop("none of the requested genes were found in ", gtf)
    absent <- setdiff(genes, unique(gr$gene_name))
    if (length(absent))
      stop("not found in the GTF: ", paste(absent, collapse = ", "))
  } else {
    # Window selection is half-open on purpose: a gene counts as overlapping if
    # any part of it is inside, so genes running off either edge still render
    # (the track clamps them to the panel).
    keep <- as.character(GenomicRanges::seqnames(gr)) == chrom &
            GenomicRanges::start(gr) <= end & GenomicRanges::end(gr) >= start
    gr <- gr[keep]
  }

  spans <- .gr_to_dt(gr[gr$type == "gene"])

  # One symbol can carry several gene records (PAR_Y copies, scaffold
  # duplicates). Fig5's track draws one line, label and TSS arrow per span row,
  # so a duplicate silently renders the gene twice. Keep the longest record and
  # say so, rather than emitting a double.
  if (anyDuplicated(spans$gene)) {
    dup <- unique(spans$gene[duplicated(spans$gene)])
    warning("multiple gene records for ", paste(dup, collapse = ", "),
            "; keeping the longest of each", call. = FALSE)
    data.table::setorder(spans[, .len := end - start], gene, -.len)
    spans <- unique(spans, by = "gene")[, .len := NULL]
  }

  ex_gr <- gr[gr$type == "exon"]
  exons <- if (length(ex_gr)) {
    merged <- unlist(GenomicRanges::reduce(
      GenomicRanges::split(ex_gr, ex_gr$gene_name)))
    data.table::data.table(
      gene   = names(merged),
      chrom  = as.character(GenomicRanges::seqnames(merged)),
      start  = GenomicRanges::start(merged),
      end    = GenomicRanges::end(merged),
      strand = as.character(GenomicRanges::strand(merged)))
  } else {
    spans[0, ]
  }
  # Exons only for genes that produced a span, so the two tables stay consistent.
  exons <- exons[gene %in% spans$gene]

  data.table::setorder(spans, chrom, start)
  data.table::setorder(exons, chrom, start)
  list(spans = spans, exons = exons)
}
