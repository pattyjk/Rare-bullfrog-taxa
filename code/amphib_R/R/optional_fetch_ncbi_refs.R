## OPTIONAL. Re-download the NCBI 16S reference sequences used for exact-V4 matching -------
## The copies used in the analysis (downloaded 2026-10-16) are already in data/raw/ncbi_refs/.
## Re-running gives a newer GenBank snapshot, so counts may differ slightly. Requires internet access.
## Query per genus: "<Genus>"[Organism] AND (16S ribosomal RNA[Title] OR 16S rRNA[Title]) AND 1000:1700[SLEN]; first 1,500 records.
source("R/00_functions.R")
fetch_genus <- function(g, outdir = file.path(RAW, "ncbi_refs"), retmax = 1500) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  term <- utils::URLencode(sprintf('"%s"[Organism] AND (16S ribosomal RNA[Title] OR 16S rRNA[Title]) AND 1000:1700[SLEN]', g), reserved = TRUE)
  x <- paste(readLines(sprintf("%sesearch.fcgi?db=nuccore&term=%s&retmax=%d&usehistory=y", base, term, retmax), warn = FALSE), collapse = "")
  n <- min(as.integer(sub(".*?<Count>(\\d+)</Count>.*", "\\1", x)), retmax)
  we <- sub(".*<WebEnv>(.*?)</WebEnv>.*", "\\1", x); qk <- sub(".*<QueryKey>(\\d+)</QueryKey>.*", "\\1", x)
  fa <- unlist(lapply(seq(0, n - 1, by = 500), function(s) { Sys.sleep(0.4)
    readLines(sprintf("%sefetch.fcgi?db=nuccore&query_key=%s&WebEnv=%s&rettype=fasta&retmode=text&retstart=%d&retmax=500", base, qk, we, s), warn = FALSE) }))
  writeLines(fa, file.path(outdir, paste0(g, ".fasta"))); message(g, ": ", sum(grepl("^>", fa)), " sequences")
}
if (sys.nframe() == 0) for (g in REF_GENERA) fetch_genus(g)
