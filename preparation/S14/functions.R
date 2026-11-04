# Data preparation from the original eQTL discovery comparison scripts.
discovery_read <- function(path, columns) {
  read_path <- path
  if (identical(readBin(path, "raw", n = 2), as.raw(c(31, 139))) && !grepl("[.]gz$", path)) {
    read_path <- tempfile(fileext = ".csv.gz")
    stopifnot(file.symlink(path, read_path))
    on.exit(unlink(read_path))
  }
  data.table::fread(read_path, select = columns, nThread = 2) %>% as_tibble()
}

discovery_sets <- function(data, id_column, study, analysis, genes) {
  data %>% transmute(molecular_trait_id = as.character(.data[[id_column]]),
                      cs_id = as.character(cs_id)) %>%
    distinct() %>%
    mutate(gene_id = stringr::str_remove(molecular_trait_id, "\\.[0-9]+$")) %>%
    left_join(genes %>% select(gene_id, gene_name, gene_type), by = "gene_id", relationship = "many-to-one") %>%
    mutate(study = study, analysis = analysis,
           protein_coding = !is.na(gene_type) & gene_type == "protein_coding")
}

discovery_counts <- function(sets) {
  sets %>% group_by(study, analysis) %>%
    summarise(credible_sets = n(), genes_with_credible_sets = n_distinct(gene_id), .groups = "drop")
}

shared_discovery_genes <- function(sets) {
  sets %>% filter(protein_coding,
    (study == "AoU" & analysis == "COMB") | study %in% c("INTERVAL", "TOPMed")) %>%
    distinct(study, gene_id) %>% count(gene_id, name = "n_studies") %>%
    filter(n_studies == 3L) %>% arrange(gene_id) %>% pull(gene_id)
}

# The original minimum enclosing ellipse is a display envelope, not a CI.
discovery_envelope <- function(data) {
  coords <- cbind(data$sample_size, data$credible_sets)
  low <- apply(coords, 2, min)
  span <- apply(coords, 2, function(x) diff(range(x)))
  scaled <- sweep(sweep(coords, 2, low, "-"), 2, span, "/")
  P <- t(scaled)
  Q <- rbind(P, rep(1, nrow(scaled)))
  u <- rep(1 / nrow(scaled), nrow(scaled))
  for (iteration in seq_len(15000)) {
    X <- sweep(Q, 2, u, "*") %*% t(Q)
    distances <- colSums(Q * (solve(X) %*% Q))
    j <- which.max(distances)
    if (distances[j] <= 3.001) break
    step <- (distances[j] - 3) / (3 * (distances[j] - 1))
    u <- (1 - step) * u
    u[j] <- u[j] + step
  }
  center <- as.numeric(P %*% u)
  A <- solve(sweep(P, 2, u, "*") %*% t(P) - tcrossprod(center)) / 2
  e <- eigen(A)
  delta <- sweep(scaled, 2, center, "-")
  distance2 <- rowSums((delta %*% A) * delta)
  stretch <- sqrt(max(distance2)) * 1.02
  angles <- seq(0, 2 * pi, length.out = 721)
  boundary <- sweep(stretch * cbind(cos(angles), sin(angles)) %*%
                      diag(1 / sqrt(e$values)) %*% t(e$vectors), 2, center, "+")
  boundary <- sweep(sweep(boundary, 2, span, "*"), 2, low, "+")
  list(boundary = tibble(sample_size = boundary[, 1], credible_sets = boundary[, 2]),
       membership = tibble(tissue = data$tissue, normalized_distance2 = distance2 / stretch^2),
       details = list(method = "Minimum enclosing ellipse in display coordinates; 2% padding",
                      tissue_count = nrow(data), all_tissues_enclosed = all(distance2 / stretch^2 < 1),
                      confidence_region = FALSE, iterations = iteration))
}
