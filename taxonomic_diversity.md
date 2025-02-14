Taxonomic Diversity
================

``` r
library("picante")
library("vegan")
library("rbiom")
library("geiger")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("stringr")
library("readr")
library("tidyverse")
library("viridis")
library("scales")
library("ggnewscale")
library("colorspace")
library("gridExtra")

sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     colorspace_2.1-0  ggnewscale_0.4.8  scales_1.2.1     
    ##  [5] viridis_0.6.2     viridisLite_0.4.1 lubridate_1.9.2   forcats_1.0.0    
    ##  [9] purrr_1.0.1       tidyr_1.3.0       tibble_3.2.1      tidyverse_2.0.0  
    ## [13] readr_2.1.4       stringr_1.5.0     dplyr_1.1.2       ggpubr_0.6.0     
    ## [17] ggplot2_3.4.2     geiger_2.0.11     phytools_1.5-1    maps_3.4.1       
    ## [21] rbiom_1.0.3       picante_1.8.2     nlme_3.1-162      vegan_2.6-4      
    ## [25] lattice_0.21-8    permute_0.9-7     ape_5.7-1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] splines_4.2.2           foreach_1.5.2           carData_3.0-5          
    ##  [4] RcppParallel_5.1.7      expm_0.999-7            yaml_2.3.7             
    ##  [7] slam_0.1-50             numDeriv_2016.8-1.1     pillar_1.9.0           
    ## [10] backports_1.4.1         glue_1.6.2              quadprog_1.5-8         
    ## [13] phangorn_2.11.1         digest_0.6.31           ggsignif_0.6.4         
    ## [16] htmltools_0.5.5         Matrix_1.5-4            pkgconfig_2.0.3        
    ## [19] broom_1.0.4             mvtnorm_1.1-3           tzdb_0.3.0             
    ## [22] optimParallel_1.0-2     timechange_0.2.0        combinat_0.0-8         
    ## [25] mgcv_1.8-42             generics_0.1.3          car_3.1-2              
    ## [28] withr_2.5.0             cli_3.6.1               mnormt_2.1.1           
    ## [31] magrittr_2.0.3          evaluate_0.20           fansi_1.0.4            
    ## [34] doParallel_1.0.17       MASS_7.3-58.3           rstatix_0.7.2          
    ## [37] tools_4.2.2             hms_1.1.3               lifecycle_1.0.3        
    ## [40] munsell_0.5.0           cluster_2.1.4           plotrix_3.8-2          
    ## [43] compiler_4.2.2          clusterGeneration_1.3.7 rlang_1.1.0            
    ## [46] grid_4.2.2              iterators_1.0.14        rstudioapi_0.14        
    ## [49] subplex_1.8             igraph_1.4.2            rmarkdown_2.25         
    ## [52] gtable_0.3.3            codetools_0.2-19        abind_1.4-5            
    ## [55] deSolve_1.35            R6_2.5.1                knitr_1.44             
    ## [58] fastmap_1.1.1           utf8_1.2.3              fastmatch_1.1-3        
    ## [61] stringi_1.7.12          parallel_4.2.2          Rcpp_1.0.10            
    ## [64] vctrs_0.6.2             scatterplot3d_0.3-43    tidyselect_1.2.0       
    ## [67] xfun_0.39               coda_0.19-4

``` r
comm_phylo <- read.tree("./input/taxIDs.tree") # phylogeny with taxonomy id as tip labels, get from the phyloplus portal
comm_file <- readsample("./input/community_file.txt") # summarized from Bracken reports
metadata <- read_csv("./input/metadata.csv", col_names = TRUE)
metadata$collection_site <- str_replace(metadata$collection_site, "-", "_")
lineage_info <- read_csv("./input/lineage_info.csv", col_names = TRUE) # extracted from NCBI taxonomy database
for (colname in c("TaxPhylum", "TaxClass", "TaxOrder", "TaxFamily", "TaxGenus", "TaxSpecies")) {
  lineage_info[[colname]] <- make.names(lineage_info[[colname]])
}

normalized_comm_file <- sweep(comm_file, 1, rowSums(comm_file), "/") * 100
metadata <- metadata %>% filter(classified_reads >= 2500000)

# calculate weighted faith's index
swenson_weighted_faith <- function(my_phylo, my_sample) {
  weighted_faith_function <- function(my_sub_sample) {
    tmp_tree <- ape::drop.tip(my_phylo, setdiff(my_phylo$tip.label, names(my_sub_sample)[which(my_sub_sample > 0)]))
    branches <- matrix(NA, nrow = nrow(tmp_tree$edge), ncol = 4)
    branches[, 1:2] <- tmp_tree$edge
    branches[, 3] <- tmp_tree$edge.length
    get_leaves <- function(x) {
      leaves_node <- geiger::tips(tmp_tree, x[2])
      return(leaves_node)
    }
    leaves <- apply(branches, MARGIN = 1, get_leaves)
    for (i in 1:length(leaves)) {
      branches[i, 4] <- mean(my_sub_sample[leaves[[i]]], na.rm = TRUE)
    }
    return(nrow(tmp_tree$edge) * ((sum(branches[, 3] * branches[, 4])) / sum(branches[, 4])))
  }
  outt <- apply(my_sample, MARGIN = 1, weighted_faith_function)
  return(outt)
}


weighted_faith_output <- swenson_weighted_faith(comm_phylo, normalized_comm_file)
outputs <- data.frame("weighted_faith" = weighted_faith_output)
outputs <- merge(outputs, metadata, by.x = "row.names", by.y = "sample_id")
outputs$region <- factor(outputs$region, levels = c(c("Mapocho river region, Chile", "Maipo river region, Chile", 
                                                      "State of Paraiba, Brazil", "State of Rio de Janeiro, Brazil")))

region_colors <- setNames(c("#8da0cb", "#66c2a5", "#fc8d62", "#e78ac3"), 
                          c("Mapocho river region, Chile", "Maipo river region, Chile", 
                            "State of Paraiba, Brazil", "State of Rio de Janeiro, Brazil"))

violinplot_region <- ggplot(outputs, aes(x = region, y = weighted_faith)) + 
  geom_violin(aes(fill = region), scale = "count", alpha = 0.3, adjust = 0.5) + 
  geom_boxplot(aes(color = region), width = 0.15, outlier.shape = NA, show.legend = FALSE) +
  geom_point(aes(shape = water_type, color = region), size = 1.5, position = position_jitter(width = 0.25)) +
  scale_fill_manual("Sampling\n  Region", values = region_colors, 
                    guide = guide_legend(order = 1, nrow = 2, byrow = TRUE),
                    labels = c("Mapocho river region, Chile" = "Mapocho river, Chile",
                               "Maipo river region, Chile" = "Maipo river, Chile",
                               "State of Paraiba, Brazil" = "Paraíba, Brazil",
                               "State of Rio de Janeiro, Brazil" = "Rio de Janeiro, Brazil")) + 
  scale_shape_manual(name = "Water\nTypes", guide = guide_legend(order = 2, nrow = 2, byrow = TRUE),
                     values = setNames(c(24, 25, 21, 22, 23), 
                                       c("Creek", "River", "Reservoir", "Pond", "Irrigation Canal"))) +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_x_discrete(labels = c("Mapocho river region, Chile" = "Mapocho", 
                              "Maipo river region, Chile" = "Maipo",
                              "State of Paraiba, Brazil" = "Paraíba", 
                              "State of Rio de Janeiro, Brazil" = "Rio de Janeiro")) +
  labs(x = "Sampling Region", y = "Weighted Faith's Index") +
  theme_pubr() +
  theme(legend.direction = "horizontal", legend.box = "vertical")

# Combine abundance data by phylum
comm_file_phylum <- normalized_comm_file
for (i_taxo in unique(lineage_info[["TaxPhylum"]])) {
  ind_list <- which(lineage_info[["TaxPhylum"]] == i_taxo)
  taxid_list <- as.character(lineage_info$CurrentTaxID[ind_list])
  if (length(taxid_list) == 1) {
    comm_file_phylum[[i_taxo]] = comm_file_phylum[, taxid_list]
  } else {
    comm_file_phylum[[i_taxo]] = rowSums(comm_file_phylum[, taxid_list])
  }
}
comm_file_phylum <- comm_file_phylum[, !str_detect(names(comm_file_phylum), "^[0-9]")]
abundant_phylum <- names(comm_file_phylum)[colSums(comm_file_phylum) / nrow(comm_file_phylum) > 1]
abundant_phylum <- abundant_phylum[abundant_phylum != "NA."]
comm_file_phylum_modified <- comm_file_phylum[, abundant_phylum]
comm_file_phylum_modified$Others <- 100.00 - rowSums(comm_file_phylum_modified)
ordered_phylum <- c(sort(unique(colnames(comm_file_phylum_modified)[colnames(comm_file_phylum_modified) != "Others"])), 
                    "Others")
comm_file_phylum_modified <- comm_file_phylum_modified %>% 
  rownames_to_column(var = "Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Phylum", values_to = "Abundance") %>% 
  mutate(Region = setNames(metadata$region, metadata$sample_id)[Sample],
         Sample = factor(Sample, levels = outputs$Row.names[order(outputs$weighted_faith, decreasing = TRUE)]),
         Phylum = factor(Phylum, levels = ordered_phylum))
composition_plot <- ggplot(comm_file_phylum_modified) + 
  geom_bar(aes(x = Sample, y = Abundance, fill = Phylum), stat = "identity", width = 1) +
  geom_segment(aes(x = Sample, xend = Sample, color = Region), y = 0, yend = -2, alpha = 0.5) +
  scale_fill_manual(values = setNames(c(viridis(5, begin = 0.4, end = 0.9, option = "B"), "#d9d9d9"), 
                                      ordered_phylum),
                    guide = guide_legend(nrow = 2)) +
  scale_color_manual(values = region_colors) +
  ylim(c(0, 100)) +
  labs(y = "Relative Abundance") +
  theme_pubr(legend = "top") +
  geom_hline(yintercept = 100, color = "white") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = 9),
        legend.title = element_text(size = 12), legend.key.size = unit(16, "points")) +
  guides(color = "none")

ggarrange(violinplot_region, composition_plot, labels = c("a", "b"), widths = c(4, 7))
```

![](taxonomic_diversity_files/figure-gfm/alpha_diversity-1.png)<!-- -->

``` r
# Slightly modify functions from the vagan package so that the UMAP ordinates can be used as input.

vectorfit_modified <- function(X, P, permutations = 0, strata = NULL, ...) {
  EPS <- sqrt(.Machine$double.eps)
  # directly make the fitting unweighted
  w <- rep(1, nrow(X))
  P <- as.matrix(P)
  if (nrow(P) != nrow(X)) {stop("input data have non-matching numbers of observations")}
  Xw <- .Call("do_wcentre", X, w)
  Pw <- .Call("do_wcentre", P, w)
  colnames(Pw) <- colnames(P)
  nc <- ncol(X)
  Q <- qr(Xw)
  H <- qr.fitted(Q, Pw)
  heads <- qr.coef(Q, Pw)
  r <- diag(cor(H, Pw) ^ 2)
  r[is.na(r)] <- 0
  heads <- decostand(heads, "norm", 2)
  heads <- t(heads)
  if (is.null(colnames(X))) {
    colnames(heads) <- paste("Dim", 1:nc, sep = "")
  } else {
    colnames(heads) <- colnames(X)
  }
  ## make permutation matrix for all variables handled in the next loop
  nr <- nrow(X)
  permat <- vegan:::getPermuteMatrix(permutations, nr, strata = strata)
  if (ncol(permat) != nr) {
    stop(gettextf("'permutations' have %d columns, but data have %d rows", ncol(permat), nr))
  }
  permutations <- nrow(permat)
  if (permutations) {
    ptest <- function(indx, ...) {
      take <- P[indx, , drop = FALSE]
      take <- .Call("do_wcentre", take, w)
      Hperm <- qr.fitted(Q, take)
      diag(cor(Hperm, take)) ^ 2
    }
    permstore <- sapply(1:permutations, function(indx, ...) ptest(permat[indx,], ...))
    ## Single variable is dropped to a vector, and otherwise
    ## permutations are the matrix columns and variables are rows
    if (!is.matrix(permstore)) {permstore <- matrix(permstore, ncol = permutations)}
    permstore <- sweep(permstore, 1, r - EPS, ">=")
    validn <- rowSums(is.finite(permstore))
    pvals <- (rowSums(permstore, na.rm = TRUE) + 1) / (validn + 1)
  } else {pvals <- NULL}
  sol <- list(arrows = heads, r = r, permutations = permutations, pvals = pvals)
  sol$control <- attr(permat, "control")
  class(sol) <- "vectorfit"
  sol
}

envfit_modified <- function(coord_matrix, env, permutations = 999, na.rm = FALSE, strata = NULL, ...) {
  vectors <- NULL
  seed <- NULL
  keep <- complete.cases(coord_matrix) & complete.cases(env)
  if (any(!keep)) {
    if (!na.rm) {stop("missing values in data: consider na.rm = TRUE")}
    coord_matrix <- coord_matrix[keep, , drop = FALSE]
    ## drop any lost levels, explicitly don't include NA as a level
    env <- droplevels(env[keep,, drop = FALSE], exclude = NA)
    na.action <- structure(seq_along(keep)[!keep], class = "omit")
  }
  ## make permutation matrix for all variables handled in the next loop
  nr <- nrow(coord_matrix)
  permat <-  vegan:::getPermuteMatrix(permutations, nr)
  if (ncol(permat) != nr) {
    stop(gettextf("'permutations' have %d columns, but data have %d rows", ncol(permat), nr))
  }
  vectors <- vectorfit_modified(coord_matrix, env, permutations, ...)
  sol <- list(vectors = vectors)
  sol
}

ordisurf_modified <- function(coord_matrix, y, choices = c(1, 2), knots = 10, family = "gaussian", col = "red", 
                              isotropic = TRUE, thinplate = TRUE, bs = "tp", fx = FALSE, add = FALSE, display = "sites", 
                              main, nlevels = 10, levels, npoints = 31, labcex = 0.6, 
                              bubble = FALSE, cex = 1, select = TRUE, method = "REML", gamma = 1, plot = TRUE,
                              lwd.cl = par("lwd"), ...) {
  weights.default <- function(object, ...) {NULL}
  if (!missing(thinplate)) {
    warning("use of 'thinplate' is deprecated and will soon be removed;\nuse 'isotropic' instead")
    isotropic <- thinplate
  }
  GRID <- npoints
  w <- rep(1, nrow(coord_matrix))
  yname <- deparse(substitute(y))
  kk <- complete.cases(coord_matrix) & !is.na(y)
  if (!all(kk)) {
    coord_matrix <- coord_matrix[kk, , drop = FALSE]
    y <- y[kk]
    w <- w[kk]
  }
  x1 <- coord_matrix[, 1]
  x2 <- coord_matrix[, 2]
  if (!(missfx <- missing(fx)) && missing(knots)) {
    warning("requested fixed d.f. splines but without specifying 'knots':\nswitching to 'fx = FALSE'")
  }
  if (length(fx) > 2L) {
    warning("length of 'fx' supplied exceeds '2': using the first two")
  }
  fx <- rep(fx, length.out = 2)
  if (!missfx) {
    if ((miss.select <- missing(select)) && any(fx)) {
      warning("'fx = TRUE' requested; using 'select = FALSE'")
      select <- FALSE
    } else if (!miss.select && isTRUE(select)) {
      stop("fixed d.f. splines ('fx = TRUE') incompatible with 'select = TRUE'")
    }
  }
  if (length(knots) > 2L)
    warning("length of 'knots' supplied exceeds '2': using the first two")
  knots <- rep(knots, length.out = 2)
  if (length(bs) > 2L)
    warning("number of basis types supplied exceeds '2': using the first two")
  bs <- rep(bs, length.out = 2)
  BS <- c("tp","ts","cr","cs","ds","ps","ad")
  want <- match(bs, BS)
  user.bs <- bs
  bs <- BS[want]
  wrong <- is.na(bs)
  if (any(wrong)) {
    stop(gettextf("supplied basis type of '%s' not supported", paste(unique(user.bs[wrong]), collapse = ", ")))
  }
  if (isTRUE(isotropic) && any(bs %in% c("cr", "cs", "ps"))) {
    stop("bases \"cr\", \"cs\", and \"ps\" not allowed in isotropic smooths")
  }
  if (knots[1] <= 0) {
    f <- formula(y ~ x1 + x2)
  } else if (knots[1] == 1) {
    f <- formula(y ~ poly(x1, 1) + poly(x2, 1))
  } else if (knots[1] == 2) {
    f <- formula(y ~ poly(x1, 2) + poly(x2, 2) + poly(x1, 1):poly(x2, 1))
  } else if (isotropic) {
    f <- formula(paste0("y ~ s(x1, x2, k = ", knots[1], ", bs = \"", bs[1], "\", fx = ", fx[1],")"))
  } else {
    if (any(bs %in% c("ad"))) {
      f <- formula(paste0("y ~ s(x1, k = ", knots[1], ", bs = \"", bs[1], "\", fx = ", fx[1], ") + s(x2, k = ",
                          knots[2], ", bs = \"", bs[2], "\", fx = ", fx[2], ")"))
    } else {
      f <- formula(paste0("y ~ te(x1, x2, k = c(", paste0(knots, collapse = ", "), "), bs = c(",
                          paste0("\"", bs, "\"", collapse = ", "), "), fx = c(",paste0(fx, collapse = ", "),"))"))
    }
  }
  mod <- mgcv::gam(f, family = family, weights = w, select = select, method = method, gamma = gamma)
  xn1 <- seq(min(x1), max(x1), len = GRID)
  xn2 <- seq(min(x2), max(x2), len = GRID)
  newd <- expand.grid(x1 = xn1, x2 = xn2)
  fit <- predict(mod, type = "response", newdata = as.data.frame(newd))
  poly <- chull(cbind(x1, x2))
  xhull1 <- x1[poly] + sign(x1[poly] - mean(x1[poly])) * diff(range(x1))/(GRID - 1)
  xhull2 <- x2[poly] + sign(x2[poly] - mean(x2[poly])) * diff(range(x2))/(GRID - 1)
  npol <- length(poly)
  np <- nrow(newd)
  inpoly <- numeric(np)
  inpoly <- .C("pnpoly", as.integer(npol), as.double(xhull1), as.double(xhull2), as.integer(np), as.double(newd[,1]),
               as.double(newd[,2]), inpoly = as.integer(inpoly))$inpoly
  is.na(fit) <- inpoly == 0
  if (plot) {
    if (!add) {
      if (bubble) {
        if (is.numeric(bubble))
          cex <- bubble
        cex <- (y -  min(y))/diff(range(y)) * (cex - 0.4) + 0.4
      }
      plot(coord_matrix, asp = 1, cex = cex, ...)
    }
    if (!missing(main) || (missing(main) && !add)) {
      if (missing(main)) {main <- yname}
      title(main = main)
    }
    if (missing(levels)) {levels <- pretty(range(fit, finite = TRUE), nlevels)}
    if (!select || (select && !isTRUE(all.equal(as.numeric(summary(mod)$edf), 0)))) {
      contour(xn1, xn2, matrix(fit, nrow = GRID), col = col, add = TRUE, levels = levels, labcex = labcex,
              drawlabels = !is.null(labcex) && labcex > 0, lwd = lwd.cl)
    }
  }
  mod$grid <- list(x = xn1, y = xn2, z = matrix(fit, nrow = GRID))
  class(mod) <- c("ordisurf", class(mod))
  mod
}
```

``` r
if (file.exists("./input/taxonomic_diversity/weighted_unifrac.csv")) {
  dist_wunifrac <- as.data.frame(read_csv("./input/taxonomic_diversity/weighted_unifrac.csv"))
  rownames(dist_wunifrac) <- dist_wunifrac[[1]]
  dist_wunifrac <- dist_wunifrac[, -1]
} else {
  dist_wunifrac <- unifrac(t(as.matrix(normalized_comm_file)), weighted = TRUE, tree = comm_phylo)
  dist_wunifrac = as.data.frame(as.matrix(dist_wunifrac))
  write.csv(dist_wunifrac, "./input/taxonomic_diversity/weighted_unifrac.csv")
}

if (!file.exists("./input/taxonomic_diversity/umap_output.csv")) {
  system("python ./input/umap_call.py --dist ./input/taxonomic_diversity/weighted_unifrac.csv --seed 2023 --output ./input/taxonomic_diversity/umap_output.csv")
}

umap_coords <- read.csv("./input/taxonomic_diversity/umap_output.csv")
rownames(umap_coords) <- umap_coords$Sample
umap_coords <- as.matrix(umap_coords[, c(2, 3)])

# The list `community_df_list` includes 5 data frames, each summarizing the abundance data at the corresponding taxonomic level.
taxonomy_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
community_df_list <- list()
for (taxo in taxonomy_levels) {
  comm_file_tmp <- normalized_comm_file
  for (i_taxo in unique(lineage_info[[paste0("Tax", taxo)]])) {
    if (i_taxo == "NA.") {next}
    ind_list <- which(lineage_info[[paste0("Tax", taxo)]] == i_taxo)
    taxid_list <- as.character(lineage_info$CurrentTaxID[ind_list])
    if (length(taxid_list) == 1) {
      comm_file_tmp[[i_taxo]] = comm_file_tmp[, taxid_list]
    } else {
      comm_file_tmp[[i_taxo]] = rowSums(comm_file_tmp[, taxid_list])
    }
  }
  comm_file_tmp <- comm_file_tmp[, !str_detect(names(comm_file_tmp), "^[0-9]")]
  community_df_list[[taxo]] <- comm_file_tmp
}

# The list `envfit_df_list` is based on `community_df_list`, and it contains envfit result table describing
# UMAP coordinates and R squared values for each of the taxon identified at different taxonomic ranks.

if (!file.exists("./input/taxonomic_diversity/envfit_df_list.rds")) {
  envfit_df_list <- list()
  for (taxo in taxonomy_levels) {
    cat("Processing taxonomic level", taxo, end = "\r")
    temp_envfit <- envfit_modified(umap_coords, community_df_list[[taxo]], permutations = 999)
    temp_df <- as.data.frame(temp_envfit$vectors$arrows) %>% 
      merge(data.frame("R2" = temp_envfit$vectors$r, "pvals" = temp_envfit$vectors$pvals), by = "row.names") %>%
      column_to_rownames("Row.names")
    envfit_df_list[[taxo]] <- temp_df
  }
  saveRDS(envfit_df_list, file = "./input/taxonomic_diversity/envfit_df_list.rds")
}
envfit_df_list <- readRDS(file = "./input/taxonomic_diversity/envfit_df_list.rds")

umap_plot <- as.data.frame(umap_coords) %>%
  merge(metadata, by.x = "row.names", by.y = "sample_id") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = region, fill = region, shape = water_type)) +
  geom_point(alpha = 0.7) +
  scale_fill_manual(name = "Sampling Region", values = region_colors, guide = "none") +
  scale_color_manual(name = "Samplign Region", values = region_colors, guide = "none") +
  scale_shape_manual(name = "Water Types", guide = "none",
                     values = setNames(c(24, 25, 21, 22, 23), 
                              c("Creek", "River", "Reservoir", "Pond", "Irrigation Canal"))) +
  labs(x = "UMAP1", y = "UMAP2")

umap_x_range <- ggplot_build(umap_plot)$layout$panel_scales_x[[1]]$range$range
umap_y_range <- ggplot_build(umap_plot)$layout$panel_scales_y[[1]]$range$range

# Plot ordisurf results
contour_stats <- function(ordisur_object) {
  temp_out <- expand.grid(x = ordisur_object$grid$x, y = ordisur_object$grid$y)
  out <- cbind(temp_out, c(ordisur_object$grid$z))
  names(out) <- c("x", "y", "z")
  return(out)
}

# The list `ordisurf_plot_list` contains 5 different lists, each of the sub-list contains ggplot objects for key
# taxa identified at the corresponding taxonomic level.
if (!file.exists("./input/taxonomic_diversity/ordisurf_plot_list.rds")) {
  ordisurf_plot_list <- list()
  message("Not all taxonomies are processed. Only include taxa where a significant linear model can explain at least 30% variance.")
  for (taxo in taxonomy_levels) {
    taxonomic_names <- rownames(envfit_df_list[[taxo]] %>% filter(R2 >= 0.30 & pvals <= 0.05))
    cat(paste("Processing:", taxo), end = "\n")
    cat(paste("Total taxa to process:", length(taxonomic_names)), end = "\n")
    ordisurf_results <- list()
    taxa_colors <- setNames(viridis(length(taxonomic_names)), taxonomic_names)
    for (taxonomic_name in taxonomic_names) {
      cat(paste(taxonomic_name, paste(rep(" ", 10), collapse = "")), end = "\r")
      ordisurf_object <- ordisurf_modified(umap_coords, community_df_list[[taxo]][[taxonomic_name]], 
                                           plot = FALSE, knots = 400)
      dev_expl <- percent(summary(ordisurf_object)$dev.expl, accuracy = 0.01)
      if (taxo == "Genus") {
        title_name <- bquote(italic(.(taxonomic_name)) * ": " ~ .(dev_expl))
      } else {
        title_name <- bquote(.(taxonomic_name) * ": " ~ .(dev_expl))
      }
      contour_values <- contour_stats(ordisurf_object)
      temp_plot <- umap_plot +
        new_scale_color() +
        geom_contour(data = contour_values, mapping = aes(x = x, y = y, z = z, color = after_stat(level)),
                     linewidth = 0.75, inherit.aes = FALSE) +
        scale_color_binned(name = "Relative\nAbundance", low = adjust_transparency(taxa_colors[[taxonomic_name]], 0.05), 
                           high = taxa_colors[[taxonomic_name]]) +
        new_scale_color() +
        geom_segment(data = envfit_df_list[[taxo]] %>% rownames_to_column("taxon") %>% filter(taxon == taxonomic_name), 
                     aes(x = mean(umap_x_range), xend = mean(umap_x_range) + UMAP_1 * 3, 
                         y = mean(umap_y_range), yend = mean(umap_y_range) + UMAP_2 * 3), 
                     color = taxa_colors[[taxonomic_name]], 
                     arrow = arrow(length = unit(0.3, "cm"), type = "closed"), linewidth = 1, inherit.aes = FALSE) +
        theme_pubr() +
        theme(legend.position = c(0.01, 1), legend.justification = c("left", "top"), legend.direction = "vertical",
              legend.background = element_rect(fill = "transparent", color = "transparent"), 
              legend.title = element_text(size = 10), legend.text = element_text(size = 9), 
              legend.key.size = unit(12, "points")) +
        guides(color = guide_legend(title.position = "top", title.hjust = 0.5)) +
        ggtitle(title_name)
      if (taxo %in% c("Genus")) {
        temp_plot <- temp_plot + theme(legend.position = "none")
      }
      ordisurf_results[[taxonomic_name]] <- temp_plot
    }
    ordisurf_plot_list[[taxo]] <- ordisurf_results
  }
  saveRDS(ordisurf_plot_list, file = "./input/taxonomic_diversity/ordisurf_plot_list.rds")
}
ordisurf_plot_list <- readRDS(file = "./input/taxonomic_diversity/ordisurf_plot_list.rds")

plot_modification <- function(x){
  x + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
}

p1 <- as.data.frame(umap_coords) %>%
  merge(metadata, by.x = "row.names", by.y = "sample_id") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, fill = region, shape = water_type), color = "black") +
  geom_point(size = 4.5, alpha = 0.7) +
  scale_fill_manual(name = "Sampling Region", values = region_colors, 
                    labels = c("State of Paraiba, Brazil" = "State of Paraíba, Brazil")) +
  scale_shape_manual(name = "Water Types", 
                     values = setNames(c(24, 25, 21, 22, 23), 
                                       c("Creek", "River", "Reservoir", "Pond", "Irrigation Canal"))) +
  labs(x = "UMAP1", y = "UMAP2") +
  guides(fill = guide_legend(title.position = "top", nrow = 2, byrow = TRUE, order = 1, override.aes = list(shape = 21)), 
         shape = guide_legend(title.position = "top", nrow = 2, byrow = TRUE, order = 2),
         color = "none") + 
  theme_pubr(legend = "top", base_size = 21)

p2 <- do.call(grid.arrange, c(lapply(ordisurf_plot_list$Phylum, plot_modification), nrow = 4))

grid.arrange(p1, ggplot() + theme_void(), p2, widths = c(3, 0.1, 1))
```

![](taxonomic_diversity_files/figure-gfm/beta_diversity-1.png)<!-- -->
