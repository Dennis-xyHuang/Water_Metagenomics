Parameter Validation
================

``` r
library("ggplot2")
library("dplyr")
library("readr")
library("ggpubr")
library("vegan")
library("gridExtra")
library("rbiom")
library("tidyr")

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
    ##  [1] tidyr_1.3.0    rbiom_1.0.3    gridExtra_2.3  vegan_2.6-4    lattice_0.21-8
    ##  [6] permute_0.9-7  ggpubr_0.6.0   readr_2.1.4    dplyr_1.1.2    ggplot2_3.4.2 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.0   xfun_0.39          slam_0.1-50        purrr_1.0.1       
    ##  [5] splines_4.2.2      carData_3.0-5      colorspace_2.1-0   vctrs_0.6.2       
    ##  [9] generics_0.1.3     htmltools_0.5.5    yaml_2.3.7         mgcv_1.8-42       
    ## [13] utf8_1.2.3         rlang_1.1.0        pillar_1.9.0       glue_1.6.2        
    ## [17] withr_2.5.0        lifecycle_1.0.3    munsell_0.5.0      ggsignif_0.6.4    
    ## [21] gtable_0.3.3       evaluate_0.20      knitr_1.44         tzdb_0.3.0        
    ## [25] fastmap_1.1.1      parallel_4.2.2     fansi_1.0.4        Rcpp_1.0.10       
    ## [29] broom_1.0.4        scales_1.2.1       backports_1.4.1    RcppParallel_5.1.7
    ## [33] abind_1.4-5        hms_1.1.3          digest_0.6.31      rstatix_0.7.2     
    ## [37] grid_4.2.2         cli_3.6.1          tools_4.2.2        magrittr_2.0.3    
    ## [41] tibble_3.2.1       cluster_2.1.4      car_3.1-2          pkgconfig_2.0.3   
    ## [45] MASS_7.3-58.3      Matrix_1.5-4       rmarkdown_2.25     rstudioapi_0.14   
    ## [49] R6_2.5.1           nlme_3.1-162       compiler_4.2.2

``` r
# Rarefaction
sample_names <- c("1C52", "6C51", "2021-Dec-MAI22", "2022-Mar-MAI26", "2022-May-MAI22", 
                  "MAP8-54H", "MAP10-90K", "MAP11-102L", "W218", "W224")
read_choices <- c("10000", "50000", "100000", "500000", "1000000", "1500000", "2000000", "2500000", "3000000", 
                  "3500000", "4000000")
random_seed_choices <- c("1996", "2002", "2008", "2014", "2023")

run_random_rarefaction <- FALSE
if (run_random_rarefaction) {
  for (sample in sample_names) {
    file_path_1 <- file.path(getwd(), "Rarefaction/Raw_data", paste0("trimmed.bacterial_", sample, "_R1.fastq.gz"))
    file_path_2 <- file.path(getwd(), "Rarefaction/Raw_data", paste0("trimmed.bacterial_", sample, "_R2.fastq.gz"))
    for (n_reads in read_choices) {
      for (seed_number in random_seed_choices) {
        out_file_path_1 <- file.path(getwd(), "Rarefaction/Rarefied_data", 
                                     paste0("rarefied_", sample, "_d", n_reads, "_s", seed_number, "_R1.fastq"))
        out_file_path_2 <- file.path(getwd(), "Rarefaction/Rarefied_data", 
                                     paste0("rarefied_", sample, "_d", n_reads, "_s", seed_number, "_R2.fastq"))
        report_path <- file.path(getwd(), "Rarefaction/Kraken_report", 
                                 paste0("rarefied_", sample, "_d", n_reads, "_s", seed_number, "_kraken_report.txt"))
        commands <- c(paste("seqtk sample", paste0("-s", seed_number), file_path_1, n_reads, ">", out_file_path_1),
                      paste("seqtk sample", paste0("-s", seed_number), file_path_2, n_reads, ">", out_file_path_2),
                      paste("pigz", out_file_path_1), paste("pigz", out_file_path_2),
                      paste("kraken2 --db /kraken/db/path --threads 16 --report", report_path, "--paired --output -",
                            "--use-names --gzip-compressed", paste0(out_file_path_1, ".gz"), 
                            paste0(out_file_path_2, ".gz")))
        for (cmd_str in commands) {
          cat(cmd_str, end = "\r")
          system(cmd_str)
        }
      }
    }
  }
}

if (file.exists("./input/parameter_validation/rarefaction_summary.csv")) {
  rarefaction_summary_df <- read_csv("./input/parameter_validation/rarefaction_summary.csv")
} else {
  rarefaction_summary_df <- data.frame("sample" = character(), "reads" = numeric(), "seed" = numeric(), 
                                       "n_spp" = numeric(), "n_genus" = numeric())
  for (sample in sample_names) {
    for (n_reads in read_choices) {
      for (seed_number in random_seed_choices) {
        report_path <- file.path(getwd(), "Rarefaction/Kraken_report", 
                                 paste0("rarefied_", sample, "_d", n_reads, "_s", seed_number, "_kraken_report.txt"))
        report_df <- read_tsv(report_path, col_names = FALSE)
        report_df_spp <- report_df %>% filter(X4 == "S")
        n_spp <- nrow(report_df_spp)
        report_df_genus <- report_df %>% filter(X4 == "G")
        n_genus <- nrow(report_df_genus)
        temp_list <- list("sample" = sample, "reads" = as.numeric(n_reads), "seed" = as.numeric(seed_number), 
                          "n_spp" = n_spp, "n_genus" = n_genus)
        rarefaction_summary_df <- rbind(rarefaction_summary_df, temp_list)
      }
    }
  }
  write_csv(rarefaction_summary_df, "./input/parameter_validation/rarefaction_summary.csv")
}

rarefaction_plot_df <- rarefaction_summary_df %>% 
  group_by(sample, reads) %>% 
  summarize(min_spp = min(n_spp), max_spp = max(n_spp), avg_spp = mean(n_spp),
            min_genus = min(n_genus), max_genus = max(n_genus), avg_genus = mean(n_genus))

p1 <- ggplot(rarefaction_plot_df) +
  geom_line(aes(x = reads, y = avg_spp, color = sample), linewidth = 0.5, show.legend = FALSE) +
  geom_ribbon(aes(x = reads, ymin = min_spp, ymax = max_spp, fill = sample), alpha = 0.3, show.legend = FALSE) +
  geom_vline(xintercept = 2500000, linetype = "dashed") +
  scale_x_continuous("Number of Rarefied Reads", 
                     breaks = c(10000, 1000000, 2000000, 3000000, 4000000),
                     labels = scales::comma) +
  scale_y_continuous("Number of Species", labels = scales::comma) +
  scale_color_brewer(type = "qual", palette = 3) +
  scale_fill_brewer(type = "qual", palette = 3) +
  guides(legend.position = "None") +
  ggtitle("a") +
  theme_pubclean() +
  theme(plot.title = element_text(size = 25), axis.text.x = element_text(size = 8))

p2 <- ggplot(rarefaction_plot_df) +
  geom_line(aes(x = reads, y = avg_genus, color = sample), linewidth = 0.5, show.legend = FALSE) +
  geom_ribbon(aes(x = reads, ymin = min_genus, ymax = max_genus, fill = sample), alpha = 0.3, show.legend = FALSE) +
  geom_vline(xintercept = 2500000, linetype = "dashed") +
  scale_x_continuous("Number of Rarefied Reads", 
                     breaks = c(10000, 1000000, 2000000, 3000000, 4000000),
                     labels = scales::comma) +
  scale_y_continuous("Number of Genera", labels = scales::comma) +
  scale_color_brewer(type = "qual", palette = 3) +
  scale_fill_brewer(type = "qual", palette = 3) +
  guides(legend.position = "None") +
  ggtitle("b") +
  theme_pubclean() +
  theme(plot.title = element_text(size = 25), axis.text.x = element_text(size = 8))

# Bracken Threshold

# Zymo standard microbial community
# Listeria monocytogenes - 12%;  taxID: 1639
# Pseudomonas aeruginosa - 12%; taxID: 287
# Bacillus subtilis - 12%; taxID: 1423
# Escherichia coli - 12%; taxID: 562
# Salmonella enterica - 12%; taxID: 28901
# Lactobacillus fermentum - 12%; taxID: 1613
# Enterococcus faecalis - 12%; taxID: 1351
# Staphylococcus aureus - 12%; taxID: 1280
# Saccharomyces cerevisiae - 2%
# Cryptococcus neoformans - 2%

control_names <- c("Contr-A", "Contr-G")
bracken_read_thresholds <- c(seq(0, 2000, 10))

run_bracken <- FALSE
if (run_bracken) {
  for (sample in c(control_names, sample_names)) {
    temp_kraken_report <- file.path(getwd(), "Bracken_Threshold/Kraken_Report", paste0(sample, "_kraken_report.txt"))
    for (bracken_read_i in bracken_read_thresholds) {
      temp_bracken_report <- file.path(getwd(), "Bracken_Threshold/Bracken_Report",
                                       paste0(sample, "_bracken_threshold_", bracken_read_i, "_report.txt"))
      cmd_str <- paste("bracken -d /bracken/db/path", "-i", temp_kraken_report, "-o", temp_bracken_report, 
                       "-r 150 -l S", "-t", bracken_read_i)
      cat(cmd_str, end = "\r")
      system(cmd_str)
    }
  }
}


ref_phylo <- read.tree("./input/parameter_validation/parameter_validation_taxIDs.tree")
lineage_df <- read_csv("./input/parameter_validation/parameter_validation_lineage_info.csv")
taxid_to_genus <- setNames(lineage_df$TaxGenus, lineage_df$CurrentTaxID)
zymo_spp_taxids = c("1639", "287", "1423", "562", "28901", "1613", "1351", "1280")
zymo_spp_names = setNames(c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", 
                            "Escherichia coli", "Salmonella enterica", "Limosilactobacillus fermentum", 
                            "Enterococcus faecalis", "Staphylococcus aureus"), 
                          zymo_spp_taxids)
zymo_colors <- setNames(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", 
                          "#fccde5", "#969696"),
                        c(zymo_spp_names, "Others"))
zymo_community_df <- as.data.frame(rbind(setNames(rep(2500000 / 8, 8), zymo_spp_taxids)))

if (file.exists("./input/parameter_validation/bracken_control_summary.csv")) {
  bracken_control_summary_df <- read_csv("./input/parameter_validation/bracken_control_summary.csv")
} else {
  bracken_control_summary_df <- data.frame(matrix(nrow = 0, ncol = 17))
  colnames(bracken_control_summary_df) <- c("sample", "type", "bracken_threshold", "total_spp", "total_genus",
                                            zymo_spp_taxids, "others", "bray-curtis", "unifrac", "FP_n10")
  for (sample in control_names) {
    temp_kraken_report <- file.path(getwd(), "Bracken_Threshold/Kraken_Report", paste0(sample, "_kraken_report.txt"))
    original_kraken_report <- read_tsv(temp_kraken_report , col_names = FALSE)
    filter_final_n10 = (original_kraken_report$X4 >= original_kraken_report$X5 * 10) & 
      (original_kraken_report$X6 == "S")
    minimizer_based_FP_n10 = original_kraken_report[filter_final_n10,] %>% pull(X7)
    for (bracken_read_i in bracken_read_thresholds) {
      temp_bracken_report <- file.path(getwd(), "Bracken_Threshold/Bracken_Report",
                                       paste0(sample, "_bracken_threshold_", bracken_read_i, "_report.txt"))
      bracken_report <- read_tsv(temp_bracken_report)
      bracken_report <- bracken_report %>% 
        mutate(new_est_reads_adj = new_est_reads / sum(bracken_report$new_est_reads) * 2500000)
      sample_df <- rbind(bracken_report$new_est_reads_adj)
      colnames(sample_df) <- bracken_report$taxonomy_id
      sample_df <- as.data.frame(sample_df)
      comparison_df <- bind_rows(zymo_community_df, sample_df)
      rownames(comparison_df) <- c("standard", "sample")
      comparison_df[is.na(comparison_df)] <- 0
      bray_curtis_distance <- as.matrix(vegdist(comparison_df))[1, 2]
      unifrac_distance <- as.matrix(unifrac(t(as.matrix(comparison_df)), weighted = TRUE, tree = ref_phylo))[1, 2]
      tmp_abund <- list()
      for (spp_id in zymo_spp_taxids) {
        filtered_df <- bracken_report %>% filter(taxonomy_id == spp_id)
        tmp_abund[[spp_id]] <- ifelse(nrow(filtered_df) == 0, 0, (filtered_df %>% pull(new_est_reads_adj))[[1]])
      }
      tmp_list = list("sample" = sample, "type" = "control", "bracken_threshold" = bracken_read_i, 
                      "total_spp" = nrow(bracken_report), 
                      "total_genus" = length(unique(taxid_to_genus[bracken_report$taxonomy_id])),
                      "1639" = tmp_abund$"1639", "287" = tmp_abund$"287", "1423" = tmp_abund$"1423",
                      "562" = tmp_abund$"562", "28901" = tmp_abund$"28901", "1613" = tmp_abund$"1613",
                      "1351" = tmp_abund$"1351", "1280" = tmp_abund$"1280", "others" = 2500000 - sum(unlist(tmp_abund)), 
                      "bray-curtis" = bray_curtis_distance, "unifrac" = unifrac_distance,
                      "FP_n10" = length(bracken_report %>% 
                                          filter(taxonomy_id %in% minimizer_based_FP_n10) %>% 
                                          pull(fraction_total_reads)))
      tmp_df <- data.frame(tmp_list, check.names = FALSE)
      bracken_control_summary_df = rbind(bracken_control_summary_df, tmp_df)
    }
  }
  write_csv(bracken_control_summary_df, "./input/parameter_validation/bracken_control_summary.csv")
}

if (file.exists("./input/parameter_validation/bracken_sample_summary.csv")) {
  bracken_sample_summary_df <- read_csv("./input/parameter_validation/bracken_sample_summary.csv")
} else {
  bracken_sample_summary_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(bracken_sample_summary_df) <- c("sample", "type", "bracken_threshold", "total_spp", "total_genus", 
                                           "FP_n10") 
  for (sample in sample_names) {
    temp_kraken_report <- file.path(getwd(), "Bracken_Threshold/Kraken_Report", paste0(sample, "_kraken_report.txt"))
    original_kraken_report <- read_tsv(temp_kraken_report , col_names = FALSE)
    filter_final_n10 = (original_kraken_report$X4 >= original_kraken_report$X5 * 10) & 
      (original_kraken_report$X6 == "S")
    minimizer_based_FP_n10 = original_kraken_report[filter_final_n10,] %>% pull(X7)
    for (bracken_read_i in bracken_read_thresholds) {
      temp_bracken_report <- file.path(getwd(), "Bracken_Threshold/Bracken_Report",
                                       paste0(sample, "_bracken_threshold_", bracken_read_i, "_report.txt"))
      bracken_report <- read_tsv(temp_bracken_report)
      tmp_list = list("sample" = sample, type = "sample", "bracken_threshold" = bracken_read_i, 
                      "total_spp" = nrow(bracken_report),
                      "total_genus" = length(unique(taxid_to_genus[bracken_report$taxonomy_id])),
                      "FP_n10" = length(bracken_report %>% 
                                          filter(taxonomy_id %in% minimizer_based_FP_n10) %>% 
                                          pull(fraction_total_reads)))
      tmp_df <- data.frame(tmp_list, check.names = FALSE)
      bracken_sample_summary_df = rbind(bracken_sample_summary_df, tmp_df)
    }
  }
  write_csv(bracken_sample_summary_df, "./input/parameter_validation/bracken_sample_summary.csv")
}


composition_df <- data.frame("sample" = character(), "bracken_threshold" = numeric(), 
                             "species" = character(), "abundance" = numeric())
for (i_row in 1:nrow(bracken_control_summary_df)) {
  tmp_df = data.frame("sample" = rep(bracken_control_summary_df[[i_row, "sample"]], 9),
                      "bracken_threshold" = rep(bracken_control_summary_df[[i_row, "bracken_threshold"]], 9),
                      "species" = c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", 
                                    "Escherichia coli", "Salmonella enterica", "Limosilactobacillus fermentum", 
                                    "Enterococcus faecalis", "Staphylococcus aureus", "Others"),
                      "abundance" = unlist(bracken_control_summary_df[i_row, c("1639", "287", "1423", "562", "28901", 
                                                                               "1613", "1351", "1280", "others")]))
  composition_df <- rbind(composition_df, tmp_df)
}

bracken_summary_df <- bind_rows(bracken_control_summary_df, bracken_sample_summary_df)

bracken_counts <- bracken_summary_df %>% 
  dplyr::select(c("sample", "type", "bracken_threshold", "total_spp", "total_genus")) %>% 
  rename(setNames(c("total_spp", "total_genus"), c("Species", "Genus"))) %>% 
  pivot_longer(cols = c("Species", "Genus"), names_to = "Rank", values_to = "Counts") %>%
  mutate(Rank = factor(Rank, levels = c("Species", "Genus")))
p3 <- ggplot(bracken_counts %>% filter(Rank == "Species")) +
  geom_line(aes(x = bracken_threshold + 1, y = Counts + 1, color = type, group = sample), linewidth = 0.75) +
  geom_vline(xintercept = 500, linetype = "dashed") +
  scale_x_continuous("Bracken Threshold", trans = "sqrt", breaks = c(10, 100, 500, 1000)) +
  scale_y_continuous("Number of Species", trans = "sqrt", breaks = c(8, 100, 2000, 8000)) +
  scale_color_manual("Sample", values = setNames(c("#d95f02", "#7570b3"), c("control", "sample")),
                     labels = c("Zymo Control", "Water Sample")) +
  ggtitle("c") +
  theme_pubr() +
  theme(legend.position = "top", legend.justification = "right", plot.title = element_text(size = 25))

p4 <- ggplot(composition_df) + 
  facet_wrap(~sample, labeller = labeller(sample = c("Contr-A" = "Control #1", "Contr-G" = "Control #2"))) +
  geom_area(aes(x = bracken_threshold, y = abundance / 2500000 * 100, fill = species), 
            color = "white", alpha = 0.8, stat = "identity", linewidth = 0.75) + 
  geom_line(data = bracken_control_summary_df, 
            aes(x = bracken_threshold, y = `bray-curtis` / 0.2 * 100, group = sample), 
            color = "#756bb1", linewidth = 1) +
  geom_line(data = bracken_control_summary_df %>% filter(bracken_threshold <= 2000), 
            aes(x = bracken_threshold, y = unifrac / 0.2 * 100, group = sample),
            color = "#de2d26", linewidth = 1) +
  geom_line(data = data.frame(x = c(Inf, Inf), y = c(Inf, Inf), Metric = c("Weighted Bray-Curtis", "Weighted UniFrac")),
            aes(x = x, y = y, color = Metric), linewidth = 1) +
  geom_vline(xintercept = 500, linetype = "dashed") +
  scale_x_continuous(name = "Bracken Threshold") + 
  scale_y_continuous(name = "Relative Abundance (%)", expand = expansion(mult = 0.02), 
                     sec.axis = sec_axis(name = "Distance to Theoretical Composition", trans = ~ . / 100 * 0.2)) +
  scale_fill_manual("Species", values = zymo_colors, 
                    breaks = c("Bacillus subtilis", "Escherichia coli", "Enterococcus faecalis", 
                               "Limosilactobacillus fermentum", "Listeria monocytogenes", "Pseudomonas aeruginosa", 
                               "Staphylococcus aureus", "Salmonella enterica", "Others"),
                    labels = c(expression(italic("B. subtilis")), expression(italic("E. coli")), 
                               expression(italic("E. faecalis")), expression(italic("L. fermentum")), 
                               expression(italic("L. monocytogenes")), expression(italic("P. aeruginosa")),
                               expression(italic("S. aureus")), expression(italic("S. enterica")), "Others")) + 
  scale_color_manual("Metrics", 
                     values = setNames(c("#756bb1", "#de2d26"), c("Weighted Bray-Curtis", "Weighted UniFrac"))) +
  ggtitle("d") +
  theme_pubr() +
  theme(legend.position = "right", legend.text.align = 0, plot.title = element_text(size = 25))

p5 <- ggplot(bracken_summary_df) +
  geom_line(aes(x = bracken_threshold, y = FP_n10, group = sample, color = type), linewidth = 0.5) +
  geom_line(aes(x = bracken_threshold, y = FP_n10, group = sample, color = type), linewidth = 0.5) +
  geom_vline(xintercept = 500, linetype = "dashed") +
  scale_x_continuous("Bracken Threshold", trans = "log1p", limits = c(0, 1000), breaks = c(1, 10, 100, 1000)) +
  scale_y_continuous("Number of False Positive Species") +
  scale_color_manual("Sample",values = setNames(c("#d95f02", "#7570b3"), c("control", "sample")),
                     labels = c("Zymo Control", "Water Sample")) +
  ggtitle("e") +
  theme_pubr() +
  theme(legend.position = "top", legend.justification = "right", plot.title = element_text(size = 25))

grid.arrange(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 2, 3), c(4, 4, 5)))
```

![](parameter_validation_files/figure-gfm/parameter_validation-1.png)<!-- -->
