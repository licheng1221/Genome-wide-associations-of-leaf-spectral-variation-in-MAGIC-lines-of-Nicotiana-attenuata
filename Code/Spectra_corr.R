## =============================================================================
## Integrated Script: Spectra Adjustment & GWAS Visualization (ESA Compliant)
## Figure 2: Variance Partitioning
## Figure 3: Spectra as phenotypes
## Author: Cheng Li, SG@UZH
## Date: Feb 2026
## =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(tidyverse)
  library(scales)
  library(spectrolab)
  library(RColorBrewer)
  library(cowplot)
})

# --------------------------
# 1. Global User Paths (EDIT THESE)
# --------------------------
# Input data
spec_rds_raw <- "data/Pheno_leafspec.rds"
cyto_tsv      <- "data/cyto.tsv"
az_rds        <- "data/R_AZ.rds"
hsc_res_rds   <- "data/Cal_HSC/results.rds"   # HSC results
hsc_seg_rds   <- "data/Cal_HSC/segments.rds"  # HSC segments

# Output directories
out_dir <- "output"
fig_dir <- "figures"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# --------------------------
# 2. Global Plotting Settings
# --------------------------
par(
  cex.main = 0.65, font.main = 1, cex.lab = 0.65, font.lab = 1,
  cex.axis = 0.50, font.axis = 1, mar = c(3.5, 3.5, 2, 1) + 0.1, 
  mgp = c(2.0, 0.6, 0), bg = "transparent"
)

my_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    text = element_text(size = 7),
    legend.background = element_rect(fill = "transparent")
  )

RGB_alpha <- function(hex_color){
  rgb_val <- col2rgb(hex_color) / 255
  rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha = 0.4)
}

# --------------------------
# 3. Data Processing & Figure 2 (Spectra Adjustment)
# --------------------------
Spec_raw <- readRDS(spec_rds_raw)
Cyto     <- fread(cyto_tsv, data.table = FALSE)
AZ       <- readRDS(az_rds)

# Clean Identifiers
AZ$Taxa <- gsub("-", "_", gsub("M", "L", AZ$Genotype_ID))
split_Line <- strsplit(AZ$Genotype_ID, "-")
AZ$genotype_group  <- sapply(split_Line, `[`, 1)
AZ$genotype_number <- sapply(split_Line, `[`, 2)

az <- AZ %>% 
  filter(genotype_group %in% c("M1", "M2")) %>%
  select(Plant_ID, Leaf_Num, Batch, Day, Ctime, Taxa, genotype_group, genotype_number) %>%
  distinct(Taxa, .keep_all = TRUE)

UTonly <- AZ %>% filter(genotype_group == "UT")

pheno_merged <- Cyto %>%
  inner_join(Spec_raw, by = c("Line" = "Taxa")) %>%
  inner_join(az, by = c("Line" = "Taxa"))

# Indices Calculation
calc_indices <- function(df, start_col_offset = 0) {
  get_w <- function(w) df[, w - 349 + start_col_offset]
  df$NDWI      <- (get_w(865) - get_w(1614)) / (get_w(865) + get_w(1614))
  df$CIre      <- get_w(783) / get_w(704) - 1
  df$CCI       <- (get_w(560) - get_w(664)) / (get_w(560) + get_w(664))
  df$ARDSI_Cab <- (get_w(750) - get_w(730)) / (get_w(770) + get_w(720))
  df$ARDSI_Cw  <- (get_w(1360) - get_w(1080)) / (get_w(1560) + get_w(1240))
  df$ARDSI_Cm  <- (get_w(2200) - get_w(1640)) / (get_w(2240) + get_w(1720))
  return(df)
}

pheno_merged <- calc_indices(pheno_merged, start_col_offset = 2) %>% relocate(NDWI:ARDSI_Cm, .before = Plant_ID)
UTonly       <- calc_indices(UTonly, start_col_offset = 0) %>% relocate(NDWI:ARDSI_Cm, .before = Genotype_ID)

# Factor Preparation
pheno_merged <- pheno_merged %>%
  mutate(BA = factor(Batch), DT = factor(Day), CT = factor(Ctime),
         REP = factor(sapply(strsplit(Line, "_"), `[`, 1)),
         PL = factor(sapply(strsplit(Line, "_"), `[`, 2)),
         MO = factor(Cytoplasm), LN = factor(Leaf_Num))

UTonly$BA <- factor(UTonly$Batch)

# Adjustment Loop
n_wvl <- 2157
n_pl  <- nrow(pheno_merged)
newy  <- matrix(0, nrow = n_pl, ncol = n_wvl)
terms_vec <- c("LN", "REP", "ba", "LN:REP", "LN:ba", "DT:CT", "MO", "Residuals")
ss_list <- vector("list", 2157)
wl_vec  <- numeric(2157)

for(i in 3:2159) {
  pheno_merged$y <- pheno_merged[, i]
  BAmeans <- tapply(UTonly[, i-2], UTonly$BA, mean, na.rm = TRUE)
  utba  <- c(mean(BAmeans, na.rm = TRUE), BAmeans[1:34])
  pheno_merged$ba <- utba[as.numeric(pheno_merged$BA)]
  
  lmy2 <- lm(y ~ LN + REP + ba + LN:(REP + ba) + DT:CT + MO, data = pheno_merged)
  newy[, i-2] <- residuals(lmy2) + mean(pheno_merged$y)
  
  a  <- anova(lmy2)
  ss <- setNames(rep(NA_real_, length(terms_vec)), terms_vec)
  rn <- rownames(a)
  present <- intersect(terms_vec, rn)
  ss[present] <- a[present, "Sum Sq"]
  if ("Residuals" %in% rn) ss["Residuals"] <- a["Residuals", "Sum Sq"]
  ss_list[[i-2]] <- ss / sum(ss, na.rm = TRUE)
  wl_vec[i-2]    <- i - 3 + 350
}

# Export Corrected Data (Inputs for Figure 3)
newy2 <- data.frame(Taxa = pheno_merged$Line, REP = pheno_merged$REP, PL = pheno_merged$PL, newy)
colnames(newy2)[4:2160] <- c(seq(350, 2500), "NDWI", "CIre", "CCI", "ARDSI_Cab", "ARDSI_Cw", "ARDSI_Cm")

pheno_corr_path   <- file.path(out_dir, "Pheno_corr.rds")
indices_corr_path <- file.path(out_dir, "Indices_corr.rds")

saveRDS(newy2[, c(1, 4:2154)], pheno_corr_path)
saveRDS(newy2[, c(1, 2155:2160)], indices_corr_path)

# Figure 2 Plotting
var_df <- do.call(rbind, lapply(seq_along(ss_list), function(k) {
  data.frame(wavelength = wl_vec[k], term = names(ss_list[[k]]), prop = as.numeric(ss_list[[k]]))
}))

group_df <- var_df %>%
  mutate(group = case_when(
    term == "LN" ~ "Technique",
    term %in% c("REP", "ba", "LN:REP", "LN:ba") ~ "Space + tech",
    term == "DT:CT" ~ "Time",
    term == "MO" ~ "Maternal cytoplasm",
    term == "Residuals" ~ "Residuals"
  )) %>%
  group_by(wavelength, group) %>%
  summarise(prop = sum(prop, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = factor(group, levels = c("Space + tech", "Technique", "Time", "Maternal cytoplasm", "Residuals")))

Fig_2 <- ggplot(group_df, aes(x = wavelength, y = prop, fill = group)) +
  geom_area() +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(500, 2500, by = 500), expand = c(0,0)) +
  scale_fill_manual(values = c("Residuals" = "grey70", "Space + tech" = "#009E73", 
                               "Technique" = "#CC79A7", "Time" = "#56B4E9", "Maternal cytoplasm" = "#D55E00")) +
  labs(x = "Wavelength (nm)", y = "Variance explained (% of total)", fill = "Effect group") +
  my_theme

ggsave(file.path(fig_dir, "Figure_2.pdf"), Fig_2, width = 18, height = 7, units = "cm", device = cairo_pdf, dpi = 600)

# --------------------------
# 4. Figure 3: HSC-PA & Spectral Analysis
# --------------------------
# Use outputs from previous steps
Spec_corr  <- readRDS(pheno_corr_path)
Ind_corr   <- readRDS(indices_corr_path)
hsc_res    <- readRDS(hsc_res_rds)
hsc_seg    <- readRDS(hsc_seg_rds)

# Panel (a): Quantiles
Spec_no_taxa  <- Spec_raw[, -c(1:52)]
Pheno_no_taxa <- Spec_corr[, -c(1:52)]
CV_Spec <- t(as.data.frame(apply(Pheno_no_taxa, 2, sd) / apply(Pheno_no_taxa, 2, mean)))

plot.new() 
dev.control("enable")
par(cex.main = 0.65, cex.lab = 0.65, cex.axis = 0.5, mgp = c(1.8, 0.5, 0))
plot_quantile(as_spectra(Spec_no_taxa), col = RGB_alpha("#56B4E9"), main = "All RILs in this study",
              total_prob = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance",
              ylim = c(0, 0.6), xaxs = "i")
plot_quantile(as_spectra(Pheno_no_taxa), col = RGB_alpha("#D55E00"), total_prob = 1, border = FALSE, add = TRUE)
plot(as_spectra(CV_Spec), col = "#D55E00", lty = 3, lwd = 1, add = TRUE)
abline(v = c(1000, 1800), lty = 3)
plot_regions(as_spectra(Spec_no_taxa), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
Fig_3a <- recordPlot()

# Panel (b): Indices
df_indices <- Ind_corr %>%
  pivot_longer(-Taxa, names_to = "Index", values_to = "Value") %>%
  filter(Index %in% c("NDWI","CIre","CCI","ARDSI_Cab","ARDSI_Cw","ARDSI_Cm")) %>%
  mutate(Index = factor(Index, levels = c("NDWI","CIre","CCI","ARDSI_Cab","ARDSI_Cw","ARDSI_Cm")))

p_list <- lapply(levels(df_indices$Index), function(idx) {
  ggplot(filter(df_indices, Index == idx), aes(x = "", y = Value)) +
    geom_violin(fill = "grey95", color = NA, trim = FALSE) +
    geom_boxplot(width = 0.2, outlier.size = 0.3) +
    labs(title = idx, x = NULL, y = NULL) +
    my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
})
Fig_3b <- plot_grid(plotlist = p_list, ncol = 3, align = "hv")

# Panel (c): Heatmap
segments_long <- hsc_seg %>%
  spread(key = level, value = segment) %>%
  gather(key = "level", value = "segment", -wavelength) %>%
  mutate(wavelength = as.numeric(wavelength), segment = as.factor(segment), level = as.factor(level))

named_colors <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(unique(segments_long$segment))), 
                         sort(unique(segments_long$segment)))

Fig_3c <- ggplot(segments_long, aes(x = wavelength, y = level, fill = segment)) +
  geom_tile(color = NA) + 
  scale_fill_manual(values = named_colors, na.value = "transparent", name = "Segment") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Results of HSC-PA", x = "Wavelength (nm)", y = "Level") +
  my_theme + theme(legend.position = "right") + guides(fill = guide_legend(ncol = 4))

# Panel (d): Segment Selection
segments_to_read <- filter(hsc_res, N_PCs == 1)$segment
example_spectra  <- data.frame(wavelength = as.numeric(colnames(Pheno_no_taxa)), 
                               reflectance = colMeans(Pheno_no_taxa))
final_used <- hsc_seg %>% filter(segment %in% segments_to_read) %>%
  mutate(wavelength = as.numeric(wavelength)) %>% inner_join(example_spectra, by = "wavelength")

saveRDS(final_used, file.path(out_dir, "final_used.RDS"))

Fig_3d <- ggplot(final_used, aes(x = wavelength, y = reflectance, color = as.factor(segment))) +
  geom_point(size = 2) + 
  scale_color_manual(values = named_colors, name = "Segment") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Segments with 1 retained PC", x = "Wavelength (nm)", y = "Reflectance") +
  my_theme + theme(legend.position = "right", legend.key.size = unit(0.4, "cm")) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2)))

# Assemble Figure 3
top_row <- plot_grid(Fig_3a, Fig_3b, labels = c("a", "b"), label_size = 10, rel_widths = c(1, 1))
final_plot_3 <- plot_grid(top_row, Fig_3c, Fig_3d, ncol = 1, labels = c("", "c", "d"), 
                          label_size = 10, rel_heights = c(1, 1, 0.8), align = "v", axis = "l")

ggsave(file.path(fig_dir, "Figure_3.pdf"), final_plot_3, width = 18, height = 20, 
       units = "cm", device = cairo_pdf, dpi = 600, bg = "white")