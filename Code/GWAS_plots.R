## =============================================================================
## Integrated Script: GWAS Visualization (ESA Compliant) - Part 2
## Figure 4: Indices GWAS results
## Figure 5: Single Wavelength (SW) GWAS results
## Figure 6: HSC-PA GWAS results
## Author: Cheng Li, SG@UZH
## Date: Feb 2026
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(spectrolab)
  library(cowplot)
  library(magick)
})

# 1. Directory & Path Definitions ----------------------------------------------
# Input/Output Directories
data_dir <- "data"
out_dir  <- "output"
fig_dir  <- "figures"

# Specific Input File Paths (Data)
indices_csv_path <- file.path(data_dir, "Indices_Top50_Full_Details.csv")
sw_csv_path      <- file.path(data_dir, "SW_Top50_Full_Details.csv")
hsc_csv_path     <- file.path(data_dir, "HSC_Top50_Full_Details.csv")

ind_rds_path     <- file.path(data_dir, "Indices_corr.rds")
geno_rds_path    <- file.path(data_dir, "geno_new.RDS")
hsc_pc_rds_path  <- file.path(data_dir, "HSC_PCs_corr.rds")
final_used_path  <- file.path(out_dir, "final_used.RDS")

# GAPIT PDF Image Paths (Now reading from data directory)
pdf_manhattan_indices <- file.path(data_dir, "GAPIT.Association.Manhattans_Symphysic.pdf")
pdf_legend_indices    <- file.path(data_dir, "GAPIT.Association.Manhattans_Symphysic_Legend.pdf")
pdf_qq_indices          <- file.path(data_dir, "GAPIT.Association.QQs_Symphysic_Ind_cw.pdf")
pdf_manhattan_sw      <- file.path(data_dir, "GAPIT.Association.Manhattan_Geno.BLINK.952.pdf")
pdf_qq_sw             <- file.path(data_dir, "GAPIT.Association.QQ.BLINK.952.pdf")
pdf_qq_seg37          <- file.path(data_dir, "GAPIT.Association.QQ.BLINK.Seg37.pdf")
pdf_qq_seg38          <- file.path(data_dir, "GAPIT.Association.QQ.BLINK.Seg38.pdf")

# 2. Global Settings & Functions -----------------------------------------------
par(
  cex.main = 0.65, cex.lab = 0.65, cex.axis = 0.50, font.main = 1,
  mar = c(3.5, 3.5, 2, 1) + 0.1, mgp = c(2.0, 0.6, 0), bg = "transparent"
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

make_pdf_panel <- function(pdf_path, page = 1, dpi = 300, trim = TRUE) {
  img <- magick::image_read_pdf(pdf_path, density = dpi)  
  stopifnot(length(img) >= page)
  im  <- img[page]
  if (trim) im <- magick::image_trim(im)                  
  cowplot::ggdraw() + cowplot::draw_image(im)
}

# 3. Figure 4: Indices-GWAS ----------------------------------------------------
# 3.1 Data Loading
Indices <- read_delim(indices_csv_path, skip = 1, col_names = FALSE, delim = ",")
colnames(Indices) <- c("Indices","Model","SNP","Chr","Pos","p_value","MAF","nobs",
                       "FDR_Adjusted_p_value","effect","threshold")

Indices <- Indices %>%
  mutate(Indices = str_replace_all(Indices, "\\(NYC\\)", "") %>% str_squish())

# 3.2 Heatmap (4b) & Boxplots (4d, 4e)
Fig_4b <- ggplot(filter(Indices, Model == "BLINK"), 
                 aes(as.factor(Indices), SNP, fill = as.numeric(FDR_Adjusted_p_value))) + 
  geom_tile(width = 0.4) +
  scale_fill_gradient(low = "#67000d", high = "#fee0d2") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Chromosome Position", x = "Indices", fill = "FDR_p-value")

Ind  <- readRDS(ind_rds_path)
geno <- readRDS(geno_rds_path)
Ind_sig_geno <- left_join(select(geno, taxa, chr1_438155938, chr1_193505583),
                          select(Ind, Taxa, ARDSI_Cw), by = c("taxa" = "Taxa"))

box_base <- list(geom_boxplot(color = "black", outlier.shape = NA),
                 geom_jitter(color = "blue", shape = 1, width = 0.2, size = 1),
                 my_theme)

Fig_4d <- ggplot(Ind_sig_geno, aes(x = factor(chr1_438155938), y = ARDSI_Cw)) + box_base + labs(x = "chr1_438155938", y = "ARDSI_Cw")
Fig_4e <- ggplot(Ind_sig_geno, aes(x = factor(chr1_193505583), y = ARDSI_Cw)) + box_base + labs(x = "chr1_193505583", y = "ARDSI_Cw")

# 3.3 Assembly
Fig_4a_plot   <- make_pdf_panel(pdf_manhattan_indices)
Fig_4a_legend <- make_pdf_panel(pdf_legend_indices)
Fig_4c         <- make_pdf_panel(pdf_qq_indices)

row1      <- plot_grid(Fig_4a_plot, Fig_4a_legend, ncol = 2, rel_widths = c(1, 0.2), labels = c("a", ""))
row_de    <- plot_grid(Fig_4d, Fig_4e, ncol = 2, labels = c("d", "e"), label_size = 10, align = "h")
right_col <- plot_grid(Fig_4c, row_de, ncol = 1, labels = c("c", ""), rel_heights = c(1.1, 1))
bottom_row<- plot_grid(Fig_4b, right_col, ncol = 2, labels = c("b", ""), rel_widths = c(1.3, 1), align = "h", axis = "tb")

final_plot_4 <- plot_grid(row1, bottom_row, ncol = 1, rel_heights = c(1, 1.2), align = "v", axis = "l")
ggsave(file.path(fig_dir, "Figure_4.pdf"), final_plot_4, width = 18, height = 22, units = "cm", device = cairo_pdf, dpi = 600, bg = "white")

# 4. Figure 5: Single Wavelength (SW) GWAS -------------------------------------
SW <- read_delim(sw_csv_path, skip = 1, col_names = FALSE, delim = ",")
colnames(SW) <- c("Wavelength","Model","SNP","Chr","Pos","p_value","MAF","nod","FDR_Adjusted_p_value","effect","threshold")

spectral_bg <- list(
  geom_rect(aes(xmin=400, xmax=700, ymin=-Inf, ymax=Inf), fill="grey90", color=NA, alpha=0.05),
  geom_rect(aes(xmin=800, xmax=1300, ymin=-Inf, ymax=Inf), fill="grey90", color=NA, alpha=0.05),
  geom_rect(aes(xmin=1550, xmax=1800, ymin=-Inf, ymax=Inf), fill="grey90", color=NA, alpha=0.05),
  geom_rect(aes(xmin=2000, xmax=2400, ymin=-Inf, ymax=Inf), fill="grey90", color=NA, alpha=0.05),
  geom_vline(xintercept = c(1000, 1800), lty = 3)
)

Fig_5b <- ggplot(filter(SW, Model == "BLINK"), aes(as.numeric(Wavelength), SNP, color = FDR_Adjusted_p_value)) + 
  spectral_bg + geom_point(shape = 15, size = 2) + xlim(400, 2500) +
  scale_color_gradient(low = "#67000d", high = "#fee0d2") + my_theme +
  labs(y = "Chromosome Position", x = "Wavelength (nm)", colour = "FDR_p-value")

Fig_5a <- make_pdf_panel(pdf_manhattan_sw)
Fig_5c <- make_pdf_panel(pdf_qq_sw)

bottom_row_5 <- plot_grid(Fig_5b, Fig_5c, ncol = 2, labels = c("b", "c"), rel_widths = c(2.5, 1), align = "h", axis = "tb")
final_plot_5 <- plot_grid(Fig_5a, bottom_row_5, ncol = 1, labels = c("a", ""), rel_heights = c(1, 1.7), align = "v", axis = "l")
ggsave(file.path(fig_dir, "Figure_5.pdf"), final_plot_5, width = 18, height = 20, units = "cm", device = cairo_pdf, dpi = 600, bg = "white")

# 5. Figure 6: HSC Segment GWAS ------------------------------------------------
Seg <- read_delim(hsc_csv_path, col_names = TRUE, delim = ",") %>%
  mutate(Segment = str_replace_all(Segment, "Seg|\\(NYC\\)", "") %>% str_squish()) %>%
  filter(!str_detect(Segment, "\\(Kansas\\)"))

Fig_6a <- ggplot(filter(Seg, Model == "BLINK"), aes(as.factor(Segment), SNP, fill = FDR_Adjusted_p_value)) + 
  geom_tile(width = 0.8) + scale_fill_gradient(low = "#67000d", high = "#fee0d2") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 5),
                   legend.title = element_text(size = 6)) +
  labs(y = "Chromosome Position", x = "Segments", fill = "FDR_p-value")

final_used <- readRDS(final_used_path)
HSC_fig <- merge(final_used, filter(Seg, Model == "BLINK"), by.x = "segment", by.y = "Segment")

Fig_6b <- ggplot(filter(HSC_fig, Model == "BLINK"), aes(as.numeric(wavelength), SNP, color = FDR_Adjusted_p_value)) + 
  spectral_bg + geom_point(shape = 15, size = 2) + xlim(400, 2500) +
  scale_color_gradient(low = "#67000d", high = "#fee0d2") + my_theme +
  theme(legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 5), legend.title = element_text(size = 6)) +
  labs(y = "Chromosome Position", x = "Wavelength (nm)", colour = "FDR_p-value")

HSC <- readRDS(hsc_pc_rds_path)
HSC_sig_geno <- left_join(select(geno, taxa, chr1_131409683), select(HSC, Taxa, Seg37, Seg38), by = c("taxa" = "Taxa"))
Fig_6e <- ggplot(HSC_sig_geno, aes(x = factor(chr1_131409683), y = Seg37)) + box_base + labs(x = "chr1_131409683", y = "Segment 37")
Fig_6f <- ggplot(HSC_sig_geno, aes(x = factor(chr1_131409683), y = Seg38)) + box_base + labs(x = "chr1_131409683", y = "Segment 38")

Fig_6c <- make_pdf_panel(pdf_qq_seg37)
Fig_6d <- make_pdf_panel(pdf_qq_seg38)

bottom_grid_6 <- plot_grid(Fig_6c, Fig_6d, Fig_6e, Fig_6f, ncol = 4, labels = c("c", "d", "e", "f"), 
                           rel_widths = c(1, 1, 1.2, 1.2), label_size = 10, align = "hv", axis = "tblr")

final_plot_6 <- plot_grid(Fig_6a, Fig_6b, bottom_grid_6, ncol = 1, labels = c("a", "b", ""), 
                          label_size = 9, rel_heights = c(2.2, 2.2, 1))

ggsave(file.path(fig_dir, "Figure_6.pdf"), final_plot_6, width = 18, height = 24, units = "cm", device = cairo_pdf, dpi = 600, bg = "white")