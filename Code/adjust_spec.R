## ================================================================
## Spectra adjustment (Cleaned & Optimized for GitHub)
## ================================================================
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(tidyverse)
  library(scales)
  library(ggplot2)
})

# --------------------------
# 1. Global Settings & Theme
# --------------------------
par(
  cex.main = 0.8, font.main = 1, cex.lab = 0.7, font.lab = 1,
  cex.axis = 0.6, font.axis = 1, mar = c(4,4,2,1) + 0.1, bg = "transparent"
)

my_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    text = element_text(size = 8),
    legend.background = element_rect(fill = "transparent")
  )

# --------------------------
# 2. User paths (EDIT THESE)
# --------------------------
spec_rds <- "data/Pheno_leafspec.rds"
cyto_tsv <- "data/cyto.tsv"
az_rds   <- "data/R_AZ.rds"
out_dir  <- "output"
fig_dir  <- "figures"

if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# --------------------------
# 3. Data Processing
# --------------------------
Spec <- readRDS(spec_rds)
Cyto <- fread(cyto_tsv, data.table = FALSE)
AZ   <- readRDS(az_rds)

# Clean Identifiers
AZ$Taxa <- gsub("-", "_", gsub("M", "L", AZ$Genotype_ID))
split_Line <- strsplit(AZ$Genotype_ID, "-")
AZ$genotype_group  <- sapply(split_Line, `[`, 1)
AZ$genotype_number <- sapply(split_Line, `[`, 2)

# Subset and Filter
az <- AZ %>% 
  filter(genotype_group %in% c("M1", "M2")) %>%
  select(Plant_ID, Leaf_Num, Batch, Day, Ctime, Taxa, genotype_group, genotype_number) %>%
  distinct(Taxa, .keep_all = TRUE)

UTonly <- AZ %>% filter(genotype_group == "UT")

# Merge Datasets
pheno <- Cyto %>%
  inner_join(Spec, by = c("Line" = "Taxa")) %>%
  inner_join(az, by = c("Line" = "Taxa"))

# --------------------------
# 4. Indices Calculation Function
# --------------------------
# Indices defined by wavelength offsets (base wavelength index depends on column start)
calc_indices <- function(df, start_col_offset = 0) {
  # Helper to get column by wavelength
  get_w <- function(w) df[, w - 349 + start_col_offset]
  
  df$NDWI      <- (get_w(865) - get_w(1614)) / (get_w(865) + get_w(1614))
  df$CIre      <- get_w(783) / get_w(704) - 1
  df$CCI       <- (get_w(560) - get_w(664)) / (get_w(560) + get_w(664))
  df$ARDSI_Cab <- (get_w(750) - get_w(730)) / (get_w(770) + get_w(720))
  df$ARDSI_Cw  <- (get_w(1360) - get_w(1080)) / (get_w(1560) + get_w(1240))
  df$ARDSI_Cm  <- (get_w(2200) - get_w(1640)) / (get_w(2240) + get_w(1720))
  return(df)
}

pheno  <- calc_indices(pheno, start_col_offset = 2) %>% relocate(NDWI:ARDSI_Cm, .before = Plant_ID)
UTonly <- calc_indices(UTonly, start_col_offset = 0) %>% relocate(NDWI:ARDSI_Cm, .before = Genotype_ID)

# Optional Export Raw
fwrite(pheno,  file.path(out_dir, "pheno_raw.txt"), sep = "\t")
fwrite(UTonly, file.path(out_dir, "UT_raw.txt"), sep = "\t")

# --------------------------
# 5. Factor Preparation
# --------------------------
pheno <- pheno %>%
  mutate(
    BA  = factor(Batch),
    DT  = factor(Day),
    CT  = factor(Ctime),
    REP = factor(sapply(strsplit(Line, "_"), `[`, 1)),
    PL  = factor(sapply(strsplit(Line, "_"), `[`, 2)),
    MO  = factor(Cytoplasm),
    LN  = factor(Leaf_Num)
  )

UTonly$BA <- factor(UTonly$Batch)

# --------------------------
# 6. Adjustment Loop
# --------------------------
n_wvl <- 2157
n_pl  <- nrow(pheno)
newy  <- matrix(0, nrow = n_pl, ncol = n_wvl)

terms_vec <- c("LN", "REP", "ba", "LN:REP", "LN:ba", "DT:CT", "MO", "Residuals")
ss_list   <- vector("list", 2157)
wl_vec    <- numeric(2157)

# Batch mean adjustment and environment correction
for(i in 3:2159) {
  pheno$y <- pheno[, i]
  BAmeans <- tapply(UTonly[, i-2], UTonly$BA, mean, na.rm = TRUE)
  
  # Map batch means (handling first batch if necessary)
  utba  <- c(mean(BAmeans, na.rm = TRUE), BAmeans[1:34])
  pheno$ba <- utba[as.numeric(pheno$BA)]
  
  # Fit Linear Model
  lmy2 <- lm(y ~ LN + REP + ba + LN:(REP + ba) + DT:CT + MO, data = pheno)
  
  # Store corrected values (Residuals + Mean)
  newy[, i-2] <- residuals(lmy2) + mean(pheno$y)
  
  # Variance Partitioning (Type I SS)
  a  <- anova(lmy2)
  ss <- setNames(rep(NA_real_, length(terms_vec)), terms_vec)
  rn <- rownames(a)
  
  present <- intersect(terms_vec, rn)
  ss[present] <- a[present, "Sum Sq"]
  if ("Residuals" %in% rn) ss["Residuals"] <- a["Residuals", "Sum Sq"]
  
  ss_list[[i-2]] <- ss / sum(ss, na.rm = TRUE)
  wl_vec[i-2]    <- i - 3 + 350
}

# --------------------------
# 7. Results Export
# --------------------------
newy2 <- data.frame(Taxa = pheno$Line, REP = pheno$REP, PL = pheno$PL, newy)
colnames(newy2)[4:2160] <- c(seq(350, 2500), "NDWI", "CIre", "CCI", "ARDSI_Cab", "ARDSI_Cw", "ARDSI_Cm")

fwrite(newy2, file.path(out_dir, "pheno.txt"), sep = "\t")
saveRDS(newy2[, c(1, 4:2154)], file.path(out_dir, "Pheno_corr.rds"))
saveRDS(newy2[, c(1, 2155:2160)], file.path(out_dir, "Indices_corr.rds"))

# --------------------------
# 8. Plotting
# --------------------------
var_df <- do.call(rbind, lapply(seq_along(ss_list), function(k) {
  data.frame(wavelength = wl_vec[k], term = names(ss_list[[k]]), prop = as.numeric(ss_list[[k]]))
}))

var_df$group <- case_when(
  var_df$term == "LN" ~ "Technique",
  var_df$term %in% c("REP", "ba", "LN:REP", "LN:ba") ~ "Space + tech",
  var_df$term == "DT:CT" ~ "Time",
  var_df$term == "MO" ~ "Maternal cytoplasm",
  var_df$term == "Residuals" ~ "Residuals"
)

group_df <- var_df %>%
  group_by(wavelength, group) %>%
  summarise(prop = sum(prop, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = factor(group, levels = c("Space + tech", "Technique", "Time", "Maternal cytoplasm", "Residuals")))

Fig_variance_part <- ggplot(group_df, aes(x = wavelength, y = prop, fill = group)) +
  geom_area() +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(500, 2500, by = 500), expand = c(0,0)) +
  scale_fill_manual(values = c(
    "Residuals" = "grey70", "Space + tech" = "#009E73", 
    "Technique" = "#CC79A7", "Time" = "#56B4E9", "Maternal cytoplasm" = "#D55E00"
  )) +
  labs(x = "Wavelength (nm)", y = "Variance explained (% of total)", fill = "Effect group") +
  my_theme

ggsave(file.path(fig_dir, "Figure_2.pdf"), Fig_variance_part, width = 18, height = 7, units = "cm", dpi = 600)