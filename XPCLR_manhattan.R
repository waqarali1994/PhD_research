
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
getwd()
# ---------- Load XP-CLR ----------
setwd("~/Documents/PhD work/Yield_project/XPCLR/xpclr_output/Tropical_temperate/")
files <- list.files(pattern = "^xpclr_chr[0-9]+_.*_vs_Tropical\\.txt$")

all_data <- lapply(files, function(f) {
  d <- fread(f)
  chr <- as.integer(gsub("xpclr_chr([0-9]+)_.*", "\\1", f))  # numeric
  pop <- gsub(".*chr[0-9]+_(.*)_vs_Tropical.txt", "\\1", f)
  mid <- (d[[3]] + d[[4]]) / 2
  data.frame(chr = chr, pop = pop, pos = mid, xpclr = d[[12]])
})

df <- rbindlist(all_data)
df <- df[!is.na(df$xpclr), ]

# ---------- Compute cumulative positions ----------
chr_lengths_pop <- df %>%
  group_by(pop, chr) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%
  arrange(pop, chr) %>%
  group_by(pop) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

df <- df %>%
  left_join(chr_lengths_pop, by = c("pop", "chr")) %>%
  mutate(cum_pos = pos + chr_start)

axis_df <- chr_lengths_pop %>%
  group_by(pop) %>%
  mutate(center = chr_start + chr_len / 2) %>%
  ungroup()

df <- df %>%
  group_by(pop) %>%
  mutate(threshold = quantile(xpclr, 0.995, na.rm = TRUE))

# ---------- Load GWAS Hits ----------
gwas_hits <- fread("newyieldgwas_z.csv")  # Columns: CHROM, POS, Trait, support
gwas_hits <- gwas_hits %>%
  rename(chr = CHROM, pos = POS) %>%
  filter(support >= 0.1) %>%
  mutate(chr = as.integer(gsub("chr", "", chr))) %>%
  mutate(SNP_label = paste0("chr", chr, ":", pos))

# Assign manual trait colors
trait_colors <- c(
  "Ear Length"        = "#4477AA",
  "Ear Weight"        = "#66CCEE",
  "Ear Width"         = "#EE6677",
  "Kernel Row Number" = "#CCBB44",
  "Hundred Kernel Weight"     = "orange",
  "Kernels Per Row"   = "#228B22",
  "Kernel Weight Per Ear"       = "#BBBBBB"
)

# ---------- Plot by population ----------
for (pop_name in unique(df$pop)) {
  df_pop <- df %>% filter(pop == pop_name)
  chr_layout <- chr_lengths_pop %>% filter(pop == pop_name)
  axis_df_pop <- axis_df %>% filter(pop == pop_name)
  
  # Background chromosome rectangles
  rects <- chr_layout %>%
    mutate(
      xmin = chr_start / 1e6,
      xmax = (chr_start + chr_len) / 1e6,
      fill = factor(chr %% 2)
    )
  
  gwas_cum <- gwas_hits %>%
    semi_join(chr_layout, by = "chr") %>%
    left_join(chr_layout, by = "chr") %>%
    mutate(cum_pos = pos + chr_start)
  
  # Optional: highlight top 3 per trait
  gwas_cum_top <- gwas_cum  # label all GWAS SNPs with support ≥ 0.2
  
  
  p <- ggplot(df_pop, aes(x = cum_pos / 1e6, y = xpclr)) +
    # Chromosome shading
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = -20, ymax = 0, fill = fill),
              inherit.aes = FALSE, alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("black", "grey")) +
    
    # Main XPCLR scatter
    geom_point(color = "grey30", size = 0.4) +
    
    # XPCLR threshold
    geom_hline(yintercept = unique(df_pop$threshold), linetype = "dashed", color = "red", linewidth = 1) +
    
    # GWAS SNP vertical lines
    geom_segment(
      data = gwas_cum,
      aes(x = cum_pos / 1e6, xend = cum_pos / 1e6,
          y = 0, yend = unique(df_pop$threshold) * 1.8, color = Trait),
      linewidth = 0.8
    ) +
    
    # SNP labels close to lines
    geom_text(
      data = gwas_cum_top,
      aes(x = cum_pos / 1e6,
          y = unique(df_pop$threshold) * 1.8,
          label = SNP_label,
          color = Trait),
      angle = 90,
      hjust = 0,
      vjust = 0.4,
      size = 3,
      fontface = "bold",
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = trait_colors) +
    scale_x_continuous(
      breaks = axis_df_pop$center / 1e6,
      labels = paste0("chr", axis_df_pop$chr)
    ) +
    labs(
      x = "",
      y = "XP-CLR Score",
      title = paste(pop_name, "vs Parviglumis")
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(family = "Courier", size = 10)
    )
  
  ggsave(
    filename = paste0("XPCLR_TraitColored_Final_", pop_name, ".png"),
    plot = p,
    width = 14, height = 6, dpi = 300
  )
}




all_plots <- list()
populations <- c("IOD", "NSS", "SS", "Tropical")

for (i in seq_along(populations)) {
  pop_name <- populations[i]
  df_pop <- df %>% filter(pop == pop_name)
  chr_layout <- chr_lengths_pop %>% filter(pop == pop_name)
  axis_df_pop <- axis_df %>% filter(pop == pop_name)
  
  # ... Your existing plot code setup for df_pop and gwas_cum ...
  
  # Only show x-axis text for the bottom plot
  x_axis_theme <- if (i == length(populations)) {
    theme(axis.text.x = element_text(angle = 0, size = 12, color = "black"))
  } else {
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  }
  
  p <- ggplot(...) +  # Your existing plot building logic
    x_axis_theme   # 👈 add conditional x-axis theme
  
  all_plots[[pop_name]] <- p
}


library(patchwork)

combined_plot <- all_plots[["IOD"]] /
  all_plots[["NSS"]] /
  all_plots[["SS"]] /
  all_plots[["Tropical"]] +
  plot_layout(ncol = 1, guides = "collect") +  # Share legend
  plot_annotation(title = "XP-CLR Selection vs GWAS SNPs across Populations")

# Save the final figure
ggsave("XPCLR_Stacked_OneXAxis.png", combined_plot, width = 14, height = 20, dpi = 300)










library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Load XP-CLR results
setwd("~/Documents/PhD work/Yield project/XPCLR/xpclr_output/")
files <- list.files(pattern = "^xpclr_chr[0-9]+_.*_vs_parviglumis\\.txt$")

all_data <- lapply(files, function(f) {
  d <- fread(f)
  chr <- as.integer(gsub("xpclr_chr([0-9]+)_.*", "\\1", f))
  pop <- gsub(".*chr[0-9]+_(.*)_vs_parviglumis.txt", "\\1", f)
  mid <- (d[[3]] + d[[4]]) / 2
  data.frame(chr = chr, pop = pop, pos = mid, xpclr = d[[12]])
})

df <- rbindlist(all_data)
df <- df[!is.na(df$xpclr), ]

# Compute positions
chr_lengths_pop <- df %>%
  group_by(pop, chr) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%
  arrange(pop, chr) %>%
  group_by(pop) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

df <- df %>%
  left_join(chr_lengths_pop, by = c("pop", "chr")) %>%
  mutate(cum_pos = pos + chr_start)

axis_df <- chr_lengths_pop %>%
  group_by(pop) %>%
  mutate(center = chr_start + chr_len / 2) %>%
  ungroup()

df <- df %>%
  group_by(pop) %>%
  mutate(threshold = quantile(xpclr, 0.995, na.rm = TRUE))

# GWAS hits
gwas_hits <- fread("newyieldgwas_z.csv")  # CHROM, POS, Trait, support
gwas_hits <- gwas_hits %>%
  rename(chr = CHROM, pos = POS) %>%
  filter(support >= 0.2) %>%
  mutate(chr = as.integer(gsub("chr", "", chr))) %>%
  mutate(SNP_label = paste0("chr", chr, ":", pos))

trait_colors <- c(
  "Ear Length" = "#4477AA", "Ear Weight" = "#66CCEE",
  "Ear Width" = "#EE6677", "Kernel Row Number" = "#CCBB44",
  "Kernel Weight" = "orange", "Kernels Per Row" = "#228B22",
  "Plot Weight" = "#BBBBBB"
)

# Plot loop
all_plots <- list()
populations <- c("IOD", "NSS", "SS", "Tropical")

for (i in seq_along(populations)) {
  pop_name <- populations[i]
  df_pop <- df %>% filter(pop == pop_name)
  chr_layout <- chr_lengths_pop %>% filter(pop == pop_name)
  axis_df_pop <- axis_df %>% filter(pop == pop_name)
  
  rects <- chr_layout %>%
    mutate(xmin = chr_start / 1e6,
           xmax = (chr_start + chr_len) / 1e6,
           fill = factor(chr %% 2))
  
  gwas_cum <- gwas_hits %>%
    semi_join(chr_layout, by = "chr") %>%
    left_join(chr_layout, by = "chr") %>%
    mutate(cum_pos = pos + chr_start)
  
  gwas_cum_top <- gwas_cum  # label all
  
  # Axis theme control
  x_axis_theme <- if (i == length(populations)) {
    theme(axis.text.x = element_text(angle = 0, size = 12, color = "black"))
  } else {
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  p <- ggplot(df_pop, aes(x = cum_pos / 1e6, y = xpclr)) +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = -20, ymax = 0, fill = fill),
              inherit.aes = FALSE, alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("black", "grey")) +
    geom_point(color = "grey30", size = 0.4) +
    geom_hline(yintercept = unique(df_pop$threshold), linetype = "dashed", color = "red") +
    geom_segment(data = gwas_cum,
                 aes(x = cum_pos / 1e6, xend = cum_pos / 1e6,
                     y = 0, yend = unique(df_pop$threshold) * 1.8, color = Trait),
                 linewidth = 0.8) +
    geom_text(data = gwas_cum_top,
              aes(x = cum_pos / 1e6,
                  y = unique(df_pop$threshold) * 1.8,
                  label = SNP_label,
                  color = Trait),
              angle = 90, hjust = 0, vjust = 0.4, size = 3, fontface = "bold",
              show.legend = FALSE) +
    scale_color_manual(values = trait_colors) +
    scale_x_continuous(
      breaks = axis_df_pop$center / 1e6,
      labels = paste0("chr", axis_df_pop$chr)
    ) +
    labs(
      x = "", y = "XP-CLR Score",
      title = paste(pop_name, "vs Parviglumis")
    ) +
    theme_classic(base_size = 16) %+replace%
    theme(
      axis.text.x = if (i == length(populations)) element_text(angle = 0, size = 12, color = "black") else element_blank(),
      axis.ticks.x = if (i == length(populations)) element_line() else element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(family = "Courier", size = 10)
    )
  
  
  all_plots[[pop_name]] <- p
}

library(patchwork)

combined_plot <- (
  all_plots[["IOD"]] /
    all_plots[["NSS"]] /
    all_plots[["SS"]] /
    all_plots[["Tropical"]]
) +
  plot_layout(ncol = 1, guides = "collect") +  # Collect all legends
  plot_annotation(title = "") &
  theme(legend.position = "bottom")  # 👈 apply bottom legend here globally

# Save
ggsave("XPCLR_Stacked_AllPopulations.png", combined_plot, width = 14, height = 20, dpi = 300)










library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(grid)

# ----- Theme for all plots -----
common_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.2),
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.background = element_blank(),
    legend.key.size = unit(4, "lines"),
    legend.spacing.x = unit(1, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "black")
  )

# ---------- Load XP-CLR ----------
setwd("~/Documents/PhD work/Yield project/XPCLR/xpclr_output/")
files <- list.files(pattern = "^xpclr_chr[0-9]+_.*_vs_parviglumis\\.txt$")

all_data <- lapply(files, function(f) {
  d <- fread(f)
  chr <- as.integer(gsub("xpclr_chr([0-9]+)_.*", "\\1", f))
  pop <- gsub(".*chr[0-9]+_(.*)_vs_parviglumis.txt", "\\1", f)
  mid <- (d[[3]] + d[[4]]) / 2
  data.frame(chr = chr, pop = pop, pos = mid, xpclr = d[[12]])
})

df <- rbindlist(all_data)
df <- df[!is.na(df$xpclr), ]

# ---------- Compute positions ----------
chr_lengths_pop <- df %>%
  group_by(pop, chr) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%
  arrange(pop, chr) %>%
  group_by(pop) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

df <- df %>%
  left_join(chr_lengths_pop, by = c("pop", "chr")) %>%
  mutate(cum_pos = pos + chr_start)

axis_df <- chr_lengths_pop %>%
  group_by(pop) %>%
  mutate(center = chr_start + chr_len / 2) %>%
  ungroup()

df <- df %>%
  group_by(pop) %>%
  mutate(threshold = quantile(xpclr, 0.995, na.rm = TRUE))

# ---------- GWAS ----------
gwas_hits <- fread("newyieldgwas_z.csv") %>%
  rename(chr = CHROM, pos = POS) %>%
  filter(support >= 0.2) %>%
  mutate(chr = as.integer(gsub("chr", "", chr))) %>%
  mutate(SNP_label = paste0("chr", chr, ":", pos))

  # Assign manual trait colors
  trait_colors <- c(
    "Ear Length"        = "#4477AA",
    "Ear Weight"        = "#66CCEE",
    "Ear Width"         = "#EE6677",
    "Kernel Row Number" = "#CCBB44",
    "Hundred Kernel Weight"     = "orange",
    "Kernels Per Row"   = "#228B22",
    "Kernel Weight Per Ear"       = "#BBBBBB"
  )
# ---------- Create plots ----------
all_plots <- list()
populations <- c("IOD", "NSS", "SS", "Tropical")

for (i in seq_along(populations)) {
  pop_name <- populations[i]
  df_pop <- df %>% filter(pop == pop_name)
  chr_layout <- chr_lengths_pop %>% filter(pop == pop_name)
  axis_df_pop <- axis_df %>% filter(pop == pop_name)
  
  rects <- chr_layout %>%
    mutate(xmin = chr_start / 1e6,
           xmax = (chr_start + chr_len) / 1e6,
           fill = factor(chr %% 2))
  
  gwas_cum <- gwas_hits %>%
    semi_join(chr_layout, by = "chr") %>%
    left_join(chr_layout, by = "chr") %>%
    mutate(cum_pos = pos + chr_start)
  
  gwas_cum_top <- gwas_cum  # label all
  
  # X-axis theme control
  x_axis_theme <- if (pop_name %in% c("SS", "Tropical")) {
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22))
  } else {
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  
  # Build plot
  p <- ggplot(df_pop, aes(x = cum_pos / 1e6, y = xpclr)) +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = -20, ymax = 0, fill = fill),
              inherit.aes = FALSE, alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("black", "grey")) +
    geom_point(color = "grey30", size = 0.4) +
    geom_hline(yintercept = unique(df_pop$threshold), linetype = "dashed", color = "red", linewidth = 1) +
    geom_segment(data = gwas_cum,
                 aes(x = cum_pos / 1e6, xend = cum_pos / 1e6,
                     y = 0, yend = unique(df_pop$threshold) * 1.8, color = Trait),
                 linewidth = 0.8) +
    geom_text(data = gwas_cum_top,
              aes(x = cum_pos / 1e6,
                  y = unique(df_pop$threshold) * 1.8,
                  label = SNP_label, color = Trait),
              angle = 90, hjust = 0, vjust = 0.4, size = 5,
              fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = trait_colors) +
    scale_x_continuous(
      breaks = axis_df_pop$center / 1e6,
      labels = paste0("chr", axis_df_pop$chr)
    ) +
    labs(
      title = paste(pop_name, "vs Parviglumis")
    ) +
    common_theme +
    x_axis_theme
  
  all_plots[[pop_name]] <- p
}

# ---------- Patchwork stack with one Y-axis and shared legend ----------
stacked_plot <- (all_plots[["IOD"]] / all_plots[["NSS"]] / all_plots[["SS"]] / all_plots[["Tropical"]]) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = 'top')

# Add shared Y-axis
final_plot <- wrap_elements(grid::textGrob("XP-CLR Score", rot = 90, gp = gpar(fontsize = 22, fontface = "bold"))) +
  stacked_plot +
  plot_layout(widths = c(0.03, 1))

# ---------- Save plot ----------
ggsave("XPCLR_AllPops_Stacked_Beautiful.png", final_plot, width = 24, height = 18, dpi = 300)










library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(grid)
library(cowplot)

# ----- Common Theme -----
common_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.2),
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.background = element_blank(),
    legend.key.size = unit(4, "lines"),
    legend.spacing.x = unit(1, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "black")
  )

# ----- Load XP-CLR -----
setwd("~/Documents/PhD work/Yield_project/XPCLR/xpclr_output/Tropical_temperate/")
files <- list.files(pattern = "^xpclr_chr[0-9]+_.*_vs_Tropical\\.txt$")

all_data <- lapply(files, function(f) {
  d <- fread(f)
  chr <- as.integer(gsub("xpclr_chr([0-9]+)_.*", "\\1", f))
  pop <- gsub(".*chr[0-9]+_(.*)_vs_Tropical.txt", "\\1", f)
  mid <- (d[[3]] + d[[4]]) / 2
  data.frame(chr = chr, pop = pop, pos = mid, xpclr = d[[12]])
})

df <- rbindlist(all_data)
df <- df[!is.na(df$xpclr) & xpclr <= 600, ]  # 🔍 Keep only XPCLR <= 600

# ----- Compute cumulative positions -----
chr_lengths_pop <- df %>%
  group_by(pop, chr) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%
  arrange(pop, chr) %>%
  group_by(pop) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

df <- df %>%
  left_join(chr_lengths_pop, by = c("pop", "chr")) %>%
  mutate(cum_pos = pos + chr_start)

axis_df <- chr_lengths_pop %>%
  group_by(pop) %>%
  mutate(center = chr_start + chr_len / 2) %>%
  ungroup()

df <- df %>%
  group_by(pop) %>%
  mutate(threshold = quantile(xpclr, 0.995, na.rm = TRUE))

# ----- GWAS Hits -----
gwas_hits <- fread("newyieldgwas_z.csv") %>%
  rename(chr = CHROM, pos = POS) %>%
  filter(support >= 0.1) %>%
  mutate(chr = as.integer(gsub("chr", "", chr))) %>%
  mutate(SNP_label = paste0("chr", chr, ":", pos))

# Trait colors
trait_colors <- c(
  "Ear Length" = "#4477AA",
  "Ear Weight" = "#66CCEE",
  "Ear Width" = "#EE6677",
  "Kernel Row Number" = "#CCBB44",
  "Hundred Kernel Weight" = "orange",
  "Kernels Per Row" = "#228B22",
  "Kernel Weight Per Ear" = "#BBBBBB"
)

# ----- Plot Loop -----
all_plots <- list()
populations <- c("IOD", "NSS", "SS", "Tropical")

for (i in seq_along(populations)) {
  pop_name <- populations[i]
  df_pop <- df %>% filter(pop == pop_name)
  chr_layout <- chr_lengths_pop %>% filter(pop == pop_name)
  axis_df_pop <- axis_df %>% filter(pop == pop_name)
  
  rects <- chr_layout %>%
    mutate(xmin = chr_start / 1e6,
           xmax = (chr_start + chr_len) / 1e6,
           fill = factor(chr %% 2))
  
  gwas_cum <- gwas_hits %>%
    semi_join(chr_layout, by = "chr") %>%
    left_join(chr_layout, by = "chr") %>%
    mutate(cum_pos = pos + chr_start)
  
  gwas_cum_top <- gwas_cum  # Label all
  
  # Axis formatting for x-axis
  x_axis_theme <- if (pop_name %in% c("SS", "Tropical")) {
    theme(axis.text.x = element_text(angle = 0, size = 22, color = "black"))
  } else {
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  p <- ggplot(df_pop, aes(x = cum_pos / 1e6, y = xpclr)) +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = -20, ymax = 0, fill = fill),
              inherit.aes = FALSE, alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("black", "grey")) +
    geom_point(color = "grey30", size = 0.4) +
    geom_hline(yintercept = unique(df_pop$threshold), linetype = "dashed", color = "red", linewidth = 1) +
    geom_segment(data = gwas_cum,
                 aes(x = cum_pos / 1e6, xend = cum_pos / 1e6,
                     y = 0, yend = unique(df_pop$threshold) * 1.8, color = Trait),
                 linewidth = 0.8) +
    geom_text(data = gwas_cum_top,
              aes(x = cum_pos / 1e6,
                  y = unique(df_pop$threshold) * 1.8,
                  label = SNP_label, color = Trait),
              angle = 90, hjust = 0, vjust = 0.4, size = 5,
              fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = trait_colors) +
    scale_x_continuous(
      breaks = axis_df_pop$center / 1e6,
      labels = paste0("chr", axis_df_pop$chr)
    ) +
    coord_cartesian(ylim = c(0, 600)) +  
    labs(title = paste(pop_name, "vs Parviglumis")) +
    common_theme + x_axis_theme
  
  all_plots[[pop_name]] <- p
}

# ----- Arrange plots in 2x2 with shared y-axis and legend -----
combined <- (all_plots[["IOD"]] | all_plots[["NSS"]]) /
  (all_plots[["SS"]]  | all_plots[["Tropical"]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# Add shared y-axis label
final_plot <- wrap_elements(grid::textGrob("XP-CLR Score", rot = 90, gp = gpar(fontsize = 22, fontface = "bold"))) +
  combined +
  plot_layout(widths = c(0.03, 1))


ggsave("XPCLR_AllPops_Stacked_Beautiful.png", final_plot, width = 24, height = 18, dpi = 300)






















library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(grid)
library(cowplot)
# ---------- Theme ----------
common_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.2),
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.background = element_blank(),
    legend.key.size = unit(4, "lines"),
    legend.spacing.x = unit(1, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "black")
  )
# ---------- Load XP-CLR ----------
setwd("~/Documents/PhD work/Yield_project/XPCLR/xpclr_output/Tropical_temperate/")
files <- list.files(pattern = "^xpclr_chr[0-9]+_.*_vs_Tropical\\.txt$")
all_data <- lapply(files, function(f) {
  d <- fread(f)
  chr <- as.integer(gsub("xpclr_chr([0-9]+)_.*", "\\1", f))
  pop <- gsub(".*chr[0-9]+_(.*)_vs_Tropical.txt", "\\1", f)
  mid <- (d$start + d$stop) / 2
  data.frame(chr = chr, pop = pop, pos = mid, xpclr = d$xpclr)
})
df <- rbindlist(all_data)
df <- df[!is.na(df$xpclr) & df$xpclr <= 600, ]  # Filter top XP-CLR
# ---------- Cumulative Genome Position ----------
chr_lengths_pop <- df %>%
  group_by(pop, chr) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%
  arrange(pop, chr) %>%
  group_by(pop) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))
df <- df %>%
  left_join(chr_lengths_pop, by = c("pop", "chr")) %>%
  mutate(cum_pos = pos + chr_start)
axis_df <- chr_lengths_pop %>%
  group_by(pop) %>%
  mutate(center = chr_start + chr_len / 2) %>%
  ungroup()
df <- df %>%
  group_by(pop) %>%
  mutate(threshold = quantile(xpclr, 0.995, na.rm = TRUE))
# ---------- GWAS Hits ----------
gwas_hits <- fread("newyieldgwas_z.csv") %>%
  rename(chr = CHROM, pos = POS) %>%
  filter(support >= 0.1) %>%
  mutate(chr = as.integer(gsub("chr", "", chr))) %>%
  mutate(SNP_label = paste0("chr", chr, ":", pos))
# ---------- Manual Trait Colors ----------
trait_colors <- c(
  "Ear Length" = "#4477AA",
  "Ear Weight" = "#66CCEE",
  "Ear Width" = "#EE6677",
  "Kernel Row Number" = "#CCBB44",
  "Hundred Kernel Weight" = "orange",
  "Kernels Per Row" = "#228B22",
  "Kernel Weight Per Ear" = "#BBBBBB"
)
# ---------- Build Plots ----------
all_plots <- list()
populations <- c("IOD", "NSS", "SS")
for (pop_name in populations) {
  df_pop <- df %>% filter(pop == pop_name)
  chr_layout <- chr_lengths_pop %>% filter(pop == pop_name)
  axis_df_pop <- axis_df %>% filter(pop == pop_name)
  rects <- chr_layout %>%
    mutate(xmin = chr_start / 1e6, xmax = (chr_start + chr_len) / 1e6, fill = factor(chr %% 2))
  gwas_cum <- gwas_hits %>%
    semi_join(chr_layout, by = "chr") %>%
    left_join(chr_layout, by = "chr") %>%
    mutate(cum_pos = pos + chr_start)
  p <- ggplot(df_pop, aes(x = cum_pos / 1e6, y = xpclr)) +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -20, ymax = 0, fill = fill),
              inherit.aes = FALSE, alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("black", "grey")) +
    geom_point(color = "grey30", size = 0.4) +
    geom_hline(yintercept = unique(df_pop$threshold), linetype = "dashed", color = "red", linewidth = 1) +
    geom_segment(data = gwas_cum,
                 aes(x = cum_pos / 1e6, xend = cum_pos / 1e6,
                     y = 0, yend = unique(df_pop$threshold) * 1.8, color = Trait),
                 linewidth = 0.8) +
    geom_text(data = gwas_cum,
              aes(x = cum_pos / 1e6, y = unique(df_pop$threshold) * 1.8,
                  label = SNP_label, color = Trait),
              angle = 90, hjust = 0, vjust = 0.4, size = 5,
              fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = trait_colors) +
    scale_x_continuous(
      breaks = axis_df_pop$center / 1e6,
      labels = paste0("chr", axis_df_pop$chr)
    ) +
    coord_cartesian(ylim = c(0, 600)) +
    labs(title = paste(pop_name, "vs Tropical")) +
    common_theme +
    theme(axis.text.x = element_text(angle = 0, size = 22, color = "black"))
  all_plots[[pop_name]] <- p
}
# ---------- Combine and Save ----------
combined_plot <- wrap_elements(grid::textGrob("XP-CLR Score", rot = 90,
                                              gp = gpar(fontsize = 22, fontface = "bold"))) +
  (all_plots[["IOD"]] / all_plots[["NSS"]] / all_plots[["SS"]]) +
  plot_layout(ncol = 2, guides = "collect", widths = c(0.03, 1)) &
  theme(legend.position = "top")
ggsave("XPCLR_Temperate_vs_Tropical.png", combined_plot, width = 24, height = 18, dpi = 300)
