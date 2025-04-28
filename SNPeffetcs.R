
library(dplyr)

# Define which columns you want to keep
columns_to_keep <- c("genotype", "year", "location", "nitrogenTreatment",
                     "kernelRowNumber", "hundredKernelMass", "kernelMassPerEar",
                     "earWeight", "earWidth", "earLength", "kernelsPerRow")

# Apply filtering and selection
phenotypic_data_filtered <- phenotypic_data %>%
  filter(nitrogenTreatment == "Medium") %>%
  select(all_of(columns_to_keep))


write.csv(phenotypic_data_filtered,"filtered_phenotype_data.csv")



# Load dataset
snp_trait_raw <- fread("newyieldgwas_z.csv")  # Replace with your actual filename

# Combine CHROM and POS into a new SNP column
snp_trait_clean <- snp_trait_raw %>%
  mutate(SNP = paste0(CHROM, "_", POS)) %>%
  select(SNP, Trait)  # Reorder: SNP first, then Trait

# View result
print(snp_trait_clean)

# Save to CSV
fwrite(snp_trait_clean, "combinedsnp.csv", row.names = FALSE)



library(data.table)
library(dplyr)

# ---- Step 1: Load your raw SNP–trait mapping file ----
snp_trait_raw <- fread("combinedsnp.csv")  # Replace with your actual filename

# ---- Step 2: Combine CHROM and POS to create SNP ID ----
snp_trait_clean <- snp_trait_raw %>%
  mutate(SNP = paste0(CHROM, "_", POS)) %>%
  select(SNP, Trait)

# ---- Step 3: Define a mapping from readable trait names to exact column names ----
trait_name_map <- c(
  "Kernel Weight" = "hundredKernelMass",
  "Plot Weight" = "kernelMassPerEar",
  "Kernel Row Number" = "kernelRowNumber",
  "Kernels Per Row" = "kernelsPerRow",
  "Ear Width" = "earWidth",
  "Ear Length" = "earLength",
  "Ear Weight" = "earWeight"
)

# ---- Step 4: Standardize trait names ----
snp_trait_clean <- snp_trait_clean %>%
  mutate(Trait = trait_name_map[Trait]) %>%
  filter(!is.na(Trait))  # Remove any unmatched trait names

# ---- Step 5: Save the cleaned SNP–Trait map ----
fwrite(snp_trait_clean, "combinedsnp.csv", row.names = FALSE)



# ---- Load Required Libraries ----
library(data.table)
library(dplyr)
library(tidyr)

# ---- Step 1: Load SNP-Trait Mapping ----
snp_trait_map <- fread("combinedsnp.csv") %>%
  filter(!is.na(Trait) & Trait != "")  # Remove empty or NA trait names

# ---- Step 2: Load Genotype and Phenotype Data ----
allele_data <- fread("transposed.csv")           # SNPs as columns, 'genotype' in first column
phenotypic_data <- fread("filtered_phenotype_data.csv")

# ---- Step 3: Merge and Create Environment ----
merged_data <- inner_join(phenotypic_data, allele_data, by = "genotype") %>%
  filter(nitrogenTreatment == "Medium") %>%
  mutate(Environment = paste(location, year))

# ---- Step 4: No global homozygous filtering ----
filtered_data <- merged_data  # Keep full dataset, filter per SNP inside loop

# ---- Step 5: Function to Compute Per-SNP Effect Sizes ----
compute_all_effects <- function(data, snp_trait_map) {
  results <- list()
  
  for (i in seq_len(nrow(snp_trait_map))) {
    snp <- snp_trait_map$SNP[i]
    trait <- snp_trait_map$Trait[i]
    
    # Skip invalid mappings
    if (is.na(snp) || is.na(trait) || trait == "" || !(trait %in% colnames(data))) next
    
    message("Processing: ", snp, " → ", trait)
    
    # 🔸 Filter to homozygotes for this SNP only
    snp_data <- data %>%
      filter(.data[[snp]] %in% c("0|0", "1|1"))
    
    # Compute mean trait per genotype per environment
    effect_df <- snp_data %>%
      group_by(Env = Environment, Genotype = .data[[snp]]) %>%
      summarise(mean_trait = mean(.data[[trait]], na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Genotype, values_from = mean_trait)
    
    # Assign values safely
    ref_val <- if ("0|0" %in% colnames(effect_df)) effect_df$`0|0` else NA_real_
    alt_val <- if ("1|1" %in% colnames(effect_df)) effect_df$`1|1` else NA_real_
    
    # Compute percent effect size
    effect_df <- effect_df %>%
      mutate(
        SNP = snp,
        Trait = trait,
        REF = ref_val,
        ALT = alt_val,
        Effect_Size = ((ALT / REF) - 1) * 100
      ) %>%
      select(SNP, Trait, Env, REF, ALT, Effect_Size)
    
    results[[paste(snp, trait, sep = "_")]] <- effect_df
  }
  
  bind_rows(results)
}

# ---- Step 6: Run and Save Final Output ----
all_effects <- compute_all_effects(filtered_data, snp_trait_map)
fwrite(all_effects, "Final_Effect_Summary_SpecificSNPs.csv")













# ---- Load Required Libraries ----
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Step 1: Preprocess Phenotype Data ----
phenotypic_data <- fread("your_original_phenotype_file.csv")  # Update if needed

columns_to_keep <- c("genotype", "year", "location", "nitrogenTreatment",
                     "kernelRowNumber", "hundredKernelMass", "kernelMassPerEar",
                     "earWeight", "earWidth", "earLength", "kernelsPerRow")

phenotypic_data_filtered <- phenotypic_data %>%
  filter(nitrogenTreatment == "Medium") %>%
  select(all_of(columns_to_keep))

fwrite(phenotypic_data_filtered, "filtered_phenotype_data.csv")


# ---- Step 2: Prepare SNP-Trait Mapping ----
snp_trait_raw <- fread("newyieldgwas_z.csv")

trait_name_map <- c(
  "Kernel Weight" = "hundredKernelMass",
  "Plot Weight" = "kernelMassPerEar",
  "Kernel Row Number" = "kernelRowNumber",
  "Kernels Per Row" = "kernelsPerRow",
  "Ear Width" = "earWidth",
  "Ear Length" = "earLength",
  "Ear Weight" = "earWeight"
)

snp_trait_clean <- snp_trait_raw %>%
  mutate(SNP = paste0(CHROM, "_", POS)) %>%
  mutate(Trait = trait_name_map[Trait]) %>%
  filter(!is.na(Trait)) %>%
  select(SNP, Trait)

fwrite(snp_trait_clean, "combinedsnp.csv")


# ---- Step 3: Load and Merge Datasets ----
allele_data <- fread("transposed.csv")
phenotypic_data <- fread("filtered_phenotype_data.csv")
snp_trait_map <- fread("combinedsnp.csv")

merged_data <- inner_join(phenotypic_data, allele_data, by = "genotype") %>%
  filter(nitrogenTreatment == "Medium") %>%
  mutate(Environment = paste(location, year))


# ---- Step 4: Compute Effect Sizes (Per-SNP, Homozygous Only) ----
compute_all_effects <- function(data, snp_trait_map) {
  results <- list()
  
  for (i in seq_len(nrow(snp_trait_map))) {
    snp <- snp_trait_map$SNP[i]
    trait <- snp_trait_map$Trait[i]
    
    if (is.na(snp) || is.na(trait) || trait == "" || !(trait %in% colnames(data))) next
    message("Processing: ", snp, " → ", trait)
    
    snp_data <- data %>% filter(.data[[snp]] %in% c("0|0", "1|1"))
    
    effect_df <- snp_data %>%
      group_by(Env = Environment, Genotype = .data[[snp]]) %>%
      summarise(mean_trait = mean(.data[[trait]], na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Genotype, values_from = mean_trait)
    
    ref_val <- if ("0|0" %in% colnames(effect_df)) effect_df$`0|0` else NA_real_
    alt_val <- if ("1|1" %in% colnames(effect_df)) effect_df$`1|1` else NA_real_
    
    effect_df <- effect_df %>%
      mutate(
        SNP = snp,
        Trait = trait,
        REF = ref_val,
        ALT = alt_val,
        Effect_Size = ((ALT / REF) - 1) * 100
      ) %>%
      select(SNP, Trait, Env, REF, ALT, Effect_Size)
    
    results[[paste(snp, trait, sep = "_")]] <- effect_df
  }
  
  bind_rows(results)
}

all_effects <- compute_all_effects(merged_data, snp_trait_map)
fwrite(all_effects, "Final_Effect_Summary_SpecificSNPs.csv")





install.packages('R.utils')

fread("mlm_yield_results_earLength.y.csv.gz", nrows = 5)


library(data.table)
library(dplyr)
library(R.utils)

# Load effect data
effect_data <- fread("Final_Effect_Summary_SpecificSNPs.csv")

# Get unique Trait–SNP pairs
trait_snp_list <- unique(effect_data[, .(Trait, SNP)])

# Function to read p-values for a trait
read_trait_pvalues <- function(trait) {
  file <- paste0("mlm_yield_results_", trait, ".y.csv.gz")
  if (!file.exists(file)) return(NULL)
  
  gwas_data <- fread(file)
  gwas_data <- gwas_data %>%
    mutate(SNP = paste0(CHR, "_", POS)) %>%
    select(SNP, P.value = `P`) %>%
    mutate(Trait = trait)
  
  return(gwas_data)
}

# Run for all 7 traits
all_pvals <- bind_rows(lapply(unique(trait_snp_list$Trait), read_trait_pvalues))

# Merge with effect data
effect_annotated <- effect_data %>%
  left_join(all_pvals, by = c("SNP", "Trait")) %>%
  mutate(
    Significance = case_when(
      is.na(P.value) ~ "",
      P.value < 0.001 ~ "***",
      P.value < 0.01  ~ "**",
      P.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

fwrite(effect_annotated, "Effect_Summary_With_Pvalues.csv")







read_trait_pvalues <- function(trait) {
  file <- paste0("mlm_yield_results_", trait, ".y.csv.gz")
  if (!file.exists(file)) return(NULL)
  
  gwas_data <- fread(file)
  
  p_col <- paste0(trait, ".y.MLM")  # e.g., earLength.y.MLM
  
  if (!p_col %in% colnames(gwas_data)) {
    warning("Column ", p_col, " not found in ", file)
    return(NULL)
  }
  
  gwas_data <- gwas_data %>%
    select(SNP, P.value = all_of(p_col)) %>%
    mutate(Trait = trait)
  
  return(gwas_data)
}



all_pvals <- bind_rows(lapply(unique(trait_snp_list$Trait), read_trait_pvalues))

effect_annotated <- effect_data %>%
  left_join(all_pvals, by = c("SNP", "Trait")) %>%
  mutate(
    Significance = case_when(
      is.na(P.value) ~ "",
      P.value < 0.001 ~ "***",
      P.value < 0.01  ~ "**",
      P.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

fwrite(effect_annotated, "Effect_Summary_With_Pvalues.csv")






library(data.table)
library(dplyr)
library(ggplot2)

# ---- Load data ----
data <- fread("Effect_Summary_With_Pvalues.csv")

# ---- Prepare data ----
data_long <- data %>%
  mutate(
    Env_Name = case_when(
      grepl("Ames.*2022", Env) ~ "Ames 2022",
      grepl("Ames.*2023", Env) ~ "Ames 2023",
      grepl("Lincoln.*2022", Env) ~ "Lincoln 2022",
      grepl("Lincoln.*2023", Env) ~ "Lincoln 2023",
      grepl("Missouri.*2022", Env) ~ "Missouri Valley 2022",
      grepl("Missouri.*2023", Env) ~ "Missouri Valley 2023",
      grepl("Crawfordsville.*2022", Env) ~ "Crawfordsville 2022",
      grepl("Crawfordsville.*2023", Env) ~ "Crawfordsville 2023",
      TRUE ~ Env
    ),
    Significance = ifelse(Significance == "", NA, Significance)
  )

env_colors <- c(
  "Ames 2022"             = "#4477AA",  # Blue
  "Ames 2023"             = "#66CCEE",  # Sky Blue
  "Missouri Valley 2022" = "#228833",  # Green
  "Missouri Valley 2023" = "#CCBB44",  # Olive/Yellow
  "Crawfordsville 2022"  = "#EE6677",  # Coral
  "Crawfordsville 2023"  = "#AA3377",  # Plum
  "Lincoln 2022"         = "#BBBBBB",  # Light Gray
  "Lincoln 2023"         = "#000000"   # Black
)

# ---- Loop over each trait ----
for (trait_name in unique(data_long$Trait)) {
  trait_data <- data_long %>% filter(Trait == trait_name)
  
  p <- ggplot(trait_data, aes(x = Effect_Size, y = factor(SNP, levels = unique(SNP)))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    
    # Solid circles, filled by environment, no border
    geom_point(
      aes(color = Env_Name),
      shape = 19,
      size = 4
    ) +
    
    # Stars above significant points
    geom_text(
      data = trait_data %>% filter(!is.na(Significance)),
      aes(label = Significance),
      vjust = -1.2,
      fontface = "bold",
      size = 3.5,
      color = "black",
      show.legend = FALSE
    ) +
    
    # Trait name inside plot at top right
    annotate(
      "text",
      x = max(trait_data$Effect_Size, na.rm = TRUE),
      y = Inf,
      label = trait_name,
      hjust = 1.1,
      vjust = 1.5,
      fontface = "bold",
      size = 5
    ) +
    
    scale_color_manual(values = env_colors, name = "Environment") +
    
    labs(
      x = "Effect Size (% Change)",
      y = "SNP"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.line.y = element_line(color = "grey70"),
      legend.position = "bottom",
      legend.title = element_blank(),  # 👈 Remove "Environment" label
      legend.text = element_text(size = 11)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 5, shape = 19))
    )
  
  # Save the plot
  ggsave(
    filename = paste0("effect_size_", trait_name, "_clean_final.png"),
    plot = p,
    width = 12, height = 7, dpi = 300
  )
}
