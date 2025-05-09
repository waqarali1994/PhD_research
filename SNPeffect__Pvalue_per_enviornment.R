
# Set working directory
setwd("~/Documents/PhD work/Yield_project/SNP effect/SNP_effect_v2/")
# Load required libraries
library(dplyr)
library(data.table)

# Load phenotype data
data <- fread('HIPS_INBREDS_2022_2023_v12 (1).csv')

# Columns to keep for analysis
columns_to_keep <- c("genotype", "year", "row", "range", "location", "nitrogenTreatment",
                     "kernelRowNumber", "hundredKernelMass", "kernelMassPerEar",
                     "earWeight", "earWidth", "earLength", "kernelsPerRow")

# Filter data for "Medium" nitrogen treatment and select relevant columns
data_filtered <- data %>%
  filter(nitrogenTreatment == "Medium") %>%
  select(all_of(columns_to_keep))

# Save filtered data to a CSV file
fwrite(data_filtered, "filtered_data.csv")

# Plot histograms for each trait
hist(data_filtered$kernelRowNumber, breaks = 40, main = "Kernel Row Number Distribution", xlab = "Kernel Row Number")
hist(data_filtered$kernelsPerRow, breaks = 40, main = "Kernels Per Row Distribution", xlab = "Kernels Per Row")
hist(data_filtered$earLength, breaks = 40, main = "Ear Length Distribution", xlab = "Ear Length")
hist(data_filtered$earWidth, breaks = 40, main = "Ear Width Distribution", xlab = "Ear Width")
hist(data_filtered$hundredKernelMass, breaks = 40, main = "Hundred Kernel Mass Distribution", xlab = "Hundred Kernel Mass")
hist(data_filtered$kernelMassPerEar, breaks = 40, main = "Kernel Mass Per Ear Distribution", xlab = "Kernel Mass Per Ear")
hist(data_filtered$earWeight, breaks = 40, main = "Ear Weight Distribution", xlab = "Ear Weight")

# Boxplot statistics for 'earLength'
boxplot.stats(data_filtered$earLength)

# Define the cutoff values for each trait
cutoffs <- list(
  earWeight = c(NA, 180),
  kernelMassPerEar = c(-1.8, 1.2),
  earLength = c(5, 21),
  earWidth = c(NA, 5.2),
  hundredKernelMass = c(NA, 40)
)

# Function to apply cutoff filters to the data
apply_cutoffs <- function(data, cutoffs) {
  filtered_data <- data
  
  for (trait in names(cutoffs)) {
    if (trait %in% colnames(filtered_data)) {  # Ensure trait exists in the dataset
      lower <- cutoffs[[trait]][1]
      upper <- cutoffs[[trait]][2]
      
      # Apply lower and upper cutoff values
      if (!is.na(lower) & !is.na(upper)) {
        filtered_data[[trait]] <- ifelse(filtered_data[[trait]] >= lower & filtered_data[[trait]] <= upper, filtered_data[[trait]], NA)
      } else if (!is.na(lower)) {
        filtered_data[[trait]] <- ifelse(filtered_data[[trait]] >= lower, filtered_data[[trait]], NA)
      } else if (!is.na(upper)) {
        filtered_data[[trait]] <- ifelse(filtered_data[[trait]] <= upper, filtered_data[[trait]], NA)
      }
    } else {
      warning(paste("Trait", trait, "not found in data. Skipping."))
    }
  }
  
  return(filtered_data)
}

# Apply the cutoff filter three times on data
filtered_data_1 <- apply_cutoffs(data_filtered, cutoffs)

# Save the filtered data to a CSV file
write.csv(filtered_data_1, "Hips_extremev1.csv", row.names = FALSE)

# Helper function to generate alphabetic identifiers
generate_alphabets <- function(n) {
  if (n <= 26) {
    return(LETTERS[1:n])
  } else {
    return(sapply(1:n, function(i) {
      paste0(LETTERS[((i - 1) %% 26) + 1], 
             ifelse(i > 26, ceiling(i / 26), ""))
    }))
  }
}

# View the first few rows of the final filtered data
print(head(filtered_data_1))

cat("Filtered data saved to 'Hips_extremev1.csv'.\n")





# Function to compute spatial corrections for multiple traits across environments
getSpatialCorrectionsEnvironment <- function(grouped_dataset, phenotypes) {
  # Initialize an empty dataframe to store results
  df.sp <- tibble(group = NULL, plotNumber = NULL, trait = NULL, value = NULL, genotype = NULL)
  
  # Get the unique groups
  groups <- unique(grouped_dataset$group)
  
  # Loop over each trait
  for (trait in phenotypes) {
    # Loop over each group
    for (currGroup in groups) {
      # Filter data for the current group and trait, ensuring no NA values
      group.df <- filter(grouped_dataset, 
                         group == currGroup & 
                           !is.na(row) & !is.na(range) & !is.na(.data[[trait]]))
      
      # Skip if no valid data for the current group and trait
      if (nrow(group.df) == 0) {
        next
      }
      
      # Calculate knots for spatial modeling, ensuring there are enough data points
      rangeKnots <- floor(max(group.df$range, na.rm = TRUE) / 2) + 1
      rowKnots <- floor(max(group.df$row, na.rm = TRUE) / 2) + 1
      
      # Check if enough data points exist for model fitting (at least two)
      if (rangeKnots < 2 | rowKnots < 2) {
        warning(paste("Not enough data for group:", currGroup, "trait:", trait))
        next
      }
      
      # Fit the SpATS model
      model <- tryCatch({
        SpATS(
          response = trait, 
          genotype = "plotNumber", 
          genotype.as.random = TRUE, 
          spatial = ~ SAP(range, row, nseg = c(rangeKnots, rowKnots)), 
          data = group.df
        )
      }, error = function(e) {
        warning(paste("Error fitting model for group:", currGroup, "trait:", trait))
        return(NULL)
      })
      
      # Skip if model fitting failed
      if (is.null(model)) {
        next
      }
      
      # Extract intercept and BLUPs (Best Linear Unbiased Predictions)
      intercept <- model$coeff["Intercept"]
      sp <- as_tibble(model$coeff, rownames = "plotNumber") %>%
        mutate(
          group = currGroup,  # Add the group column
          trait = trait,  # Add the trait name
          plotNumber = as.character(plotNumber),  # Keep plotNumber as character
          value = value + intercept  # Add intercept to coefficients
        ) %>%
        left_join(group.df %>% select(plotNumber, genotype) %>% distinct(), by = "plotNumber") %>%
        select(group, plotNumber, genotype, trait, value)
      
      # Append results to the master dataframe
      df.sp <- bind_rows(df.sp, sp)
    }
  }
  
  # Return the final dataframe with spatially corrected values (BLUPs)
  return(df.sp)
}

# Define the traits to correct
phenotypes <- c("earLength", "earWidth", "kernelRowNumber", "kernelMassPerEar", "hundredKernelMass", "earWeight", "kernelsPerRow")

# Run the spatial correction function
spatiallyCorrectedBLUPs <- getSpatialCorrectionsEnvironment(
  grouped_dataset = grouped_dataset,  # Correct argument passed
  phenotypes = phenotypes
)

# Check the corrected results
print(head(spatiallyCorrectedBLUPs))

# Merge the corrected BLUPs with the original dataset
final_dataset <- grouped_dataset %>%
  left_join(
    spatiallyCorrectedBLUPs %>%
      pivot_wider(names_from = trait, values_from = value), 
    by = c("group", "plotNumber", "genotype")  # Match correct column names
  )

# Save the final dataset to a CSV file
write.csv(final_dataset, "spatially_corrected_final_dataset.csv", row.names = FALSE)

# Check the final dataset
print(head(final_dataset))




# Load required libraries
library(data.table)
library(dplyr)

# Load the spatially corrected dataset
df <- fread("spatially_corrected_final_dataset.csv")

# List of traits to analyze (ensure these match your column names)
traits <- c('earLength.y', 'earWidth.y', 'kernelsPerRow.y', 'kernelRowNumber.y', 
            'hundredKernelMass.y', 'earWeight.y', 'kernelMassPerEar.y')

# Get unique environments (groups)
groups <- unique(df$group)

# Loop through each environment
for (currGroup in groups) {
  
  # Subset data for the current group
  group_df <- df[df$group == currGroup, ]
  
  # Initialize a dataframe to collect BLUEs for all traits in this environment
  blues_all <- data.frame()
  
  for (trait in traits) {
    
    # Check if the trait exists and has enough non-missing values
    if (!(trait %in% colnames(group_df))) {
      warning(paste("Trait", trait, "not found in data. Skipping."))
      next
    }
    if (sum(!is.na(group_df[[trait]])) < 2) {
      warning(paste("Trait", trait, "in", currGroup, "has too few values. Skipping."))
      next
    }
    
    # Fit a linear model: genotype as fixed effect per environment
    model <- lm(
      as.formula(paste(trait, "~ genotype")),
      data = group_df
    )
    
    # Extract coefficients (including intercept)
    blues <- coef(model)
    
    # Separate intercept and genotype coefficients
    intercept <- blues["(Intercept)"]
    genotype_blues <- blues[names(blues) != "(Intercept)"]
    
    # Add intercept to get genotype BLUEs
    blues_with_intercept <- genotype_blues + intercept
    
    # Create a data.frame for this trait
    blues_df <- data.frame(
      genotype = names(blues_with_intercept),
      BLUE = as.vector(blues_with_intercept)
    )
    
    # Rename the BLUE column to reflect just the trait (no environment)
    colnames(blues_df)[2] <- trait
    
    # Merge into the full result
    if (nrow(blues_all) == 0) {
      blues_all <- blues_df
    } else {
      blues_all <- merge(blues_all, blues_df, by = "genotype", all = TRUE)
    }
  }
  
  # Save results for current environment
  output_file <- paste0(currGroup, "_blues.csv")
  write.csv(blues_all, output_file, row.names = FALSE)
  print(paste("✅ BLUEs saved:", output_file))
}






      
# Set working directory
setwd('~/Documents/PhD work/Yield_project/SNP effect/SNP_effect_v2/')
      
      # Load necessary libraries
      library(data.table)
      library(ggplot2)
      
      # Traits to process
      traits <- c('earLength.y', 'earWidth.y', 'kernelsPerRow.y', 'kernelRowNumber.y', 
                  'hundredKernelMass.y', 'earWeight.y', 'kernelMassPerEar.y')
      
      # Create folders if they don't exist
      if (!dir.exists("trait_histograms")) dir.create("trait_histograms")
      if (!dir.exists("cleaned_data")) dir.create("cleaned_data")
      
      # Define cutoffs per environment
      cutoffs_by_environment <- list(
        "2022_Ames_Medium" = list(
          earLength = c(10.5, 14.2), earWidth = c(3.4, 4.1), kernelRowNumber = c(12.5, 16.2),
          hundredKernelMass = c(21, 27.7), earWeight = c(NA, 100), kernelMassPerEar = c(NA, 72)
        ),
        "2022_Crawfordsville_Medium" = list(
          earLength = c(9.8, 14), earWidth = c(3.25, 4.14), kernelRowNumber = c(12, 16.5),
          hundredKernelMass = c(NA, 28), earWeight = c(NA, 98), kernelMassPerEar = c(33, 76)
        ),
        "2023_Ames_Medium" = list(
          earLength = c(11, 15), earWidth = c(3.7, 4.45), kernelRowNumber = c(NA, 16.5),
          hundredKernelMass = c(NA, 30), earWeight = c(45, 230), kernelMassPerEar = c(NA, 95)
        ),
        "2023_Crawfordsville_Medium" = list(
          earLength = c(11, 15.5), earWidth = c(3.3, 4.3), kernelRowNumber = c(NA, 16),
          hundredKernelMass = c(19, 29), earWeight = c(50, 235), kernelMassPerEar = c(38, 84)
        ),
        "2022_Lincoln_Medium" = list(
          earLength = c(9.2, 13.5), earWidth = c(2.7, 3.5), kernelsPerRow = c(NA, 24),
          kernelRowNumber = c(10.5, 14.2), hundredKernelMass = c(13, 21),
          earWeight = c(22, 57), kernelMassPerEar = c(NA, 47)
        ),
        "2023_Lincoln_Medium" = list(
          earLength = c(9.5, 13.2), earWidth = c(2.9, 3.8), kernelsPerRow = c(11, 24),
          kernelRowNumber = c(11, 14.4), hundredKernelMass = c(15.5, 24.5),
          earWeight = c(25, 66), kernelMassPerEar = c(12, 55)
        ),
        "2022_Missouri_Valley_Medium" = list(
          earLength = c(9.8, NA), earWidth = c(3.1, 3.9), kernelsPerRow = c(12, 27),
          kernelRowNumber = c(11.8, 15.5), hundredKernelMass = c(15, 24),
          earWeight = c(NA, 90), kernelMassPerEar = c(21, 70)
        ),
        "2023_Missouri_Valley_Medium" = list(
          earLength = c(11.8, 16.5), earWidth = c(3.5, 4.2), kernelsPerRow = c(17, 30),
          kernelRowNumber = c(11.8, 16.7), hundredKernelMass = c(18, 26.5),
          earWeight = c(NA, 120), kernelMassPerEar = c(NA, 99)
        )
      )
      
      
      # Loop through each _blues.csv file
      for (file in list.files(pattern = "_blues.csv$")) {
        
        group_df <- fread(file)
        
        # Extract standardized environment name
        currGroup <- gsub("_blues.csv$", "", file)
        currGroup <- gsub(" ", "_", currGroup)  # convert spaces to underscores to match cutoff keys
        
        # Print current environment
        message("Processing: ", currGroup)
        
        # Get cutoffs for this environment
        cutoffs <- cutoffs_by_environment[[currGroup]]
        if (is.null(cutoffs)) {
          message("⚠️ No cutoffs defined for: ", currGroup)
          next
        }
        
        for (trait in traits) {
          
          if (!(trait %in% colnames(group_df))) next
          if (sum(!is.na(group_df[[trait]])) < 2) next
          
          base_trait <- sub(".y", "", trait)
          if (!(base_trait %in% names(cutoffs))) {
            message("Skipping trait (no cutoff): ", base_trait, " in ", currGroup)
            next
          }
          
          # Get cutoffs
          lower_cutoff <- cutoffs[[base_trait]][1]
          upper_cutoff <- cutoffs[[base_trait]][2]
          
          # Plot histogram BEFORE cutoff
          breaks_before <- pretty(range(group_df[[trait]], na.rm = TRUE), n = 40)
          hist_before <- ggplot(group_df, aes(x = .data[[trait]])) +
            geom_histogram(breaks = breaks_before, fill = "blue", color = "black", alpha = 0.7) +
            ggtitle(paste("Histogram of", trait, "Before Cutoff -", currGroup)) +
            theme_classic()
          
          ggsave(paste0("trait_histograms/", currGroup, "_", trait, "_before_cutoff.png"), plot = hist_before)
          
          # ✂️ Apply cutoff
          group_df[[trait]] <- ifelse(
            (!is.na(lower_cutoff) & group_df[[trait]] < lower_cutoff) |
              (!is.na(upper_cutoff) & group_df[[trait]] > upper_cutoff),
            NA, group_df[[trait]]
          )
          
          # Plot histogram AFTER cutoff
          breaks_after <- pretty(range(group_df[[trait]], na.rm = TRUE), n = 40)
          hist_after <- ggplot(group_df, aes(x = .data[[trait]])) +
            geom_histogram(breaks = breaks_after, fill = "green", color = "black", alpha = 0.7) +
            ggtitle(paste("Histogram of", trait, "After Cutoff -", currGroup)) +
            theme_classic()
          
          ggsave(paste0("trait_histograms/", currGroup, "_", trait, "_after_cutoff.png"), plot = hist_after)
        }
        
        # Save cleaned data
        cleaned_file_name <- paste0("cleaned_data/", gsub("[^A-Za-z0-9]", "_", currGroup), "_cleaned.csv")
        write.csv(group_df, cleaned_file_name, row.names = FALSE)
        print(paste("Cleaned data saved:", cleaned_file_name))
      }
      
      
      
      
      
      
      
      # Set to your cleaned data directory
      setwd('~/Documents/PhD work/Yield_project/SNP effect/SNP_effect_v2/cleaned_data/')
      
      # List all *_cleaned.csv files
      files <- list.files(pattern = "_cleaned.csv$")
      
      # Loop through each file
      for (file in files) {
        
        # Read the file
        df <- fread(file)
        
        # Only process if 'genotype' column exists
        if (!("genotype" %in% colnames(df))) {
          message("⚠️ No 'genotype' column in ", file)
          next
        }
        
        # Remove 'genotype' prefix (e.g., "genotype2369" → "2369")
        df$genotype <- gsub("^genotype", "", df$genotype)
        
        # Optionally: replace with "G" prefix (e.g., "genotype2369" → "G2369")
        # df$genotype <- gsub("^genotype", "G", df$genotype)
        
        # Overwrite the original file
        write.csv(df, file, row.names = FALSE)
        message("enotype column cleaned in: ", file)
      }
      
      
      
      
      
      
      
df1 <- fread('filtered_phenotype_data.csv')
df2 <- fread('2022_Ames_Medium_cleaned.csv')      

common_ids <- intersect(df1$genotype, df2$genotype)
length(common_ids)  # how many are common
