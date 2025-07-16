# Load required libraries
library(terra)
library(dplyr)
library(stringr)

# Define study area shapefile
study_area <- vect("./base_files/btsl_90m_polygon.shp")

# Define scenarios and species
scenarios <- c("NCNC", "NC45", "NC85", "BAUNC", "BAU45", "BAU85") 
species_list <- c("BLBR", "MAAM", "PLCI", "RANA", "URAM")

# Initialize summary table
summary_table <- data.frame()

# Initialize a full combined table
percent_changes_all <- data.frame()

# Loop through scenarios and species
for (scenario in scenarios) {
  for (species in species_list) {
    
    # Define initial raster path (ts2010)
    init_path <- file.path("./Circuitscape/outputs/", scenario, paste0("curmap_", scenario, "_Resistance.", species, ".it1.ts2010.tif"))
    r_init <- crop(rast(init_path), study_area) %>% mask(study_area)
    init_sum <- global(r_init, "sum", na.rm = TRUE)[[1]]
    
    # Initialize vector to store changes for each iteration
    changes <- numeric()
    
    for (i in 1:40) {
      final_path <- file.path("./Circuitscape/outputs/", scenario, paste0("curmap_", scenario, "_Resistance.", species, ".it", i, ".ts2110.tif"))
      r_final <- crop(rast(final_path), study_area) %>% mask(study_area)
      final_sum <- global(r_final, "sum", na.rm = TRUE)[[1]]
      percent_change <- (final_sum - init_sum) / init_sum * 100
      changes[i] <- percent_change
      
      # Append to full combined table
      percent_changes_all <- rbind(
        percent_changes_all,
        data.frame(
          scenario = scenario,
          species = species,
          iteration = i,
          initial_sum = init_sum,
          final_sum = final_sum,
          percent_change = percent_change
        )
      )
    }
    
    # Compute summary stats
    summary_table <- rbind(
      summary_table,
      data.frame(
        scenario = scenario,
        species = species,
        mean_percent_change = mean(changes),
        se_percent_change = sd(changes) / sqrt(length(changes))
      )
    )
  }
}

# Save summary results
write.csv(summary_table, file = "./Circuitscape/outputs/percent_change_summary.csv", row.names = FALSE)

# Save full combined table
write.csv(percent_changes_all, file = "./habitat/percent_changes_all.csv", row.names = FALSE)
