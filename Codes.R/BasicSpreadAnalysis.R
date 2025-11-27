

rm(list = ls())
Sys.setLanguage(lang="en")

# Load packages
source("./R.Codes/02 Claude code.R")

pacman::p_load(ggplot2, sf, rgbif, raster, terra, sp,
               rnaturalearth, xls, xlsx, readxl, writexl,
               dplyr, ggspatial, purrr, units, dbscan,
               geosphere, concaveman, patchwork, tidyr)

# Load functions
#source("./R.Codes/functions.R")



cat("\n=== STEP 1: Loading and cleaning data ===\n")

# Load raw data
old.occ <- readRDS("./Database/cleaned_occ1.rds")
df_aquatic <- readRDS("./Database/cleaned.occ.dispersal.50.rds")
terr <- readRDS("./Database/cleaned.occ.terrestrial.dispersal.rds")

# Clean old.occ
old.occ <- old.occ %>%
  filter(!is.na(decimalLatitude),
         !is.na(decimalLongitude),
         !is.na(year)) %>%  
  filter(year >= 1950) %>%
  filter(decimalLongitude > -10) %>%
  filter(is.na(coordinateUncertaintyInMeters) | 
           coordinateUncertaintyInMeters <= 5000) %>%
  filter(!(decimalLongitude > 1 & decimalLongitude < 4.5 &
             decimalLatitude > 38.5 & decimalLatitude < 40.5)) %>%
  filter(decimalLatitude >= -90, decimalLatitude <= 90,
         decimalLongitude >= -180, decimalLongitude <= 180) %>%
  filter(!(decimalLatitude >= 35.85 & decimalLatitude <= 35.92 &
             decimalLongitude >= -5.35 & decimalLongitude <= -5.28),
         !(decimalLatitude >= 35.26 & decimalLatitude <= 35.32 &
             decimalLongitude >= -2.96 & decimalLongitude <= -2.90))

# Remove problematic species
remove_species <- c("Lampropeltis getula", "Tamias sibiricus", "Nyctereutes procyonoides",
                    "Pueraria montana", "Pycnonotus cafer", "Humulus scandens", 
                    "Gambusia affinis", "Asclepias syriaca", "Gunnera tinctoria",
                    "Procambarus virginalis", "Lithobates catesbeianus", "Elodea nuttallii",
                    "Cipangopaludina chinensis", "Prosopis juliflora", "Corvus splendens",
                    "Acridotheres cristatellus")

df_aquatic <- df_aquatic %>% filter(!species %in% remove_species)
terr <- terr %>% filter(!species %in% remove_species)

df_aquatic$habitat <- "Aquatic-Mix"
terr$habitat <- "Terrestrial"

df_aquatic <- st_as_sf(df_aquatic, crs = 3857)
df_aquatic <- st_transform(df_aquatic, 4326)

terr <- st_as_sf(terr, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)


colnames(terr)[colnames(terr) == "dispersal_type"] <- "expansion_type"
terr$expansion_type <- as.character(terr$expansion_type)

final.df <- bind_rows(df_aquatic, terr)

essential_cols <- intersect(c("species", "year", "habitat", "expansion_type", "n_records", "geometry"),
                            names(final.df))
final.df <- final.df[, essential_cols]

final.df <- final.df %>%
  mutate(expansion_type = case_when(
    expansion_type %in% c("origin", "natural_expansion", "primary_diffusion") ~ "Natural",
    expansion_type %in% c("jump_introduction", "blocked_by_dam", "reservoir_mediated", "jump") ~ "Jump",
    TRUE ~ "Unknown"
  ))

cat(sprintf("Combined dataset: %d species\n", n_distinct(final.df$species)))  # 51 species here

# STEP 2: FIX DATA QUALITY ISSUES

cat("\n=== STEP 2: Fixing data quality issues ===\n")

# A. Fix n_records (if all NA, set to 1)
if (all(is.na(final.df$n_records))) {
  final.df$n_records <- 1
  cat("Set n_records = 1 for all rows (assuming 1 record per row)\n")
}

# B. Create temporal bins
final.df <- final.df %>%
  mutate(
    decade = floor(year / 10) * 10,
    period_5yr = floor(year / 5) * 5,
    period_3yr = floor(year / 3) * 3
  )

# C. Add time since introduction
final.df <- final.df %>%
  group_by(species) %>%
  mutate(
    intro_year = min(year, na.rm = TRUE),
    years_since_intro = year - intro_year,
    establishment_phase = case_when(
      years_since_intro < 10 ~ "Early",
      years_since_intro < 30 ~ "Expansion",
      TRUE ~ "Established"
    )
  ) %>%
  ungroup()

# D. Add sampling effort metrics
final.df <- final.df %>%
  group_by(species, year) %>%
  mutate(
    records_this_year = n(),
    sampling_weight = 1 / records_this_year
  ) %>%
  ungroup()

# E. Project to metric CRS for analysis
df_proj <- st_transform(final.df, crs = 3857)

cat(sprintf("Final dataset: %d records, %d species\n", nrow(df_proj), n_distinct(df_proj$species))) # same 51 species :)


# STEP 3: DATA VALIDATION AND SPECIES FILTERING


cat("\n=== STEP 3: Validating species data quality ===\n")
# It seems we have some NA for some year * species

validation_checks <- df_proj %>%
  st_drop_geometry() %>% filter(year !="") %>% 
  group_by(species) %>%
  summarise(
    n_total = n(),
    n_years = n_distinct(year),
    year_span = max(year) - min(year),
    first_year = min(year),
    last_year = max(year),
    n_natural = sum(expansion_type == "Natural"),
    n_jump = sum(expansion_type == "Jump"),
    pct_natural = round(mean(expansion_type == "Natural") * 100, 1),
    pct_jump = round(mean(expansion_type == "Jump") * 100, 1),
    habitat = first(habitat),
    .groups = "drop"
  ) %>%
  mutate(
    has_enough_data = n_years >= 5 & year_span >= 5 & n_total >= 20, # This is super important
    has_natural = n_natural >= 10,
    has_jumps = n_jump >= 5
  )

print(validation_checks)

# Export validation table
write_xlsx(validation_checks, "species_validation_checks.xlsx")

# Filter species
species_to_analyze <- validation_checks %>%
  filter(has_enough_data) %>%
  pull(species) 

cat(sprintf("\nSpecies passing quality checks: %d out of %d\n", 
            length(species_to_analyze), nrow(validation_checks)))  # now we have 41/51

# Filter data
df_proj <- df_proj %>% filter(species %in% species_to_analyze)


# STEP 4: CALCULATE ADAPTIVE DBSCAN PARAMETERS

cat("\n=== STEP 4: Calculating adaptive clustering parameters ===\n")

calculate_adaptive_eps <- function(species_data) {
  coords <- st_coordinates(species_data)
  
  if (nrow(coords) >= 10) {
    # Calculate k-nearest neighbor distances
    knn_dist <- dbscan::kNNdist(coords, k = 5)
    
    # Use 90th percentile of kNN distances
    eps_meters <- quantile(knn_dist, 0.90)
    eps_km <- eps_meters / 1000
    
    # Constrain between 25-100 km
    eps_km <- max(min(eps_km, 100), 25)
  } else {
    eps_km <- 50  # Default
  }
  
  return(eps_km)
}

# Calculate eps for each species
species_eps <- data.frame()
for (sp in species_to_analyze) {
  sp_data <- df_proj %>% filter(species == sp)
  eps_km <- calculate_adaptive_eps(sp_data)
  species_eps <- rbind(species_eps, data.frame(species = sp, eps_km = eps_km))
}

print(species_eps)
write_xlsx(species_eps, "species_adaptive_eps.xlsx")

# =============================================================================
# STEP 5: IDENTIFY AND FLAG SPATIAL OUTLIERS
# =============================================================================

cat("\n=== STEP 5: Flagging spatial outliers ===\n")

df_proj$is_spatial_outlier <- FALSE

for (sp in species_to_analyze) {
  sp_data <- df_proj %>% filter(species == sp)
  coords <- st_coordinates(sp_data)
  
  if (nrow(coords) > 10) {
    eps_km <- species_eps$eps_km[species_eps$species == sp]
    db <- dbscan::dbscan(coords, eps = eps_km * 1000, minPts = 5)
    df_proj$is_spatial_outlier[df_proj$species == sp] <- (db$cluster == 0)
  }
}

outlier_summary <- df_proj %>%
  st_drop_geometry() %>%
  group_by(species) %>%
  summarise(
    n_outliers = sum(is_spatial_outlier),
    pct_outliers = round(mean(is_spatial_outlier) * 100, 1),
    .groups = "drop"
  )

print(outlier_summary)

#let me take a close look - i am not sure: 
s <- species_to_analyze[10] 
plot_data <- df_proj %>%   filter(species == s)
ggplot() +
  geom_sf(data = plot_data, 
          aes(color = is_spatial_outlier), 
          size = 2.5, alpha = 0.7) +
  scale_color_manual(
    name = "Style",
    values = c("FALSE" = "blue", "TRUE" = "red"),
    labels = c("FALSE" = "Cluster Point", "TRUE" = "Outlier")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# STEP 6: DEFINE EXPANSION GROUPS FOR ANALYSIS

cat("\n=== STEP 6: Creating analysis groups ===\n")
unique(df_proj$species) # 41 sp

expansion_groups <- list(
  Natural_only = df_proj %>% filter(expansion_type == "Natural"),
  Jump_only = df_proj %>% filter(expansion_type == "Jump"),
  All_combined = df_proj %>% filter(expansion_type %in% c("Natural", "Jump"))
)

for (grp in names(expansion_groups)) {
  cat(sprintf("%s: %d records, %d species\n", 
              grp, 
              nrow(expansion_groups[[grp]]),
              n_distinct(expansion_groups[[grp]]$species)))
}

# =============================================================================
# STEP 7: RUN SPREAD ANALYSIS FOR EACH GROUP
# =============================================================================

cat("\n=== STEP 7: Running spread analysis ===\n")

results_by_group <- list()

for (grp in names(expansion_groups)) {
  cat(sprintf("\n>>> Processing group: %s <<<\n", grp))
  df_subset <- expansion_groups[[grp]]
  
  species_metrics_grp <- data.frame()
  
  for (sp in species_to_analyze) {
    species_data <- df_subset %>% filter(species == sp)
    
    if (nrow(species_data) < 5) {
      cat(sprintf("  Skipping %s - insufficient data (n = %d)\n", sp, nrow(species_data)))
      next
    }
    
    # Get adaptive eps for this species
    eps_km <- species_eps$eps_km[species_eps$species == sp]
    
    # Run analysis
    species_metrics <- tryCatch({
      analyze_species_spread_improved(
        species_data = species_data,
        species_name = sp,
        apply_filtering = TRUE,
        eps_km = eps_km,
        time_resolution = "5yr"  # Use 5-year periods
      )
    }, error = function(e) {
      cat(sprintf("  ERROR for %s: %s\n", sp, e$message))
      return(NULL)
    })
    
    if (!is.null(species_metrics)) {
      species_metrics$analysis_group <- grp
      species_metrics$eps_km_used <- eps_km
      species_metrics_grp <- rbind(species_metrics_grp, species_metrics)
      cat(sprintf("  ✓ Added %d time periods for %s\n", nrow(species_metrics), sp))
    } else {
      cat(sprintf("  ✗ No valid metrics generated for %s\n", sp))
    }
  }
  
  results_by_group[[grp]] <- species_metrics_grp
}
