

analyze_species_spread_improved <- function(species_data, 
                                           species_name, 
                                           apply_filtering = FALSE,
                                           eps_km = 50,
                                           time_resolution = c("5yr", "10yr", "3yr", "annual")) {
  
  time_resolution <- match.arg(time_resolution)
  
  cat(paste("Analyzing", species_name, "...\n"))
  
 
  if (apply_filtering && nrow(species_data) > 3) {
    coords <- st_coordinates(species_data)
    if (nrow(coords) >= 4) {
      # Use adaptive eps for filtering
      db_result <- dbscan::dbscan(coords, eps = eps_km * 1000, minPts = 3)
      keep_indices <- db_result$cluster > 0
      species_data <- species_data[keep_indices, ]
    }
  }
  
  if (nrow(species_data) == 0) {
    cat(paste("  -> No valid records found for", species_name, "after filtering.\n"))
    return(NULL)
  }
  

  species_data <- species_data %>%
    mutate(time_bin = case_when(
      time_resolution == "annual" ~ year,
      time_resolution == "3yr" ~ floor(year / 3) * 3,
      time_resolution == "5yr" ~ floor(year / 5) * 5,
      time_resolution == "10yr" ~ floor(year / 10) * 10
    ))
  
  time_bins_available <- sort(unique(species_data$time_bin))
  
  if (length(time_bins_available) < 2) {
    cat(paste("  -> Not enough time periods for", species_name, "\n"))
    return(NULL)
  }
  

  period_metrics_list <- list()
  
  for (period in time_bins_available) {
    period_data <- species_data %>% filter(time_bin == !!period)
    if (nrow(period_data) == 0) next
    
    # Calculate metrics
    range_metrics <- calculate_range_metrics(period_data)
    occupancy_metrics <- calculate_occupancy_metrics(period_data, grid_size_km = 10)
    
    # Get habitat from data
    habitat_val <- if ("habitat" %in% names(period_data)) {
      first(period_data$habitat)
    } else {
      NA
    }
    
    period_metrics <- data.frame(
      species = species_name,
      time_period = period,
      time_resolution = time_resolution,
      habitat = habitat_val,
      n_records = nrow(period_data),
      cluster_hull_area_km2 = range_metrics$cluster_hull_area,
      max_distance_km = range_metrics$max_distance,
      range_span_ew_deg = range_metrics$range_span_ew,
      range_span_ns_deg = range_metrics$range_span_ns,
      n_occupied_cells_10km = occupancy_metrics$n_occupied_cells,
      occupied_area_km2 = occupancy_metrics$total_area_occupied,
      centroid_lat = range_metrics$centroid_lat,
      centroid_lon = range_metrics$centroid_lon
    )
    
    period_metrics_list[[as.character(period)]] <- period_metrics
  }
  
  if (length(period_metrics_list) == 0) {
    cat(paste("  -> No valid metrics for", species_name, "\n"))
    return(NULL)
  }
  
  final_metrics_df <- bind_rows(period_metrics_list)
  

  if (nrow(final_metrics_df) >= 1) {
    tmp_df <- final_metrics_df %>% rename(year = time_period)
    tmp_df <- calculate_clustered_invasion_front(species_data, tmp_df, cluster_eps_km = eps_km)
    final_metrics_df <- tmp_df %>% rename(time_period = year)
  } else {
    final_metrics_df$invasion_front_km <- NA
  }
  
  # ---------------------------------------------------------------------
  # 5. CALCULATE TEMPORAL DYNAMICS
  # ---------------------------------------------------------------------
  if (nrow(final_metrics_df) >= 2) {
    # Rename for temporal metrics calculation
    tmp_df <- final_metrics_df %>% rename(year = time_period)
    tmp_df <- calculate_yearly_spread_metrics(tmp_df)
    tmp_df$dominant_direction <- calculate_overall_direction(tmp_df)
    final_metrics_df <- tmp_df %>% rename(time_period = year)
  } else {
    final_metrics_df$spread_rate_km_year <- NA
    final_metrics_df$centroid_shift_km_year <- NA
    final_metrics_df$area_change_km2_year <- NA
    final_metrics_df$new_cells_occupied <- NA
    final_metrics_df$dominant_direction <- NA
  }
  
  return(final_metrics_df)
}


calculate_clustered_invasion_front <- function(species_data, 
                                              yearly_metrics_df, 
                                              cluster_eps_km = 50) {
  
  if (nrow(yearly_metrics_df) < 1 || nrow(species_data) < 3) {
    yearly_metrics_df$invasion_front_km <- NA
    return(yearly_metrics_df)
  }
  
  # --- Step 1: Project data and identify persistent clusters ---
  # We use the time_bin column later, so let's add it to the projected data
  species_proj <- st_transform(species_data, 3035)
  species_proj$time_bin <- species_data$time_bin # Carry over the time_bin
  
  coords_mat <- st_coordinates(species_proj)
  
  db_result <- dbscan::dbscan(coords_mat, eps = cluster_eps_km * 1000, minPts = 3)
  species_proj$cluster <- db_result$cluster # Add cluster ID to projected data
  
  clusters <- unique(species_proj$cluster[species_proj$cluster > 0])
  
  if (length(clusters) == 0) {
    warning("DBSCAN did not identify any persistent clusters.")
    yearly_metrics_df$invasion_front_km <- NA
    return(yearly_metrics_df)
  }
  
  # --- Step 2: Calculate max cluster diameter for each time period ---
  # (This replaces the old "find origin" step)
  
  invasion_fronts <- numeric(nrow(yearly_metrics_df))
  
  for (i in 1:nrow(yearly_metrics_df)) {
    # 'year' here is actually the time_period (e.g., 1990, 1995)
    current_period <- yearly_metrics_df$year[i] 
    
    period_data_proj <- species_proj %>% filter(time_bin == current_period)
    
    active_clusters <- unique(period_data_proj$cluster[period_data_proj$cluster > 0])
    
    cluster_diameters <- c() 
    
    if (length(active_clusters) > 0) {
      for (cid in active_clusters) {
        
        cluster_period_data <- period_data_proj %>% filter(cluster == cid)
        
        # Need at least 2 points to calculate a pairwise distance
        if (nrow(cluster_period_data) >= 2) {
          
          # Get coordinates (which are already projected in meters)
          cluster_coords <- st_coordinates(cluster_period_data)
          
          # Calculate all pairwise distances and find the max
          dist_matrix <- dist(cluster_coords)
          max_dist_m <- max(dist_matrix, na.rm = TRUE)
          
          cluster_diameters <- c(cluster_diameters, max_dist_m / 1000) # Convert to km
          
        } else {
          cluster_diameters <- c(cluster_diameters, 0)
        }
      }
    }
    
    if (length(cluster_diameters) > 0) {
      invasion_fronts[i] <- max(cluster_diameters, na.rm = TRUE)
    } else {
      invasion_fronts[i] <- NA 
    }
  } 
  
  yearly_metrics_df$invasion_front_km <- invasion_fronts
  return(yearly_metrics_df)
}



calculate_range_metrics <- function(coords_sf) {
  if (nrow(coords_sf) == 0) {
    return(list(
      cluster_hull_area = NA,
      max_distance = NA,
      range_span_ew = NA,
      range_span_ns = NA,
      centroid_lat = NA,
      centroid_lon = NA
    ))
  }
  
  if (nrow(coords_sf) == 1) {
    coords_4326 <- st_transform(coords_sf, 4326)
    coords_mat <- st_coordinates(coords_4326)
    return(list(
      cluster_hull_area = 0,
      max_distance = 0,
      range_span_ew = 0,
      range_span_ns = 0,
      centroid_lat = coords_mat[1,2],
      centroid_lon = coords_mat[1,1]
    ))
  }
  
  # Project to metric CRS
  coords_proj <- tryCatch({
    st_transform(coords_sf, 3035)
  }, error = function(e) {
    warning("Projection failed, using original coordinates")
    coords_sf
  })
  
  coords_mat <- st_coordinates(coords_proj)[,1:2]
  
  hull_area <- NA
  clusters <- dbscan::dbscan(coords_mat, eps = 50000, minPts = 3)
  
  df_clustered <- coords_proj %>% mutate(cluster = clusters$cluster)
  cluster_ids <- unique(df_clustered$cluster[df_clustered$cluster > 0])
  
  if (length(cluster_ids) > 0) {
    cluster_hulls <- list()
    for (cid in cluster_ids) {
      cluster_points <- df_clustered %>% filter(cluster == cid)
      coords_cluster <- unique(round(st_coordinates(cluster_points)[,1:2], 1))
      if (nrow(coords_cluster) >= 3) {
        tryCatch({
          hull_coords <- concaveman::concaveman(coords_cluster,
                                               concavity = 5,
                                               length_threshold = 0)
          hull_poly <- st_polygon(list(hull_coords)) %>% 
            st_sfc(crs = st_crs(coords_proj))
          cluster_hulls[[length(cluster_hulls)+1]] <- hull_poly
        }, error = function(e) {
          NULL
        })
      }
    }
    if (length(cluster_hulls) > 0) {
      all_hulls <- do.call(c, cluster_hulls)
      cluster_hulls_union <- st_union(all_hulls)
      hull_area <- as.numeric(st_area(cluster_hulls_union)) / 1e6  # km²
    }
  }
  
  max_dist <- tryCatch({
    if (nrow(coords_mat) > 1) {
      dist_matrix <- dist(coords_mat)
      max(dist_matrix, na.rm = TRUE) / 1000
    } else {
      0
    }
  }, error = function(e) NA)
  
  # Range spans in degrees
  coords_4326 <- st_transform(coords_sf, 4326) %>% st_coordinates()
  range_ew <- max(coords_4326[,1], na.rm=TRUE) - min(coords_4326[,1], na.rm=TRUE)
  range_ns <- max(coords_4326[,2], na.rm=TRUE) - min(coords_4326[,2], na.rm=TRUE)
  
  # Centroid
  centroid_coords <- tryCatch({
    centroid <- st_centroid(st_union(coords_proj))
    centroid_lonlat <- st_transform(centroid, 4326)
    st_coordinates(centroid_lonlat)
  }, error = function(e) matrix(c(NA,NA), nrow=1))
  
  return(list(
    cluster_hull_area = hull_area,
    max_distance = max_dist,
    range_span_ew = range_ew,
    range_span_ns = range_ns,
    centroid_lat = centroid_coords[1,2],
    centroid_lon = centroid_coords[1,1]
  ))
}

calculate_occupancy_metrics <- function(coords_sf, grid_size_km = 10) {
  if (nrow(coords_sf) == 0) {
    return(list(
      n_occupied_cells = NA,
      total_area_occupied = NA
    ))
  }
  
  coords_proj <- tryCatch({
    st_transform(coords_sf, 3035)
  }, error = function(e) coords_sf)
  
  bbox <- tryCatch(st_bbox(coords_proj), error = function(e) NULL)
  if (is.null(bbox)) {
    return(list(n_occupied_cells = NA, total_area_occupied = NA))
  }
  
  grid_size_m <- grid_size_km * 1000
  bbox[1] <- bbox[1] - grid_size_m
  bbox[2] <- bbox[2] - grid_size_m
  bbox[3] <- bbox[3] + grid_size_m
  bbox[4] <- bbox[4] + grid_size_m
  
  n_occupied <- tryCatch({
    grid <- st_make_grid(st_as_sfc(bbox), cellsize = grid_size_m, square = TRUE)
    occupied <- st_intersects(grid, coords_proj)
    sum(lengths(occupied) > 0)
  }, error = function(e) NA)
  
  total_area <- ifelse(is.na(n_occupied), NA, n_occupied * (grid_size_km^2))
  
  return(list(
    n_occupied_cells = n_occupied,
    total_area_occupied = total_area
  ))
}



calculate_yearly_spread_metrics <- function(yearly_metrics_df) {
  if (nrow(yearly_metrics_df) < 2) {
    yearly_metrics_df$spread_rate_km_year <- NA
    yearly_metrics_df$centroid_shift_km_year <- NA
    yearly_metrics_df$area_change_km2_year <- NA
    yearly_metrics_df$new_cells_occupied <- NA
    return(yearly_metrics_df)
  }
  
  yearly_metrics_df <- yearly_metrics_df %>% arrange(year)
  n_years <- nrow(yearly_metrics_df)
  
  spread_rates <- rep(NA, n_years)
  centroid_shifts <- rep(NA, n_years)
  area_changes <- rep(NA, n_years)
  new_cells <- rep(NA, n_years)
  
  for (i in 2:n_years) {
    current_year <- yearly_metrics_df$year[i]
    prev_year <- yearly_metrics_df$year[i-1]
    year_diff <- current_year - prev_year
    
    # Centroid shift
    if (all(!is.na(yearly_metrics_df[i, c("centroid_lat","centroid_lon")])) &&
        all(!is.na(yearly_metrics_df[i-1, c("centroid_lat","centroid_lon")]))) {
      lat1 <- yearly_metrics_df$centroid_lat[i-1] * pi/180
      lat2 <- yearly_metrics_df$centroid_lat[i] * pi/180
      dlat <- (yearly_metrics_df$centroid_lat[i] - yearly_metrics_df$centroid_lat[i-1]) * pi/180
      dlon <- (yearly_metrics_df$centroid_lon[i] - yearly_metrics_df$centroid_lon[i-1]) * pi/180
      
      a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
      c <- 2 * atan2(sqrt(a), sqrt(1-a))
      centroid_dist <- 6371 * c
      
      centroid_shifts[i] <- centroid_dist / year_diff
    }
    
    # Range expansion rate
    if (!is.na(yearly_metrics_df$max_distance_km[i]) &&
        !is.na(yearly_metrics_df$max_distance_km[i-1])) {
      spread_rates[i] <- (yearly_metrics_df$max_distance_km[i] -
                           yearly_metrics_df$max_distance_km[i-1]) / year_diff
    }
    
    # Area change
    if (!is.na(yearly_metrics_df$cluster_hull_area_km2[i]) &&
        !is.na(yearly_metrics_df$cluster_hull_area_km2[i-1])) {
      area_changes[i] <- (yearly_metrics_df$cluster_hull_area_km2[i] -
                           yearly_metrics_df$cluster_hull_area_km2[i-1]) / year_diff
    }
    
    # New cells
    if (!is.na(yearly_metrics_df$n_occupied_cells_10km[i]) &&
        !is.na(yearly_metrics_df$n_occupied_cells_10km[i-1])) {
      new_cells[i] <- (yearly_metrics_df$n_occupied_cells_10km[i] -
                        yearly_metrics_df$n_occupied_cells_10km[i-1]) / year_diff
    }
  }
  
  yearly_metrics_df$spread_rate_km_year <- spread_rates
  yearly_metrics_df$centroid_shift_km_year <- centroid_shifts
  yearly_metrics_df$area_change_km2_year <- area_changes
  yearly_metrics_df$new_cells_occupied <- new_cells
  
  return(yearly_metrics_df)
}



calculate_overall_direction <- function(yearly_metrics_df) {
  if (nrow(yearly_metrics_df) < 2) return(NA)
  
  valid_centroids <- yearly_metrics_df[!is.na(yearly_metrics_df$centroid_lat) & 
                                        !is.na(yearly_metrics_df$centroid_lon), ]
  if (nrow(valid_centroids) < 2) return(NA)
  
  first_row <- valid_centroids[1, ]
  last_row <- valid_centroids[nrow(valid_centroids), ]
  
  total_lon_change <- last_row$centroid_lon - first_row$centroid_lon
  total_lat_change <- last_row$centroid_lat - first_row$centroid_lat
  
  if (abs(total_lon_change) > abs(total_lat_change)) {
    ifelse(total_lon_change > 0, "Eastward", "Westward")
  } else {
    ifelse(total_lat_change > 0, "Northward", "Southward")
  }
}


# =============================================================================
# HELPER FUNCTION: Clean and filter species
# =============================================================================

clean.filter.species <- function(data) {
  remove1 <- c("Pueraria montana", "Pycnonotus cafer", "Humulus scandens", 
               "Gambusia affinis", "Cipangopaludina chinensis", 
               "Prosopis juliflora", "Corvus splendens", 
               "Acridotheres cristatellus")
  
  data1 <- data[!data$species %in% remove1, ]
  
  # Add taxonomy if available
  if (file.exists("Database/Species.Taxonomy.xlsx")) {
    taxonomy <- read_xlsx("Database/Species.Taxonomy.xlsx")
    colnames(taxonomy)[1] <- "species"
    data1 <- data1 %>% 
      left_join(taxonomy[, c("species", "Group")], by = "species")
    
    # Manual fixes for missing groups
    data1$Group[data1$species == "Faxonius limosus"] <- "Crustaceans"
    data1$Group[data1$species == "Trachemys scripta"] <- "Reptiles"
    data1$Group[data1$species == "Vespa velutina"] <- "Insects"
    data1$Group[data1$species == "Cenchrus setaceus"] <- "Vascular plants"
    data1$Group[data1$species == "Pontederia crassipes"] <- "Vascular plants"
    data1$Group[data1$species == "Mustela vison"] <- "Mammals"
  }
  
  return(data1)
}


# =============================================================================
# LOGISTIC MODEL FUNCTIONS (unchanged)
# =============================================================================

logistic_model <- function(t, N0, K, r, tmin) {
  K / (1 + ((K - N0) / N0) * exp(-r * (t - tmin)))
}

fit_logistic_species <- function(df) {
  df <- df %>%
    arrange(decade) %>%
    mutate(cum_cells = cumsum(n_occupied_cells_10km))
  
  tmin <- min(df$decade)
  
  out <- tryCatch({
    fit <- nlsLM(
      cum_cells ~ logistic_model(decade, N0, K, r, tmin),
      data = df,
      start = list(N0 = max(1, df$cum_cells[1]),
                   K = max(df$cum_cells) * 1.5,
                   r = 0.2),
      lower = c(0, 0, 1e-6),
      upper = c(1e4, max(df$cum_cells) * 10, 5),
      control = nls.lm.control(maxiter = 500)
    )
    
    pars <- coef(fit)
    N0 <- pars["N0"]; K <- pars["K"]; r <- pars["r"]
    
    # Derived metrics
    t_inf <- tmin + (1/r) * log((K - N0)/N0)
    N_inf <- K/2
    
    t_sat <- tmin + (1/r) * log(((K - N0)/N0) * (1/0.01))
    N_sat <- 0.99 * K
    
    n_max <- (r * K) / (4 * N0)
    
    alpha_intro <- 0.01
    alpha_lag <- 0.1
    
    t_intro <- tmin - (1/r) * log((2 - alpha_intro + 2*sqrt(1 - alpha_intro)) / 
                                   (alpha_intro * ((K/N0) - 1)))
    t_lag <- tmin - (1/r) * log((2 - alpha_lag + 2*sqrt(1 - alpha_lag)) / 
                                 (alpha_lag * ((K/N0) - 1)))
    lag_duration <- t_lag - t_intro
    
    df_pred <- data.frame(decade = seq(min(df$decade), max(df$decade), 1))
    df_pred$predicted <- logistic_model(df_pred$decade, N0, K, r, tmin)
    
    p <- ggplot(df, aes(x = decade, y = cum_cells)) +
      geom_point(size = 3, color = "red") +
      geom_line(data = df_pred, aes(x = decade, y = predicted), 
                color = "blue", size = 1) +
      geom_vline(xintercept = t_inf, linetype = "dashed", color = "darkgreen") +
      geom_hline(yintercept = N_inf, linetype = "dotted", color = "darkgreen") +
      geom_vline(xintercept = t_sat, linetype = "dashed", color = "purple") +
      geom_hline(yintercept = N_sat, linetype = "dotted", color = "purple") +
      annotate("text", x = t_inf, y = N_inf, label = "Inflection", 
               vjust = -1, color = "darkgreen") +
      annotate("text", x = t_sat, y = N_sat, label = "Saturation", 
               vjust = -1, color = "purple") +
      labs(title = paste("Logistic fit:", unique(df$species)),
           x = "Year", y = "Cumulative occupied cells") +
      theme_bw()
    
    return(list(
      params = data.frame(
        species = unique(df$species),
        N0 = N0, K = K, r = r,
        t_inf = t_inf, N_inf = N_inf,
        t_sat = t_sat, N_sat = N_sat,
        n_max = n_max,
        t_intro = t_intro, t_lag = t_lag,
        lag_duration = lag_duration,
        converged = TRUE,
        stringsAsFactors = FALSE
      ),
      plot = p
    ))
  }, error = function(e) {
    cat("Error for species", unique(df$species), ":", e$message, "\n")
    list(
      params = data.frame(
        species = unique(df$species),
        N0 = NA, K = NA, r = NA,
        t_inf = NA, N_inf = NA,
        t_sat = NA, N_sat = NA,
        n_max = NA, t_intro = NA,
        t_lag = NA, lag_duration = NA,
        converged = FALSE,
        stringsAsFactors = FALSE
      ),
      plot = NULL
    )
  })
  return(out)
}

round_to <- function(x, base) base * round(x / base)
