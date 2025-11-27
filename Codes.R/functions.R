# =======================================================================
# FUNCTIONS for SPREAD ANALYSIS (DBSCAN + concaveman version)
# =======================================================================

library(sf)
library(dplyr)
library(dbscan)
library(concaveman)
library(geosphere)
library(units)

# -----------------------------------------------------------------------
# 1. GEOGRAPHIC RANGE METRICS
# -----------------------------------------------------------------------
calculate_range_metrics <- function(coords_sf) {
  # Handle empty or single point data
  if (nrow(coords_sf) == 0) {
    return(list(
      cluster_hull_area = NA,
      max_distance      = NA,
      range_span_ew     = NA,
      range_span_ns     = NA,
      centroid_lat      = NA,
      centroid_lon      = NA
    ))
  }
  
  if (nrow(coords_sf) == 1) {
    coords_4326 <- st_coordinates(coords_sf)
    return(list(
      cluster_hull_area = 0,
      max_distance      = 0,
      range_span_ew     = 0,
      range_span_ns     = 0,
      centroid_lat      = coords_4326[1,2],
      centroid_lon      = coords_4326[1,1]
    ))
  }
  
  # Project to metric CRS for geometry calculations
  coords_proj <- tryCatch({
    st_transform(coords_sf, 3035)   # ETRS89 / LAEA Europe
  }, error = function(e) {
    warning("Projection failed, using original coordinates")
    coords_sf
  })
  
  coords_mat <- st_coordinates(coords_proj)[,1:2]
  
  # === DBSCAN + concaveman cluster hull area ===
  hull_area <- NA
  clusters <- dbscan::dbscan(coords_mat, eps = 50000, minPts = 3)  # eps = 50 km
  
  df_clustered <- coords_proj %>% mutate(cluster = clusters$cluster)
  cluster_ids  <- unique(df_clustered$cluster[df_clustered$cluster > 0])
  
  if (length(cluster_ids) > 0) {
    cluster_hulls <- list()
    for (cid in cluster_ids) {
      cluster_points <- df_clustered %>% filter(cluster == cid)
      coords_cluster <- unique(round(st_coordinates(cluster_points)[,1:2], 1))
      if (nrow(coords_cluster) >= 3) {
        hull_coords <- concaveman::concaveman(coords_cluster,
                                              concavity = 5,
                                              length_threshold = 0)
        hull_poly <- st_polygon(list(hull_coords)) |> st_sfc(crs = st_crs(coords_proj))
        cluster_hulls[[length(cluster_hulls)+1]] <- hull_poly
      }
    }
    if (length(cluster_hulls) > 0) {
      all_hulls <- do.call(c, cluster_hulls)
      cluster_hulls_union <- st_union(all_hulls)
      hull_area <- as.numeric(st_area(cluster_hulls_union)) / 1e6  # km²
    }
  }
  
  # Max pairwise distance (km)
  max_dist <- tryCatch({
    dist_matrix <- dist(coords_mat)
    max(dist_matrix, na.rm = TRUE) / 1000
  }, error = function(e) NA)
  
  # Range spans in degrees (lon/lat must be taken from unprojected coords)
  coords_4326 <- st_transform(coords_sf, 4326) |> st_coordinates()
  range_ew <- max(coords_4326[,1], na.rm=TRUE) - min(coords_4326[,1], na.rm=TRUE)
  range_ns <- max(coords_4326[,2], na.rm=TRUE) - min(coords_4326[,2], na.rm=TRUE)
  
  # Centroid (in projected space, then reproject back to lon/lat)
  centroid_coords <- tryCatch({
    centroid <- st_centroid(st_union(coords_proj))
    centroid_lonlat <- st_transform(centroid, 4326)
    st_coordinates(centroid_lonlat)
  }, error = function(e) matrix(c(NA,NA), nrow=1))
  
  return(list(
    cluster_hull_area = hull_area,
    max_distance      = max_dist,
    range_span_ew     = range_ew,
    range_span_ns     = range_ns,
    centroid_lat      = centroid_coords[1,2],
    centroid_lon      = centroid_coords[1,1]
  ))
}


# -----------------------------------------------------------------------
# 2. OCCUPANCY METRICS
# -----------------------------------------------------------------------
calculate_occupancy_metrics <- function(coords_sf, grid_size_km = 10) {
  if (nrow(coords_sf) == 0) {
    return(list(
      n_occupied_cells   = NA,
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
    n_occupied_cells   = n_occupied,
    total_area_occupied = total_area
  ))
}


# -----------------------------------------------------------------------
# 3. YEAR-BY-YEAR SPREAD DYNAMICS
# -----------------------------------------------------------------------
calculate_yearly_spread_metrics <- function(yearly_metrics_df) {
  if (nrow(yearly_metrics_df) < 2) {
    yearly_metrics_df$spread_rate_km_year      <- NA
    yearly_metrics_df$centroid_shift_km_year   <- NA
    yearly_metrics_df$area_change_km2_year     <- NA
    yearly_metrics_df$new_cells_occupied       <- NA
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
    prev_year    <- yearly_metrics_df$year[i-1]
    year_diff    <- current_year - prev_year
    
    # Centroid shift (using Euclidean distance since centroids are in lon/lat degrees)
    if (all(!is.na(yearly_metrics_df[i, c("centroid_lat","centroid_lon")]))) {
      # Convert degrees to approximate km using Haversine-like formula
      lat1 <- yearly_metrics_df$centroid_lat[i-1] * pi/180
      lat2 <- yearly_metrics_df$centroid_lat[i] * pi/180
      dlat <- (yearly_metrics_df$centroid_lat[i] - yearly_metrics_df$centroid_lat[i-1]) * pi/180
      dlon <- (yearly_metrics_df$centroid_lon[i] - yearly_metrics_df$centroid_lon[i-1]) * pi/180
      
      a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
      c <- 2 * atan2(sqrt(a), sqrt(1-a))
      centroid_dist <- 6371 * c  # Earth radius in km
      
      centroid_shifts[i] <- centroid_dist / year_diff
    }
    
    # Range expansion rate
    if (!is.na(yearly_metrics_df$max_distance_km[i]) &&
        !is.na(yearly_metrics_df$max_distance_km[i-1])) {
      spread_rates[i] <- (yearly_metrics_df$max_distance_km[i] -
                            yearly_metrics_df$max_distance_km[i-1]) / year_diff
    }
    
    # Area change (uses cluster hull area now)
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
  
  yearly_metrics_df$spread_rate_km_year    <- spread_rates
  yearly_metrics_df$centroid_shift_km_year <- centroid_shifts
  yearly_metrics_df$area_change_km2_year   <- area_changes
  yearly_metrics_df$new_cells_occupied     <- new_cells
  
  yearly_metrics_df
}


# -----------------------------------------------------------------------
# 4. OVERALL SPECIES DIRECTION
# -----------------------------------------------------------------------
calculate_overall_direction <- function(yearly_metrics_df) {
  if (nrow(yearly_metrics_df) < 2) return(NA)
  
  valid_centroids <- yearly_metrics_df[!is.na(yearly_metrics_df$centroid_lat) & 
                                         !is.na(yearly_metrics_df$centroid_lon), ]
  if (nrow(valid_centroids) < 2) return(NA)
  
  first_row <- valid_centroids[1, ]
  last_row  <- valid_centroids[nrow(valid_centroids), ]
  
  total_lon_change <- last_row$centroid_lon - first_row$centroid_lon
  total_lat_change <- last_row$centroid_lat - first_row$centroid_lat
  
  if (abs(total_lon_change) > abs(total_lat_change)) {
    ifelse(total_lon_change > 0, "Eastward", "Westward")
  } else {
    ifelse(total_lat_change > 0, "Northward", "Southward")
  }
}


# -----------------------------------------------------------------------
# 5. INVASION FRONT ANALYSIS
# -----------------------------------------------------------------------
calculate_invasion_front <- function(species_data, yearly_metrics_df) {
  if (nrow(yearly_metrics_df) < 2) {
    yearly_metrics_df$invasion_front_km <- NA
    return(yearly_metrics_df)
  }
  
  # First occurrence = origin
  first_year_data <- species_data %>%
    filter(year == min(year, na.rm = TRUE))
  if (nrow(first_year_data) == 0) {
    yearly_metrics_df$invasion_front_km <- NA
    return(yearly_metrics_df)
  }
  
  origin_coords <- tryCatch({
    origin_centroid <- st_centroid(st_union(first_year_data))
    st_coordinates(origin_centroid)
  }, error = function(e) {
    coords <- st_coordinates(first_year_data)
    c(mean(coords[,1], na.rm = TRUE), mean(coords[,2], na.rm = TRUE))
  })
  
  invasion_fronts <- numeric(nrow(yearly_metrics_df))
  for (i in 1:nrow(yearly_metrics_df)) {
    year_data <- species_data %>% filter(year == yearly_metrics_df$year[i])
    if (nrow(year_data) > 0 && !any(is.na(origin_coords))) {
      year_coords <- st_coordinates(year_data)
      distances <- apply(year_coords, 1, function(pt) {
        # Manual Haversine calculation since coords are in lon/lat degrees
        lat1 <- origin_coords[2] * pi/180
        lat2 <- pt[2] * pi/180
        dlat <- (pt[2] - origin_coords[2]) * pi/180
        dlon <- (pt[1] - origin_coords[1]) * pi/180
        
        a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
        c <- 2 * atan2(sqrt(a), sqrt(1-a))
        6371 * c  # Distance in km
      })
      invasion_fronts[i] <- max(distances, na.rm = TRUE)
    } else {
      invasion_fronts[i] <- NA
    }
  }
  yearly_metrics_df$invasion_front_km <- invasion_fronts
  yearly_metrics_df
}


# -----------------------------------------------------------------------
# 6. MASTER WRAPPER
# -----------------------------------------------------------------------
analyze_species_spread <- function(species_data, species_name, apply_filtering = FALSE) {
  cat(paste("Analyzing", species_name, "...\n"))
  
  # DBSCAN filtering (outlier cleanup)
  if (apply_filtering && nrow(species_data) > 3) {
    coords <- st_coordinates(species_data)
    if (nrow(coords) >= 4) {
      db_result <- dbscan(coords, eps = 0.5, minPts = 2)
      keep_indices <- db_result$cluster > 0 | (db_result$cluster == 0 & 
                                                 apply(coords, 1, function(x) {
                                                   dists <- sqrt(rowSums((coords - matrix(x, nrow(coords), 2, byrow=TRUE))^2)) / 1000
                                                   sum(dists <= 50000 & dists > 0) >= 1   # 50 km threshold
                                                 }))
      species_data <- species_data[keep_indices, ]
    }
  }
  if (nrow(species_data) == 0) return(NULL)
  
  # Create decade bins
  species_data <- species_data %>% mutate(decade = floor(year/10)*10) # lets see the numbers of years.. just for try
  # species_data <- species_data %>% mutate(decade = floor(year/5)*5)
  decades_available <- sort(unique(species_data$decade))
  
  decadal_metrics <- data.frame()
  for (dec in decades_available) {
    dec_data <- species_data %>% filter(decade == !!dec)
    if (nrow(dec_data) == 0) next
    
    range_metrics    <- calculate_range_metrics(dec_data)
    occupancy_metrics <- calculate_occupancy_metrics(dec_data, grid_size_km = 10)
    
    dec_metrics <- data.frame(
      species                = species_name,
      decade                 = dec,
      n_records              = nrow(dec_data),
      cluster_hull_area_km2  = range_metrics$cluster_hull_area,
      max_distance_km        = range_metrics$max_distance,
      range_span_ew_deg      = range_metrics$range_span_ew,
      range_span_ns_deg      = range_metrics$range_span_ns,
      n_occupied_cells_10km  = occupancy_metrics$n_occupied_cells,
      occupied_area_km2      = occupancy_metrics$total_area_occupied,
      centroid_lat           = range_metrics$centroid_lat,
      centroid_lon           = range_metrics$centroid_lon
    )
    decadal_metrics <- rbind(decadal_metrics, dec_metrics)
  }
  
  # Invasion front & temporal dyn
  if (nrow(decadal_metrics) >= 1) {
    tmp <- decadal_metrics %>% rename(year = decade)
    tmp <- calculate_invasion_front(species_data, tmp)
    decadal_metrics <- tmp %>% rename(decade = year)
  } else {
    decadal_metrics$invasion_front_km <- NA
  }
  
  if (nrow(decadal_metrics) >= 2) {
    tmp <- decadal_metrics %>% rename(year = decade)
    tmp <- calculate_yearly_spread_metrics(tmp)
    decadal_metrics <- tmp %>% rename(decade = year)
    decadal_metrics$dominant_direction <- calculate_overall_direction(
      decadal_metrics %>% rename(year = decade))
  } else {
    decadal_metrics$spread_rate_km_year   <- NA
    decadal_metrics$centroid_shift_km_year<- NA
    decadal_metrics$area_change_km2_year  <- NA
    decadal_metrics$new_cells_occupied    <- NA
    decadal_metrics$dominant_direction    <- NA
  }
  
  decadal_metrics
}


clean.filter.species <- function(data){
  remove1 <- c( "Pueraria montana", "Pycnonotus cafer","Humulus scandens","Gambusia affinis",
                "Cipangopaludina chinensis","Prosopis juliflora","Corvus splendens","Acridotheres cristatellus")
  # remove <- c( "Corvus splendens", "Faxonius limosus", "Prosopis juliflora","Gambusia affinis",
  #             "Humulus scandens",  "Pycnonotus cafer" ,"Asclepias syriaca","Gunnera tinctoria")
  
  data1 <- data[!data$species %in% remove1, ] 
  
  
  taxonomy <- read_xlsx("Database/Species.Taxonomy.xlsx")
  colnames(taxonomy)[1] <- "species"
  data1 <- data1 %>% left_join(taxonomy[,c("species","Group")], by = "species")
  data1$Group[data1$species =="Faxonius limosus"] <- "Crustaceans"
  data1$Group[data1$species =="Trachemys scripta"] <- "Reptiles"
  data1$Group[data1$species =="Vespa velutina"] <- "Insects"
  data1$Group[data1$species =="Cenchrus setaceus"] <- "Vascular plants"
  data1$Group[data1$species =="Pontederia crassipes"] <- "Vascular plants"
  data1$Group[data1$species =="Mustela vison"] <- "Mammals"
  data1$Group[data1$species =="Faxonius limosus"] <- "Crustaceans"
  
  return(data1)
}


# Logistic model function
logistic_model <- function(t, N0, K, r, tmin) {
  K / (1 + ((K - N0) / N0) * exp(-r * (t - tmin)))
}

# Fitting function
# fit_logistic_species <- function(df) {
  df <- df %>%
    arrange(year) %>% filter(year > 0) %>% 
    mutate(cum_cells = cumsum(n_cells_occupied))
  
  tmin <- min(df$year)
  
  out <- tryCatch({
    fit <- nlsLM(
      cum_cells ~ logistic_model(year, N0, K, r, tmin),
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
    
    # Marginal abundance peak
    n_max <- (r * K) / (4 * N0)
    
    # Intro and lag times (using Soto et al. Eq. 3)
    alpha_intro <- 0.01
    alpha_lag <- 0.1
    
    t_intro <- tmin - (1/r) * log((2 - alpha_intro + 2*sqrt(1 - alpha_intro)) / 
                                    (alpha_intro * ((K/N0) - 1)))
    t_lag <- tmin - (1/r) * log((2 - alpha_lag + 2*sqrt(1 - alpha_lag)) / 
                                  (alpha_lag * ((K/N0) - 1)))
    lag_duration <- t_lag - t_intro
    
    # Predictions for plot
    df_pred <- data.frame(year = seq(min(df$year), max(df$year), 1))
    df_pred$predicted <- logistic_model(df_pred$year, N0, K, r, tmin)

   df_pred$species =  unique(df$species)
df_pred = df_pred %>% left_join(df[,c(1,4)], by ='species')
df_pred = df_pred[!duplicated(df_pred[,c(1:3)]) , ]

p <- ggplot(df, aes(x = year, y = cum_cells)) + ylim(0,max(df_pred$predicted)) + 

  geom_point(aes(color = Group), size = 3, alpha = 0.4) +
  geom_line(data = df_pred, aes(x = year, y = predicted, color = Group), size = 1) +
  geom_vline(xintercept = t_inf, linetype = "dashed", color = "black") +
  #geom_hline(yintercept = N_inf, linetype = "dotted", color = "black") +
  geom_vline(xintercept = t_sat, linetype = "dotted", color = "black") +
 # geom_hline(yintercept = N_sat, linetype = "dotted", color = "black") +
  annotate("text", x = t_inf, y = N_inf, label = "Inflection", vjust = -1, color = "black") +
  annotate("text", x = t_sat, y = N_sat, label = "Saturation", vjust = -1, color = "black") +
  scale_color_manual(values = group_colors) + 
  labs(title = paste(unique(df$species)),
       x = "Year", y = "Cumulative occupied cells") +
  theme_bw() + theme(legend.position = 'none') 
p
    
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


assign_grid_cells <- function(sf_data, cell_size = 10000) {

  sf_projected <- st_transform(sf_data, crs = 3035)
  coords <- st_coordinates(sf_projected)
  
  grid_x <- floor(coords[, "X"] / cell_size)
  grid_y <- floor(coords[, "Y"] / cell_size)
  

  cell_id <- paste(grid_x, grid_y, sep = "_")
  
  return(cell_id)
}



#### Improved version:----------
# -----------------------------------------------------------------------
# 5. INVASION FRONT ANALYSIS (CLUSTERED APPROACH)
# -----------------------------------------------------------------------
calculate_clustered_invasion_front <- function(species_data, 
                                               yearly_metrics_df, 
                                               cluster_eps_km = 50) {
  
  # Ensure there's enough data to work with
  if (nrow(yearly_metrics_df) < 1 || nrow(species_data) < 3) {
    yearly_metrics_df$invasion_front_km <- NA
    return(yearly_metrics_df)
  }
  
  # --- Step 1: Identify persistent clusters across all data using DBSCAN ---
  # Project data to a metric CRS for distance calculations
  species_proj <- st_transform(species_data, 3035)
  coords_mat <- st_coordinates(species_proj)
  
  # Run DBSCAN. `eps` is the search radius in meters.
  db_result <- dbscan::dbscan(coords_mat, eps = cluster_eps_km * 1000, minPts = 3)
  species_data$cluster <- db_result$cluster
  
  # Filter out noise points (cluster 0)
  clusters <- unique(species_data$cluster[species_data$cluster > 0])
  
  if (length(clusters) == 0) {
    warning("DBSCAN did not identify any persistent clusters. Returning NA for invasion front.")
    yearly_metrics_df$invasion_front_km <- NA
    return(yearly_metrics_df)
  }
  
  # --- Step 2: Find the origin (first appearance centroid) for EACH cluster ---
  cluster_origins <- list()
  for (cid in clusters) {
    cluster_points <- species_data %>% filter(cluster == cid)
    
    # Find the first year this cluster appears
    first_year <- min(cluster_points$year, na.rm = TRUE)
    
    # Get the points from that first year to define the origin
    origin_points <- cluster_points %>% filter(year == first_year)
    
    # Calculate the centroid of these first points (lon/lat)
    origin_centroid <- st_centroid(st_union(origin_points))
    cluster_origins[[as.character(cid)]] <- st_coordinates(origin_centroid)
  }
  
  # --- Step 3: Calculate max distance from origin for each cluster, for each year ---
  invasion_fronts <- numeric(nrow(yearly_metrics_df))
  
  for (i in 1:nrow(yearly_metrics_df)) {
    current_year <- yearly_metrics_df$year[i]
    year_data <- species_data %>% filter(year == current_year)
    
    # Store the front size for each cluster active in this year
    fronts_this_year <- c()
    
    active_clusters_this_year <- unique(year_data$cluster[year_data$cluster > 0])
    
    if (length(active_clusters_this_year) > 0) {
      for (cid in active_clusters_this_year) {
        
        cid_char <- as.character(cid)
        
        # Ensure the cluster has a defined origin
        if (cid_char %in% names(cluster_origins)) {
          origin_coords <- cluster_origins[[cid_char]]
          
          # Get points for this cluster in this year
          year_cluster_data <- year_data %>% filter(cluster == cid)
          year_coords <- st_coordinates(year_cluster_data)
          
          # Calculate Haversine distance from this cluster's origin to all its points this year
          distances <- apply(year_coords, 1, function(pt) {
            lat1 <- origin_coords[2] * pi/180
            lat2 <- pt[2] * pi/180
            dlat <- (pt[2] - origin_coords[2]) * pi/180
            dlon <- (pt[1] - origin_coords[1]) * pi/180
            
            a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
            c <- 2 * atan2(sqrt(a), sqrt(1-a))
            6371 * c # Distance in km
          })
          
          fronts_this_year <- c(fronts_this_year, max(distances, na.rm = TRUE))
        }
      }
    }
    
    # The final metric is the MAXIMUM front found across all active clusters
    if (length(fronts_this_year) > 0) {
      invasion_fronts[i] <- max(fronts_this_year, na.rm = TRUE)
    } else {
      invasion_fronts[i] <- NA
    }
  }
  
  yearly_metrics_df$invasion_front_km <- invasion_fronts
  return(yearly_metrics_df)
}


# -----------------------------------------------------------------------
# 6. MASTER WRAPPER (MODIFIED)
# -----------------------------------------------------------------------
 # =======================================================================
# MASTER WRAPPER FUNCTION (FINAL VERSION)
# =======================================================================
analyze_species_spread <- function(species_data, species_name, apply_filtering = FALSE) {
  cat(paste("Analyzing", species_name, "...\n"))

  # ---------------------------------------------------------------------
  # 1. PRE-PROCESSING AND FILTERING (Optional)
  # ---------------------------------------------------------------------
  # DBSCAN filtering for outlier cleanup if requested
  if (apply_filtering && nrow(species_data) > 3) {
    coords <- st_coordinates(species_data)
    if (nrow(coords) >= 4) {
      db_result <- dbscan::dbscan(coords, eps = 0.5, minPts = 2) # Note: eps is in degrees here
      keep_indices <- db_result$cluster > 0
      species_data <- species_data[keep_indices, ]
    }
  }
  
  # Exit if no data remains after filtering
  if (nrow(species_data) == 0) {
    cat(paste("  -> No valid records found for", species_name, "after filtering. Skipping.\n"))
    return(NULL)
  }

  # ---------------------------------------------------------------------
  # 2. CALCULATE METRICS FOR EACH TIME PERIOD (DECADE/QUINQUENNIUM)
  # ---------------------------------------------------------------------
  # Create time period bins (e.g., 5-year periods)
  species_data <- species_data %>% mutate(time_period = floor(year / 5) * 5)
 # species_data <- species_data %>% mutate(decade = floor(year/10)*10) 
  time_periods_available <- sort(unique(species_data$time_period))

  period_metrics_list <- list()
  for (period in time_periods_available) {
    period_data <- species_data %>% filter(time_period == !!period)
    if (nrow(period_data) == 0) next

    # Calculate geographic and occupancy metrics for the current time period
    range_metrics <- calculate_range_metrics(period_data)
    occupancy_metrics <- calculate_occupancy_metrics(period_data, grid_size_km = 10)

    # Combine metrics into a single row
    period_metrics <- data.frame(
      species = species_name,
      year = period, # Using 'year' as the standard time column name
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
  
  # Combine all periodic metrics into a single data frame
  if (length(period_metrics_list) == 0) {
    cat(paste("  -> Not enough periodic data to analyze for", species_name, ". Skipping.\n"))
    return(NULL)
  }
  final_metrics_df <- bind_rows(period_metrics_list)

  # ---------------------------------------------------------------------
  # 3. CALCULATE DYNAMIC METRICS ACROSS TIME
  # ---------------------------------------------------------------------
  # A. INVASION FRONT (using the new clustered method)
  if (nrow(final_metrics_df) >= 1) {
    final_metrics_df <- calculate_clustered_invasion_front(species_data, final_metrics_df, cluster_eps_km = 50)
  } else {
    final_metrics_df$invasion_front_km <- NA
  }

  # B. SPREAD RATE, CENTROID SHIFT, AND DIRECTION
  if (nrow(final_metrics_df) >= 2) {
    # Calculate year-over-year changes
    final_metrics_df <- calculate_yearly_spread_metrics(final_metrics_df)
    
    # Calculate the single dominant direction over the entire period
    final_metrics_df$dominant_direction <- calculate_overall_direction(final_metrics_df)
  } else {
    # Add empty columns if there's not enough data for temporal analysis
    final_metrics_df$spread_rate_km_year <- NA
    final_metrics_df$centroid_shift_km_year <- NA
    final_metrics_df$area_change_km2_year <- NA
    final_metrics_df$new_cells_occupied <- NA
    final_metrics_df$dominant_direction <- NA
  }
  
  # Rename 'year' column back to 'decade' or 'time_period' for clarity
  final_metrics_df <- final_metrics_df %>%
    rename(time_period = year)

  return(final_metrics_df)
}
