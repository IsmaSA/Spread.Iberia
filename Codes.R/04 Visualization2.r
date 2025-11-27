setwd("C:/Users/IsmaSA/Desktop/My projects/Spread")

### Visualization 2.0
rm(list =ls())

source("./R.Codes/functions.R")
pacman::p_load(ggplot2,directlabels,sf,rgbif,raster,terra,sp,rnaturalearth, xl,xlsx, scales, readxl,writexl, rnaturalearth,dplyr,ggspatial,
purrr,units,alphahull,dbscan,geosphere,alphahull,concaveman,patchwork, mapdeck, colourvalues)
library(colourvalues)

sp <- read_xlsx("./Database/combined_summary_all_groups.xlsx")
length( unique(sp$species)) -1 
table(sp[sp$analysis_group=="All_combined",]$Group)

sp <- sp[sp$species!="Cardiospermum grandiflorum", ]



iberia <- read_xlsx("./Database/ListNNS.Iberia.v.2.xlsx", sheet=2)
iberia <- iberia[iberia$Country %in% c("Spain", "Portugal"), ]

setdiff(sp$species, iberia$canonicalName)
###   Doubtfull species, they may be not established: 
# "Mustela vison"     ----> It is but with the name "Neogale vison"
# "Hakea sericea"     ---> We keep it based on https://doi.org/10.1080/23818107.2024.2318761       
# "Acridotheres cristatellus" --> Remove
# "Bipalium kewense"       ---> Keep  
# "Obama nungara"       ---> Keep     
# "Corvus splendens"     ---> Remove     
# "Prosopis juliflora"     ---> Remove   
# "Cipangopaludina chinensis" ---> Remove
# "Gambusia affinis"      ---> Remove   
# "Humulus scandens"       ---> Remove    
# "Pycnonotus cafer"    ---> Remove     
# "Pueraria montana"    ---> Remove  

#remove1 <- c( "Pueraria montana", "Pycnonotus cafer","Humulus scandens","Gambusia affinis",
 #             "Cipangopaludina chinensis","Prosopis juliflora","Corvus splendens","Acridotheres cristatellus")
# remove <- c( "Corvus splendens", "Faxonius limosus", "Prosopis juliflora","Gambusia affinis",
#             "Humulus scandens",  "Pycnonotus cafer" ,"Asclepias syriaca","Gunnera tinctoria")


#sp <- sp[!sp$species %in% remove1, ] 
#sp1 <- sp1[!sp1$species %in% remove1, ] 

#data <- rbind(sp,sp1)
#data <- data[data$species !="Asclepias syriaca",]
#data <- data[data$species !="Gunnera tinctoria",]
#data <- data[data$species !="Procambarus virginalis",] # Few records
#data <- data[data$species !="Lithobates catesbeianus",] # No established
#data <- data[data$species !="Elodea nuttallii",] # Few records 

data <- sp

options(scipen = 999)
unique(data$species) # this shoudl contains only 39 species


# Figure 2 a =================================================================
summary1 <- readRDS("./Database/cleaned_occ1.rds") 

summary1 <- summary1[summary1$species %in% unique(sp$species), ]
a <- summary1 %>% group_by(species) %>% summarise(n =n())
median(a$n)

df <- summary1 
df_clean <- df %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(year)) %>%
  filter(year >= 1950) %>% 
  filter(df$decimalLongitude > -10) %>% 
  filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 5000)  %>%
  filter(!(decimalLongitude > 1 & decimalLongitude < 4.5 &
           decimalLatitude > 38.5 & decimalLatitude < 40.5))  %>% filter(year < 2025) 

df_clean <- df_clean[df_clean$species %in% data$species, ]

fran <- read_xlsx("Database/Database_Spain.xlsx") %>%
  filter(!(decimalLongitude > 1 & decimalLongitude < 4.5 &
             decimalLatitude > 38.5 & decimalLatitude < 40.5))
fran <- fran[names(df_clean)]
missing_cols <- setdiff(names(df_clean), names(fran))
for (col in missing_cols) {
  fran[[col]] <- NA
}

fran$species[fran$species=="Trachemys scripta elegans"] <- "Trachemys scripta"
fran$species[fran$species=="Trachemys scripta scripta"] <- "Trachemys scripta"
fran$species[fran$species=="Vespa velutina nigrithorax"] <- "Vespa velutina"
fran$species[fran$species=="Neogale vison"] <- "Mustela vison"
fran$species[fran$species=="Eichhornia crassipes"] <- "Pontederia crassipes"
fran$species[fran$species=="Ameirus melas"] <- "Ameiurus melas"
fran$species[fran$species=="Ludwigia peploides subsp, montevidensis"] <- "Ludwigia peploides"
fran$species[fran$species=="Pennisetum setaceum"] <- "Cenchrus setaceus"

df_clean <- rbind(df_clean, fran)
df_clean <- df_clean[!df_clean$species %in% c("Lampropeltis getula","Tamias sibiricus","Nyctereutes procyonoides"),]
df_clean <- df_clean %>% mutate(
  decimalLatitude = as.numeric(decimalLatitude),
  decimalLongitude = as.numeric(decimalLongitude)) %>% filter(
    !is.na(decimalLatitude), 
    !is.na(decimalLongitude),
    decimalLatitude >= -90 & decimalLatitude <= 90,
    decimalLongitude >= -180 & decimalLongitude <= 180  )

df_clean <- df_clean %>% # to remove ceuta and melilla approxxxx
  filter(
    !(decimalLatitude >= 35.85 & decimalLatitude <= 35.92 &
        decimalLongitude >= -5.35 & decimalLongitude <= -5.28) & 
      !( decimalLatitude >= 35.26 & decimalLatitude <= 35.32 & decimalLongitude >= -2.96 & decimalLongitude <= -2.90) )



# Number of species at square 0.1
points_sf <- st_as_sf(df_clean, coords = c("decimalLongitude","decimalLatitude"), crs = 4326)
bbox <- st_bbox(points_sf)
hex_grid <- st_make_grid(st_as_sfc(bbox), cellsize = 0.1, square = FALSE) |> 
  st_sf(geometry = _)

  
  df_clean1 <- df_clean
  pts <- st_as_sf(df_clean1, coords = c("decimalLongitude","decimalLatitude"), crs = 4326)
  join <- st_join(hex_grid, pts)
  
  richness <- join |>
    group_by(geometry) |>
    summarise(n_species = n_distinct(species)) |>
    ungroup()
yellow_to_red <- c( "#FFD700", "#FFA500", "#FF6347", "#FF4500", "#DC143C", "#8B0000")


richness$n_species <- as.numeric(richness$n_species)
richness$n_species[richness$n_species == 0] <- NA
pal <- colorRampPalette(c("#FFD700", "#8B0000"))
richness$col <- pal(100)[as.numeric(cut(richness$n_species, breaks = 100))]
richness$n_species_squared <- richness$n_species^2

library(mapdeck)
set_token('pk.eyJ1IjoiaXNtYXNhIiwiYSI6ImNtZTJmejE3cTB0Ym8ybHNhaTJyNWt0bHYifQ.eFeXSonlPcSvrWKJp5gsUw')

  
p1 <- mapdeck(
    style = "https://basemaps.cartocdn.com/gl/dark-matter-nolabels-gl-style/style.json",
    zoom = 5.3, pitch = 35, location = c(-3.5, 40)
  ) %>%
add_polygon(
  data = richness[richness$n_species > 1,],
  fill_colour = "n_species",
  palette = c("ylorrd"),  # Your desired yellow to red palette
  update_view = FALSE,
elevation = "n_species_squared",
  elevation_scale = 500,  # Increase for more dramatic height
  extruded = TRUE, 
  legend = TRUE,
  tooltip = "n_species",
  elevation_function = "sum",
  radius = 1500,
  legend_options = list(
    title = "Number of species",
    title_font_size = 18,
    title_font_color = "#EEE",
    text_font_size = 14,
    text_font_color = "#EEE",
    orientation = "vertical"
  )
)
p1      
htmlwidgets::saveWidget(p1, "mapdeck_map.html", selfcontained = TRUE)


#  figure 2 b ----------------------------------
group_colors <- c(
  "Algae"            = "#dede00",  # teal green
  "Birds"            = "#f20a0e",  # red
  "Crustaceans"      = "#984ea3",  # strong blue
  "Fishes"           = "#377eb8",  # purple
  "Insects"          = "#ff7f00",  # orange
  "Mammals"          = "#af7451",  # warm brown
  "Molluscs"         = "#dede00",  # medium green
  "Platyhelminthes"  = "#999999",  # grey
  "Reptiles"         = "#f781bf",  # pink
  "Amphibians" ="#dede00",
  "Plants"  = "#4daf4a"   # yellow
)
taxonomic_colors<- group_colors

sp <- read_xlsx("./Database/combined_summary_all_groups.xlsx")
sp <- sp[sp$analysis_group=="All_combined", ]
sp$species[sp$species=="Mustela vison"] <- "Neogale vison"
sp <- sp[sp$species!="Cardiospermum grandiflorum", ]

sp$last_period[sp$species=="Acridotheres tristis"] <- 2025
sp$last_period[sp$species=="Misgurnus anguillicaudatus"] <- 2025
sp$last_period[sp$species=="Ludwigia grandiflora"] <- 2025
sp$last_period[sp$species=="Oxyura jamaicensis"] <- 2025
sp$last_period[sp$species=="Myriophyllum heterophyllum"] <- 2025
sp$last_period[sp$species=="Threskiornis aethiopicus"] <- 2025
sp$last_period[sp$species=="Misgurnus anguillicaudatus"] <- 2025
sp$last_period[sp$species=="Hydrocotyle ranunculoides"] <- 2025


sp1 <- sp %>% arrange(-first_period)
sp1$species <- factor(sp1$species, levels = sp1$species)

ggplot(sp1, aes(y = species)) +
  geom_segment(aes(x = first_period, xend = last_period, 
                   y = species, yend = species, color = Group), size = 1) +
  geom_point(aes(x = first_period, color = Group), size = 3) +
  geom_point(aes(x = last_period, color = Group), size = 3) +
  labs(x = "Year", y = "") +
  scale_x_continuous(breaks = seq(1950, 2020, 10)) + 
  scale_color_manual(values = group_colors) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "italic"),
    legend.position = c(0.05, 0.1),      # bottom-left inside the plot
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),    
    legend.key.size = unit(0.9, "cm")         
  )

ggsave("Figure 1b.svg",
 width = 7, height = 7, dpi = 300)



##### FIGURE 3 ======================
rm(list=ls())

sp <- read_xlsx("./Database/combined_summary_all_groups.xlsx")
sp$species[sp$species=="Mustela vison"] <- "Neogale vison"
sp <- sp[sp$species!="Cardiospermum grandiflorum", ]
sp <- sp[sp$species!="Alternanthera philoxeroides"| sp$analysis_group =="All_combined", ]

#sp1<- clean.filter.species(sp)

# max_area_km2
# Maximum convex hull area (km²) occupied by the species across decades.
# max_distance_km
# Maximum pairwise distance (km) between occurrence points across decades.
# mean_spread_rate
# Mean annual spread rate of the effective range radius (km/year).
# mean_area_change
# Mean annual change in convex hull area (km²/year).
# Max cells occupied
# Max cells occupied per species from all sp


sp1 <- sp[order(sp$species, decreasing = FALSE), ]
names(sp1)

cols <- c("N.cells", "max_area_km2", "max_invasion_front_km", "mean_spread_rate")  

# Here I need to combine with the fiiiiiiiiiiiiiiiiiiother data

plots <- list()
species_order <- sp1 %>% #filter(analysis_group=="All_combined") %>% 
  group_by(species) %>%
  summarise(n.cells.occupied = sum(N.cells, na.rm = TRUE), .groups = "drop") %>%
  arrange(n.cells.occupied) %>% pull(species)

c = cols[1]

custom_labels <- c(
  "max_area_km2" = "Maximum area occupied (km²)",
  "max_invasion_front_km" = "Maximum radial distance (km)", 
  "mean_spread_rate" = "Spread rate (km/year)",
  "N.cells" = "Number of cells occupied (10×10 km)")

sp1 <- sp1 %>%
  mutate(across(c(max_area_km2, max_invasion_front_km, mean_spread_rate, N.cells), as.numeric))

plots <- list()
var<- cols[4]
for (var in cols) {
  sp1_plot <- sp1 %>%
    mutate(species = factor(species, levels = species_order))
  
  if (var %in% c(cols[1], cols[3])) {
    # Variables to invert (plotted as negative)
    p <- ggplot(sp1_plot, aes(x = species, 
                              y = -pmax(.data[[var]], 0.1), 
                              fill = Group,
                              alpha = analysis_group)) +
      geom_col(data = filter(sp1_plot, analysis_group == "Natural_only"), 
               width = 0.8, alpha = 0.5) +
      geom_col(data = filter(sp1_plot, analysis_group == "All_combined"), 
               width = 0.2, alpha = 1.0, position = position_nudge(x = 0.1)) +
      scale_fill_manual(values = taxonomic_colors) +
      scale_alpha_manual(values = c("Natural_only" = 0.5, "All_combined" = 1.0),
                         name = "Spread Type") +  
      labs(x = "", 
           y = custom_labels[var],
           fill = "Taxonomic Group") +
      coord_flip() + 
      theme_bw(base_size = 12) + 
      theme(
        legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        panel.grid.major.x = element_line(color = "grey", size = 0.5), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank()
      ) +
      scale_y_continuous(breaks = pretty_breaks(n = 9),
                         labels = function(x) abs(x)) +  
      scale_x_discrete(position = "top")
    
  } else {
    # Normal direction variables
    p <- ggplot(sp1_plot, aes(x = species, 
                              y = pmax(.data[[var]], 0.1), 
                              fill = Group,
                              alpha = analysis_group)) +
      geom_col(data = filter(sp1_plot, analysis_group == "Natural_only"), 
               width = 0.7, alpha = 0.5) +
      geom_col(data = filter(sp1_plot, analysis_group == "All_combined"), 
               width = 0.3, alpha = 1.0, position = position_nudge(x = 0.1)) +
      scale_fill_manual(values = taxonomic_colors) +
      scale_alpha_manual(values = c("Natural_only" = 0.5, "All_combined" = 1.0),
                         name = "Spread Type") +
      labs(x = "", 
           y = custom_labels[var],
           fill = "Taxonomic Group") +
      coord_flip() + 
      theme_bw(base_size = 12) + 
      theme(
        legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        panel.grid.major.x = element_line(color = "grey", size = 0.5), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank()
      ) +
      scale_y_continuous(breaks = pretty_breaks(n = 9))
  }
  
  plots[[var]] <- p
}
library(gridExtra)
a=grid.arrange(grobs = plots, ncol = 2)

names(plots)

ggsave(
  filename = "03 summary.svg",  
  plot = a,            
  width = 15,                       
  height = 16,                      
  dpi = 300 )




#############  FIGURE 4 ---------------------------------------------------------------------

# dominant_direction
# mean_centroid_shift

rm(list=ls())
sp <- read_xlsx("./Database/combined_summary_all_groups.xlsx")
sp$species[sp$species=="Mustela vison"] <- "Neogale vison"
sp <- sp[sp$species!="Cardiospermum grandiflorum", ]
sp <- sp[sp$species!="Alternanthera philoxeroides"| sp$analysis_group =="All_combined", ]


summary1 <- readRDS("./Database/cleaned_occ1.rds") 
df_clean <- summary1 %>%  filter(
    species %in% sp$species, !is.na(decimalLatitude),
    !is.na(decimalLongitude),  !is.na(year),
    year >= 1950,  decimalLongitude > -10,
    coordinateUncertaintyInMeters <= 5000, !(decimalLongitude > 1 & decimalLongitude < 4.5 &
      decimalLatitude > 38.5 & decimalLatitude < 40.5) )

str(df_clean)

calculate_centroids <- function(df) {
  df %>%
    group_by(species) %>%
    summarise(
      min_year = min(year, na.rm = TRUE),
      max_year = max(year, na.rm = TRUE)
    ) %>%
    rowwise() %>%
    mutate(
      past_range = list(min_year:(min_year+9)),
      last_range = list((max_year-9):max_year)
    ) %>%
    left_join(df, by = "species") %>%
    group_by(species) %>%
    summarise(
      past_centroid_lat = mean(decimalLatitude[year %in% unlist(first(past_range))], na.rm = TRUE),
      past_centroid_lon = mean(decimalLongitude[year %in% unlist(first(past_range))], na.rm = TRUE),
      last_centroid_lat = mean(decimalLatitude[year %in% unlist(first(last_range))], na.rm = TRUE),
      last_centroid_lon = mean(decimalLongitude[year %in% unlist(first(last_range))], na.rm = TRUE),
      .groups = "drop"
    )
}

centroids_df <- calculate_centroids(df_clean)

head(centroids_df)

taxonomy <- read_xlsx("Database/Species.Taxonomy.xlsx")
colnames(taxonomy)[1] <- "species"

centroids_df <- centroids_df %>%  left_join(sp[,c(1,18)], by = "species")

iberia <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(admin %in% c("Spain", "Portugal"))

centroids_df <- centroids_df %>%
  mutate(
    past_centroid_lat = as.numeric(unlist(past_centroid_lat)),
    past_centroid_lon = as.numeric(unlist(past_centroid_lon)),
    last_centroid_lat = as.numeric(unlist(last_centroid_lat)),
    last_centroid_lon = as.numeric(unlist(last_centroid_lon)),
    bearing = as.numeric(unlist(bearing))
  ) %>%
  distinct()  

centroids_long <- centroids_df %>%
  pivot_longer(
    cols = c(past_centroid_lat, past_centroid_lon, last_centroid_lat, last_centroid_lon),
    names_to = c("period", "coordinate"),
    names_pattern = "(past|last)_centroid_(lat|lon)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = coordinate,
    values_from = value
  ) %>%
  rename(latitude = lat, longitude = lon) %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    genus = word(species, 1),
    epithet = word(species, 2),
    species_short = paste0(substr(genus, 1, 1), ". ", epithet),
    species_label = paste0("italic('", species_short, "')")
  )

group_colors


library(ggrepel)
p1<- ggplot(data = iberia) +
  geom_sf(fill = "antiquewhite") +
  geom_segment(data = centroids_df,
               aes(x = past_centroid_lon, y = past_centroid_lat,
                   xend = last_centroid_lon, yend = last_centroid_lat,
                   color = Group),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7) +
  # Plot centroid points
  geom_point(data = centroids_long,
             aes(x = longitude, y = latitude, color = Group, shape = period),
             size = 3, alpha = 0.8) +
  geom_text_repel(
    data = centroids_long %>% filter(period == "last" & !is.na(longitude) & !is.na(latitude)),
    aes(x = longitude, y = latitude, label = species_label, color = Group),
    parse = TRUE,           # <= important: tells ggplot to interpret italics
    size = 3, max.overlaps = 30
  )  +scale_color_manual(values=group_colors) +
  scale_shape_manual(values = c(past = 17, last = 19)) + # triangle = past, circle = last
  coord_sf(xlim = c(-10, 5), ylim = c(36, 44), expand = FALSE) +
  theme_void() +
  labs(color = "Taxonomic Group",
       shape = "Period") +
  theme(legend.position = "none")
p1

centroids_df <- centroids_df %>%
  mutate(
    distance_km = distHaversine(
      cbind(past_centroid_lon, past_centroid_lat),
      cbind(last_centroid_lon, last_centroid_lat)
    ) / 1000,  
    
    bearing = bearing(
      cbind(past_centroid_lon, past_centroid_lat),
      cbind(last_centroid_lon, last_centroid_lat)),
    
    direction = case_when(
      bearing >= -22.5 & bearing < 22.5 ~ "N",
      bearing >= 22.5 & bearing < 67.5 ~ "NE", 
      bearing >= 67.5 & bearing < 112.5 ~ "E",
      bearing >= 112.5 & bearing < 157.5 ~ "SE",
      bearing >= 157.5 | bearing < -157.5 ~ "S",
      bearing >= -157.5 & bearing < -112.5 ~ "SW",
      bearing >= -112.5 & bearing < -67.5 ~ "W",
      bearing >= -67.5 & bearing < -22.5 ~ "NW"
    )
  )


library(dplyr)
library(ggplot2)

plot_df <- centroids_df %>%
  filter(!is.na(distance_km)) %>%
  mutate(
    label_expr = paste0("italic('", species, "')")  # full name in italics
  )

# Max for label placement
xmax <- max(plot_df$distance_km, na.rm = TRUE)

p2<- ggplot(plot_df,
       aes(x = distance_km,
           y = reorder(label_expr, distance_km),
           fill = Group)) +
  geom_col(width = 0.75) +
      scale_fill_manual(values = group_colors) +
  geom_text(aes(label = direction),
            hjust = -0.2, size = 3.3, color = "black") +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  coord_cartesian(xlim = c(0, xmax * 1.12)) +
  labs(
    x = "Distance (km)",
    y = "Species",
    fill = "Taxonomic group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 10))
  )

j<- p1 /p2 
ggsave(filename = "04 centroids.svg", plot = j, width = 10, height = 14, dpi = 300)
table(centroids_df$direction)


# This is the same but coloured by the direction  --
direction_colors <- c(
  N  = "#1F78B4",
  NE = "#33A8C7",
  E  = "#66C2A5",
  SE = "#B2DF8A",
  S  = "#FDB863",
  SW = "#E66101",
  W  = "#7570B3",
  NW = "#CC79A7")


p1 <- ggplot(data = iberia) +
  geom_sf(fill = "antiquewhite") +
  geom_segment(
    data = centroids_df,
    aes(
      x = past_centroid_lon, y = past_centroid_lat,
      xend = last_centroid_lon, yend = last_centroid_lat,
      color = direction
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    alpha = 0.7
  ) +
  geom_point(
    data = centroids_long,
    aes(x = longitude, y = latitude, color = direction, shape = period),
    size = 3, alpha = 0.8
  ) +
  geom_text_repel(
    data = centroids_long %>%
      filter(period == "last" & !is.na(longitude) & !is.na(latitude)),
    aes(x = longitude, y = latitude, label = species_label, color = direction),
    parse = TRUE,
    size = 3, max.overlaps = 30
  ) +
  scale_color_manual(values = direction_colors) +
  scale_shape_manual(values = c(past = 17, last = 19)) +
  coord_sf(xlim = c(-10, 5), ylim = c(36, 44), expand = FALSE) +
  theme_void() +
  labs(color = "Direction", shape = "Period") +
  theme(legend.position = "right")

p1


p2 <- ggplot(plot_df,
       aes(x = distance_km,
           y = reorder(label_expr, distance_km),
           fill = direction)) +  # <-- use direction column here
  geom_col(width = 0.75, col='black', size= 0.2) +
  scale_fill_manual(values = direction_colors) +  # <-- your named vector of colors
  geom_text(aes(label = direction),
            hjust = -0.2, size = 3.3, color = "black") +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  coord_cartesian(xlim = c(0, xmax * 1.12)) +
  labs(
    x = "Distance (km)",
    y = "Species",
    fill = "Direction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 10))
  )
p2

j<- p1 /p2 
ggsave(filename = "Figure S3 centroids.svg", plot = j, width = 10, height = 14, dpi = 300)


###############  FIGURE  5  -------------------------------------------
rm(list=ls())

sp <- read_xlsx("./Database/detailed_metrics_All_combined.xlsx")
sp <- sp[sp$analysis_group=="All_combined", ]
sp$species[sp$species=="Mustela vison"] <- "Neogale vison"
sp <- sp[sp$species!="Cardiospermum grandiflorum", ]

all_species_metrics = sp
library(scales)

# Broussonetia papyrifera
highlighted_species <- c("Procambarus clarkii", 'Ludwigia grandiflora','Myocastor coypus')

highlighted_species <- all_species_metrics %>% arrange(desc(.[["species"]])) %>%  slice(251:299)
highlighted_species <-highlighted_species$species
all_species_metrics <- all_species_metrics[order(all_species_metrics$species, decreasing = FALSE), ]

all_species_metrics <- all_species_metrics %>%
  mutate(is_highlighted = species %in% highlighted_species)

highlighted_data <- all_species_metrics %>% filter(is_highlighted) 
background_data <- all_species_metrics %>% filter(!is_highlighted) 

highlight_species <- unique(highlighted_data$species)


last_points <- highlighted_data %>%
  group_by(species) %>%
  filter(time_period == max(time_period)) %>%
  ungroup() %>%
  mutate(label = paste0("italic('", species, "')"))
last_points <- last_points %>%
  mutate(
    binomial = str_extract(species, "^[A-Za-z]+\\s+[A-Za-z-]+"),   # "Genus epithet" or NA
    binomial = if_else(is.na(binomial), str_extract(species, "^[A-Za-z]+"), binomial),
    
    genus = word(binomial, 1),
    epithet = word(binomial, 2),
    
    genus_init = substr(genus, 1, 1),
    epithet = ifelse(is.na(epithet), "", epithet),
    
    short_plain = ifelse(epithet == "", paste0(genus_init, "."), paste0(genus_init, ". ", epithet)),
    
    short_expr = ifelse(
      epithet == "",
      paste0("italic(", genus_init, ".)"),
      paste0("italic(", genus_init, ".~", epithet, ")")
    )
  )
p <- ggplot() +
  geom_rect(aes(xmin = 2014, xmax = 2025, ymin = -Inf, ymax = Inf),
            fill = "turquoise", alpha = 0.2) +
  geom_line(data = background_data, 
            aes(x = time_period, y = spread_rate_km_year, group = species),
            size = 0.7, color = "grey80", alpha = 0.6) +
  
  geom_line(data = highlighted_data,
            aes(x = time_period, y = spread_rate_km_year, group = species, color = Group),
            size = 1.2) +
  
  geom_text_repel(data = last_points,
                  aes(x = time_period, y = spread_rate_km_year, 
                      label = short_expr, color = Group),  # <- use short_label
                  parse = TRUE,    # no need for parsing if not adding expressions
                  nudge_x = 1,
                  direction = "y",
                  hjust = 0,
                  segment.color = NA,
                  size = 3.5,
                  inherit.aes = FALSE) +
  labs(x="Year", y = "Spread rate (km/year)") + 
  
  geom_vline(xintercept = 2014, color = "grey60", linetype = "dashed", size = 0.8) +
  
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 1),
    breaks = scales::pretty_breaks(n = 7)
  ) +  
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 7),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")
p

# the other
all_species_metrics <- read_xlsx("./Database/species_spread_metrics_detailed_All.xlsx")
all_species_metrics<- clean.filter.species(all_species_metrics)
all_species_metrics <- all_species_metrics[all_species_metrics$species %in% data$species,]

library(scales)

# Broussonetia papyrifera
highlighted_species <- c("Trachemys scripta", "Procambarus clarkii","Ondatra zibethicus",
"Myriophyllum heterophyllum","Ludwigia peploides","Ludwigia grandiflora","Gambusia holbrooki",
"Cardiospermum grandiflorum", "Acacia saligna"      , "Acridotheres tristis")

#highlighted_species <- all_species_metrics %>% arrange(desc(.[["species"]])) %>%  slice(171:195)
#highlighted_species <-highlighted_species$species
all_species_metrics <- all_species_metrics[order(all_species_metrics$species, decreasing = FALSE), ]

all_species_metrics <- all_species_metrics %>%
  mutate(is_highlighted = species %in% highlighted_species)

highlighted_data <- all_species_metrics %>% filter(is_highlighted) 
background_data <- all_species_metrics %>% filter(!is_highlighted) 

highlight_species <- unique(highlighted_data$species)


last_points <- highlighted_data %>%
  group_by(species) %>%
  filter(decade == max(decade)) %>%
  ungroup() %>%
  mutate(label = paste0("italic('", species, "')"))
last_points <- last_points %>%
  mutate(
    binomial = str_extract(species, "^[A-Za-z]+\\s+[A-Za-z-]+"),   # "Genus epithet" or NA
    binomial = if_else(is.na(binomial), str_extract(species, "^[A-Za-z]+"), binomial),
    
    genus = word(binomial, 1),
    epithet = word(binomial, 2),
    
    genus_init = substr(genus, 1, 1),
    epithet = ifelse(is.na(epithet), "", epithet),
    
    short_plain = ifelse(epithet == "", paste0(genus_init, "."), paste0(genus_init, ". ", epithet)),
    
    short_expr = ifelse(
      epithet == "",
      paste0("italic(", genus_init, ".)"),
      paste0("italic(", genus_init, ".~", epithet, ")")
    )
  )
p1 <- ggplot() +
  geom_rect(aes(xmin = 2010, xmax = 2020, ymin = -Inf, ymax = Inf),
            fill = "turquoise", alpha = 0.2) +
  geom_line(data = background_data, 
            aes(x = decade, y = spread_rate_km_year, group = species),
            size = 0.7, color = "grey80", alpha = 0.6) +
  
  geom_line(data = highlighted_data,
            aes(x = decade, y = spread_rate_km_year, group = species, color = Group),
            size = 1.2) +
  
  geom_text_repel(data = last_points,
                  aes(x = decade, y = spread_rate_km_year, 
                      label = short_expr, color = Group),  # <- use short_label
                  parse = TRUE,    # no need for parsing if not adding expressions
                  nudge_x = 1,
                  direction = "y",
                  hjust = 0,
                  segment.color = NA,
                  size = 3.5,
                  inherit.aes = FALSE) +
  labs(x="Year", y = "Spread rate (km/year)") + 
  
  geom_vline(xintercept = 2014, color = "grey60", linetype = "dashed", size = 0.8) +
  
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 1),
    breaks = scales::pretty_breaks(n = 7)
  ) +  
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 7),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")
p1

pp <- p +p1

ggsave(
  filename = "05 Spread.Rate.svg",  
  plot = pp,            
  width = 14,                       
  height = 7,                      
  dpi = 300 )



# This is other style for the plot:            ATENTIONNNNNNNNNNNNNN

sp <- read_xlsx("./Database/detailed_metrics_All_combined.xlsx")
sp$species[sp$species=="Mustela vison"] <- "Neogale vison"
sp <- sp[sp$species!="Cardiospermum grandiflorum", ]
s <- read_xlsx("./Database/detailed_metrics_Natural_only.xlsx")
s$species[s$species=="Mustela vison"] <- "Neogale vison"
s <- s[s$species!="Cardiospermum grandiflorum", ]


names(a)
names(s)
a= sp[,c(1,2,15,20,22)]
spp= s[,c(1,2,15,20,22)]
sppp= rbind(a,spp)

sp1 <- read_xlsx("./Database/combined_summary_all_groups.xlsx")
sp1$species[sp1$species=="Mustela vison"] <- "Neogale vison"
sp1 <- sp1[sp1$species!="Cardiospermum grandiflorum", ]


names(sp1)
b= sp1[,c(1,2,12:13,18)]

super.df = sppp %>% left_join(b, by = c('species','analysis_group'))
super.df = super.df[super.df$spread_rate_km_year !=0, ]
super.df = super.df[!is.na(super.df$spread_rate_km_year), ]
names(super.df)
super.df$analysis_group



opacity_levels_plot4 <- c(
  "All_combined" = 1,
  "Natural_only" = 0.6)

species_summary <- super.df %>%
  group_by(Group.x, species, analysis_group) %>%
  summarise(
    mean_spread = mean(spread_rate_km_year, na.rm = TRUE),
    min_spread = min(spread_rate_km_year, na.rm = TRUE),
    max_spread = max(spread_rate_km_year, na.rm = TRUE),
    .groups = "drop" )

species_summary <- species_summary %>%
  mutate(
    species = reorder(species, mean_spread, FUN = max),
    analysis_group = factor(analysis_group, levels = c("Natural_only", "All_combined")))

plot4 <- ggplot(
  data = species_summary, aes(
    y = species,
    color = Group.x,
    fill = Group.x,
    alpha = analysis_group)
) + geom_errorbarh(
    aes(xmin = min_spread, xmax = max_spread),
    height = 0.25,
    linewidth = 0.7,
    position = position_dodge(width = 0.5)
  ) +geom_point(
    aes(x = mean_spread),
    shape = 21,              
    size = 3.8,
    stroke = 0.6,
    color = "black",
    position = position_dodge(width = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_alpha_manual(
    name = "Expansion type:",
    values = opacity_levels_plot4,
    labels = c("Natural only", "Natural + jump")
  ) +
  labs(
    x = "Spread rate (km/year)",
    y = "" ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(size = 9, face = "italic"),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12)  )

plot4

ggsave("./Plots/05 Spread rates.svg", 
       plot = plot4, 
width = 9,height =10, dpi = 300)


###############  FIGURE  6  -------------------------------------------
rm(list=ls())

library(readxl)
library(dplyr)
library(ggplot2)
library(minpack.lm)
library(openxlsx)

all_params <- read_xlsx("./Database/logistic_parameters.xlsx")
summary_stats <- all_params %>%
  filter(converged == TRUE) %>%
  summarise(
    n_species = n(),
    mean_r = mean(r, na.rm = TRUE),
    median_r = median(r, na.rm = TRUE),
    mean_K = mean(K, na.rm = TRUE),
    median_K = median(K, na.rm = TRUE),
    mean_lag_duration = mean(lag_duration, na.rm = TRUE),
    median_lag_duration = median(lag_duration, na.rm = TRUE))
summary_stats
#n_species    mean_r  median_r   mean_K median_K mean_lag_duration median_lag_duration
#       40    0.236    0.194      1685.     453.              12.3
mean(all_params$lag_duration)
sd(all_params$lag_duration)
range(all_params$lag_duration)


if (sum(all_params$converged, na.rm = TRUE) > 0) {
  top_spreaders <- all_params %>%
    filter(converged == TRUE) %>%
    arrange(desc(r)) %>%
    head(5)
  
  cat("\nTop 5 fastest spreaders (highest r):\n")
  print(top_spreaders[, c("species", "r", "K", "lag_duration")])
}

library(gridExtra)

plot_list <- lapply(all_results, function(res) res$plot)
plot_list <- plot_list[!sapply(plot_list, is.null)]  # remove failed fits

plot_list_modified <- lapply(seq_along(all_results), function(i) {
  res <- all_results[[i]]
  if (is.null(res$plot) | is.null(res$params)) return(NULL)
  
  sp_name <- res$params$species[1]
  if (sp_name == "Mustela vison") sp_name <- "Neogale vison"
  
  t_sat   <- res$params$t_sat[1]
  N_sat   <- res$params$N_sat[1]
  t_lag <- res$params$t_lag[1]

  p <- res$plot +
    #geom_vline(xintercept = t_sat, linetype = "dashed", color = "black", size = 0.3) +
    #geom_vline(xintercept = t_lag, linetype = "dotted", color = "black", size = 0.3) +
    #annotate("text", x = t_lag, y = N_sat, label = "Saturation",
    #         color = "black", size = 2, vjust = -1) +
    theme_bw() +
    labs(x = "Year", y = "Cumulative occupied cells",
         title = paste0(sp_name)) +  
    theme( legend.position = 'none',
      plot.title = element_text(size = 8, hjust = 0.5, face = "italic"),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.text = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.2, colour = "grey80")
    ) +
    scale_x_continuous(breaks = seq(1950, 2020, by = 20),
                       limits = c(1950, 2020),
                       expand = expansion(mult = c(0.02, 0.02)))
  return(p)
})

plot_list_modified <- plot_list_modified[!sapply(plot_list_modified, is.null)]
combined_plot <- wrap_plots(plot_list_modified, ncol = 5)
combined_plot[[1]]


ggsave("./Plots/06 logistic_fits.svg", 
       plot = combined_plot, 
       width = 3 * 3.5,      # 4 inches per column
       height = 6 * 2,     # 3 inches per row
       dpi = 300)

# to check 
plot_list_modified[[27]] + ylim(0,50)
plot_list_modified[[28]] + ylim(0,50)
plot_list_modified[[31]] + ylim(0,75)

#  Figure S2  ----------------------------------


# use sp1 from figure 3

sp1 <- sp1[sp1$analysis_group=="All_combined",]
sp1
custom_labels
cols 

c<-cols[3]
plots<- list()

for(c in cols){
  print(c)
  # Sort the dataframe by Group column instead
  sp1_sorted <- sp1 %>% arrange(Group) 
  library(scales)
  
  p <- ggplot(sp1_sorted, aes(x = Group, 
                              y = pmax(.data[[c]], 0.1),  # Replace zeros with small value
                              fill = Group, color = Group)) +   
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # boxplots, no outlier dots (to avoid huge spikes)
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = group_colors) +   scale_y_log10(
      breaks = 10^(0:6),
      labels = scales::label_number(accuracy = 1, big.mark = ",")
    )+
    scale_color_manual(values = group_colors) + 
    labs(x = "", y = custom_labels[c])+
    theme_bw(base_size = 12) + theme(legend.position = "none") +
   # scale_y_continuous(breaks = pretty_breaks(n = 9) ) +
    #scale_y_continuous(
    # trans = "log10",
    #breaks = c( 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000,10000000),
    #labels = c( "0.01", "0.1", "1", "10", "100", "1,000", "10,000", "100,000", "1,000,000","10000000")
    #) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
  plots[[c]] <- p
  # ggsave(filename = paste0("Figure2_", c, ".svg"), plot = p,
  #       width = 8, height = 6, units = "in")
}
names(plots)
combined_plots <- (
    plots[["N.cells"]] +
  plots[["max_area_km2"]] + 
    plots[["max_invasion_front_km"]] + 
    plots[["mean_spread_rate"]] 
) + plot_layout(ncol = 2, nrow = 2)

ggsave(
  filename = "S02 summary.svg",  
  plot = combined_plots,            
  width = 15,                       
  height = 12,                      
  dpi = 300 )













