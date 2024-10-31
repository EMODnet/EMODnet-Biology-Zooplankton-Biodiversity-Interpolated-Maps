library(tidyverse)
library(ggOceanMaps)
library(sf)
library(sp)
library(geosphere)
library(ggspatial)

# Load data
load("data/derived_data/data_diversity_all.rda")

# Define file paths for the shapefiles and basin names
shapefilesDir <- "data/raw_data/shapefiles/sharkweb_shapefiles/"          # Directory containing shapefiles
basin_shapefile <- "Havsomr_SVAR_2016_3b_CP1252.shp"                      # Main shapefile for sea basins
basin_names <- "sea_basin_utf8.txt"                                       # Text file with basin names

# Read shapefiles and list of basin names
basins <- st_read(file.path(shapefilesDir, basin_shapefile))              # Load basin shapefiles using sf package
basin_names <- read_delim(file.path(shapefilesDir, basin_names),          # Read basin names as a UTF-8 encoded table
                          delim = ";", 
                          col_names = TRUE, 
                          locale = locale(encoding = "UTF-8"))

# Create basemap
baltic_sea_map <- basemap(
  limits = c(min(diversity_all$decimalLongitude) - 1, max(diversity_all$decimalLongitude) + 1, min(diversity_all$decimalLatitude) - 1, max(diversity_all$decimalLatitude) + 1),
  land.col = "#eeeac4",
  land.border.col = "black",
  rotate = TRUE,
  bathymetry = FALSE
)
# Define WGS84 coordinate reference system
proWG = CRS("+proj=longlat +datum=WGS84")

# Select unique coordinates from diversity data
cooridinates <- diversity_all %>%
  dplyr::select(decimalLongitude, decimalLatitude) %>%  # Select longitude and latitude columns
  distinct()  # Keep only unique coordinate pairs

# Create spatial points data frame for clustering
xy <- SpatialPointsDataFrame(
  matrix(c(cooridinates$decimalLongitude, cooridinates$decimalLatitude), ncol = 2),  # Create matrix of coordinates
  data.frame(ID = seq(1:nrow(cooridinates))),  # Data frame with ID for each coordinate
  proj4string = proWG  # Assign WGS84 coordinate system
)

# Calculate distance matrix between points
mdist <- distm(xy)  # Compute distance matrix in meters using geosphere package

# Perform hierarchical clustering on the distance matrix
hc = hclust(as.dist(mdist), method = "complete")

# Define the maximum distance for clustering nearby sample together
d <- 20000  # Set distance threshold (in meters) for clustering nearby stations

# Apply hierarchical clustering cut based on the distance threshold
xy$clust <- cutree(hc, h = d) 

# Add cluster assignments to the coordinates data frame
cooridinates <- cooridinates %>%
  mutate("stationCluster" = xy$clust)

# Merge cluster information with the main diversity dataset
diversity_all = diversity_all %>%
  left_join(cooridinates, by = c("decimalLongitude", "decimalLatitude"))

# Calculate mean
mean_diversity <- diversity_all %>%
  mutate(position = paste(decimalLatitude, decimalLongitude)) %>%
  group_by(stationCluster) %>%
  summarise(shannon = mean(shannon, na.rm = TRUE),
            uniqueTaxa = mean(uniqueTaxa, na.rm = TRUE),
            decimalLongitude = mean(decimalLongitude, na.rm = TRUE),
            decimalLatitude = mean(decimalLatitude, na.rm = TRUE))


# Set CRS of basin layer
basins <- st_set_crs(basins, 3006)

# Aggregate basins by the 17 sea basins
all_basins <- basins %>%
  group_by(BASIN_NR) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Change CRS
all_basins <- st_transform(all_basins, 4326)

# Add geometry information to data
data_diversity_coordinates <- diversity_all %>%
  mutate(lon = decimalLongitude,
         lat = decimalLatitude)

# Gather all unique positions
cords = data_diversity_coordinates %>%
  distinct(decimalLongitude, decimalLatitude, lat, lon)

# Convert data points to sf
points_sf <- st_as_sf(cords, coords = c("lon", "lat"), crs = st_crs(all_basins))

# Assign basin number by position
data_diversity_st <- st_join(points_sf, all_basins)

# Add sea basin name and create translate list
data_diversity_st <- data_diversity_st %>%
  as.data.frame() %>%
  left_join(basin_names) %>%
  select(-geometry, -BASIN_NR)

# Join with basin names
diversity_all <- diversity_all %>%
  left_join(data_diversity_st)

# Plot points on map
map <- baltic_sea_map +
  geom_spatial_point(
    data = mean_diversity,
    aes(x = decimalLongitude, y = decimalLatitude, size = uniqueTaxa),
    pch = 21,
    # size = 2,
    fill = "red",
    colour = "black"
  ) +
  scale_size_continuous(range = c(1, 10), name = "Species richness") +
  ggtitle("Average Zooplankton richness 2000-2022")

# Plot time series of abundance at each location
time_series <- ggplot(diversity_all, aes(x = eventDate, y = shannon)) +
  geom_point(aes(color = factor(location_sea_basin_en))) +
  scale_color_discrete() +
  labs(x = "Date", y = "Shannon index", title = "Zooplankton diversity", color = "Sea basin") +
  theme_minimal()

# Save time series
ggsave(
  plot = time_series,
  path = "product/plots/",
  filename = "times_series_all.png",
  device = "png",
  dpi = 300,
  width = 7,
  height = 4,
  bg = "white"
)

# Save map
ggsave(
  plot = map,
  path = "product/maps/",
  filename = "map_all.png",
  device = "png",
  dpi = 300,
  width = 7,
  height = 7,
  bg = "white"
)
