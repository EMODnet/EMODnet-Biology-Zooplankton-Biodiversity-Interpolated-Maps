library(tidyverse)
library(finch)
library(vegan)
library(ggOceanMaps)
library(sf)


shapefilesDir <- "data/raw_data/shapefiles/sharkweb_shapefiles/"
basin_shapefile <- "Havsomr_SVAR_2016_3b_CP1252.shp"
basin_names <- "sea_basin_utf8.txt"

# Read shapefiles and list of basin names
basins <- st_read(file.path(shapefilesDir, basin_shapefile))
basin_names <- read_delim(file.path(shapefilesDir, basin_names), 
                          delim = ";", 
                          col_names = TRUE, 
                          locale = locale(encoding = "UTF-8"))



zooplankton_dwca <- dwca_read("data/raw_data/dwca-shark-zooplankton-nat-v1.8.zip")

event <- read_tsv(zooplankton_dwca$files$txt_files[1])
occurrence <- read_tsv(zooplankton_dwca$files$txt_files[2])
emof <- read_tsv(zooplankton_dwca$files$txt_files[3])

event_sel <- event %>%
  select(id, eventID, parentEventID, eventDate, verbatimLocality, minimumDepthInMeters, maximumDepthInMeters, decimalLatitude, decimalLongitude)

occurrence_sel <- occurrence %>%
  select(id, occurrenceID, measurementType, measurementValue, measurementUnit)

all <- occurrence_sel %>%
  left_join(event_sel) %>%
  left_join(emof)

abundance <- all %>%
  filter(measurementType == "Abundance") %>%
  group_by(id) %>%
  mutate(sample_sum = sum(as.numeric(measurementValue))) %>%
  ungroup()

# Create basemap
baltic_sea_map <- basemap(
  limits = c(min(abundance$decimalLongitude) - 1, max(abundance$decimalLongitude) + 1, min(abundance$decimalLatitude) - 1, max(abundance$decimalLatitude) + 1),
  land.col = "#eeeac4",
  land.border.col = "black",
  rotate = TRUE,
  bathymetry = FALSE
)

# Plot points on map
map <- baltic_sea_map +
  geom_spatial_point(
    data = abundance,
    aes(x = decimalLongitude, y = decimalLatitude, size = sample_sum),
    pch = 21,
    # size = 2,
    fill = "red",
    colour = "black"
  ) +
  scale_size_continuous(range = c(1, 10), name = "Abundance (ind/m3)")


# Calculate alpha diversity
data_diversity <- abundance %>%  
  group_by(id) %>%
  summarise(uniqueTaxa = length(unique(scientificName)),
            shannon = diversity(as.numeric(measurementValue))) %>%
  left_join(event)

# Calculate mean
mean_diversity <- data_diversity %>%
  group_by(verbatimLocality) %>%
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
data_diversity_coordinates <- data_diversity %>%
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
data_diversity <- data_diversity %>%
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
  ggtitle("Average Zooplankton richness 1979-2022")


# Plot time series of abundance at each location
time_series <- ggplot(data_diversity, aes(x = eventDate, y = shannon)) +
  geom_point(aes(color = factor(location_sea_basin_en))) +
  scale_color_discrete() +
  labs(x = "Date", y = "Shannon index", title = "Zooplankton diversity", color = "Sea basin") +
  theme_minimal()

# Save map
ggsave(
  plot = time_series,
  path = "output/",
  filename = "times_series.png",
  device = "png",
  dpi = 300,
  width = 7,
  height = 4
)

# Save map
ggsave(
  plot = map,
  path = "output/",
  filename = "map.png",
  device = "png",
  dpi = 300,
  width = 7,
  height = 7
)
