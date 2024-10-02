library(JuliaCall)
julia_setup(JULIA_HOME = "/Applications/Julia-1.10.app/Contents/Resources/julia/bin")
julia_install_package("DIVAnd")
julia_command("using DIVAnd")
library(logger)
library(tidyverse)
library(ggOceanMaps)
library(ggspatial)
library(rnaturalearth)
library(sf)
library(oce)

load("data/derived_data/data_diversity.rda")
load("data/derived_data/data_diversity_all.rda")

# Get land polygons from Natural Earth at 1:10m scale
land_polygons <- ne_countries(scale = "medium", returnclass = "sf")

subset <- data_diversity %>%
  filter(year(data_diversity$eventDate) == 2010) %>%  group_by(verbatimLocality) %>%
  summarise(shannon = mean(shannon, na.rm = TRUE),
            uniqueTaxa = mean(uniqueTaxa, na.rm = TRUE),
            decimalLongitude = mean(decimalLongitude, na.rm = TRUE),
            decimalLatitude = mean(decimalLatitude, na.rm = TRUE))

subset_all <- diversity_all %>%
  filter(year(eventDate) == 2010) 
# %>%  #group_by(verbatimLocality) %>%
#   summarise(shannon = mean(shannon, na.rm = TRUE),
#             uniqueTaxa = mean(uniqueTaxa, na.rm = TRUE),
#             decimalLongitude = mean(decimalLongitude, na.rm = TRUE),
#             decimalLatitude = mean(decimalLatitude, na.rm = TRUE))

# Create basemap
baltic_sea_map <- basemap(
  limits = c(min(subset_all$decimalLongitude) - 1, max(subset_all$decimalLongitude) + 1, min(subset_all$decimalLatitude) - 1, max(subset_all$decimalLatitude) + 1),
  land.col = "#eeeac4",
  land.border.col = "black",
  rotate = TRUE,
  bathymetry = FALSE
)

# Plot points on map
map <- baltic_sea_map +
  geom_spatial_point(
    data = subset_all,
    aes(x = decimalLongitude, y = decimalLatitude, size = shannon),
    pch = 21,
    # size = 2,
    fill = "red",
    colour = "black"
  ) +
  scale_size_continuous(range = c(1, 10), name = "Shannon index") +
  ggtitle("Average Zooplankton shannon 2010")

x <- subset_all$decimalLongitude
y <- subset_all$decimalLatitude
f <- subset_all$shannon

# Set the interpolation grid
NX <- 500
NY <- 550
xx <- seq(9.5, 30, length.out=NX)
yy <- seq(53.5, 66, length.out=NY)
julia_assign("xx", xx)
julia_assign("yy", yy)
julia_command("xi, yi = ndgrid(xx, yy)");
xi = julia_eval("xi")
yi = julia_eval("yi")

# Load land polygons from Natural Earth
land_polygons <- ne_countries(scale = "medium", returnclass = "sf")

# Create a grid of spatial points from xx and yy
grid_points <- expand.grid(x = xx, y = yy)

# Convert grid to sf object with the correct CRS
grid_points_sf <- st_as_sf(grid_points, coords = c("x", "y"), crs = st_crs(land_polygons))

# Check which grid points intersect the land polygons
land_mask <- st_intersects(grid_points_sf, land_polygons, sparse = FALSE)

# 'land_mask' is a logical vector indicating which points are on land (TRUE means it's on land)
land_mask <- apply(land_mask, 1, any)  # Collapse the logical matrix to a single vector

# Apply the additional condition: Latitude > 63 and Longitude < 15
additional_mask <- with(grid_points, y > 63 & x < 15)

# Combine the two masks: land and the additional conditions
combined_mask <- land_mask | additional_mask  # TRUE if it's either on land or meets the condition

# Reshape the combined_mask vector into a matrix with dimensions NX by NY
mask <- matrix(!combined_mask, nrow = NX, ncol = NY, byrow = FALSE)  # Points over water and outside the condition are TRUE

# Metrics: pm is the inverse of the resolution along the 1st dimension
pm <- matrix(1, NX, NY) / (xi[2,1]-xi[1,1])
pn <- matrix(1, NX, NY) / (xi[2,1]-xi[1,1])

# Analysis parameters
# Correlation length
len <- 10

# Obs. error variance normalized by the background error variance
epsilon2 <- 1.0

# From R variable to Julia variable
# (maybe possible to do it differently)
julia_assign("pm", pm)
julia_assign("pn", pn)
julia_assign("xi", xi)
julia_assign("yi", yi)
julia_assign("x", x)
julia_assign("y", y)
julia_assign("f", f)
julia_assign("len", len)
julia_assign("mask", mask)
julia_assign("epsilon2", epsilon2)

# DIVAnd execution
julia_command("fi, s = DIVAndrun(mask, (pm, pn) ,(xi, yi), (x, y), f, len, epsilon2; alphabc=2);") 

# From Julia variable to R variable
fi = julia_eval("fi")


# Load land polygons from Natural Earth
land_polygons <- ne_countries(scale = "medium", returnclass = "sf")

# Convert your interpolation result into a data frame for ggplot2
interp_grid <- expand.grid(lon = xx, lat = yy)
interp_grid$fi <- as.vector(fi)

# Define color palette and limits
pal <- oce.colorsTemperature(100)  # Generate the color palette with 100 color breaks
zlim <- range(fi, na.rm = TRUE)

# Plot the land polygons with ggplot2
interpolated_map <- ggplot() +
  
  # Add the interpolated field as a raster layer
  geom_raster(data = interp_grid, aes(x = lon, y = lat, fill = fi)) +
  
  # Add land polygons from rnaturalearth
  geom_sf(data = land_polygons, fill = "#eeeac4", color = "black") +
  
  # Set color scale for interpolation
  scale_fill_gradientn(colors = pal, limits = zlim, na.value = "transparent", name = "Shannon index") +
  
  # Adjust the plot limits based on your longitude and latitude ranges
  coord_sf(xlim = range(xx), ylim = range(yy), expand = FALSE) +
  
  # Add labels and titles
  labs(title = "Zooplankton diversity",
       x = "Longitude", y = "Latitude") +
  
  # Use a minimal theme
  theme_minimal() +
  geom_spatial_point(
    data = subset_all,
    aes(x = decimalLongitude, y = decimalLatitude),
    pch = 21,
    # size = 2,
    fill = "white",
    colour = "black"
  ) 

ggsave("product/interpolated_map4.png",
       interpolated_map,
       device = "png",
       bg = "white")
