library(JuliaCall)
library(logger)
library(tidyverse)
library(ggOceanMaps)
library(ggspatial)
library(rnaturalearth)
library(sf)
library(oce)

# Set up Julia
julia_setup(JULIA_HOME = "/Applications/Julia-1.10.app/Contents/Resources/julia/bin")
julia_install_package("DIVAnd")
julia_command("using DIVAnd")

# Load data
load("data/derived_data/data_diversity.rda")
load("data/derived_data/data_diversity_all.rda")

# Get land polygons from Natural Earth at 1:10m scale
land_polygons <- ne_countries(scale = "medium", returnclass = "sf")

# Subset samples from year 2000 and forward
subset_all <- diversity_all %>%
  filter(year(eventDate) >= 2000)

# Function to classify seasons
get_season <- function(date) {
  month <- month(date)
  
  if (month %in% c(3, 4, 5)) {
    return("Spring")
  } else if (month %in% c(6, 7, 8)) {
    return("Summer")
  } else if (month %in% c(9, 10, 11)) {
    return("Autumn")
  } else {
    return("Winter")
  }
}

# Apply the function to the dataset
subset_all <- subset_all %>%
  mutate(season = sapply(eventDate, get_season))

# Function to convert season to a numeric value with consideration for Winter overlap
get_season_numeric <- function(date) {
  month <- month(date)
  
  if (month %in% c(3, 4, 5)) {
    return(1)  # Spring
  } else if (month %in% c(6, 7, 8)) {
    return(2)  # Summer
  } else if (month %in% c(9, 10, 11)) {
    return(3)  # Autumn
  } else {
    if (month == 12) {
      return(4)  # Winter (December, stay in current year)
    } else {
      return(0)  # Winter (January, February, go to previous year)
    }
  }
}

# Create a combined numeric time variable with proper handling of Winter
subset_all <- subset_all %>%
  mutate(year = as.numeric(format(eventDate, "%Y")),        # Extract and convert year to numeric
         season_num = sapply(eventDate, get_season_numeric), # Convert season to numeric, handle Winter overlap
         time_numeric = ifelse(season_num == 0, (year - 1) + 0.75, 
                               year + (season_num - 1) / 4)) # Correct year for January-February Winter

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

# Save plot
ggsave("product/sample_map.png",
       map,
       device = "png",
       bg = "white")

x <- subset_all$decimalLongitude
y <- subset_all$decimalLatitude
t <- subset_all$time_numeric
f <- subset_all$shannon

# Set the interpolation grid
NX <- 100
NY <- 110
# NT <- length(unique(t))

xx <- seq(9.5, 30, length.out=NX)
yy <- seq(53.5, 66, length.out=NY)
tt <- sort(unique(subset_all$time_numeric))  # Time grid in days (or months, etc.)

NT <- length(tt)

julia_assign("xx", xx)
julia_assign("yy", yy)
julia_assign("tt", tt)

julia_command("xi, yi, ti = ndgrid(xx, yy, tt)")
xi = julia_eval("xi")
yi = julia_eval("yi")
ti <- julia_eval("ti")

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

mask_3d <- array(rep(mask, times = NT), dim = c(NX, NY, NT))

# Calculate the resolution (spacing) for longitude, latitude, and time
dx <- xi[2,1,1] - xi[1,1,1]  # Longitude resolution
dy <- yi[1,2,1] - yi[1,1,1]  # Latitude resolution
dt <- ti[1,1,2] - ti[1,1,1]  # Time resolution (assuming ti is the correct grid)

# Metrics: pm and pn for spatial dimensions, pt for time
pm <- matrix(1, NX, NY) / dx  # Longitude metric
pn <- matrix(1, NX, NY) / dy  # Latitude metric
pt <- 1.0 / dt  # Time metric

# Convert pm, pn, pt into 3D arrays to match the 3D grid
pm_3d <- array(pm, dim = c(NX, NY, NT))
pn_3d <- array(pn, dim = c(NX, NY, NT))
pt_3d <- array(pt, dim = c(NX, NY, NT))

# Assign metrics to Julia
julia_assign("pm", pm_3d)
julia_assign("pn", pn_3d)
julia_assign("pt", pt_3d)

# Analysis parameters
# Correlation length
len <- 5
len_time <- .25

# Obs. error variance normalized by the background error variance
epsilon2 <- 1.0

# From R variable to Julia variable
# (maybe possible to do it differently)
julia_assign("xi", xi)
julia_assign("yi", yi)
julia_assign("x", x)
julia_assign("y", y)
julia_assign("t", t)
julia_assign("f", f)
julia_assign("len", len)
julia_assign("len_time", len_time)
julia_assign("mask", mask_3d)
julia_assign("epsilon2", epsilon2)

# Run the DIVAnd interpolation in 3D (longitude, latitude, time)
julia_command("fi, s = DIVAndrun(mask, (pm, pn, pt), (xi, yi, ti), (x, y, t), f, (len, len, len_time), epsilon2);")

# From Julia variable to R variable
fi = julia_eval("fi")

# Clip the interpolated values to be non-negative
fi <- pmax(fi, 0)

save(xx, yy, tt, xi, yi, ti, fi, file = "data/derived_data/diva3d.rda")

# Select a specific time slice (e.g., the first time step)
time_index <- 66  # Or select any other index to visualize different time steps
fi_slice <- fi[,, time_index]  # Extract the 2D slice of the interpolated values for the chosen time step

# Create a mapping of time_step to year and season based on tt
time_labels <- data.frame(
  time_step = 1:NT,
  year = floor(tt),
  season = factor(
    (tt - floor(tt)) * 4 + 1,
    labels = c("Spring", "Summer", "Autumn", "Winter")
  )
)

# Create a combined year and season numeric column (time_numeric)
time_labels$time_numeric <- tt
time_labels$year_season <- paste(time_labels$year, time_labels$season)

# Extract title
title <- time_labels$year_season[time_index]

# Convert the slice to a data frame for plotting
interp_grid <- expand.grid(lon = xx, lat = yy)
interp_grid$fi <- as.vector(fi_slice)

# Define color palette and limits
pal <- oce.colorsTemperature(100)  # Generate the color palette with 100 color breaks
zlim <- range(fi, na.rm = TRUE)

# Plot the interpolated map for a specific time slice
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
  labs(title = paste0("Zooplankton diversity: ", title),
       x = "Longitude", y = "Latitude") +
  
  # Use a minimal theme
  theme_minimal() 

ggsave("product/interpolated_map.png",
       interpolated_map,
       device = "png",
       bg = "white")

# Load land polygons from Natural Earth
land_polygons <- ne_countries(scale = "medium", returnclass = "sf")

# Convert your 3D interpolation result (fi) into a long format for ggplot2
interp_df <- expand.grid(lon = xx, lat = yy, time_step = 1:NT)
interp_df$fi <- as.vector(fi)  # Flatten the 3D array into a 1D vector


# Merge the time labels into the interpolation data frame
interp_df <- merge(interp_df, time_labels, by = "time_step")

# Define color palette and limits
pal <- oce.colorsTemperature(100)  # Generate the color palette with 100 color breaks
zlim <- range(fi, na.rm = TRUE)

# Create the animated plot using gganimate
animated_map <- ggplot() +
  
  # Add the interpolated field as a raster layer
  geom_raster(data = interp_df, aes(x = lon, y = lat, fill = fi)) +
  
  # Add land polygons from rnaturalearth
  geom_sf(data = land_polygons, fill = "#eeeac4", color = "black") +
  
  # Set color scale for interpolation
  scale_fill_gradientn(
    colors = pal, limits = zlim, na.value = "transparent", name = "Shannon index"
  ) +
  
  # Adjust the plot limits based on your longitude and latitude ranges
  coord_sf(xlim = range(xx), ylim = range(yy), expand = FALSE) +
  
  # Add labels and dynamic title with year and season
  labs(
    title = 'Zooplankton diversity: {time_labels$year[which.min(abs(time_labels$time_numeric - frame_time))]} {time_labels$season[which.min(abs(time_labels$time_numeric - frame_time))]}',
    x = "Longitude", y = "Latitude"
  ) +
  
  # Use a minimal theme
  theme_minimal() +
  
  # Animate over the time_numeric variable using transition_time
  transition_time(time_numeric) +
  ease_aes('linear')

# Animate with adjusted fps and duration
animation <- animate(
  animated_map,
  fps = 5,          # Decrease fps to slow down animation playback
  duration = 60,    # Increase duration to slow down the animation
  renderer = ffmpeg_renderer(format = "webm")
)

# Save the animation as a WebM video
anim_save("product/zooplankton_diversity_animation.webm", animation)
