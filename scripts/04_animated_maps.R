library(tidyverse)
library(ggOceanMaps)
library(ggspatial)
library(oce)
library(rnaturalearth)
library(gganimate)

# Load data
load("data/derived_data/diva3d.rda")

# Get land polygons from Natural Earth
land_polygons <- ne_countries(scale = "medium", returnclass = "sf")

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
anim_save("product/animations/zooplankton_diversity_animation.webm", animation)
