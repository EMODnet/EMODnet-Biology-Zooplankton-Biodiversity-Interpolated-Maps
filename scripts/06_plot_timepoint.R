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

# Select a specific time slice (e.g., the first time step)
time_index <- 66  # Or select any other index to visualize different time steps
e_slice <- e[,, time_index]  # Extract the 2D slice of the error values for the chosen time step
e_slice[e_slice > 21.65] <- NA # Mask values above threshold

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
interp_grid$fi <- as.vector(e_slice)

# Define color palette and limits
pal <- oce.colorsTemperature(100)  # Generate the color palette with 100 color breaks
zlim <- range(e, na.rm = TRUE)

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
