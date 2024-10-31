library(JuliaCall)
library(tidyverse)
library(rnaturalearth)
library(sf)

# Set up Julia
julia_setup(JULIA_HOME = "/Applications/Julia-1.10.app/Contents/Resources/julia/bin")
julia_install_package("DIVAnd")
julia_command("using DIVAnd")

# Load data
# load("data/derived_data/data_diversity.rda")
load("data/derived_data/data_diversity_all.rda")

# Get land polygons from Natural Earth
land_polygons <- ne_countries(scale = "medium", returnclass = "sf")

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
diversity_all <- diversity_all %>%
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
diversity_all <- diversity_all %>%
  mutate(year = as.numeric(format(eventDate, "%Y")),        # Extract and convert year to numeric
         season_num = sapply(eventDate, get_season_numeric), # Convert season to numeric, handle Winter overlap
         time_numeric = ifelse(season_num == 0, (year - 1) + 0.75, 
                               year + (season_num - 1) / 4)) # Correct year for January-February Winter

# Assign the dimensions
x <- diversity_all$decimalLongitude
y <- diversity_all$decimalLatitude
t <- diversity_all$time_numeric
f <- diversity_all$shannon

# Set the interpolation grid
NX <- 100
NY <- 110
NT <- length(unique(t))

xx <- seq(9.5, 30, length.out=NX)
yy <- seq(53.5, 66, length.out=NY)
tt <- sort(unique(diversity_all$time_numeric))

# Function to convert decimal year to days since 1970 and round
decimal_year_to_days <- function(decimal_year) {
  year <- floor(decimal_year)
  fraction <- decimal_year - year
  
  # Calculate the date corresponding to the decimal year
  date <- as.Date(paste0(year, "-01-01")) + round(fraction * 365)  # Rounding to nearest day
  
  # Calculate days since January 1, 1970
  days_since_epoch <- as.numeric(difftime(date, as.Date("1970-01-01"), units = "days"))
  
  return(days_since_epoch)
}

# Apply the function to the tt vector
tt_days <- sapply(tt, decimal_year_to_days)

# Assign values in julia
julia_assign("xx", xx)
julia_assign("yy", yy)
julia_assign("tt", tt)
# julia_assign("tt_days", tt_days)

julia_command("xi, yi, ti = ndgrid(xx, yy, tt)")
xi = julia_eval("xi")
yi = julia_eval("yi")
ti <- julia_eval("ti")

# # Same operations again, but expressed as days
# julia_command("xi, yi, ti_days = ndgrid(xx, yy, tt_days)")
# ti_days <- julia_eval("ti_days")

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

# Save the DIVAnd output
save(xx, yy, tt, xi, yi, ti, fi, tt_days, NT, file = "data/derived_data/diva3d.rda")
