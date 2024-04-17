library(JuliaCall)
# julia_setup(JULIA_HOME = path.expand("~/.juliaup/bin/"))
julia_command("using DIVAnd")
library(logger)
library(tidyverse)
library(ggOceanMaps)
library(ggspatial)


load("data/derived_data/data_diversity.rda")

subset <- data_diversity %>%
  filter(year(data_diversity$eventDate) == 2010) %>%  group_by(verbatimLocality) %>%
  summarise(shannon = mean(shannon, na.rm = TRUE),
            uniqueTaxa = mean(uniqueTaxa, na.rm = TRUE),
            decimalLongitude = mean(decimalLongitude, na.rm = TRUE),
            decimalLatitude = mean(decimalLatitude, na.rm = TRUE))

# Create basemap
baltic_sea_map <- basemap(
  limits = c(min(subset$decimalLongitude) - 1, max(subset$decimalLongitude) + 1, min(subset$decimalLatitude) - 1, max(subset$decimalLatitude) + 1),
  land.col = "#eeeac4",
  land.border.col = "black",
  rotate = TRUE,
  bathymetry = FALSE
)


# Plot points on map
map <- baltic_sea_map +
  geom_spatial_point(
    data = subset,
    aes(x = decimalLongitude, y = decimalLatitude, size = shannon),
    pch = 21,
    # size = 2,
    fill = "red",
    colour = "black"
  ) +
  scale_size_continuous(range = c(1, 10), name = "Shannon index") +
  ggtitle("Average Zooplankton shannon 2010")


x <- subset$decimalLongitude
y <- subset$decimalLatitude
f <- subset$shannon

# Set the interpolation grid
NX <- 100
NY <- 110
xx <- seq(10, 24, length.out=NX)
yy <- seq(55, 66, length.out=NY)
julia_assign("xx", xx)
julia_assign("yy", yy)
julia_command("xi, yi = ndgrid(xx, yy)");
xi = julia_eval("xi")
yi = julia_eval("yi")

# Mask: all points are valid points
mask <- matrix(TRUE, NX, NY)

# Metrics: pm is the inverse of the resolution along the 1st dimension
pm <- matrix(1, NX, NY) / (xi[2,1]-xi[1,1])
pn <- matrix(1, NX, NY) / (xi[2,1]-xi[1,1])

# Analysis parameters
# Correlation length
len <- 0.1

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
