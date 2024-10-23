library(RNetCDF)
library(JuliaCall)

load("data/derived_data/diva3d.rda")

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
tt <- sapply(tt, decimal_year_to_days)

julia_assign("xx", xx)
julia_assign("yy", yy)
julia_assign("tt", tt)

julia_command("xi, yi, ti = ndgrid(xx, yy, tt)")

xi = julia_eval("xi")
yi = julia_eval("yi")
ti <- julia_eval("ti")

# Create nc file
nc <- create.nc(file.path("product","alpha_diversity.nc")) 

## Define dimensions of the NetCDF

# Define lon dimension
dim.def.nc(nc, dimname = "lon", dimlength = length(xx))

# Define lon variable
var.def.nc(nc, varname = "lon", vartype = "NC_DOUBLE", dimensions = "lon")

# Add attributes
att.put.nc(nc, variable = "lon", name = "units", type = "NC_CHAR", value = "degrees_east")
att.put.nc(nc, variable = "lon", name = "standard_name", type = "NC_CHAR", value = "longitude")
att.put.nc(nc, variable = "lon", name = "long_name", type = "NC_CHAR", value = "Longitude")

# Put data
var.put.nc(nc, variable = "lon", data = xx)

# Check
paste("longitude")
var.get.nc(nc, variable = "lon")

### Latitude

# Define lat dimension
dim.def.nc(nc, dimname = "lat", dimlength = length(yy))

# Define lat variable
var.def.nc(nc, varname = "lat", vartype = "NC_DOUBLE", dimensions = "lat")

# Add attributes
att.put.nc(nc, variable = "lat", name = "units", type = "NC_CHAR", value = "degrees_north")
att.put.nc(nc, variable = "lat", name = "standard_name", type = "NC_CHAR", value = "latitude")
att.put.nc(nc, variable = "lat", name = "long_name", type = "NC_CHAR", value = "Latitude")

# Put data
var.put.nc(nc, variable = "lat", data = yy)

# Check
paste("latitudes")
var.get.nc(nc, variable = "lat")

### Time

# Define time dimension
dim.def.nc(nc, dimname = "time", dimlength = length(tt)) 

# Define time variable
var.def.nc(nc, varname = "time", vartype = "NC_DOUBLE", dimensions = "time")

# Add attributes
att.put.nc(nc, variable = "time", name = "standard_name", type = "NC_CHAR", value = "time")
att.put.nc(nc, variable = "time", name = "long_name", type = "NC_CHAR", value = "Time")
att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "days since 1970-01-01 00:00:00")
att.put.nc(nc, variable = "time", name = "calendar", type = "NC_CHAR", value = "gregorian")

# Put data
var.put.nc(nc, variable = "time", data = tt)

# Check
paste("timepoints")
var.get.nc(nc, variable = "time")

### Define non-dimensional crs variable 
var.def.nc(nc, varname = "crs", vartype = "NC_CHAR", dimensions = NA)

# Add attributes
att.put.nc(nc, variable = "crs", name = "long_name", type = "NC_CHAR", value = "Coordinate Reference System")
att.put.nc(nc, variable = "crs", name = "geographic_crs_name", type = "NC_CHAR", value = "WGS 84")
att.put.nc(nc, variable = "crs", name = "grid_mapping_name", type = "NC_CHAR", value = "latitude_longitude")
att.put.nc(nc, variable = "crs", name = "reference_ellipsoid_name", type = "NC_CHAR", value = "WGS 84")
att.put.nc(nc, variable = "crs", name = "horizontal_datum_name", type = "NC_CHAR", value = "WGS 84")
att.put.nc(nc, variable = "crs", name = "prime_meridian_name", type = "NC_CHAR", value = "Greenwich")
att.put.nc(nc, variable = "crs", name = "longitude_of_prime_meridian", type = "NC_DOUBLE", value = 0.)
att.put.nc(nc, variable = "crs", name = "semi_major_axis", type = "NC_DOUBLE", value = 6378137.)
att.put.nc(nc, variable = "crs", name = "semi_minor_axis", type = "NC_DOUBLE", value = 6356752.314245179)
att.put.nc(nc, variable = "crs", name = "inverse_flattening", type = "NC_DOUBLE", value = 298.257223563)
att.put.nc(nc, variable = "crs", name = "spatial_ref", type = "NC_CHAR", value = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]')
att.put.nc(nc, variable = "crs", name = "GeoTransform", type = "NC_CHAR", value = '-180 0.08333333333333333 0 90 0 -0.08333333333333333 ')

## Define variables and store the NC file

# Create the shannon variable defined by the four dimensions
var.def.nc(nc, varname = "shannon", vartype = "NC_DOUBLE", dimensions = c("lon", "lat", "time"))

# Add attributes
att.put.nc(nc, variable = "shannon", name = "_FillValue", type = "NC_DOUBLE", value = -99999)
att.put.nc(nc, variable = "shannon", name = "long_name", type = "NC_CHAR", value = "Shannon index")

var.put.nc(nc, variable = "shannon", data = fi) 

## Global attributes and save nc file

attributes <- list(
  title = "Zooplankton alpha diversity in the greater Baltic Sea area",
  summary = "This dataset compiles monthly Shannon diversity index (H') and species richness (n taxa) of phytoplankton data in the greater Baltic Sea area between years 2000-2021, calculated from abundance data originating from the Swedish National phytoplankton monitoring programme. Samples have been rarified before calculating alpha diversity measures, and nearby spatial data points (stations) have been clustered together",                       
  Conventions = "CF-1.8",
  naming_authority = "emodnet-biology.eu",
  history = "https://www.vliz.be/imis?dasid=XXXX",
  source = "https://www.vliz.be/imis?dasid=XXXX",
  license = "CC-BY",
  standard_name_vocabulary = "CF Standard Name Table v1.8",
  date_created = as.character(Sys.Date()),
  creator_name = "Anders Torstensson",
  creator_email = "anders.torstensson@smhi.se",
  creator_url = "www.smhi.se",
  institution = "Swedish Meteorological and Hydrological Institute",
  project = "EMODnet-Biology",
  publisher_name = "EMODnet-Biology",                 
  publisher_email = "bio@emodnet.eu",                
  publisher_url = "www.emodnet-biology.eu",                  
  geospatial_lat_min = min(yy),
  geospatial_lat_max = max(yy),
  geospatial_lon_min = min(xx),
  geospatial_lon_max = max(xx),
  creator_institution = "Swedish Meteorological and Hydrological Institute (SMHI)",            
  publisher_institution = "Swedish Meteorological and Hydrological Institute (SMHI)",        
  geospatial_lat_units = "degrees_north",           
  geospatial_lon_units = "degrees_east",           
  comment = "Uses attributes recommended by http://cfconventions.org",
  license = "CC-BY", 
  publisher_name = "EMODnet Biology Data Management Team",
  citation = "Torstensson A (2024) Phytoplankton alpha diversity in the greater Baltic Sea area. Integrated data products created under the European Marine Observation  Data Network (EMODnet) Biology project Phase IV (EMFF/2019/1.3.1.9/Lot  6/SI2.837974), funded by the by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund",
  acknowledgement = "European Marine Observation Data Network (EMODnet) Biology project (EMFF/2019/1.3.1.9/Lot 6/SI2.837974), funded by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund"
)

# Define function that detects if the data type should be character of 
# integer and add to global attributes
add_global_attributes <- function(nc, attributes){
  
  stopifnot(is.list(attributes))
  
  for(i in 1:length(attributes)){
    if(is.character(attributes[[i]])){
      type <- "NC_CHAR"
    }else if(is.numeric(attributes[[i]])){
      type <- "NC_DOUBLE"
    }
    att.put.nc(nc, variable = "NC_GLOBAL", name = names(attributes[i]), type = type, value = attributes[[i]])
  }
  sync.nc(nc)
}

# Add attributes
add_global_attributes(nc, attributes)

# Close nc file
close.nc(nc)

paste0("NetCDF data containing shannon index and species richness (n taxa) has been stored in ", derivedDir ,"/alpha_diversity.nc")

