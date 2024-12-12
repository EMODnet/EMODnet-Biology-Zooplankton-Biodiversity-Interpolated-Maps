library(RNetCDF)

# Define path to NetCDF output
net_cdf_file <- file.path("product", "netcdf" , "zooplankton_diversity.nc")

# Load DIVAnd data
load("data/derived_data/diva3d.rda")

# Create nc file
nc <- create.nc(net_cdf_file) 

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
dim.def.nc(nc, dimname = "time", dimlength = length(tt_days)) 

# Define time variable
var.def.nc(nc, varname = "time", vartype = "NC_DOUBLE", dimensions = "time")

# Add attributes
att.put.nc(nc, variable = "time", name = "standard_name", type = "NC_CHAR", value = "time")
att.put.nc(nc, variable = "time", name = "long_name", type = "NC_CHAR", value = "Time")
att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "days since 1970-01-01 00:00:00")
att.put.nc(nc, variable = "time", name = "calendar", type = "NC_CHAR", value = "gregorian")

# Put data
var.put.nc(nc, variable = "time", data = tt_days)

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
  title = "Spatio-Temporal Interpolation of Zooplankton Alpha Diversity in the Greater Baltic Sea Region",
  summary = "",                       
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
  citation = "Torstensson A (2024) Spatio-Temporal Interpolation of Zooplankton Alpha Diversity in the Greater Baltic Sea Region. Integrated data products created under the European Marine Observation Data Network (EMODnet) Biology project Phase V (CINEA/EMFAF/2022/3.5.2/SI2.895681), funded by the by the European Union under Regulation (EU) No 2021/1139 of the European Parliament and of the Council of 7 July 2021 on the European Maritime and Fisheries Fund.",
  acknowledgement = "European Marine Observation Data Network (EMODnet) Biology project (CINEA/EMFAF/2022/3.5.2/SI2.895681), funded by the European Union under Regulation (EU) No 2021/1139 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund"
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

cat("NetCDF data containing shannon index has been stored in", net_cdf_file, "\n")
