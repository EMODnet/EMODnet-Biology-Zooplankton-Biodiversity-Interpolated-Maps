# Load necessary libraries for data manipulation, biodiversity analysis, and shapefile handling
library(tidyverse)  # For data manipulation and wrangling
library(finch)      # For handling Darwin Core Archive (DwC-A) biodiversity data
library(vegan)      # For community ecology analysis, including species diversity calculations
library(worrms)     # For taxonomic information retrieval from the World Register of Marine Species (WoRMS)

### Finnish Data - Zooplankton Monitoring - v.1.0 does not contain any sample dates, use different version below
# Source: https://obis.org/dataset/7f29807d-c940-4136-9ccd-3baa1e7e9bab
# Load zooplankton data from a Darwin Core Archive (DwC-A) zip file
# fin_dwca <- dwca_read("data/raw_data/dwca-syke-zplank-v1.0.zip")
# 
# # Read individual text files from the DwC-A archive: event, occurrence, and measurement or fact data (EMOF)
# event_fin <- read_tsv(fin_dwca$files$txt_files[1])             # Event data (e.g., sampling locations and dates)
# occurrence_fin <- read_tsv(fin_dwca$files$txt_files[3])         # Occurrence data (species presence/absence)
# emof_fin <- read_tsv(fin_dwca$files$txt_files[2])               # EMOF data (e.g., measurements like abundance)
# 
# # Extract unique species identifiers (Life Science Identifiers - LSID) from occurrence data
# lsid_fin <- unique(occurrence_fin$scientificNameID)
# 
# # Extract AphiaIDs from LSIDs, assuming AphiaID is the numeric part at the end of each LSID
# aphia_id_fin <- as.numeric(str_extract(lsid_fin, "\\d+$"))
# 
# # Initialize an empty data frame to store species taxonomic records from WoRMS (World Register of Marine Species)
# worms_records_fin <- data.frame()
# 
# # Loop through each AphiaID to retrieve corresponding species record from WoRMS
# for (id in seq_along(aphia_id_fin)) {
#   worms_record_fin <- wm_record(aphia_id_fin[id])  # Retrieve taxonomic record for each AphiaID
#   worms_records_fin <- bind_rows(worms_records_fin, worms_record_fin)  # Append each record to worms_records_est data frame
# }
# 
# # Create a taxonomic table linking scientificNameID and AphiaID, and merge with WoRMS records
# tax_table_fin <- data.frame(scientificNameID = lsid_fin,
#                             AphiaID = aphia_id_fin) %>%
#   left_join(worms_records_fin, by = "AphiaID")  # Join WoRMS records by AphiaID to add taxonomic details
# 
# # Select key columns from the event data for later merging
# event_fin_sel <- event_fin %>%
#   select(id, eventID, eventDate, decimalLatitude, decimalLongitude)
# 
# # Merge, filter, and process the data to retain relevant records
# all_fin <- occurrence_fin %>%
#   left_join(emof_fin) %>%  # Join EMOF data to occurrence data by matching columns
#   left_join(event_fin_sel) %>%  # Join selected event data columns
#   left_join(tax_table_fin, by = "scientificNameID") %>%  # Join with taxonomic table to add taxonomic details
#   filter(measurementType == "Abundance of biological entity specified elsewhere per unit volume of the water body") %>%  # Select only abundance measurements ("arvukus" means abundance)
#   # filter(minimumDepthInMeters == 0) %>%  # Select samples taken at the surface (0 meters)
#   filter(!phylum.y == "Ciliophora") %>%
#   filter(year(eventDate) >= 2000)  # Include only data from the year 2000 onwards
# 
# # Create a species abundance matrix for the selected samples
# fin_species_matrix <- all_fin %>%
#   mutate(sample_name = eventID) %>%  # Create a sample name based on the unique eventID
#   select(scientificNameID, sample_name, measurementValue) %>%  # Select relevant columns for the matrix
#   
#   # Group by species and sample name to calculate total abundance per sample
#   group_by(scientificNameID, sample_name) %>%
#   summarise(abundance = sum(as.numeric(measurementValue), na.rm = TRUE), .groups = 'drop') %>%  # Sum abundances, handle NA values
#   
#   # Transform data to wide format: each column represents a sample, each row a species
#   pivot_wider(names_from = sample_name, values_from = abundance, values_fill = list(abundance = 0))

# Load CSV files related to events and occurrences of species
occurrence_fin <- read_csv("data/raw_data/d6a04582-1fd2-4443-82ab-95fbccd3df74/Occurrence.csv")  # Occurrence data
emof_fin <- read_csv("data/raw_data/d6a04582-1fd2-4443-82ab-95fbccd3df74/MeasurementOrFact.csv")  # Measurement data (EMOF)

# Filter the occurrence and measurement data
all_fin <- occurrence_fin %>%
  left_join(emof_fin) %>%                                                   # Join occurrence and measurement data by common columns
  filter(measurementtype == "Abundance of biological entity specified elsewhere per unit volume of the water body") %>%  # Select records with abundance measurements
  filter(!phylum == "Ciliophora") %>%                                       # Exclude ciliates (Phylum Ciliophora)
  filter(minimumdepthinmeters == 0) %>%                                     # Select samples taken at the surface (0 meters)
  filter(year(eventdate) >= 2000 & year(eventdate) <= 2022)                 # Include only data between year 2000 to 2022

# Create a species abundance matrix
fin_species_matrix <- all_fin %>%
  mutate(sample_name = eventid) %>%                                         # Create a sample name column based on event ID
  select(scientificnameid, sample_name, measurementvalue) %>%               # Select species ID, sample name, and measurement value columns
  
  # Group by species ID and sample name to calculate total abundance per sample
  group_by(scientificnameid, sample_name) %>%
  summarise(abundance = sum(as.numeric(measurementvalue), na.rm = TRUE), .groups = 'drop') %>%  # Sum abundances, handle NA values
  
  # Transform data to wide format: each column represents a sample, each row a species
  pivot_wider(names_from = sample_name, values_from = abundance, values_fill = list(abundance = 0)) %>%
  rename(scientificNameID = scientificnameid)                               # Rename scientific name ID for consistency

### Swedish Data - Zooplankton Monitoring
# Source: https://www.gbif.se/ipt/resource?r=shark-zooplankton-nat

# Load zooplankton data from a Darwin Core Archive (DwC-A) zip file
swe_dwca <- dwca_read("data/raw_data/dwca-shark-zooplankton-nat-v1.9.zip")

# Read individual text files from the DwC-A archive: event, occurrence, and measurement or fact data (EMOF)
event <- read_tsv(swe_dwca$files$txt_files[1])             # Event data (e.g., sampling locations and dates)
occurrence <- read_tsv(swe_dwca$files$txt_files[2])         # Occurrence data (species presence/absence)
emof <- read_tsv(swe_dwca$files$txt_files[3])               # EMOF data (e.g., measurements like abundance)

# Select relevant columns from the event data
event_sel <- event %>%
  select(id, eventID, parentEventID, eventDate, verbatimLocality, minimumDepthInMeters, maximumDepthInMeters, decimalLatitude, decimalLongitude)

# Select relevant columns from the occurrence data
occurrence_sel <- occurrence %>%
  select(id, occurrenceID, measurementType, measurementValue, measurementUnit)

# Combine the occurrence, event, and measurement data into one dataset
all_swe <- occurrence_sel %>%
  left_join(event_sel) %>%                                           # Join event data to occurrence data by `id`
  left_join(emof)                                                    # Join EMOF data to the combined dataset

# Filter to select surface zooplankton abundance, excluding specific phyla and microplankton data
abundance_swe <- all_swe %>%
  filter(measurementType == "Abundance") %>%                         # Filter for records of type "Abundance"
  filter(!phylum %in% c("Ciliophora", "Myzozoa", "Radiozoa", "Foraminifera")) %>%  # Exclude specified microplankton phyla
  filter(minimumDepthInMeters == 0) %>%                              # Select samples taken at the surface (0 meters)
  filter(year(eventDate) >= 2000 & year(eventDate) <= 2022)          # Include only data between year 2000 to 2022

# Create a species abundance matrix for the selected samples
swe_species_matrix <- abundance_swe %>%
  mutate(sample_name = id) %>%                                       # Create a sample name based on the unique ID
  select(scientificNameID, sample_name, measurementValue) %>%        # Select relevant columns for the matrix
  
  # Group by species and sample name to calculate total abundance per sample
  group_by(scientificNameID, sample_name) %>%
  summarise(abundance = sum(as.numeric(measurementValue), na.rm = TRUE), .groups = 'drop') %>%  # Sum abundances, handle NA values
  
  # Transform data to wide format: each column represents a sample, each row a species
  pivot_wider(names_from = sample_name, values_from = abundance, values_fill = list(abundance = 0))


### Poland Data - Zooplankton Monitoring
# Source: https://ipt.vliz.be/eurobis/resource?r=pl_zo_monitoring

# Load zooplankton data for Poland from a Darwin Core Archive (DwC-A) zip file
pl_dwca <- dwca_read("data/raw_data/dwca-pl_zo_monitoring-v1.0.zip")

# Read occurrence and EMOF (measurement or fact) data from the DwC-A files
occurrence_pl <- read_tsv(pl_dwca$files$txt_files[2])  # Occurrence data (species presence/absence)
emof_pl <- read_tsv(pl_dwca$files$txt_files[1])        # EMOF data (e.g., measurements like abundance)

# Extract unique species identifiers (Life Science Identifiers - LSID) from occurrence data
lsid <- unique(occurrence_pl$scientificNameID)

# Extract AphiaIDs from LSIDs for species information, assuming AphiaID is the numeric part at the end of each LSID
aphia_id <- as.numeric(str_extract(lsid, "\\d+$"))

# Initialize an empty data frame to store species taxonomic records from WoRMS (World Register of Marine Species)
worms_records <- data.frame()

# Loop through each AphiaID to retrieve corresponding species record from WoRMS
for (id in seq_along(aphia_id)) {
  worms_record <- wm_record(aphia_id[id])  # Retrieve taxonomic record for each AphiaID
  worms_records <- bind_rows(worms_records, worms_record)  # Append each record to worms_records data frame
}

# Create a taxonomic table linking scientificNameID and AphiaID, and merge with WoRMS records
tax_table_pl <- data.frame(scientificNameID = lsid,
                           AphiaID = aphia_id) %>%
  left_join(worms_records, by = "AphiaID")  # Join WoRMS records by AphiaID to add taxonomic details

# Filter and merge data
all_pl <- occurrence_pl %>%
  left_join(emof_pl) %>%  # Join EMOF data to occurrence data by matching columns
  left_join(tax_table_pl, by = "scientificNameID") %>%  # Join with taxonomic table to add taxonomic details
  filter(measurementType == "abundance") %>%  # Select only abundance measurements
  filter(minimumDepthInMeters == 0) %>%  # Select samples taken at the surface (0 meters)
  filter(year(eventDate) >= 2000)  # Include only data from the year 2000 onwards

# Create a species abundance matrix for the selected samples
pl_species_matrix <- all_pl %>%
  mutate(sample_name = eventID) %>%  # Create a sample name based on the unique eventID
  select(scientificNameID, sample_name, measurementValue) %>%  # Select relevant columns for the matrix
  
  # Group by species and sample name to calculate total abundance per sample
  group_by(scientificNameID, sample_name) %>%
  summarise(abundance = sum(as.numeric(measurementValue), na.rm = TRUE), .groups = 'drop') %>%  # Sum abundances, handle NA values
  
  # Transform data to wide format: each column represents a sample, each row a species
  pivot_wider(names_from = sample_name, values_from = abundance, values_fill = list(abundance = 0))

### Estonian Data - Zooplankton Monitoring
# Source: https://ipt.vliz.be/eurobis/resource?r=zooplankton_data_estonian_territorial_waters_1993-2016

# Load Estonian zooplankton data from a Darwin Core Archive (DwC-A) zip file
est_dwca <- dwca_read("data/raw_data/dwca-zooplankton_data_estonian_territorial_waters_1993-2016-v1.0.zip")

# Read event, occurrence, and EMOF (measurement or fact) data from the DwC-A files
event_est <- read_tsv(est_dwca$files$txt_files[1])  # Event data (sampling events)
occurrence_est <- read_tsv(est_dwca$files$txt_files[3]) %>% 
  select(-decimalLatitude, -decimalLongitude, -eventDate)  # Occurrence data, excluding specific columns
emof_est <- read_tsv(est_dwca$files$txt_files[2])  # EMOF data

# Extract unique species identifiers (Life Science Identifiers - LSID) from occurrence data
lsid_est <- unique(occurrence_est$scientificNameID)

# Extract AphiaIDs from LSIDs, assuming AphiaID is the numeric part at the end of each LSID
aphia_id_est <- as.numeric(str_extract(lsid_est, "\\d+$"))

# Initialize an empty data frame to store species taxonomic records from WoRMS (World Register of Marine Species)
worms_records_est <- data.frame()

# Loop through each AphiaID to retrieve corresponding species record from WoRMS
for (id in seq_along(aphia_id_est)) {
  worms_record_est <- wm_record(aphia_id_est[id])  # Retrieve taxonomic record for each AphiaID
  worms_records_est <- bind_rows(worms_records_est, worms_record_est)  # Append each record to worms_records_est data frame
}

# Create a taxonomic table linking scientificNameID and AphiaID, and merge with WoRMS records
tax_table_est <- data.frame(scientificNameID = lsid_est,
                            AphiaID = aphia_id_est) %>%
  left_join(worms_records_est, by = "AphiaID")  # Join WoRMS records by AphiaID to add taxonomic details

# Select key columns from the event data for later merging
event_est_sel <- event_est %>%
  select(id, eventID, eventDate, decimalLatitude, decimalLongitude)

# Merge, filter, and process the data to retain relevant records
all_est <- occurrence_est %>%
  left_join(emof_est) %>%  # Join EMOF data to occurrence data by matching columns
  left_join(event_est_sel) %>%  # Join selected event data columns
  left_join(tax_table_est, by = "scientificNameID") %>%  # Join with taxonomic table to add taxonomic details
  filter(measurementType == "arvukus") %>%  # Select only abundance measurements ("arvukus" means abundance)
  filter(minimumDepthInMeters == 0) %>%  # Select samples taken at the surface (0 meters)
  mutate(measurementValue = as.numeric(measurementValue)) %>%  # Convert measurement values to numeric
  filter(!is.na(measurementValue)) %>%  # Filter out any NA values in measurement values
  filter(year(eventDate) >= 2000)  # Include only data from the year 2000 onwards

# Create a species abundance matrix for the selected samples
est_species_matrix <- all_est %>%
  mutate(sample_name = eventID) %>%  # Create a sample name based on the unique eventID
  select(scientificNameID, sample_name, measurementValue) %>%  # Select relevant columns for the matrix
  
  # Group by species and sample name to calculate total abundance per sample
  group_by(scientificNameID, sample_name) %>%
  summarise(abundance = sum(as.numeric(measurementValue), na.rm = TRUE), .groups = 'drop') %>%  # Sum abundances, handle NA values
  
  # Transform data to wide format: each column represents a sample, each row a species
  pivot_wider(names_from = sample_name, values_from = abundance, values_fill = list(abundance = 0))

### Danish Data - Zooplankton Monitoring
# Source: https://ipt.vliz.be/eurobis/resource?r=odam_zooplankton_1985

# Load Danish zooplankton data from a Darwin Core Archive (DwC-A) zip file
dk_dwca <- dwca_read("data/raw_data/dwca-odam_zooplankton_1985-v1.0.zip")

# Read occurrence and EMOF (measurement or fact) data from the DwC-A files
occurrence_dk <- read_tsv(dk_dwca$files$txt_files[2])  # Occurrence data
emof_dk <- read_tsv(dk_dwca$files$txt_files[1])  # EMOF data

# Extract unique species identifiers (Life Science Identifiers - LSID) from occurrence data
lsid_dk <- unique(occurrence_dk$scientificNameID)

# Extract AphiaIDs from LSIDs, assuming AphiaID is the numeric part at the end of each LSID
aphia_id_dk <- as.numeric(str_extract(lsid_dk, "\\d+$"))

# Initialize an empty data frame to store species taxonomic records from WoRMS (World Register of Marine Species)
worms_records_dk <- data.frame()

# Loop through each AphiaID to retrieve corresponding species record from WoRMS
for (id in seq_along(aphia_id_dk)) {
  worms_record_dk <- wm_record(aphia_id_dk[id])  # Retrieve taxonomic record for each AphiaID
  worms_records_dk <- bind_rows(worms_records_dk, worms_record_dk)  # Append each record to worms_records_dk data frame
}

# Create a taxonomic table linking scientificNameID and AphiaID, and merge with WoRMS records
tax_table_dk <- data.frame(scientificNameID = lsid_dk,
                           AphiaID = aphia_id_dk) %>%
  left_join(worms_records_dk, by = "AphiaID")  # Join WoRMS records by AphiaID to add taxonomic details

# Merge and process the data to retain relevant records, converting abundance units
all_dk <- occurrence_dk %>%
  left_join(emof_dk) %>%  # Join EMOF data to occurrence data by matching columns
  left_join(tax_table_dk, by = "scientificNameID") %>%  # Join with taxonomic table to add taxonomic details
  filter(measurementType == "WaterAbund_BE007117....l.") %>%  # Select only abundance measurements (specific measurement type)
  filter(year(eventDate) >= 2000) %>%  # Include only data from the year 2000 onwards
  mutate(measurementValue = measurementValue * 1000)  # Convert abundance from "#/l" to "#/m3" (multiply by 1000)

# Create a species abundance matrix for the selected samples
dk_species_matrix <- all_dk %>%
  mutate(sample_name = eventID) %>%  # Create a sample name based on the unique eventID
  select(scientificNameID, sample_name, measurementValue) %>%  # Select relevant columns for the matrix
  
  # Group by species and sample name to calculate total abundance per sample
  group_by(scientificNameID, sample_name) %>%
  summarise(abundance = sum(as.numeric(measurementValue), na.rm = TRUE), .groups = 'drop') %>%  # Sum abundances, handle NA values
  
  # Transform data to wide format: each column represents a sample, each row a species
  pivot_wider(names_from = sample_name, values_from = abundance, values_fill = list(abundance = 0))

### Calculate diversity metrics

# Combine all species matrices into one
species_matrix <- fin_species_matrix %>%
  full_join(swe_species_matrix) %>%  # Join Swedish species matrix
  # full_join(pl_species_matrix) %>%    # Join Polish species matrix, skip - full time series (2000-2022) is not covered
  # full_join(est_species_matrix) %>%   # Join Estonian species matrix, skip - full time series (2000-2022) is not covered
  # full_join(dk_species_matrix) %>%    # Join Danish species matrix, skip - full time series (2000-2022) is not covered
  mutate(across(everything(), ~replace_na(., 0)))  # Replace NA values with 0 for counts

# Gather the combined data into long format
long_species_matrix <- species_matrix %>%
  pivot_longer(
    cols = -scientificNameID,  # Exclude the scientificNameID column from pivoting
    names_to = "sample",  # Create a new column for sample names
    values_to = "count"   # Create a new column for counts
  )

# Pivot back to a wide format to have species as columns
pivoted_species_matrix <- long_species_matrix %>%
  pivot_wider(
    names_from = scientificNameID,  # Species names become columns
    values_from = count,              # Fill with counts
    values_fill = 0                    # Replace NAs with 0
  )

# Prepare the counts for further analysis
species_counts <- pivoted_species_matrix[,-1]  # Remove the identifier column

# Determine the minimum sample size across all samples
sample_sizes <- rowSums(species_counts)  # Calculate the total counts per sample
min_sample_size <- min(sample_sizes)  # Find the minimum sample size

# Set a threshold for minimum sample size
min_threshold <- 1000

# Filter out samples below the minimum threshold from species_matrix
clean_species_matrix <- species_counts[rowSums(species_counts, na.rm = TRUE) >= min_threshold,]
pivoted_species_matrix <- pivoted_species_matrix[rowSums(species_counts, na.rm = TRUE) >= min_threshold,]

# Ensure counts are integers for further analysis
clean_species_matrix <- round(clean_species_matrix)

# Re-determine the minimum sample size after filtering
sample_sizes <- rowSums(clean_species_matrix)
min_sample_size <- min(sample_sizes)

# Rarefy the species counts to the minimum sample size
rarefied_counts <- rrarefy(as.matrix(clean_species_matrix), min_sample_size)

# Calculate Shannon diversity index
shannon <- diversity(rarefied_counts, index = "shannon")

# Calculate species richness (number of unique species per sample)
uniqueTaxa <- rowSums(rarefied_counts > 0)

# Combine diversity metrics with species identifiers
diversity_metrics <- as_tibble(cbind(pivoted_species_matrix[1], shannon, uniqueTaxa))

# Extract metadata for each dataset and ensure distinct values
swe_metadata <- abundance_swe %>%
  select(eventID, eventDate, decimalLatitude, decimalLongitude) %>%
  distinct()

syke_metadata <- all_fin %>%
  select(eventid, eventdate, decimallatitude, decimallongitude) %>%
  rename(eventID = eventid,
         eventDate = eventdate,
         decimalLatitude = decimallatitude,
         decimalLongitude = decimallongitude) %>%
  distinct()

est_metadata <- all_est %>%
  select(eventID, eventDate, decimalLatitude, decimalLongitude) %>%
  distinct()

pl_metadata <- all_pl %>%
  select(eventID, eventDate, decimalLatitude, decimalLongitude) %>%
  distinct()

dk_metadata <- all_dk %>%
  select(eventID, eventDate, decimalLatitude, decimalLongitude) %>%
  distinct()

# Combine all metadata into one data frame
metadata <- rbind(swe_metadata, syke_metadata, est_metadata, pl_metadata, dk_metadata)

# Combine diversity metrics with metadata based on sample/event IDs
diversity_all <- diversity_metrics %>%
  left_join(metadata, by = c("sample" = "eventID")) %>%
  rename(id = sample)  # Rename sample column to id for clarity

# Save the combined diversity metrics and metadata to a .rda file
save(diversity_all, file= "data/derived_data/data_diversity_all.rda")
