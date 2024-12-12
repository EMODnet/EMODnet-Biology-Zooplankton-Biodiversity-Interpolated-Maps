# Spatio-Temporal Interpolation of Zooplankton Alpha Diversity in the Greater Baltic Sea Region

## Introduction

Biodiversity loss due to human activities is an increasing threat for marine ecosystems and the services we obtain from them. As biodiversity is directly related to the resilience of ecosystems to temporary disturbance, biodiversity monitoring is a vital task for areas subjected to conservation goals. Environmental factors often control the community composition and biodiversity of marine plankton, such as the pronounced salinity gradient in the Baltic Sea (e.g. Hu et al. 2016). Time series data of biodiversity can therefore provide an indication of changes in community composition due to environmental stressors, such as climate change or eutrophication.

As many biodiversity estimates are biased by sampling effort, caution must be taken when interpreting alpha diversity from microscopy counts. By rarefaction and evenness estimation, these biases can be reduced, but not ignored.

## Directory structure

```
EMODnet-Biology-Zooplankton-Biodiversity-Interpolated-Maps/
├── analysis
├── data/
│   ├── derived_data/
│   └── raw_data/
├── docs/
├── product/
│   ├── animations/
│   ├── maps/
│   ├── netcdf/
│   └── plots/
└── scripts/
```

* **analysis** - Markdown or Jupyter notebooks
* **data** - Raw and derived data
* **docs** - Rendered reports
* **product** - Output product files
* **scripts** - Reusable code

## Data series

This data product use the following datasets:

- Stockholm University, Umeå University, Swedish Meteorological and Hydrological Institute, Gothenburg University, Swedish Agency for Marine and Water Management and Swedish Environmental Protection Agency (2022). SHARK - National marine environmental monitoring of zooplankton in Sweden since 1979 https://doi.org/10.15468/edkb5r

- Finnish Environment Institute SYKE; (2018); Finnish Baltic Sea zooplankton monitoring

```
https://obis.org/dataset/7f29807d-c940-4136-9ccd-3baa1e7e9bab
https://obis.org/dataset/f1ad2d12-8f7e-4b56-848c-ffa9c6103eb1
```

## Data product

This data product provides a spatial-temporal interpolation of zooplankton alpha diversity in the greater Baltic Sea area. Using observed Shannon diversity indices from environmental data, the product employs the DIVAnd (Data-Interpolating Variational Analysis in n-Dimensions) to create a three-dimensional interpolation across longitude, latitude, and seasonal time dimensions.

The output includes interpolated diversity values, an associated error map, and metadata on the spatial-temporal grid resolution. It is tailored for ecological and oceanographic studies, aiding in the visualization and analysis of diversity trends over time and space.

[Zooplankton_Diversity_Animation](https://github.com/user-attachments/assets/4ad56d24-6bb2-4328-a3cf-65f6922e88c7)

---

### Data Input

#### 1. Observation Data
- **Observation Data**: The zooplankton data are derived from Swedish and Finnish national monitoring programs. Data sources include environmental surveys and sampling programs that collect zooplankton samples using standard plankton nets and methodologies.
- **Ecological Metrics**: Shannon diversity index values are calculated from these samples after species abundance data are aggregated at each sampling station.
- **Time Range**: Observations are included starting from 2007, when data collection became more regular, ensuring sufficient coverage for accurate spatial-temporal interpolation.

#### 2. Normalization Using Rarefaction
- To standardize data across varying sample sizes, rarefaction is applied to species abundance data. This ensures comparability of diversity metrics by subsampling to a uniform number of individuals (e.g., 100 individuals per sample). Rarefaction minimizes biases due to differences in sampling effort or organism density.
- The rarefaction algorithm iteratively subsamples without replacement, calculating diversity metrics for the standardized dataset.

#### 3. Diversity Calculation
- **Shannon Diversity Index**: The Shannon diversity index (`H`) is calculated as:
  \[
  H = -\sum_{i=1}^S p_i \ln(p_i)
  \]
  where \(S\) is the total number of species and \(p_i\) is the proportional abundance of species \(i\).
- **Input Data**: The index is derived from the rarefied species abundance data, ensuring that all samples have equal weight regardless of original size.
- **Interpretation**: Higher values indicate greater species richness and evenness, reflecting a more diverse zooplankton community.

#### 4. Land Mask
- Natural Earth polygons used to mask grid points that fall on land.
- Additional masking for points north of 63° latitude and west of 15° longitude.

#### 5. Grid Settings
- **Longitude**: 9.5°E to 30°E (`xx` grid, 100 steps).
- **Latitude**: 53.5°N to 66°N (`yy` grid, 110 steps).
- **Temporal Dimension**: Numeric seasons calculated as unique time steps (`tt` grid).

---

### Data Processing Workflow

#### 1. Season Classification
- Observations were assigned to **Spring**, **Summer**, **Autumn**, or **Winter** based on event date (`eventDate`).
- A numeric representation accounts for Winter overlap across years.

#### 2. Spatial-Temporal Grid Construction
- A regular 3D grid was created, combining longitude (`xx`), latitude (`yy`), and numeric time (`tt`).
- Land and additional condition masking were applied to exclude invalid grid points.

#### 3. Interpolation with DIVAnd
- Interpolation used a correlation length of 5 degrees for space and 0.25 units for time.
- Observational error variance was normalized by background error variance by epsilon = 1.0.

#### 4. Error Map
- An error map was computed using the `DIVAnd_errormap` function, employing a "cheap" approximation for computational efficiency (Beckers et al. 2014).

#### 5. Post-Processing
- Interpolated values were clipped to non-negative values.
- High-error regions were set to `NA` using a relative threshold error value.

---

### Data Outputs

The output from the spatial-temporal interpolation of zooplankton alpha diversity in the greater Baltic Sea area includes the following components:

1. **Interpolated 3D Grid (`fi`)**:
   - This array contains Shannon diversity index values, representing the diversity of zooplankton across the spatial-temporal grid where values are below a error threshold. It integrates longitude, latitude, and time dimensions.

2. **Error Map (`e`)**:
   - An array indicating the error estimates associated with each interpolated grid point. This helps assess the reliability of the interpolated values.

3. **Grid Metadata**:
   - Includes the coordinates of the grid (longitude, latitude, and time) as well as spatial-temporal resolution metrics (grid spacing in longitude, latitude, and time).
   - Also includes a mask to identify valid grid points.

4. **Saved Data**:
   - Additionally, a NetCDF file (`alpha_diversity.nc`) is created. This file stores the grid shannon data along with relevant attributes such as geographic coordinates, time reference, and metadata for the Shannon diversity index.

These outputs provide a detailed, spatially-explicit representation of zooplankton diversity over time, suitable for further analysis or visualization.

---

### Limitations
- Spatial and temporal gaps in observations may introduce interpolation uncertainty.
- Error thresholds and masking criteria may require adjustment for specific applications.

---

This product combines advanced spatial interpolation techniques with ecological modeling, enabling a deeper understanding of biodiversity dynamics in the Baltic Sea region.

## More information:

### References

- Beckers, J., A. Barth, C. Troupin, and A. Alvera-Azcárate, 2014: Approximate and Efficient Methods to Assess Error Fields in Spatial Gridding with Data Interpolating Variational Analysis (DIVA). J. Atmos. Oceanic Technol., 31, 515–530, https://doi.org/10.1175/JTECH-D-13-00130.1.

### Code and methodology

{{link_code}}

### Citation and download link

This product should be cited as:

{{product_citation}}

Available to download in:

{{link_download}}

### Authors

Anders Torstensson, Swedish Meteorological and Hydrological Institute
