# Buonanduci_etal_2023_EcolLett
Data and code accompanying the manuscript 'Consistent spatial scaling of high-severity wildfire can inform expected future patterns of burn severity' by Buonanduci, Donato, Halofsky, Kennedy, and Harvey. Manuscript provisionally accepted for publication in Ecology Letters.


[![CC BY 4.0][cc-by-shield]][cc-by]

This information is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by]. Any user of these data ("User" hereafter) is required to cite it appropriately in any publication that results from its use. These data may be actively used by others for ongoing research, so coordination may be necessary to prevent duplicate publication. The User is urged to contact the authors of these data for questions about methodology or results.  The User is encouraged to consider collaboration or co-authorship with authors where appropriate. Misinterpretation of data may occur if used out of context of the original study. Substantial efforts are made to ensure accuracy of the data and documentation, however complete accuracy of data sets cannot be guaranteed. All data are made available as is. Data may be updated periodically and it is the responsibility of the User to check for new versions of the data. The authors and the repository where these data were obtained shall not be liable for damages resulting from any use or misinterpretation of the data.

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg


For reproducibility, the following files are made available:

Data files

- fire_metrics.csv
- fire_patch_distributions.csv
- fire_DTS_distributions.csv

Geospatial data files

- fire_perims.gpkg
- region_NR.gpkg
- region_PNW.gpkg
- states.gpkg
- OR4310612253920020714_HS.tif
- OR4310612253920020714_DTS.tif
- OR4317212252420090913_HS.tif
- OR4317212252420090913_DTS.tif
- WY4437710988020190902_HS.tif
- WY4437710988020190902_DTS.tif

R scripts

- install_quantregGrowth.R
- figure_study_region.R
- figure_patch_metrics.R
- figure_scaling.R
- figure_scaling_year.R
- figure_scaling_regimes.R
- figure_scaling_regions.R
- figure_scaling_time_period.R
- cross_val_region.R
- cross_val_year_time_period.R
- fit_truncated_ln.R
- kfold_fun.R
- predict.gcrq.inputform.R
- predict.gcrq.lpmatrix.R


**Note**: Version `1.4-0` of the `quantregGrowth` package must be installed by the user to ensure all R scripts run as intended.



### fire_metrics.csv
This file contains high-severity patch size and structure metrics for all fire events included in the analysis. The following columns are included:

- **Fire_ID**: unique fire identifier, as assigned by MTBS
- **Fire_Name**: fire name, as assigned by MTBS
- **Year**: year in which fire occurred
- **prp_forest**: proportion of burned area that is forested, as classified by LANDFIRE Environmental Site Potential
- **Fire_Regime**: historical fire regime; either "Low" (frequent and low severity), "Mixed" (moderately frequent and mixed severity), or "High" (infrequent and high severity), as classified by LANDFIRE Fire Regime Group
- **Region**: geographic region; either "Northern Rockies" or "Pacific Northwest"
- **Time_Period**: time period in which fire occurred; either "Early" (1985 - 2000) or "Late" (2001 - 2020)
- **fire_area**: total area burned (ha)
- **log_fire_area**: log10-transformed total area burned
- **HS_area**: total area burned at high severity (ha)
- **HS_prop**: proportion of area burned at high severity
- **HS_forest_area**: total forested area burned at high severity (ha)
- **HS_forest_prop**: proportion of area that was both forested and burned at high severity
- **patch_area_mean**: arithmetic mean of high-severity patch sizes (ha)
- **patch_area_AW_mean**: area-weighted mean of high-severity patch sizes (ha)
- **log_patch_area_AW_mean**: log10-transformed area-weighted mean of high-severity patch sizes
- **beta**: beta parameter for high-severity patch size distribution
- **psi**: psi parameter for high-severity patch size distribution
- **n_patches**: number of high-severity patches
- **n_patches_gte1**: number of high-severity patches equal to or exceeding 1 ha in size
- **total_core**: total core area (ha)
- **log_total_core**: log10-transformed total core area
- **SDC**: stand-replacing decay coefficient; parameter for distance-to-seed distribution
- **log_SDC**: log10-transformed SDC parameter



### fire_patch_distributions.csv
This file contains high-severity patch size distributions for all fire events. The following columns are included:

- **Fire_ID**: unique fire identifier, as assigned by MTBS
- **patch_area_ha**: size of high-severity patch within fire event (ha)



### fire_DTS_distributions.csv
This file contains distance-to-seed distributions for all fire events. The following columns are included:

- **Fire_ID**: unique fire identifier, as assigned by MTBS
- **DTS_dist_m**: distance to potential seed source (m)
- **DTS_cells**: number of 30-m high-severity cells exceeding each distance-to-seed threshold
- **DTS_cells_prop**: proportion of high-severity cells exceeding each distance-to-seed threshold



### fire_perims.gpkg
This file contains polygon perimeters for all fire events. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). The following columns are included:

- **Fire_ID**: unique fire identifier, as assigned by MTBS
- **Fire_Name**: fire name, as assigned by MTBS
- **Year**: year in which fire occurred
- **Fire_Regime**: historical fire regime; either "Low" (frequent and low severity), "Mixed" (moderately frequent and mixed severity), or "High" (infrequent and high severity), as classified by LANDFIRE Fire Regime Group



### region_NR.gpkg
This file contains a polygon perimeter for the Northern Rockies region. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). 



### region_PNW.gpkg
This file contains a polygon perimeter for the Northern Rockies region. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). 



### states.gpkg
This file contains US state outlines, for plotting purposes. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). 



### OR4310612253920020714_HS.tif
This file contains a classified burn severity raster for the 2002 "Big Bend" fire in Oregon (Fire_ID = OR4310612253920020714). Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). The value field is categorical, equaling either 1 (low/moderate burn severity) or 2 (high severity).



### OR4310612253920020714_DTS.tif
This file contains a distance-to-seed raster for the 2002 "Big Bend" fire in Oregon (Fire_ID = OR4310612253920020714). Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). The value field is continuous and indicates the distance to potential seed source (m) for each cell that was forested and burned at high severity.



### OR4317212252420090913_HS.tif
This file contains a classified burn severity raster for the 2009 "Boze" fire in Oregon (Fire_ID = OR4317212252420090913). Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). The value field is categorical, equaling either 1 (low/moderate burn severity) or 2 (high severity).



### OR4317212252420090913_DTS.tif
This file contains a distance-to-seed raster for the 2009 "Boze" fire in Oregon (Fire_ID = OR4317212252420090913). Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). The value field is continuous and indicates the distance to potential seed source (m) for each cell that was forested and burned at high severity.



### WY4437710988020190902_HS.tif
This file contains a classified burn severity raster for the 2019 "Fishhawk" fire in Wyoming (Fire_ID = WY4437710988020190902). Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). The value field is categorical, equaling either 1 (low/moderate burn severity) or 2 (high severity).



### WY4437710988020190902_DTS.tif
This file contains a distance-to-seed raster for the 2019 "Fishhawk" fire in Wyoming (Fire_ID = WY4437710988020190902). Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). The value field is continuous and indicates the distance to potential seed source (m) for each cell that was forested and burned at high severity.



### install_quantregGrowth.R
Code for installing the appropriate version of the `quantregGrowth` package. Version `1.4-0` of the `quantregGrowth` package must be installed by the user to ensure all R scripts within this repository run as intended.



### figure_study_region.R
Code for producing Figure 1 in the manuscript.



### figure_patch_metrics.R
Code for producing components of Figure 2 in the manuscript.



### figure_scaling.R
Code for fitting quantile regression models and producing Figure 3 in the manuscript.



### figure_scaling_year.R
Code for evaluating effect of year and producing Figure 4 in the manuscript.



### figure_scaling_regimes.R
Code for evaluating effect of fire regime and producing Figure S1 in the supporting information file.



### figure_scaling_regions.R
Code for evaluating effect of region and producing Figure S2 in the supporting information file.



### figure_scaling_time_period.R
Code for evaluating effect of time period and producing Figure S3 in the supporting information file.



### cross_val_region.R
Code for cross-validation analysis evaluating the effect of geographic region.



### cross_val_year_time_period.R
Code for cross-validation analysis evaluating the effects of year and time period.



### fit_truncated_ln.R
Function for fitting truncated lognormal distributions to empirical patch size distributions.



### kfold_fun.R
Helper function for cross-validation analysis.



### predict.gcrq.inputform.R
Helper function for working with gcrq objects (quantile regression models).



### predict.gcrq.lpmatrix.R
Helper function fo extracting linear predictor matrices from gcrq objects (quantile regression models).


