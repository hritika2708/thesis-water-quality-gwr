**River Water Quality: Linear vs Geographically Weighted Regression**

R scripts to explore and model river-water quality (COD as the main response) using Linear Regression (LM) and Geographically Weighted Regression (GWR), with optional missing-data imputation, spatial diagnostics, and exportable summaries/plots.

**Requirements**

R ≥ 4.2
Install once:

install.packages(c(
  "sf","tidyverse","ggplot2","corrplot","naniar","VIM",
  "sp","spdep","spgwr","GWmodel",
  "mice","miceadds","MuMIn","car","broom","readr","tidyr","dplyr"
))

**Data**

Input: river-quality/Monthly_River_Quality_Data_SDCC.geojson  
CRS (spatial work): EPSG:2157 (Irish Transverse Mercator / ITM)  

**Outputs**

Plots (if enabled) → plots/  
Pairwise missingness PDF → pairwise_marginplots.pdf  
Optional CSVs from evaluate_gwr_model_csv()
