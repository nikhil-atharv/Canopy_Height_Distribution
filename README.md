# Canopy_Height_Distribution
Prediction of Canopy Height Distribution using Optical (Sentinel 2), SAR ( C- band &amp; L - band), GEDI L2A (Relative Height Percentile) using Machine Learning Algorithm in Google Earth Engine

**Feature Engineering**
1. Preprocessing Optical Data (Sentinel-2)
   Filter Sentinel-2 images by date, ROI, and low cloud cover (<10%).
   Apply scaling factor.
   Compute vegetation & spectral indices: NDVI, NDRI, NSMI, CIRE, SRRE, MTCI, PVI.
   Land Cover classification from Dynamic World
   
3. Embeddings: Add Satellite embeddings (Alpha Earth).
   The Addition of Google's Alpha Earth Foundation Satellite Embeddings has greatly influenced the Model's Performace (MAE, RMSE, R square)
   
4. Radar Data Integration
   Sentinel-1 (C-band): VH & VV polarizations.
   Compute RVI, NDI, and RI indices.
   Apply focal median smoothing.

   ALOS PALSAR (L-band):
   Compute RI, RVI, and NDI indices.

5. Terrain Variables (SRTM DEM)
   Derive elevation, slope, and aspect from NASA SRTM.

6. GEDI L2A Relative Height Percentile (Training Data)
   Filter GEDI footprints by quality flags.
   Extract rh98 canopy height.

**Model Building**
Training (70%)
Validation (30%)
Model Training (Random Forest Regression)
Train Random Forest (500 trees) using merged features.
Predict canopy height across ROI.

**Model Evaluation**
1. Validate predictions with GEDI footprints.
2. Correlation (R)
3. R²
4. RMSE
5. MAE
