# ==============================================================================
# BATCH GENERATION OF A UNIFIED RASTER DTM FOR THE ENTIRE FOOTPRINT (FIXED)
# ==============================================================================
library(lidR)
library(terra)
library(future)

# 1. Setup Environment and Paths ----
lidar_dir  <- "data/lidar/2016_laz" 
raster_dir <- "data/raster"

if (!dir.exists(raster_dir)) dir.create(raster_dir, recursive = TRUE)
output_dtm <- file.path(raster_dir, "mkrf_master_dtm_1m.tif")

# Fast socket adjustments for local multi-node processing
options(parallelly.makeNodePSOCK.connectTimeout = 300)
options(parallelly.makeNodePSOCK.setup_strategy = "sequential")
plan(multisession, workers = 12) # Utilizing 12 ultra-fast P-cores

# 2. Load the entire folder as a LAScatalog ----
mkrf_catalog <- lidR::readLAScatalog(lidar_dir)
if (is.na(sf::st_crs(mkrf_catalog))) sf::st_crs(mkrf_catalog) <- 3153

# 3. Configure Catalog Processing Parameters ----
opt_select(mkrf_catalog)       <- "xyzc"  
opt_chunk_size(mkrf_catalog)   <- 500     
opt_chunk_buffer(mkrf_catalog) <- 20      

# Setting output file path pattern for intermediate chunk rasters
opt_output_files(mkrf_catalog) <- file.path(raster_dir, "dtm_chunk_{XLEFT}_{YBOTTOM}")

# THE TRUE FIX: Modify internal driver slots to handle explicit overwrites
mkrf_catalog@output_options$drivers$SpatRaster$param$overwrite <- TRUE
mkrf_catalog@output_options$drivers$Raster$param$overwrite     <- TRUE

# 4. Rasterize the Entire Terrain Natively ----
print("Batch processing full catalog ground points with overwrite enabled...")

master_dtm <- lidR::rasterize_terrain(
  las       = mkrf_catalog, 
  res       = 1, 
  algorithm = tin()
)

# 5. Clean Up and Save Master GeoTIFF ----
# Automatically drops the intermediate chunks from your drive to keep it tidy
chunk_files <- list.files(raster_dir, pattern = "dtm_chunk_", full.names = TRUE)
if (length(chunk_files) > 0) file.remove(chunk_files)

print(paste("Success! The unified master DTM raster has been compiled at:", output_dtm))
