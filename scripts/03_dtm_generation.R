# ==============================================================================
# LOW-RAM ONDISK BATCH GENERATION OF A UNIFIED RASTER DTM
# ==============================================================================
library(lidR)
library(terra)
library(future)

# 1. Setup Environment and Paths ----
lidar_dir  <- "data/lidar/2016_laz" 
raster_dir <- "data/raster"
chunk_dir  <- file.path(raster_dir, "chunks") # Dedicated subfolder for tiles

if (!dir.exists(chunk_dir)) dir.create(chunk_dir, recursive = TRUE)
output_dtm <- file.path(raster_dir, "2016_dtm_30cm.tif")

# Fast socket adjustments for local multi-node processing
plan(multisession, workers = 8) 

# 2. Load the entire folder as a LAScatalog ----
mkrf_catalog <- lidR::readLAScatalog(lidar_dir)
if (is.na(sf::st_crs(mkrf_catalog))) sf::st_crs(mkrf_catalog) <- 3153

# 3. Configure Catalog Processing Parameters ----
opt_select(mkrf_catalog)       <- "xyzc"  
opt_chunk_size(mkrf_catalog)   <- 0       # Process natively by file (faster, less overhead)
opt_chunk_buffer(mkrf_catalog) <- 30      # 30m buffer for tin() stability

# CRITICAL FIX 1: Write chunks straight to disk to save your RAM from crashing!
opt_output_files(mkrf_catalog) <- file.path(chunk_dir, "dtm_chunk_{ID}")

# 4. Rasterize the Entire Terrain Natively ----
print("Starting batch processing of ground points directly to disk...")
lidR::rasterize_terrain(
  las       = mkrf_catalog, 
  res       = 0.3, 
  algorithm = tin()
)

# 5. Compile and Save Master GeoTIFF via VRT (Virtual Raster) ----
print("Stitching files via low-memory Virtual Raster (VRT)...")
chunk_files <- list.files(chunk_dir, pattern = "dtm_chunk_", full.names = TRUE)

# Create a virtual layout map of all tiles (takes 1 second and uses no RAM)
vrt_layout <- terra::vrt(chunk_files)

print("Writing final unified master DTM to your hard drive...")
# Write the virtual raster out to a single physical GeoTIFF
terra::writeRaster(vrt_layout, output_dtm, overwrite = TRUE)

# Clean up the intermediate chunk subfolder
unlink(chunk_dir, recursive = TRUE)

print(paste("Success! The unified master DTM raster has been compiled at:", output_dtm))
