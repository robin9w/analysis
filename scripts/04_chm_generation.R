# ==============================================================================
# LOW-RAM ONDISK BATCH GENERATION OF A CANOPY HEIGHT MODEL (CHM)
# ==============================================================================
library(lidR)
library(terra)
library(future)

# 1. Setup Environment and Paths ----
lidar_dir  <- "data/lidar/2016_laz"
raster_dir <- "data/raster"
chm_chunk_dir <- file.path(raster_dir, "chm_chunks")

if (!dir.exists(chm_chunk_dir)) dir.create(chm_chunk_dir, recursive = TRUE)
output_chm <- file.path(raster_dir, "2016_chm_30cm.tif")

# Fast socket adjustments for local multi-node processing
plan(multisession, workers = 8)

# 2. Load Your Spatial Datasets ----
print("Loading master DTM and LAScatalog...")
mkrf_dtm     <- terra::rast(file.path(raster_dir, "2016_dtm_30cm.tif"))
mkrf_catalog <- lidR::readLAScatalog(lidar_dir)
if (is.na(sf::st_crs(mkrf_catalog))) sf::st_crs(mkrf_catalog) <- 3153

# 3. Configure Catalog Parameters for CHM ----
opt_select(mkrf_catalog)       <- "xyz" # Only need 3D coordinates for canopy structure
opt_chunk_size(mkrf_catalog)   <- 0     # Process natively by file (low overhead)
opt_chunk_buffer(mkrf_catalog) <- 0     # No buffer needed for pit-free point-to-raster

# CRITICAL RAM SAFETY: Save intermediate CHM tiles directly to disk
opt_output_files(mkrf_catalog) <- file.path(chm_chunk_dir, "chm_chunk_{ID}")

# 4. Generate the Canopy Height Model ----
print("Rasterizing canopy heights natively on 8 cores...")
# Using p2r (Point-to-Raster) algorithm at 30cm resolution matching your project scale
lidR::rasterize_canopy(
  las       = mkrf_catalog,
  res       = 0.3,
  algorithm = p2r(subcircle = 0.15) # Fills tiny gaps between needles
)

# 5. Compile Master CHM via VRT ----
print("Compiling CHM layout map via Virtual Raster...")
chunk_files <- list.files(chm_chunk_dir, pattern = "chm_chunk_", full.names = TRUE)
vrt_layout  <- terra::vrt(chunk_files)

# CRITICAL FIX: Align extents to match the master DTM boundary before raster math
print("Aligning surface layout bounds with terrain baseline...")
vrt_cropped <- terra::crop(vrt_layout, mkrf_dtm)

# Normalize the canopy height: CHM = Digital Surface Model (Cropped) - Digital Terrain Model
print("Subtracting terrain baseline to compute true canopy heights...")
final_chm <- vrt_cropped - mkrf_dtm

print("Writing final master CHM GeoTIFF to disk...")
terra::writeRaster(final_chm, output_chm, overwrite = TRUE)

# Clean up temporary tile subfolder
unlink(chm_chunk_dir, recursive = TRUE)
print(paste("Success! Your Canopy Height Model is saved at:", output_chm))
