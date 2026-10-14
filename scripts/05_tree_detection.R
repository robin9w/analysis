# ==============================================================================
# INDIVIDUAL TREE DETECTION (ITD) & CANOPY DELINEATION - MKRF 2016 (t = 0)
# Title: 2016_itd
# This script runs inside the 'analysis/' folder using relative paths.
# ==============================================================================

library(lidR)
library(terra)
library(sf)
library(future)
library(ForestTools)

# 1. PARALLEL SETUP: Configure strict 8-core VDI utilization limit [1]
plan(multisession, workers = 8)

# 2. FILE PATHS: Relative to the analysis directory
laz_dir   <- "data/lidar/2016_laz/"
chm_path  <- "data/raster/2016_chm_30cm.tif"
vri_path  <- "data/vector/2016_mkrf_vri.gpkg"
gpkg_out  <- "data/vector/2016_itd_results.gpkg"
checkpoint_dir <- "data/vector/itd_checkpoints/"

# Create the checkpoint directory if it doesn't exist yet
if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)

# ==============================================================================
# 🌲 STEP 1: Load and Configure the Point Cloud Catalog
# ==============================================================================
message("Loading point cloud catalog...")
ctg <- readLAScatalog(laz_dir)

# Configure memory-safe streaming parameters for the VDI
opt_chunk_size(ctg)   <- 300   # Process in 300m x 300m chunks
opt_chunk_buffer(ctg) <- 15    # 15m buffer prevents edge artifacts between chunks
opt_select(ctg)       <- "xyz"  # Drop intensity/returns arrays to clear RAM headroom
opt_progress(ctg)     <- TRUE  # Enables dynamic grid rendering during execution

# Set up the native lidR driver to write vector shapes as GeoPackages (.gpkg)
ctg@output_options$drivers$sf$extension <- ".gpkg"

# Configure real-time save-to-disk checkpointing paths using coordinate templates
opt_output_files(ctg) <- paste0(checkpoint_dir, "ttops_{XLEFT}_{YBOTTOM}")

# ==============================================================================
# 🗺️ STEP 2: Extract Young Stands & Mask the Catalog Search Space
# ==============================================================================
message("Loading and parsing MKRF VRI multi-layer structure...")
vri_raw <- vect(vri_path, layer = "mkrf_fc_2016")

# Build query tracking both Layer 1 (Overstory) and Layer 2 (Understory cohort)
young_condition <- (vri_raw$TCL_1_TSP_1_AGE < 20 & !is.nan(vri_raw$TCL_1_TSP_1_AGE)) | 
  (vri_raw$TCL_2_TSP_1_AGE < 20 & !is.nan(vri_raw$TCL_2_TSP_1_AGE))

vri_young_terra <- vri_raw[young_condition, ]
vri_young       <- st_as_sf(vri_young_terra)

# Ensure CRS alignment between the point cloud catalog and VRI geometries
if (st_crs(ctg) != st_crs(vri_young)) {
  vri_young <- st_transform(vri_young, st_crs(ctg))
}

# Combine young stands into a single spatial footprint mask
message("Dissolving young stand footprints into a spatial mask...")
vri_footprint <- st_union(vri_young)

# Intersect the catalog with the boundary footprint
ctg_subset <- catalog_intersect(ctg, vri_footprint)

# ==============================================================================
# 🚀 STEP 3: Run Checkpointed Tree Detection 
# ==============================================================================
message("Locating tree tops across young stands with real-time checkpointing...")
locate_trees(ctg_subset, lmf(ws = 3, hmin = 2))

# ==============================================================================
# 📥 STEP 4: Collapse Checkpoint Layers & Run Watershed Delineation
# ==============================================================================
message("Collapsing written checkpoint files back into a unified dataset...")
# Scans for either .gpkg or .shp formats to ensure robust recovery execution
checkpoint_files <- list.files(checkpoint_dir, pattern = "\\.(gpkg|shp)$", full.names = TRUE)

if (length(checkpoint_files) == 0) {
  stop("🚨 Critical Error: No checkpoint files found in your checkpoint directory!")
}

ttops_list       <- lapply(checkpoint_files, st_read, quiet = TRUE)
ttops_all        <- do.call(rbind, ttops_list)

# Explicitly ensure the VRI dataset matches the exact CRS of the loaded checkpoints
message("Verifying coordinate reference systems alignment...")
if (st_crs(ttops_all) != st_crs(vri_young)) {
  vri_young <- st_transform(vri_young, st_crs(ttops_all))
}

message("Cleaning buffer redundancies and filtering to true young stand footprints...")
ttops <- st_intersection(ttops_all, st_geometry(vri_young))

# FIX: Reloading the raster freshly right here creates a healthy, stable pointer 
# that completely prevents any multi-threading memory address corruption bugs.
message("Loading 30cm CHM proxy and clearing stale pointer references...")
chm <- rast(chm_path)

message("Running marker-controlled watershed segmentation on young stands...")
crown_raster <- ForestTools::mcws(ttops, chm, minHeight = 2) 

message("Vectorizing pixel groups into dripline polygons... (This step will take time)")
crown_polygons <- as.polygons(crown_raster, values = TRUE, dissolve = TRUE)
crown_sf       <- st_as_sf(crown_polygons)
names(crown_sf) <- "tree_id"
crown_sf       <- crown_sf[crown_sf$tree_id > 0, ]

# ==============================================================================
# 💾 STEP 5: Export Derived Data Layers to Final GeoPackage
# ==============================================================================
message("Writing consolidated vector layers to GeoPackage...")
st_write(crown_sf, dsn = gpkg_out, layer = "dripline_polygons", delete_layer = TRUE)
st_write(st_as_sf(ttops), dsn = gpkg_out, layer = "tree_tips", delete_layer = TRUE)

message("--- PROCESS COMPLETE ---")
message("Driplines and Tree Tips successfully stored in: ", gpkg_out)
