# ==============================================================================
# CORRECTED BATCH POINT-CLOUD SEGMENTATION PIPELINE (Li 2012)
# ==============================================================================
library(lidR)
library(sf)
library(dplyr)
library(future)

# 1. Setup Environment and Paths ----
gpkg_path  <- "data/vector/mkrf_vri.gpkg"
lidar_dir  <- "data/lidar/2016_laz" 
output_dir <- "data/lidar/segmented" 

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Establish parallel multisession (8 workers)
plan(multisession, workers = 8)

# 2. Load and Validate Young Stand Footprint ----
young_stands <- sf::st_read(gpkg_path, layer = "mkrf_fc_2016") %>%
  dplyr::filter(TCL_1_TSP_1_AGE > 0 & TCL_1_TSP_1_AGE <= 20) %>%
  sf::st_make_valid()

# 3. Initialize and Configure the LAScatalog ----
mkrf_catalog <- lidR::readLAScatalog(lidar_dir)

# Ensure spatial consistency with BC Albers (EPSG:3153)
if (is.na(sf::st_crs(mkrf_catalog))) sf::st_crs(mkrf_catalog) <- 3153
young_stands <- sf::st_transform(young_stands, sf::st_crs(mkrf_catalog))

# Catalog Engine Constraints for Point Cloud Manipulation
opt_select(mkrf_catalog)       <- "*"    # Li2012 requires all geometric attributes
opt_chunk_size(mkrf_catalog)   <- 250    # Controlled memory foot-print for dense young stands
opt_chunk_buffer(mkrf_catalog) <- 20     # Buffering stops tree splitting artifacts on edges

# Target file naming template for segmented tiles
opt_output_files(mkrf_catalog) <- file.path(output_dir, "segmented_chunk_{XLEFT}_{YBOTTOM}")

# 4. FIX: Mask the Catalog Spatially ----
# Dropping `opt_filter` command strings. Instead, intersect the catalog grid
# boundaries directly with your young stands geometry. This filters out non-target
# tiles before sending the collection to the execution block.
mkrf_catalog <- lidR::clip_roi(mkrf_catalog, young_stands)

# 5. Execute Li 2012 Algorithm via Catalog Engine ----
# Tuned Parameters for Young Stands:
# hmin = 2                : Ignores low brush/understory under 2 metres.
# speed_up = 5            : Limits search radius to 5m to avoid n^2 processing traps.
# uniqueness = "bitmerge" : Stitche tree points seamlessly on chunk lines.
segmented_catalog <- lidR::segment_trees(
  las = mkrf_catalog,
  algorithm = li2012(dt1 = 1.5, dt2 = 2.0, R = 2, Zu = 15, hmin = 2, speed_up = 5),
  uniqueness = "bitmerge"
)

print("Segmentation Pipeline Complete! Output files saved to data/lidar/segmented/")
