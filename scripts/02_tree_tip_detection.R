library(lidR)
library(sf)
library(dplyr)

# 1. Load your raw tile and get the clean 3D tree tips
las <- lidR::readLAS("data/lidar/2016_laz/AC10.copc.laz")
tree_tips <- lidR::locate_trees(las, lmf(ws = 3, hmin = 2))

# 2. Convert to an sf vector object and isolate the heights relative to the ground
# (Assuming your point cloud Z represents absolute elevation, we calculate localized height)
# If your point cloud is already normalized (ground = 0), use Z directly as height.
tree_tips_sf <- sf::st_as_sf(tree_tips) %>%
  mutate(
    # ALLOMETRIC ESTIMATE: e.g., radius = 10% of tree height + a baseline of 0.5 metres
    # Adjust the multiplier (0.10) and offset (0.5) based on your field knowledge of MKRF stands
    tree_height = Z - Z_ground_placeholder, # Or just use Z if your cloud is normalized
    dripline_radius = (tree_height * 0.10) + 0.5
  )

# 3. Create the circular buffers dynamically based on each individual tree's calculated radius
# 'by_variable = TRUE' forces sf to read the unique dripline_radius column for each row
synthetic_crowns <- sf::st_buffer(tree_tips_sf, dist = tree_tips_sf$dripline_radius)

# 4. Save out as a clean, fast GeoPackage vector layer
sf::st_write(synthetic_crowns, "data/vector/mkrf_synthetic_crowns.gpkg", delete_layer = TRUE)
