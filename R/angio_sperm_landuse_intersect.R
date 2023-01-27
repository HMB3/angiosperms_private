#########################################################################################################################
################################### -------- PHYLO LANDUSE ------ #######################################################
#########################################################################################################################


# \
# 
# This code prepares all the data and code needed to plot the phylo grids ::
#   
#   
#   \


## Set env
rm(list = ls())
options(warn=0)


## Function to load or install packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, 
                     repos = "https://cran.csiro.au/")
  sapply(pkg, require, character.only = TRUE)
}


## Load packages
#devtools::install_github("HMB3/habitatIntersect")
library(habitatIntersect)
data('sdmgen_packages')
ipak(sdmgen_packages)


## The functions expect these folders,
main_dir       <- paste0(getwd(), "/")
out_dir        <- './output/'


km_conversion      <- 1000
hectare_conversion <- 10000
sqkm_conversion    <- 1000000


message('Use GEOS geometry for sf operations to speed up intersections')
sf_use_s2(FALSE)


ensure_multipolygons <- function(X) {
  
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  st_write(X, tmp1)
  ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- st_read(tmp2)
  st_sf(st_drop_geometry(X), geom = st_geometry(Y))
  
}


# STEP 1 :: Read in Phylo grids  ----


## 10km grids
angio.phylo.grids.10km <- raster::stack(
  list.files('./data/angio_genera_canape_10km/', 
             pattern =".tif", full.names = TRUE))


## read in the land use layer - pre-aggregated
AUS_land_use <- 
  st_read('./data/land_use/Land_Use_AGIND_repiar.gpkg') %>% 
  st_cast(., "POLYGON")


AUS_angio_redundancy <-
  
  st_read('./data/angio_genera_canape_10km/aus_angio_genera_redun.shp') %>% 
  st_set_crs(st_crs(8058)) %>% st_transform(., st_crs(8058)) %>% 
  st_cast(., "POLYGON") %>% 
  filter(KEY == "REDUNDANCY_ALL")


AUS_AGIND_areas <- AUS_land_use %>%  
  
  st_set_crs(st_crs(8058)) %>% st_transform(., st_crs(8058)) %>% 
  
  ## This creates slivers, which need to be removed...
  dplyr::mutate(AGIND_Sq_km = st_area(geometry)/sqkm_conversion,
                AGIND_Sq_km = drop_units(AGIND_Sq_km)) %>% 
  
  dplyr::select(AGIND,
                AGIND_Sq_km,
                geometry)

AUS_AGIND_areas_df <- AUS_AGIND_areas %>% as_tibble() %>% 
  
  select(-geometry) %>% 
  group_by(AGIND) %>% 
  summarise(AGIND_Sq_km = sum(AGIND_Sq_km))


## create a loop that does the intersect of every phylo grid with the main
## land use categories.
angio_canape_NEO_feat <-
  
  terra::as.polygons(terra::rast(angio.phylo.grids.10km[["aus_genera_canape_10km__NEO"]])) %>% 
  st_as_sf() %>% st_set_crs(st_crs(8058)) %>% st_transform(., st_crs(8058)) %>% 
  st_make_valid() %>% st_buffer(., 0.0) %>% 
  rename(CANAPE_NEO = aus_genera_canape_10km__NEO) %>% 
  st_cast(., "POLYGON") %>%
  filter(CANAPE_NEO == 1 )




angio_canape_PALAEO_feat <-
  
  terra::as.polygons(terra::rast(angio.phylo.grids.10km[["aus_genera_canape_10km__PALAEO"]])) %>% 
  st_as_sf() %>% st_set_crs(st_crs(8058)) %>% st_transform(., st_crs(8058)) %>% 
  st_make_valid() %>% st_buffer(., 0.0) %>% 
  rename(CANAPE_PALAEO = aus_genera_canape_10km__PALAEO) %>% st_cast(., "POLYGON") %>%  
  
  # ## This creates slivers, which need to be removed...
  # dplyr::mutate(PALAEO_Sq_km = st_area(geometry)/km_conversion,
  #               PALAEO_Sq_km = drop_units(PALAEO_Sq_km)) %>% 
  # 
  dplyr::select(CANAPE_PALAEO,
                # PALAEO_Sq_km,
                geometry) %>%
  filter(CANAPE_PALAEO == 1 )


angio_canape_MIXED_feat <-
  
  terra::as.polygons(terra::rast(angio.phylo.grids.10km[["aus_genera_canape_10km__MIXED"]])) %>% 
  st_as_sf() %>% st_set_crs(st_crs(8058)) %>% st_transform(., st_crs(8058)) %>% 
  st_make_valid() %>% st_buffer(., 0.0) %>% 
  rename(CANAPE_MIXED = aus_genera_canape_10km__MIXED) %>% st_cast(., "POLYGON")  %>% ## This creates slivers, which need to be removed...
  
  # dplyr::mutate(PALAEO_Sq_km = st_area(geometry)/sqkm_conversion,
  #               PALAEO_Sq_km = drop_units(PALAEO_Sq_km)) %>% 
  
  dplyr::select(CANAPE_MIXED,
                # PALAEO_Sq_km,
                geometry) %>%
  filter(CANAPE_MIXED == 1 )


angio_PD_feat <-
  
  terra::as.polygons(terra::rast(angio.phylo.grids.10km[["aus_genera_canape_10km__PD"]])) %>% 
  st_as_sf() %>% st_set_crs(st_crs(8058)) %>% st_transform(., st_crs(8058)) %>% 
  st_make_valid() %>% st_buffer(., 0.0)   %>% 
  rename(PD = aus_genera_canape_10km__PD) %>% st_cast(., "POLYGON") %>% 
  
  ## This creates slivers, which need to be removed...
  # dplyr::mutate(PD_Sq_km = st_area(geometry)/sqkm_conversion,
  #               PD_Sq_km = drop_units(PD_Sq_km)) %>% 
  
  dplyr::select(PD,
                # PD_Sq_km,
                geometry) %>% 
  filter(PD == 1 )


angio_PD_P_feat <-
  
  terra::as.polygons(terra::rast(angio.phylo.grids.10km[["aus_genera_canape_10km__PD_P"]])) %>% 
  st_as_sf() %>% st_set_crs(st_crs(8058))     %>% st_transform(., st_crs(8058)) %>% 
  st_make_valid() %>% st_buffer(., 0.0)       %>% 
  rename(PD_P = aus_genera_canape_10km__PD_P) %>% st_cast(., "POLYGON") %>% 
  
  ## This creates slivers, which need to be removed...
  # dplyr::mutate(PD_P_Sq_km = st_area(geometry)/sqkm_conversion,
  #               PD_P_Sq_km = drop_units(PD_P_Sq_km)) %>% 
  
  dplyr::select(PD_P,
                # PD_P_Sq_km,
                geometry) %>% 
  filter(PD_P == 1 )



# STEP 2 :: Intersect data  ----


AUS_AGIND_areas_canape_NEO <- st_intersection(AUS_AGIND_areas,
                                              angio_canape_NEO_feat) %>% 
  
  dplyr::mutate(NEO_Sq_km = st_area(geometry)/sqkm_conversion,
                NEO_Sq_km = drop_units(NEO_Sq_km))


AUS_AGIND_areas_canape_MIXED <- st_intersection(AUS_AGIND_areas,
                                                angio_canape_MIXED_feat) %>% 
  
  dplyr::mutate(MIXED_Sq_km = st_area(geometry)/sqkm_conversion,
                MIXED_Sq_km = drop_units(MIXED_Sq_km))


AUS_AGIND_areas_canape_PALAEO <- st_intersection(AUS_AGIND_areas,
                                                 angio_canape_PALAEO_feat) %>% 
  
  dplyr::mutate(PALAEO_Sq_km = st_area(geometry)/sqkm_conversion,
                PALAEO_Sq_km = drop_units(PALAEO_Sq_km))


AUS_AGIND_areas_canape_NEO_df <- AUS_AGIND_areas_canape_NEO %>% as_tibble() %>% 
  
  dplyr::select(-geometry, -AGIND_Sq_km) %>% 
  group_by(AGIND) %>% 
  summarise(NEO_Sq_km      = sum(NEO_Sq_km)) %>% 
  
  left_join(., AUS_AGIND_areas_df, by = "AGIND") %>% 
  
  mutate(NEO_AG_percent = (NEO_Sq_km/AGIND_Sq_km) * 100 %>% round(.))



AUS_AGIND_areas_canape_MIXED_df <- AUS_AGIND_areas_canape_MIXED %>% as_tibble() %>% 
  
  dplyr::select(-geometry, -AGIND_Sq_km) %>% 
  group_by(AGIND) %>% 
  summarise(MIXED_Sq_km = sum(MIXED_Sq_km)) %>% 
  
  left_join(., AUS_AGIND_areas_df, by = "AGIND") %>% 
  
  mutate(MIXED_AG_percent = (MIXED_Sq_km/AGIND_Sq_km) * 100 %>% round(.))


AUS_AGIND_areas_canape_PALAEO_df <- AUS_AGIND_areas_canape_PALAEO %>% as_tibble() %>% 
  
  dplyr::select(-geometry, -AGIND_Sq_km) %>% 
  group_by(AGIND) %>% 
  summarise(PALAEO_Sq_km = sum(PALAEO_Sq_km)) %>% 
  
  left_join(., AUS_AGIND_areas_df, by = "AGIND") %>% 
  
  mutate(PALAEO_AG_percent = (PALAEO_Sq_km/AGIND_Sq_km) * 100 %>% round(.))




## 
AUS_AGIND_areas_angio_redundancy <- st_intersection(AUS_AGIND_areas,
                                                    AUS_angio_redundancy) %>% 
  
  dplyr::mutate(Redun_Sq_km = st_area(geometry)/sqkm_conversion,
                Redun_Sq_km = drop_units(Redun_Sq_km))


AUS_AGIND_areas_angio_redundancy_df <- AUS_AGIND_areas_angio_redundancy %>% 
  as_tibble() %>% 
  group_by(AGIND) %>% 
  summarise(Redun_Sq_km     = sum(Redun_Sq_km),
            Sampling_effort = median(VALUE)) %>% arrange(Sampling_effort)


phylo_agricult_area_combo_df <- 
  
  dam_volume  <- list(AUS_AGIND_areas_canape_NEO_df, 
                      AUS_AGIND_areas_canape_MIXED_df,
                      AUS_AGIND_areas_canape_PALAEO_df) %>% 
  
  reduce(full_join, by = c('AGIND')) %>% 
  select(-AGIND_Sq_km.y, -AGIND_Sq_km.x) %>% 
  select(AGIND, AGIND_Sq_km, everything())



# STEP 3 :: Save data  ----
phylo_layer_list <- c('AUS_land_use',
                      'AUS_AGIND_areas',
                      'AUS_angio_redundancy',
                      'angio_canape_NEO_feat',
                      'angio_canape_PALAEO_feat',
                      'angio_canape_MIXED_feat',
                      'angio_PD_feat',
                      'angio_PD_P_feat',
                      'AUS_AGIND_areas_canape_NEO',
                      'AUS_AGIND_areas_canape_PALAEO',
                      'AUS_AGIND_areas_canape_MIXED')

phylo_database_loc <- paste0(out_dir, 'AUS_Angio_genus_canape_Agriculture.gpkg')


if(file.exists(phylo_database_loc) == TRUE) {
  message('re-create geo-package')
  file.remove(phylo_database_loc)
  
}


for(layer in phylo_layer_list) {
  
  ## layer <- water_source_layer_list[3]
  File_to_Write <- get(layer)
  
  #message('writing ', layer, ' to geo database')
  # arc.write(path     = DPEW_Water_Source_gdb_loc,
  #           data     = File_to_Write,
  #           validate = TRUE)
  
  message('writing ', layer, ' to geo-package')
  st_write(File_to_Write,
           dsn        = phylo_database_loc,
           layer      = layer,
           quiet      = TRUE,
           append     = TRUE)
  
}


save.image('angio_genius_canape_agricult.RData')

## Define CANAPE cells
## Blue = neo   = new diversity (e.g. Plants with Asian origins)
## Red  = paleo = old diversity (e.g. plants with gondwanan origins, suhc as nothofagus)
## Purple = mixed = mixed diversity (e.g. plants with mixed orgins, )


## What does the ST_intersect mean WRT to proportion of diversity on private land?
## Of the total Neo and paleo significant area across Aus -
## grazing has X% of the neo
## grazing has X% of the paleo
## grazing has X% of the mixed


## Then do the redundancy analysis
## Polys, points, etc.
# dam_layers_sf = list() 


## Create list of dam layers ----
# for(i in 1:length(dam_layer_list$name)) {
#   
#   message('reading in farm layer for ', i)
#   dam_layers_sf[[dam_layer_list$name[i]]] <- read_sf(dsn   = farm_dams_dir, 
#                                                      layer = dam_layer_list$name[i])
#   
# }


## These layers are removed...maybe not what we want...
# dam_subset      <- dam_layers_sf[!grepl('hydroareas_off_SO3_above_EXCL', names(dam_layers_sf))] 
# dam_subset_sort <- names(dam_subset) %>% sort()
# EMU_subset_sort <- sub("\\_.*", "", dam_subset_sort) %>% unique()
# gc()




