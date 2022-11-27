#########################################################################################################################
################################### -------- PHYLO LANDUSE ------ #########################################################
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




# STEP 1 :: Read in Phylo grids  ----


## 10km grids
angio.phylo.grids.10km <- raster::stack(
  list.files('./data/phylo_diversity/angio_genera_canape_10km/', 
             pattern =".tif", full.names = TRUE))


## read in the land use layer
AUS_land_use <- 
  st_read('./data/land_use/land_use_shapefiles/NLUM_ALUMV8_250m_2015_16_alb_join.shp') %>% 
  st_cast(., "POLYGON")


AUS_Agricult_land_use <- AUS_land_use <-

group_by(AGIND) %>% 
  
  ## This creates slivers, which need to be removed...
  summarize() %>% st_make_valid() %>% st_buffer(., 0.0) %>% 
  nngeo::st_remove_holes() %>% 
  
  dplyr::mutate(Hectares  = st_area(geometry)/hectare_conversion,
                Hectares  = drop_units(Hectares),
                Sq_km     = st_area(geometry)/km_conversion,
                Sq_km     = drop_units(Sq_km)) %>% 
  
  dplyr::select(CHR_Water_Source,
                Hectares,
                Sq_km,
                geometry)

## create a loop that does the intersect of every phylo grid with the main
## land use categories.


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





