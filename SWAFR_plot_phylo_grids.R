## ENVIRONMENT SETTINGS =============================================================


# \
# 
# This code prepares all the data and code needed for the analysis of inverts habitat after the 2019-2020 fires ::
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
    install.packages(new.pkg, dependencies = TRUE, repos = "https://cran.csiro.au/")
  sapply(pkg, require, character.only = TRUE)
}


## Load packages
#devtools::install_github("HMB3/habitatIntersect")
library(habitatIntersect)
data('sdmgen_packages')
ipak(sdmgen_packages)


## The functions expect these folders,
main_dir             <- paste0(getwd(), "/")
SWAFR_out_dir        <- './output/SWAFR/'


# STEP 1 :: Get spp lists ----


## To Do :
## 1). Check the errors that Finlay may have found
## 2). Clean up the folders
source('./R/raster_scatter_functions.R')


## Create taxa lists here
# Insects_ALA_1 <- read_csv('./data/ALA/Insects/recordsd-2022-05-12.csv')
# Insects_ALA_2 <- read_csv('./data/ALA/Insects/records-2022-05-12_part2.csv')
swafr_taxa <- c('Acacia',   'Adenanthos', 'Banksia', 'Calytrix', 'Daviesia', 'Epacrids', 
                'Eucalypt', 'Persoonia',  'Thysanotus')





# STEP 2 :: Create raster grids ----


## 15km grids
All.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/All_genera/',  pattern ="_15km_", full.names = TRUE))

All.clim.15km <- 
  list.files('./data/SWAFR_data/Results/All_genera/',  pattern ="_15km_", full.names = TRUE)

All.phylo.15km  <-
  list.files('./data/SWAFR_data/Results/All_genera/',  pattern ="_trimmed_", full.names = TRUE) %>% 
  .[!. %in% All.clim.15km]

All.phylo.grids.15km <- raster::stack(All.phylo.15km)


Acacia.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Acacia/',      pattern = "_climate_", full.names = TRUE))

Acacia.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Acacia/',      pattern = "_trimmed_", full.names = TRUE))


Adenanthos.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Adenanthos/',  pattern = "_15km_",    full.names = TRUE))

Adenanthos.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Adenanthos/',  pattern = "_trimmed_", full.names = TRUE))


Banksia.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Banksia/',     pattern = "_15km_",    full.names = TRUE))

Banksia.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Banksia/',     pattern = "_trimmed_", full.names = TRUE))


Calytrix.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Calytrix/',    pattern = "_15km_",    full.names = TRUE))

Calytrix.clim.15km <- list.files('./data/SWAFR_data/Results/Calytrix/',    pattern = "_15km_",    full.names = TRUE)

Calytrix.phylo.15km  <- 
  list.files('./data/SWAFR_data/Results/Calytrix/',    pattern = "_trimmed_", full.names = TRUE) %>% 
  .[!. %in% Calytrix.clim.15km ]


Calytrix.phylo.grids.15km <- raster::stack(Calytrix.phylo.15km)


Daviesia.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Daviesia/',    pattern = "_15km_", full.names = TRUE))

Daviesia.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Daviesia/',    pattern = "_trimmed_", full.names = TRUE)) 


Epacrids.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Epacrids/',    pattern ="_15km_", full.names = TRUE))

Epacrids.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Epacrids/',    pattern ="_trimmed_", full.names = TRUE)) 


Eucalypt.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Eucalypt/',    pattern ="_15km_", full.names = TRUE))

Eucalypt.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Eucalypt/',    pattern ="_trimmed_", full.names = TRUE))


Persoonia.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Persoonia/',   pattern ="_15km_", full.names = TRUE))

Persoonia.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Persoonia/',   pattern ="_trimmed_", full.names = TRUE))


Thysanotus.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Thysanotus/',  pattern ="_15km_", full.names = TRUE))

Thysanotus.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Thysanotus/',  pattern ="_trimmed_", full.names = TRUE))


## Combine the climate rasters for each taxa together
## Re-sample


## Rename the grids because I stuffed up in Biodiverse
names(All.clim.grids.15km)         <- gsub("_SW_clipped_trimmed__15km_",  "",  names(All.clim.grids.15km))
names(All.phylo.grids.15km)        <- gsub("SW_clipped_trimmed_",         "",  names(All.phylo.grids.15km))

names(Acacia.clim.grids.15km)      <- gsub("SWAFR_climate_",              "",  names(Acacia.clim.grids.15km))
names(Acacia.phylo.grids.15km)     <- gsub("SWAFR_epsg_3577_trimmed_",    "",  names(Acacia.phylo.grids.15km))

names(Adenanthos.clim.grids.15km)  <- gsub("_SWAFR__15km__",              "_", names(Adenanthos.clim.grids.15km))
names(Adenanthos.phylo.grids.15km) <- gsub("SWAFR_epsg_3577_trimmed_",    "",  names(Adenanthos.phylo.grids.15km))


names(Banksia.clim.grids.15km)     <- gsub("_SWAFR__15km__",              "_", names(Banksia.clim.grids.15km))
names(Banksia.phylo.grids.15km)    <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Banksia.phylo.grids.15km))


names(Calytrix.clim.grids.15km)    <- gsub("_SWAFR__15km__",              "_", names(Calytrix.clim.grids.15km))
names(Calytrix.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Calytrix.phylo.grids.15km))


names(Daviesia.clim.grids.15km)    <- gsub("_SWAFR__15km__",              "_", names(Daviesia.clim.grids.15km))
names(Daviesia.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Daviesia.phylo.grids.15km))


names(Epacrids.clim.grids.15km)    <- gsub("_SWAFR__15km__",              "_", names(Epacrids.clim.grids.15km))
names(Epacrids.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Epacrids.phylo.grids.15km))


names(Eucalypt.clim.grids.15km)    <- gsub("_SWAFR__15km__",              "_", names(Eucalypt.clim.grids.15km))
names(Eucalypt.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Eucalypt.phylo.grids.15km))


names(Persoonia.clim.grids.15km)   <- gsub("_SWAFR__15km__",              "_", names(Persoonia.clim.grids.15km))
names(Persoonia.phylo.grids.15km)  <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Persoonia.phylo.grids.15km))


names(Thysanotus.clim.grids.15km)  <- gsub("_SWAFR__15km__",              "_", names(Thysanotus.clim.grids.15km))
names(Thysanotus.phylo.grids.15km) <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Thysanotus.phylo.grids.15km))



## Now remove the individual rasters
# rm(list = ls(pattern = '_'))





# STEP 3 :: create plots ----


## Create all the possible combinations of the two lists
## Create a loop to do this...
## All genera 
All_mean_grids <- expand.grid(names(All.clim.grids.15km), 
                              names(All.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Acacia
Acacia_mean_grids <- expand.grid(names(Acacia.clim.grids.15km), 
                                 names(Acacia.phylo.grids.15km)) %>% 
  
  ## Just get the 'mean mean' layers
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Adenanthos
Adenanthos_mean_grids <- expand.grid(names(Adenanthos.clim.grids.15km), 
                                     names(Adenanthos.phylo.grids.15km)) %>% 
  
  ## Just get the 'mean mean' layers
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Banksia 
Banksia_mean_grids <- expand.grid(names(Banksia.clim.grids.15km), 
                                  names(Banksia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Calytrix 
Calytrix_mean_grids <- expand.grid(names(Calytrix.clim.grids.15km), 
                                   names(Calytrix.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Daviesia 
Calytrix_mean_grids <- expand.grid(names(Daviesia.clim.grids.15km), 
                                   names(Daviesia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Epacrids 
Calytrix_mean_grids <- expand.grid(names(Epacrids.clim.grids.15km), 
                                   names(Epacrids.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Eucalypt 
Calytrix_mean_grids <- expand.grid(names(Eucalypt.clim.grids.15km), 
                                   names(Eucalypt.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Persoonia 
Persoonia_mean_grids <- expand.grid(names(Persoonia.clim.grids.15km), 
                                    names(Persoonia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Thysanotus 
Thysanotus_mean_grids <- expand.grid(names(Thysanotus.clim.grids.15km), 
                                     names(Thysanotus.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


## Update function toon a list of taxa too
raster_combo_scatters(plot_list            = all_mean_grids,
                      climate_raster_stack = All.clim.grids.15km,
                      context_raser_stack  = All.phylo.grids.15km,
                      out_dir              = paste0(swafr_out_dir, 'All_genera/'))


raster_combo_scatters(plot_list            = Banksia_mean_grids,
                      climate_raster_stack = Banksia.clim.grids.15km,
                      context_raser_stack  = Banksia.phylo.grids.15km,
                      out_dir              = paste0(swafr_out_dir, 'Banksia/'))


## END =============================================================




