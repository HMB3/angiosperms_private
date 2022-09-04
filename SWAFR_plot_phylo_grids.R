#########################################################################################################################
################################### -------- PHYLO PLOTS ------ #########################################################
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


names(Calytrix.clim.grids.15km)    <- gsub("_SWAFR_epsg_3577_trimmed__15km__",  "_", names(Calytrix.clim.grids.15km))
names(Calytrix.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_",          "", names(Calytrix.phylo.grids.15km))


names(Daviesia.clim.grids.15km)    <- gsub("_SWAFR__15km__",              "_", names(Daviesia.clim.grids.15km))
names(Daviesia.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Daviesia.phylo.grids.15km))


names(Epacrids.clim.grids.15km)    <- gsub("_SWAFR__15km__",              "_", names(Epacrids.clim.grids.15km))
names(Epacrids.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Epacrids.phylo.grids.15km))


names(Eucalypt.clim.grids.15km)    <- gsub("_SWAFR__15km__",              "_", names(Eucalypt.clim.grids.15km))
names(Eucalypt.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Eucalypt.phylo.grids.15km))


names(Persoonia.clim.grids.15km)   <- gsub("_SWAFR__15km__",              "_", names(Persoonia.clim.grids.15km))
names(Persoonia.phylo.grids.15km)  <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Persoonia.phylo.grids.15km))


names(Thysanotus.clim.grids.15km)  <- gsub("SWAFR_15km__",              "", names(Thysanotus.clim.grids.15km))
names(Thysanotus.phylo.grids.15km) <- gsub("SWAFR_epsg_3577_trimmed_",     "", names(Thysanotus.phylo.grids.15km))



## Some of these rasters are wrong...
# rm(list = ls(pattern = '_'))
plot(All.clim.grids.15km[['Australian_genera_Annual_mean_temp_mean_MEAN']])




# STEP 3 :: create grid lists ----


## Create all the possible combinations of the two lists
## Create a loop to do this...
## All genera 
mean_match     <- c("mean_MEAN", "temp", "precip")
All_mean_grids <- expand.grid(names(All.clim.grids.15km), 
                              names(All.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 


All_max_grids <- expand.grid(names(All.clim.grids.15km), 
                             names(All.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 



## Acacia
Acacia_mean_grids <- expand.grid(names(Acacia.clim.grids.15km), 
                                 names(Acacia.phylo.grids.15km)) %>% 
  
  ## Just get the 'mean mean' layers
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip"))


Acacia_max_grids <- expand.grid(names(Acacia.clim.grids.15km), 
                                names(Acacia.phylo.grids.15km)) %>% 
  
  ## Just get the 'mean mean' layers
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip"))


## Adenanthos
Adenanthos_mean_grids <- expand.grid(names(Adenanthos.clim.grids.15km), 
                                     names(Adenanthos.phylo.grids.15km)) %>% 
  
  ## Just get the 'mean mean' layers
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 


Adenanthos_max_grids <- expand.grid(names(Adenanthos.clim.grids.15km), 
                                    names(Adenanthos.phylo.grids.15km)) %>% 
  
  ## Just get the 'mean mean' layers
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 



## Banksia 
Banksia_mean_grids <- expand.grid(names(Banksia.clim.grids.15km), 
                                  names(Banksia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 


Banksia_max_grids <- expand.grid(names(Banksia.clim.grids.15km), 
                                 names(Banksia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 


## Calytrix 
Calytrix_mean_grids <- expand.grid(names(Calytrix.clim.grids.15km), 
                                   names(Calytrix.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>%
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip"))


Calytrix_max_grids <- expand.grid(names(Calytrix.clim.grids.15km), 
                                  names(Calytrix.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>%
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 



## Daviesia 
Daviesia_mean_grids <- expand.grid(names(Daviesia.clim.grids.15km), 
                                   names(Daviesia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>%
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip"))


Daviesia_max_grids <- expand.grid(names(Daviesia.clim.grids.15km), 
                                  names(Daviesia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>%
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 



## Epacrids 
Epacrids_mean_grids <- expand.grid(names(Epacrids.clim.grids.15km), 
                                   names(Epacrids.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>%
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip"))


Epacrids_max_grids <- expand.grid(names(Epacrids.clim.grids.15km), 
                                  names(Epacrids.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>%
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip"))


## Eucalypt 
Eucalypt_mean_grids <- expand.grid(names(Eucalypt.clim.grids.15km), 
                                   names(Eucalypt.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip"))


Eucalypt_max_grids <- expand.grid(names(Eucalypt.clim.grids.15km), 
                                  names(Eucalypt.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 


## Persoonia 
Persoonia_mean_grids <- expand.grid(names(Persoonia.clim.grids.15km), 
                                    names(Persoonia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip"))


Persoonia_max_grids <- expand.grid(names(Persoonia.clim.grids.15km), 
                                   names(Persoonia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 


## Thysanotus 
Thysanotus_mean_grids <- expand.grid(names(Thysanotus.clim.grids.15km), 
                                     names(Thysanotus.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("mean_MEAN", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 


Thysanotus_max_grids <- expand.grid(names(Thysanotus.clim.grids.15km), 
                                    names(Thysanotus.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  
  dplyr::filter(grepl("max_MAX", raster1)) %>% 
  dplyr::filter(str_detect(raster1, "temp|precip")) 




# STEP 4 :: create grid plots ----


## Update function toon a list of taxa too
# for(taxa in swafr_taxa)
options(scipen=10000)

raster_combo_scatters(plot_list            = All_mean_grids,
                      climate_raster_stack = All.clim.grids.15km,
                      context_raser_stack  = All.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'All_genera/'))


raster_combo_scatters(plot_list            = Acacia_mean_grids,
                      climate_raster_stack = Acacia.clim.grids.15km,
                      context_raser_stack  = Acacia.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Acacia/'))


raster_combo_scatters(plot_list            = Acacia_max_grids,
                      climate_raster_stack = Acacia.clim.grids.15km,
                      context_raser_stack  = Acacia.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Acacia/'))


raster_combo_scatters(plot_list            = Adenanthos_mean_grids,
                      climate_raster_stack = Adenanthos.clim.grids.15km,
                      context_raser_stack  = Adenanthos.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Adenanthos/'))


raster_combo_scatters(plot_list            = Adenanthos_max_grids,
                      climate_raster_stack = Adenanthos.clim.grids.15km,
                      context_raser_stack  = Adenanthos.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Adenanthos/'))


raster_combo_scatters(plot_list            = Banksia_mean_grids,
                      climate_raster_stack = Banksia.clim.grids.15km,
                      context_raser_stack  = Banksia.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Banksia/'))


raster_combo_scatters(plot_list            = Banksia_max_grids,
                      climate_raster_stack = Banksia.clim.grids.15km,
                      context_raser_stack  = Banksia.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Banksia/'))


raster_combo_scatters(plot_list            = Calytrix_mean_grids,
                      climate_raster_stack = Calytrix.clim.grids.15km,
                      context_raser_stack  = Calytrix.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Calytrix/'))


raster_combo_scatters(plot_list            = Calytrix_max_grids,
                      climate_raster_stack = Calytrix.clim.grids.15km,
                      context_raser_stack  = Calytrix.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Calytrix/'))


raster_combo_scatters(plot_list            = Daviesia_mean_grids,
                      climate_raster_stack = Daviesia.clim.grids.15km,
                      context_raser_stack  = Daviesia.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Daviesia/'))


raster_combo_scatters(plot_list            = Daviesia_max_grids,
                      climate_raster_stack = Daviesia.clim.grids.15km,
                      context_raser_stack  = Daviesia.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Daviesia/'))


raster_combo_scatters(plot_list            = Epacrids_mean_grids,
                      climate_raster_stack = Epacrids.clim.grids.15km,
                      context_raser_stack  = Epacrids.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Epacrids/'))


raster_combo_scatters(plot_list            = Epacrids_max_grids,
                      climate_raster_stack = Epacrids.clim.grids.15km,
                      context_raser_stack  = Epacrids.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Epacrids/'))


raster_combo_scatters(plot_list            = Eucalypt_mean_grids,
                      climate_raster_stack = Eucalypt.clim.grids.15km,
                      context_raser_stack  = Eucalypt.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Eucalypt/'))


raster_combo_scatters(plot_list            = Eucalypt_max_grids,
                      climate_raster_stack = Eucalypt.clim.grids.15km,
                      context_raser_stack  = Eucalypt.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Eucalypt/'))


raster_combo_scatters(plot_list            = Persoonia_mean_grids,
                      climate_raster_stack = Persoonia.clim.grids.15km,
                      context_raser_stack  = Persoonia.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Persoonia/'))


raster_combo_scatters(plot_list            = Persoonia_max_grids,
                      climate_raster_stack = Persoonia.clim.grids.15km,
                      context_raser_stack  = Persoonia.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Persoonia/'))


raster_combo_scatters(plot_list            = Thysanotus_mean_grids,
                      climate_raster_stack = Thysanotus.clim.grids.15km,
                      context_raser_stack  = Thysanotus.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Thysanotus/'))


raster_combo_scatters(plot_list            = Thysanotus_max_grids,
                      climate_raster_stack = Thysanotus.clim.grids.15km,
                      context_raser_stack  = Thysanotus.phylo.grids.15km,
                      col_layer            = 'PHYLO_RPD1',
                      
                      mar                  = 1,
                      xsize                = 20, 
                      ysize                = 20, 
                      lab_size             = 10,
                      out_dir              = paste0(SWAFR_out_dir, 'Thysanotus/'))




#########################################################################################################################
################################### ------------ TBC ---------- #########################################################
#########################################################################################################################