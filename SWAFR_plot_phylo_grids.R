

## ENVIRONMENT SETTINGS =============================================================


# \
# 
# This code prepares all the data and code needed to plot grids in 
#   
#   
#   \


## Set env
rm(list = ls())
#if (!Sys.getenv("JAVA_TOOL_OPTIONS")) {
if (all(Sys.getenv("JAVA_HOME")=="")) {
  Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_321')
}

if (all(Sys.getenv("JAVA_TOOL_OPTIONS")=="")) {
  options(java.parameters = "-Xmx64G")
}

options(warn=0)

## Function to load or install packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos = "https://cran.csiro.au/")
  sapply(pkg, require, character.only = TRUE)
}

'%!in%' <- function(x,y)!('%in%'(x,y))


## Load packages
#devtools::install_github("HMB3/habitatIntersect")
library(habitatIntersect)
data('sdmgen_packages')
ipak(sdmgen_packages) 


## The functions expect these folders,
main_dir             <- paste0(getwd(), "/")
tempdir              <- './TEMP/'
swafr_data           <- './data/SWAFR_data/'
swafr_results        <- './data/SWAFR_data/Results/'
aus_results          <- './data/SWAFR/'

swafr_out_dir        <- './output/SWAFR/'
aus_out_dir          <- './output/AUS/'


## Try and set the raster temp directory to a location not on the partition, to save space
rasterOptions(memfrac = 0.9,
              tmpdir  = tempdir)

terraOptions(memfrac = 0.9, 
             tempdir = tempdir) 




# STEP 1 :: Get spp lists ----


## To Do :
## 1). Check the errors that Finlay may have found
## 2). Clean up the folders


## get target taxa
swafr_taxa <- c('Acacia',   'Adenanthos', 'Banksia', 'Calytrix', 'Daviesia', 'Epacrids', 
                'Eucalypt', 'Persoonia',  'Thysanotus')





# STEP 2 :: Combine taxa occurrence data ----


## 15km grids
All.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/All_genera/',  pattern ="_15km_", full.names = TRUE))

All.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/All_genera/', pattern ="_trimmed_", full.names = TRUE)) %>% 
  .[!. %in% list.files('./data/SWAFR_data/Results/All_genera/',  pattern ="_15km_", full.names = TRUE)]


Acacia.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Acacia/',  pattern ="_climate_", full.names = TRUE))

Acacia.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Acacia/', pattern ="_trimmed_", full.names = TRUE))


Adenanthos.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Adenanthos/',  pattern ="_15km_", full.names = TRUE))

Adenanthos.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Adenanthos/', pattern ="_trimmed_", full.names = TRUE))


Calytrix.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Calytrix/',  pattern ="_15km_", full.names = TRUE))

Calytrix.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Calytrix/', pattern ="_trimmed_", full.names = TRUE)) %>% 
  .[!. %in% list.files('./data/SWAFR_data/Results/Calytrix/',  pattern ="_15km_", full.names = TRUE) ]


Daviesia.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Daviesia/',  pattern ="_15km_", full.names = TRUE))

Daviesia.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Daviesia/', pattern ="_trimmed_", full.names = TRUE)) 


Epacrids.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Epacrids/',  pattern ="_15km_", full.names = TRUE))

Epacrids.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Epacrids/', pattern ="_trimmed_", full.names = TRUE)) 


Eucalypt.clim.grids.15km <- raster::stack(
  list.files('./data/SWAFR_data/Results/Eucalypt/',  pattern ="_15km_", full.names = TRUE))

Eucalypt.phylo.grids.15km  <- raster::stack(
  list.files('./data/SWAFR_data/Results/Eucalypt/', pattern ="_trimmed_", full.names = TRUE)) 


## Rename the grids because I stuffed up in Biodiverse
names(Acacia.clim.grids.15km)      <- gsub("SWAFR_climate_",           "",  names(Acacia.clim.grids.15km))
names(Acacia.phylo.grids.15km)     <- gsub("SWAFR_epsg_3577_trimmed_", "",  names(Acacia.phylo.grids.15km))

names(Acacia.clim.grids.15km)      <- gsub("SWAFR_climate_",           "",  names(Acacia.clim.grids.15km))
names(Acacia.phylo.grids.15km)     <- gsub("SWAFR_epsg_3577_trimmed_", "",  names(Acacia.phylo.grids.15km))

names(Adenanthos.clim.grids.15km)  <- gsub("_SWAFR__15km__",           "_", names(Adenanthos.clim.grids.15km))
names(Adenanthos.phylo.grids.15km) <- gsub("SWAFR_epsg_3577_trimmed_", "",  names(Adenanthos.phylo.grids.15km))

names(Calytrix.clim.grids.15km)    <- gsub("_SWAFR__15km__",           "_", names(Calytrix.clim.grids.15km))
names(Calytrix.phylo.grids.15km)   <- gsub("SWAFR_epsg_3577_trimmed_", "",  names(Calytrix.phylo.grids.15km))



## Now remove the individual rasters
# rm(list = ls(pattern = '_'))





# STEP 3 :: extract environmental values ----

## The below shold really all be one big table
## Where we use different columns for the analysis :: species, genus or family
# Insects_ALA_1 <- read_csv('./data/ALA/Insects/recordsd-2022-05-12.csv')
# Insects_ALA_2 <- read_csv('./data/ALA/Insects/records-2022-05-12_part2.csv')



# STEP 4 :: pLOT ----
# Create all the possible combinations of the two lists
acacia_mean_grids <- expand.grid(names(Acacia.clim.grids.15km), 
                                 names(Acacia.phylo.grids.15km)) %>% 
  
  rename(raster1 = Var1,
         raster2 = Var2) %>% 
  dplyr::filter(grepl("mean_MEAN", raster1))


raster_combo_scatters(plot_list            = acacia_mean_grids,
                      climate_raster_stack = Acacia.clim.grids.15km,
                      context_raser_stack  = Acacia.phylo.grids.15km,
                      out_dir              = paste0(swafr_out_dir, 'Acacia/'))





## END =============================================================




