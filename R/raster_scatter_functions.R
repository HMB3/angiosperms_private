#########################################################################################################################
################################### FUNCTIONS FOR PHYLO PLOTS ---- ############################################
#########################################################################################################################



'%!in%' <- function(x,y)!('%in%'(x,y))

# function to create a scatterplot from two rasters, with one point per grid cell
# (this will return an error if the two rasters do not have identical spatial attributes)
# a: file path to raster plotted on x axis
# b: file path to raster plotted on y axis
# returns a ggplot that can be further modified as needed
raster_scatter <- function(raster1, 
                           raster2){
  
  require(raster)
  require(tidyverse)
  require(tools)
  
  stack(raster1, raster2) %>%
    
    raster::rasterToPoints() %>%
    as.data.frame()  %>% na.omit() %>% 
    
    ggplot(aes_string(names(raster1), names(raster2))) +
    geom_point() 
  
}


## Create a list of histogram plots----
raster_combo_scatters <- function(plot_list, 
                                  climate_raster_stack,
                                  context_raser_stack,
                                  out_dir) {
  
  
  rows   <- nrow(plot_list)
  length <- 1:rows
  
  ## Pipe the list into Lapply
  for(i in length) {
    
    ## Pipe the list into lapply
    ## row <- length[1]
    
    ## i <- length[10]
    message("Creating scatterplot for combo ", i, '/', rows)
    
    clim_ras <- plot_list[i,][["raster1"]] %>% as.character()
    cont_ras <- plot_list[i,][["raster2"]] %>% as.character()
    
    ## Now pipe the table into the function
    # example: generate basic plot, dress it up a bit with extra ggplot code, and save
    p <- raster_scatter(raster1 = climate_raster_stack[[clim_ras]],
                        raster2 = context_raser_stack[[cont_ras]]) +
      geom_smooth() + # adding a trend line - how about R2
      theme_classic() + # changing the plot style
      labs(x = clim_ras, # adding axis labels
           y = cont_ras)
    
    png(paste0(out_dir, clim_ras, '_VS_',  cont_ras, '.png'),
        16, 10, units = 'in', res = 400)
    plot(p)
    dev.off()
    
  } 
}


