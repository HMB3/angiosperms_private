#########################################################################################################################
################################### ----- FUNCTIONS FOR PHYLO PLOTS ------ ##############################################
#########################################################################################################################


'%!in%' <- function(x,y)!('%in%'(x,y))


## Function to run a CANAPE classification using PRank rasters from Biodiverse
#  skip the super class
#  x = PE original tree 
#  y = PE alternate tree
#  z = RPE
canape_signif_fun <- function(x, y, z){
  
  #  simplify the logic below
  #  any cell with a value across any input
  x[is.na(x) & (!is.na(y) | !is.na(z))] = 0.5
  y[is.na(y) & (!is.na(x) | !is.na(z))] = 0.5
  z[is.na(z) & (!is.na(y) | !is.na(y))] = 0.5
  
  sigxy = x > 0.95 | y > 0.95
  sigzh = z > 0.975
  sigzl = z < 0.025
  
  res = x * 0
  res[sigxy] = 3          # mixed
  res[sigxy & sigzl] = 1  # neo
  res[sigxy & sigzh] = 2  # palaeo
  
  levels(res) = c("non-sig", "neo", "palaeo", "mixed")
  res  
}

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



# function to create a scatterplot from two rasters, with one point per grid cell
# (this will return an error if the two rasters do not have identical spatial attributes)
# a: file path to raster plotted on x axis
# b: file path to raster plotted on y axis
# returns a ggplot that can be further modified as needed
raster_scatter_cols <- function(raster1, 
                                raster2, 
                                rastercol,
                                
                                out_dir,
                                mar,
                                xsize, 
                                ysize, 
                                lab_size) {
  
  require(raster)
  require(tidyverse)
  require(tools)
  
  d <- stack(raster1, raster2, rastercol) %>%
    setNames(c("a", "b", "c")) %>%
    rasterToPoints() %>%
    as.data.frame() %>%
    mutate(group = case_when(c < quantile(rastercol, .25) ~ "low",
                             c > quantile(rastercol, .75) ~ "high",
                             TRUE ~ "mid") %>%
             factor(levels = c("low", "mid", "high")))
  
  d %>%
    ggplot(aes(a, b, color = group)) +
    geom_point() +
    scale_color_manual(values = c("red", "khaki", "blue")) +
    theme_classic() +
    theme(legend.position = "bottom") 
    

}


## Create a list of histogram plots----
raster_combo_scatters <- function(plot_list, 
                                  climate_raster_stack,
                                  context_raser_stack,
                                  col_layer,
                                  mar,
                                  xsize, 
                                  ysize, 
                                  lab_size,
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
    col_ras  <- plot_list[["raster2"]] %>% as.character() %>% .[grepl(col_layer, .)] %>% unique()
    
    ## Now pipe the table into the function
    # example: generate basic plot, dress it up a bit with extra ggplot code, and save
    p <- raster_scatter_cols(raster1   = climate_raster_stack[[clim_ras]],
                             raster2   = context_raser_stack[[cont_ras]],
                             rastercol = context_raser_stack[[col_ras]],
                             
                             mar,
                             xsize, 
                             ysize, 
                             lab_size,
                             out_dir) +
      
      labs(x     = clim_ras, # adding axis labels
           y     = cont_ras,
           color = col_ras) +
      
      theme(plot.margin     = unit(c(mar, mar, mar, mar), "cm"),
            # plot.title      = element_text(vjust = 5, size = tsize, face = "bold"),
            axis.text.x     = element_text(size = xsize),
            axis.title.x    = element_text(size = xsize, face = "bold"),
            
            axis.text.y     = element_text(size = ysize),
            axis.title.y    = element_text(size = ysize, face = "bold"),
            legend.text     = element_text(size = lab_size, hjust = 0.5, color = "black"))
    
    
    png(paste0(out_dir, clim_ras, '_VS_',  cont_ras, '.png'),
        16, 10, units = 'in', res = 600)
    plot(p)
    dev.off()
    
    gc()
    
  } 
}





#########################################################################################################################
####################################################  TBC  ##############################################################
#########################################################################################################################